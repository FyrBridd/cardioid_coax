"""Exterior Helmholtz BEM: pressure radiated by vibrating patches of a rigid closed surface.

Burton–Miller formulation with P1 (piecewise-linear, continuous) elements,
dense assembly by bempp-cl and a direct LU solve.

Time convention: e^{-iωt} (physics), so outgoing waves go as e^{+ikr} and
∂p/∂n = iωρ v_n. Converting to the engineering e^{+jωt} convention is the
caller's job; this module does not do it.

The mesh must be closed with outward normals (counter-clockwise triangles seen
from outside); v_n is positive outward.
"""
from __future__ import annotations

import warnings
from collections.abc import Callable
from dataclasses import dataclass

import bempp_cl.api as bempp
import numpy as np
from bempp_cl.api.operators.boundary import helmholtz, sparse
from bempp_cl.api.operators.potential import helmholtz as helmholtz_potential
from scipy.linalg import lu_factor, lu_solve

from speaker.constants import C, RHO

BEM_MEMORY_BUDGET = 22e9  # bytes


@dataclass(frozen=True)
class SurfaceMesh:
    vertices: np.ndarray  # (nv, 3), m
    triangles: np.ndarray  # (nt, 3), vertex indices
    tags: np.ndarray  # (nt,), 0 = rigid, i >= 1 = source tag_names[i-1]
    tag_names: tuple[str, ...]


def estimate_bytes(n_dofs: int) -> float:
    """Peak memory of one frequency's solve: a few dense complex n×n arrays."""
    return 64.0 * n_dofs**2


_numba_warned = False


def _select_backend() -> None:
    global _numba_warned
    if bempp.CPU_OPENCL_DRIVER_FOUND:
        bempp.DEFAULT_DEVICE_INTERFACE = "opencl"
        return
    bempp.DEFAULT_DEVICE_INTERFACE = "numba"
    if not _numba_warned:
        warnings.warn("no OpenCL CPU platform found; BEM falls back to numba (about 9x slower)", stacklevel=3)
        _numba_warned = True


def _neumann_data(space, mass, tags: np.ndarray, n_src: int) -> np.ndarray:
    """P1 node velocities per source, each normalised to unit volume velocity. (n_dofs, n_src)."""
    l2g = space.local2global.ravel()
    counts = np.bincount(l2g, minlength=space.global_dof_count)
    g = np.empty((space.global_dof_count, n_src))
    for s in range(n_src):
        on = np.repeat((tags == s + 1).astype(float), 3)
        g[:, s] = np.bincount(l2g, weights=on, minlength=space.global_dof_count) / counts
    q = np.ones(space.global_dof_count) @ (mass @ g)  # Q = ∫ v dA
    if np.any(q <= 0):
        raise ValueError(f"source tags with no area: {[i + 1 for i in np.flatnonzero(q <= 0)]}")
    return g / q


def _assemble(space, k: float, mass) -> tuple[np.ndarray, Callable[[np.ndarray], np.ndarray]]:
    """Burton–Miller system at wavenumber k: dense A, and a function applying B to Neumann data.

    A = 0.5*I - K + eta*W and B = -V - eta*(Kp + 0.5*I), with eta = i/k (bempp sign
    convention, P1/P1/P1). B is never formed; its terms are applied one at a time.
    """
    eta = 1j / k
    mass_coo = mass.tocoo()
    a = helmholtz.double_layer(space, space, space, k).weak_form().to_dense()
    np.negative(a, out=a)
    a[mass_coo.row, mass_coo.col] += 0.5 * mass_coo.data
    tmp = helmholtz.hypersingular(space, space, space, k).weak_form().to_dense()
    tmp *= eta
    a += tmp
    del tmp

    def apply_b(neumann: np.ndarray) -> np.ndarray:
        rhs = -(helmholtz.single_layer(space, space, space, k).weak_form().to_dense() @ neumann)
        rhs -= eta * (helmholtz.adjoint_double_layer(space, space, space, k).weak_form().to_dense() @ neumann)
        rhs -= eta * (mass @ neumann) * 0.5
        return rhs

    return a, apply_b


def solve_band(
    mesh: SurfaceMesh,
    freqs: np.ndarray,
    points: np.ndarray,
    *,
    memory_budget: float = BEM_MEMORY_BUDGET,
) -> np.ndarray:
    """Pressure at points (3, np) per unit volume velocity of each source, others rigid.

    Returns complex (nf, nsrc, np), e^{-iωt} convention.
    """
    used, tri = np.unique(mesh.triangles, return_inverse=True)
    n_dofs = len(used)
    need = estimate_bytes(n_dofs)
    if need > memory_budget:
        raise ValueError(
            f"BEM solve with {n_dofs} DOFs needs about {need / 1e9:.1f} GB of memory, "
            f"over the {memory_budget / 1e9:.1f} GB budget"
        )
    tags = np.asarray(mesh.tags)
    n_src = len(mesh.tag_names)
    if tags.shape != (len(mesh.triangles),) or tags.min() < 0 or tags.max() > n_src:
        raise ValueError("tags must be one integer in [0, len(tag_names)] per triangle")
    points = np.asarray(points, dtype=float)
    freqs = np.asarray(freqs, dtype=float)

    _select_backend()
    out = np.empty((len(freqs), n_src, points.shape[1]), dtype=complex)
    with warnings.catch_warnings():
        # pocl's compiler chatter on first kernel build; harmless.
        warnings.filterwarnings("ignore", message="Non-empty compiler output")
        grid = bempp.Grid(
            np.ascontiguousarray(mesh.vertices[used].T, dtype=float),
            np.ascontiguousarray(tri.reshape(-1, 3).T, dtype=np.uint32),
        )
        space = bempp.function_space(grid, "P", 1)
        ident = sparse.identity(space, space, space)
        mass = ident.weak_form().to_sparse()  # P1 mass matrix
        g = _neumann_data(space, mass, tags, n_src)

        for i, f in enumerate(freqs):
            w = 2 * np.pi * f
            k = w / C
            neumann = 1j * RHO * w * g  # ∂p/∂n = iωρ v_n
            a, apply_b = _assemble(space, k, mass)
            rhs = apply_b(neumann)

            # a is C-ordered, so a.T is Fortran-ordered and factors in place; trans=1 solves a x = rhs.
            lu = lu_factor(a.T, overwrite_a=True, check_finite=False)
            del a
            phi = lu_solve(lu, rhs, trans=1, check_finite=False)
            del lu

            dlp = helmholtz_potential.double_layer(space, points, k)
            slp = helmholtz_potential.single_layer(space, points, k)
            for s in range(n_src):
                p = dlp * bempp.GridFunction(space, coefficients=phi[:, s]) - slp * bempp.GridFunction(
                    space, coefficients=neumann[:, s]
                )
                out[i, s] = p[0]
    return out
