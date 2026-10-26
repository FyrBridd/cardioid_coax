"""Banded BEM runs over the observation sphere, with a keyed on-disk result cache.

Each frequency is solved on the mesh of the smallest band edge at or above it, with
element size h = min(h_max, C / (edge * per_wavelength)). Fields are evaluated at
reference_point + r * dirs over a fixed observation set (horizontal and vertical
rings, listening-window extras and a sphere grid for power averages).

Results are stored in the engineering e^{+jωt} convention: the complex conjugate of
solve_band's e^{-iωt} output.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from speaker.constants import C
from speaker.drivers import Driver
from speaker.layout import compute_layout
from speaker.params import Params

from .mesh import cabinet_mesh, puck_mesh
from .solver import SurfaceMesh, solve_band

BAND_EDGES = (250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0, 12000.0)  # Hz
OBS_R = 2.0  # m, observation radius for scores
AXIS_R = 1.0  # m, distance of the extra on-axis point evaluated by puck_field
PUCK_F_RANGE = (500.0, 12000.0)  # Hz
PUCK_H_MAX = 0.03  # m
PUCK_H_COAX = 0.010  # m, coax-region size at and below PUCK_HI_BAND
PUCK_HI_BAND = 4000.0  # Hz, bands above this use h_coax = h
LISTENER_AZ = 23.0  # deg, horizontal angle of the listening seat
LW_EXTRA_EL = 10.0  # deg, listening-window extras at (LISTENER_AZ, ±el)
_REPO = Path(__file__).resolve().parents[3]


def freq_grid(f_lo: float = 50.0, f_hi: float = 20000.0, ppo: int = 12) -> np.ndarray:
    """Log-spaced frequencies, ppo per octave, from f_lo up to at most f_hi."""
    n = int(np.floor(ppo * np.log2(f_hi / f_lo) + 1e-9)) + 1
    return f_lo * 2.0 ** (np.arange(n) / ppo)


@dataclass(frozen=True)
class Observation:
    dirs: np.ndarray  # (nd, 3) unit vectors, de-duplicated; rings and grids index into it
    r: float  # m
    h_angles: np.ndarray  # deg, -180..175 step 5: (sin θ, -cos θ, 0)
    h_idx: np.ndarray
    v_angles: np.ndarray  # deg, -180..175 step 5: (0, -cos φ, sin φ), positive up
    v_idx: np.ndarray
    lw_extra_idx: tuple[int, int]  # (az 23°, el +10°), (az 23°, el -10°)
    at23_idx: int  # horizontal +23° (the listener), not on the 5° ring nor in the sphere grid
    sphere_idx: np.ndarray
    sphere_w: np.ndarray  # solid-angle weights over sphere_idx, summing to 1


def _sphere_dir(az: float, el: float) -> np.ndarray:
    a, e = np.radians(az), np.radians(el)
    return np.array([np.sin(a) * np.cos(e), -np.cos(a) * np.cos(e), np.sin(e)])


def observation_set(r: float = OBS_R) -> Observation:
    """Rings, the listener direction, listening-window extras and a 10° sphere grid.

    Coincident directions are stored once and shared by index.
    """
    dirs: list[np.ndarray] = []
    index: dict[tuple[int, ...], int] = {}

    def add(v: np.ndarray) -> int:
        key = tuple(np.rint(v * 1e9).astype(np.int64))
        if key not in index:
            index[key] = len(dirs)
            dirs.append(v)
        return index[key]

    ring = np.arange(-180, 180, 5)
    h_idx = np.array([add(_sphere_dir(a, 0.0)) for a in ring])
    v_idx = np.array([add(np.array([0.0, -np.cos(t), np.sin(t)])) for t in np.radians(ring)])
    lw = (add(_sphere_dir(LISTENER_AZ, LW_EXTRA_EL)), add(_sphere_dir(LISTENER_AZ, -LW_EXTRA_EL)))
    at23 = add(_sphere_dir(LISTENER_AZ, 0.0))

    # Rings of 36 points at el -80..80 each own the band el ± 5°, an area ∝ cos(el);
    # the poles own the caps above ±85°. The areas sum to 4π exactly.
    step = np.radians(10.0)
    sphere, w = [], []
    for e in np.arange(-80, 81, 10):
        band = 2 * np.pi * (np.sin(np.radians(e) + step / 2) - np.sin(np.radians(e) - step / 2)) / 36
        for a in np.arange(0, 360, 10):
            sphere.append(add(_sphere_dir(a, e)))
            w.append(band)
    cap = 2 * np.pi * (1 - np.sin(np.pi / 2 - step / 2))
    for e in (-90.0, 90.0):
        sphere.append(add(_sphere_dir(0.0, e)))
        w.append(cap)
    w = np.array(w)
    return Observation(
        dirs=np.array(dirs), r=float(r),
        h_angles=ring.astype(float), h_idx=h_idx, v_angles=ring.astype(float), v_idx=v_idx,
        lw_extra_idx=lw, at23_idx=at23, sphere_idx=np.array(sphere), sphere_w=w / w.sum(),
    )


@dataclass
class SourceField:
    """Radiated field of each source over an observation set.

    p: pressure per unit volume velocity (Pa per m³/s), (nf, nsrc, nd) complex,
    e^{+jωt} convention (conjugate of solve_band).
    meta["bands"]: per band edge, the mesh size (nodes, triangles) and wall times.
    """

    f: np.ndarray
    names: tuple[str, ...]
    obs: Observation
    p: np.ndarray
    meta: dict


def cabinet_field(
    p: Params,
    coax: Driver,
    bass: Driver,
    f: np.ndarray,
    *,
    params_text: str,
    per_wavelength: float = 6.0,
    h_max: float = 0.05,
    cache_dir: Path | None = Path("out/cache/bem"),
) -> SourceField:
    """Field of every cabinet source at OBS_R about the coax centre on the front face."""
    obs = observation_set(OBS_R)
    cx, cz = compute_layout(p, coax, bass).coax_centre
    points = (np.array([cx, 0.0, cz]) + obs.r * obs.dirs).T
    pr, names, meta = _banded_field(
        "cabinet_field",
        lambda edge, h: cabinet_mesh(p, coax, bass, h=h),
        f, points, (params_text, repr(coax), repr(bass)),
        per_wavelength=per_wavelength, h_max=h_max, cache_dir=cache_dir,
    )
    return SourceField(np.asarray(f, dtype=float), names, obs, pr, meta)


def puck_field(
    p: Params,
    coax: Driver,
    f: np.ndarray,
    *,
    params_text: str,
    per_wavelength: float = 6.0,
    cache_dir: Path | None = Path("out/cache/bem"),
) -> SourceField:
    """Field of the coax alone on a puck, 500 Hz to 12 kHz, about the coax centre on the front face.

    meta["axis_1m"] is the (nf, nsrc) pressure 1 m on axis, in the same convention as p.
    """
    f = np.asarray(f, dtype=float)
    lo, hi = PUCK_F_RANGE
    if f.size and (f.min() < lo or f.max() > hi):
        raise ValueError(
            f"puck_field covers {lo:g} Hz to {hi:g} Hz, got {f.min():g} Hz to {f.max():g} Hz"
        )
    obs = observation_set(OBS_R)
    centre = np.array([p.cabinet.width / 2, 0.0, p.cabinet.coax_centre_z])  # puck_mesh's frame
    points = np.concatenate([obs.r * obs.dirs, [[0.0, -AXIS_R, 0.0]]]).T + centre[:, None]

    def build(edge: float, h: float) -> SurfaceMesh:
        # Above 4 kHz the coax region needs no extra refinement; below it, 10 mm capped at h.
        h_coax = h if edge > PUCK_HI_BAND else min(h, PUCK_H_COAX)
        return puck_mesh(p, coax, h=h, h_coax=h_coax)

    pr, names, meta = _banded_field(
        "puck_field", build, f, points, (params_text, repr(coax)),
        per_wavelength=per_wavelength, h_max=PUCK_H_MAX, cache_dir=cache_dir,
    )
    return SourceField(f, names, obs, pr[..., :-1], {**meta, "axis_1m": pr[..., -1]})


def _banded_field(
    kind: str,
    build: Callable[[float, float], SurfaceMesh],
    f: np.ndarray,
    points: np.ndarray,
    key_parts: Sequence[str],
    *,
    per_wavelength: float,
    h_max: float,
    cache_dir: Path | None,
) -> tuple[np.ndarray, tuple[str, ...], dict]:
    """Solve f band by band on build(edge, h) meshes; e^{+jωt} pressure at points (3, np).

    Returns (p (nf, nsrc, np), source names, meta), from the cache when the key matches.
    """
    f = np.asarray(f, dtype=float)
    points = np.ascontiguousarray(points, dtype=float)
    if f.size == 0:
        raise ValueError("no frequencies to solve")
    band = np.searchsorted(BAND_EDGES, f)  # smallest edge >= f
    if band.max() >= len(BAND_EDGES):
        raise ValueError(f"{f.max():g} Hz is above the last band edge, {BAND_EDGES[-1]:g} Hz")

    path = None
    if cache_dir is not None:
        key = _cache_key(kind, key_parts, f, points, per_wavelength, h_max)
        path = Path(cache_dir) / f"{kind}_{key}.npz"
        if path.exists():
            with np.load(path, allow_pickle=False) as z:
                return z["p"], tuple(str(n) for n in z["names"]), json.loads(str(z["meta"]))

    out, names, bands = None, None, []
    for b in np.unique(band):
        edge = BAND_EDGES[b]
        h = min(h_max, C / (edge * per_wavelength))
        sel = np.flatnonzero(band == b)
        t0 = time.perf_counter()
        mesh = build(edge, h)
        t1 = time.perf_counter()
        res = solve_band(mesh, f[sel], points)
        t2 = time.perf_counter()
        if out is None:
            names = mesh.tag_names
            out = np.empty((len(f), len(names), points.shape[1]), dtype=complex)
        elif mesh.tag_names != names:
            raise ValueError(f"band {edge:g} Hz mesh has sources {mesh.tag_names}, expected {names}")
        out[sel] = np.conj(res)  # e^{-iωt} -> e^{+jωt}
        bands.append({
            "edge": edge, "h": h, "freqs": len(sel),
            "nodes": len(mesh.vertices), "triangles": len(mesh.triangles),
            "mesh_s": t1 - t0, "solve_s": t2 - t1,
        })
    meta = {"bands": bands}

    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
        with open(tmp, "wb") as fh:
            np.savez(fh, p=out, names=np.array(names), meta=np.array(json.dumps(meta)))
        os.replace(tmp, path)
    return out, names, meta


def _cache_key(
    kind: str, key_parts: Sequence[str], f: np.ndarray, points: np.ndarray, per_wavelength: float, h_max: float
) -> str:
    """sha256 over everything the result depends on, including the code version."""
    h = hashlib.sha256()
    for part in (kind, *key_parts, repr(float(per_wavelength)), repr(float(h_max)), _code_version()):
        h.update(part.encode())
        h.update(b"\0")
    h.update(f.astype(np.float64).tobytes())
    h.update(points.astype(np.float64).tobytes())
    return h.hexdigest()


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=_REPO, check=True, capture_output=True, text=True).stdout


def _code_version() -> str:
    """git HEAD, plus a hash of the uncommitted changes under src/ (tracked diff and untracked files)."""
    head = _git("rev-parse", "HEAD").strip()
    if not _git("status", "--porcelain", "--", "src/").strip():
        return head
    h = hashlib.sha256(_git("diff", "HEAD", "--", "src/").encode())
    for name in sorted(_git("ls-files", "--others", "--exclude-standard", "-z", "--", "src/").split("\0")):
        if name:
            h.update(name.encode() + b"\0" + (_REPO / name).read_bytes())
    return f"{head}+{h.hexdigest()}"
