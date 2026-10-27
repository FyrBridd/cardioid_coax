"""Closed-form radiation solutions used to validate the BEM solver.

Time convention: e^{-iωt} (physics), so outgoing waves go as e^{+ikr}. This
matches `speaker.bem.solver`. Every function returns complex pressure (Pa) per
unit volume velocity (m^3/s) of the vibrating surface.
"""
from __future__ import annotations

import numpy as np
from scipy.special import eval_legendre, spherical_jn, spherical_yn

from speaker.constants import C, RHO


def _k(f: float) -> float:
    return 2 * np.pi * f / C


def pulsating_sphere(f: float, a: float, r: float) -> complex:
    """Sphere of radius a pulsating uniformly; field at radius r >= a."""
    k = _k(f)
    v = 1.0 / (4 * np.pi * a**2)  # Q = 4πa²v = 1
    return complex(RHO * C * (-1j * k * a) / (1 - 1j * k * a) * a / r * np.exp(1j * k * (r - a)) * v)


def _h1(n: np.ndarray, x: float, derivative: bool = False) -> np.ndarray:
    """Spherical Hankel function of the first kind, h_n = j_n + i y_n."""
    return spherical_jn(n, x, derivative) + 1j * spherical_yn(n, x, derivative)


def cap_on_sphere(f: float, a: float, alpha: float, r: float, theta: float, n_terms: int = 60) -> complex:
    """Rigid sphere of radius a with a polar cap of half-angle alpha vibrating radially (Morse).

    theta is measured from the cap axis; the field point is at radius r >= a.
    """
    k = _k(f)
    v = 1.0 / (2 * np.pi * a**2 * (1 - np.cos(alpha)))  # Q = 2πa²(1-cos α)v = 1
    n = np.arange(n_terms)
    x = np.cos(alpha)
    c = np.empty(n_terms)
    c[0] = (1 - x) / 2
    c[1:] = (eval_legendre(n[1:] - 1, x) - eval_legendre(n[1:] + 1, x)) / 2
    with np.errstate(over="ignore", invalid="ignore"):
        ratio = _h1(n, k * r) / _h1(n, k * a, derivative=True)
    # High orders at small ka overflow y_n; those terms are vanishingly small.
    ratio = np.where(np.isfinite(ratio), ratio, 0.0)
    return complex(1j * RHO * C * v * np.sum(c * ratio * eval_legendre(n, np.cos(theta))))


def baffled_piston_on_axis(f: float, a: float, r: float) -> complex:
    """Circular piston of radius a in an infinite rigid baffle; on-axis field at distance r."""
    k = _k(f)
    v = 1.0 / (np.pi * a**2)  # Q = πa²v = 1
    return complex(RHO * C * v * (np.exp(1j * k * r) - np.exp(1j * k * np.sqrt(r**2 + a**2))))
