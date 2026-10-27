"""Sealed-box bass model: electro-mechanical circuit, half-space radiation."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from speaker.constants import C, RHO
from speaker.drivers import Driver

P_REF = 20e-6


@dataclass(frozen=True)
class Alignment:
    fc: float
    qtc: float
    f3: float


def _check_f(f) -> np.ndarray:
    f = np.atleast_1d(np.asarray(f, dtype=float))
    if np.any(f <= 0):
        raise ValueError("frequencies must be > 0")
    return f


def sealed_alignment(d: Driver, vb: float, series_r: float = 0.0) -> Alignment:
    qes = d.qes * (d.re + series_r) / d.re
    qts = qes * d.qms / (qes + d.qms)
    ratio = math.sqrt(1.0 + d.vas / vb)
    fc = d.fs * ratio
    qtc = qts * ratio
    # -3 dB of 2nd-order high-pass: x^4 / ((1-x^2)^2 + x^2/Q^2) = 0.5
    def g(x: float) -> float:
        return x**4 / ((1 - x**2) ** 2 + x**2 / qtc**2) - 0.5

    x3 = brentq(g, 1e-3, 100.0)
    return Alignment(fc=fc, qtc=qtc, f3=x3 * fc)


def cone_velocity(
    d: Driver, vb: float, f: np.ndarray, v_rms: float, series_r: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """(ω, RMS cone velocity in m/s, e^{+jωt}) of one driver in a sealed box vb at v_rms."""
    w = 2 * np.pi * f
    cab = vb / (RHO * C**2)
    zm = 1j * w * d.mms + d.rms + 1.0 / (1j * w * d.cms) + d.sd**2 / (1j * w * cab)
    ze = d.re + series_r + 1j * w * d.le
    u = d.bl * v_rms / (ze * zm + d.bl**2)
    return w, u


def sealed_spl(d, vb, f, v_rms, r=1.0, n_drivers=1, series_r=0.0) -> np.ndarray:
    f = _check_f(f)
    w, u = cone_velocity(d, vb, f, v_rms, series_r)
    p = n_drivers * RHO * w * d.sd * np.abs(u) / (2 * np.pi * r)
    return 20 * np.log10(p / P_REF)


def sealed_excursion(d, vb, f, v_rms, series_r=0.0) -> np.ndarray:
    f = _check_f(f)
    w, u = cone_velocity(d, vb, f, v_rms, series_r)
    return math.sqrt(2) * np.abs(u) / w


def max_spl(d, vb, f, v_rms_limit, r=1.0, n_drivers=1, series_r=0.0) -> np.ndarray:
    f = _check_f(f)
    w = 2 * np.pi * f
    p_x = n_drivers * RHO * w**2 * d.sd * d.xmax / (2 * np.pi * r) / math.sqrt(2)
    out = 20 * np.log10(p_x / P_REF)
    if v_rms_limit is not None:
        out = np.minimum(out, sealed_spl(d, vb, f, v_rms_limit, r, n_drivers, series_r))
    return out
