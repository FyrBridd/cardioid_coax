"""Lumped passive-cardioid coax model.

The cone drives a rear chamber that is vented through resistive side slots.
Cone and vent volume velocities are combined as point sources, with the vents
displaced behind the cone, to give a free-field far-field pattern.

This is a LOW-FREQUENCY ESTIMATE ONLY. It ignores the baffle and the cabinet,
so above about 400 Hz the front/back figure it gives is pessimistic. Plan 2's
BEM replaces it there.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.optimize import minimize_scalar

from speaker.constants import C, RHO
from speaker.drivers import Driver


@dataclass(frozen=True)
class CardioidParams:
    chamber_volume: float  # m^3
    vent_area: float  # m^2, total over all slots
    vent_length: float  # m
    vent_resistance: float  # specific resistance, Pa*s/m (0 or inf allowed)
    vent_positions: tuple[tuple[float, float, float], ...]  # m, x forward
    n_slots: int = 12


@dataclass(frozen=True)
class CardioidResult:
    f: np.ndarray
    u_cone: np.ndarray
    u_vent: np.ndarray
    excursion: np.ndarray  # peak cone displacement, m


def _check_f(f) -> np.ndarray:
    f = np.atleast_1d(np.asarray(f, dtype=float))
    if np.any(f <= 0):
        raise ValueError("frequencies must be > 0")
    return f


def simulate(d: Driver, p: CardioidParams, f: np.ndarray, v_rms: float) -> CardioidResult:
    f = _check_f(f)
    w = 2 * np.pi * f
    cab = p.chamber_volume / (RHO * C**2)
    a = math.sqrt(p.vent_area / p.n_slots / math.pi)
    mv = RHO * (p.vent_length + 1.7 * a) / p.vent_area
    if math.isinf(p.vent_resistance):
        yv = np.zeros_like(w, dtype=complex)
    else:
        rv = p.vent_resistance / p.vent_area
        yv = 1.0 / (rv + 1j * w * mv)
    zbox = 1.0 / (1j * w * cab + yv)
    ze = d.re + 1j * w * d.le
    zm = 1j * w * d.mms + d.rms + 1.0 / (1j * w * d.cms) + d.sd**2 * zbox + d.bl**2 / ze
    u = d.bl * v_rms / (ze * zm)
    uc = d.sd * u
    uv = -uc * zbox * yv
    excursion = math.sqrt(2) * np.abs(u) / w
    return CardioidResult(f, uc, uv, excursion)


def polar(
    res: CardioidResult,
    p: CardioidParams,
    angles_deg: np.ndarray,
    plane: Literal["horizontal", "vertical"] = "horizontal",
) -> np.ndarray:
    th = np.radians(np.asarray(angles_deg, dtype=float))
    if plane == "horizontal":
        rhat = np.stack([np.cos(th), np.sin(th), np.zeros_like(th)], axis=1)
    elif plane == "vertical":
        rhat = np.stack([np.cos(th), np.zeros_like(th), np.sin(th)], axis=1)
    else:
        raise ValueError("plane must be 'horizontal' or 'vertical'")
    pos = np.asarray(p.vent_positions, dtype=float).reshape(-1, 3)
    n = len(pos)
    k = 2 * np.pi * res.f / C  # (nf,)
    proj = rhat @ pos.T  # (na, nv)
    phase = np.exp(1j * k[:, None, None] * proj[None, :, :])  # (nf, na, nv)
    total = res.u_cone[:, None] + res.u_vent[:, None] / n * phase.sum(axis=2)
    mag = np.abs(total)
    # Normalise to the 0 degree level at each frequency.
    front = np.abs(res.u_cone + res.u_vent / n * np.exp(1j * k[:, None] * (pos[:, 0])[None, :]).sum(axis=1))
    with np.errstate(divide="ignore"):
        return 20 * np.log10(mag / front[:, None])


def _level(res: CardioidResult, p: CardioidParams, deg: float) -> np.ndarray:
    th = math.radians(deg)
    rhat = np.array([math.cos(th), math.sin(th), 0.0])
    pos = np.asarray(p.vent_positions, dtype=float).reshape(-1, 3)
    k = 2 * np.pi * res.f / C
    ph = np.exp(1j * k[:, None] * (pos @ rhat)[None, :])
    return np.abs(res.u_cone + res.u_vent / len(pos) * ph.sum(axis=1))


def front_back(res: CardioidResult, p: CardioidParams) -> np.ndarray:
    return 20 * np.log10(_level(res, p, 0.0) / _level(res, p, 180.0))


def tune_resistance(d: Driver, p: CardioidParams, f_lo: float = 150.0, f_hi: float = 300.0) -> float:
    f = np.geomspace(f_lo, f_hi, 40)

    def score(r: float) -> float:
        q = CardioidParams(p.chamber_volume, p.vent_area, p.vent_length, r,
                           p.vent_positions, p.n_slots)
        return float(np.min(front_back(simulate(d, q, f, 2.83), q)))

    grid = np.geomspace(10.0, 5000.0, 60)
    scores = [score(r) for r in grid]
    i = int(np.argmax(scores))
    lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    res = minimize_scalar(lambda r: -score(r), bounds=(lo, hi), method="bounded")
    return float(res.x) if -res.fun >= scores[i] else float(grid[i])
