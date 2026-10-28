"""Per-variant tuning of the vent specific resistance and the coax chamber's net volume.

Neither changes a BEM field (the fillers are inside the box and the cache is keyed on the exterior mesh),
so each variant picks both from a grid, maximising the worst cardioid margin over the felt spread:

    margin(R, V) = min over s in {0.7, 1.0, 1.3} of min(hard_s - 10, aim_s - 6)

hard and aim are scores.py's cardioid_hard and cardioid_aim: the min over 150-300 Hz and 300-500 Hz of
SPL(H 0°) - SPL(H 180°) of the coax section.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from speaker.bem.run import SourceField
from speaker.drivers import Driver
from speaker.layout import compute_layout, cardioid_params
from speaker.params import Params
from speaker.scores import _spl
from speaker.system import coax_section, vent_scaled_section, way_responses

R_GRID = np.geomspace(80.0, 1000.0, 15)  # Pa*s/m
V_STEP = 0.5e-3  # m^3
V_MIN = 4e-3  # m^3
CARD_BAND = (150.0, 500.0)
HARD_BAND = (150.0, 300.0)
AIM_BAND = (300.0, 500.0)
HARD_MIN, AIM_MIN = 10.0, 6.0  # dB, the pass lines
SPREAD = (0.7, 1.0, 1.3)  # felt-factor scales on the vent resistance


@dataclass
class VentTune:
    resistance: float  # Pa*s/m
    net_volume: float  # m^3
    hard: float  # dB, at the felt factor 1.0
    aim: float
    margin: float  # dB, worst over the spread
    spread: dict[str, tuple[float, float]]  # "0.7" and "1.3": (hard, aim)
    base: dict[str, float]  # hard, aim, margin at the params' own R and V


def with_vents(p: Params, resistance: float, net_volume: float) -> Params:
    """p with the vent specific resistance and coax chamber net volume replaced."""
    return replace(p, vents=replace(p.vents, specific_resistance=resistance),
                   cabinet=replace(p.cabinet, coax_chamber_net_volume=net_volume))


def feasible_volumes(p: Params, coax: Driver, bass: Driver) -> np.ndarray:
    """V_MIN up to the empty chamber in V_STEP steps, plus the params' own value, where compute_layout succeeds."""
    empty = compute_layout(p, coax, bass).coax_chamber_empty
    n = int(np.floor((empty - V_MIN) / V_STEP + 1e-9)) + 1
    cand = {*(V_MIN + V_STEP * np.arange(max(n, 0))), p.cabinet.coax_chamber_net_volume}
    out = []
    for v in sorted(cand):
        try:
            compute_layout(with_vents(p, p.vents.specific_resistance, v), coax, bass)
        except ValueError:
            continue
        out.append(v)
    return np.array(out)


def _rejections(section: np.ndarray, cab: SourceField, f: np.ndarray) -> tuple[float, float]:
    lv = _spl(section)
    h = {int(a): int(i) for a, i in zip(cab.obs.h_angles, cab.obs.h_idx)}
    d = lv[:, h[0]] - lv[:, h[-180]]  # 180 deg is the -180 ring entry, as scores.py
    hard = float(d[(f >= HARD_BAND[0]) & (f <= HARD_BAND[1])].min())
    aim = float(d[(f >= AIM_BAND[0]) & (f <= AIM_BAND[1])].min())
    return hard, aim


def _margin(hard: float, aim: float) -> float:
    return min(hard - HARD_MIN, aim - AIM_MIN)


def tune_vents(p: Params, coax: Driver, bass: Driver, cab: SourceField, *,
               lf_recess: tuple[np.ndarray, np.ndarray] | None = None) -> VentTune:
    """The R in R_GRID plus the params' own, and the V in feasible_volumes, with the largest worst-case margin."""
    lo, hi = CARD_BAND
    if cab.f[0] > lo * (1 + 1e-9) or cab.f[-1] < hi * (1 - 1e-9):
        raise ValueError(f"cabinet field covers {cab.f[0]:g}-{cab.f[-1]:g} Hz, need {lo:g}-{hi:g} Hz")
    f = cab.f[(cab.f >= lo * (1 - 1e-9)) & (cab.f <= hi * (1 + 1e-9))]
    resistances = np.unique(np.append(R_GRID, p.vents.specific_resistance))

    def evaluate(r: float, v: float):
        q = with_vents(p, r, v)
        lay = compute_layout(q, coax, bass)
        ways = way_responses(q, coax, bass, lay, cab, None, f, trusted_limit=float(f[-1]), lf_recess=lf_recess)
        cp = cardioid_params(q, lay, coax)
        res = {s: _rejections(coax_section(ways) if s == 1.0 else vent_scaled_section(ways, cab, coax, cp, s),
                              cab, f) for s in SPREAD}
        return res, min(_margin(*x) for x in res.values())

    best = None
    for v in feasible_volumes(p, coax, bass):
        for r in resistances:
            res, m = evaluate(float(r), float(v))
            if best is None or m > best[0]:
                best = (m, float(r), float(v), res)
    m, r, v, res = best
    base_res, base_m = evaluate(p.vents.specific_resistance, p.cabinet.coax_chamber_net_volume)
    return VentTune(
        resistance=r, net_volume=v, hard=res[1.0][0], aim=res[1.0][1], margin=m,
        spread={"0.7": res[0.7], "1.3": res[1.3]},
        base={"hard": base_res[1.0][0], "aim": base_res[1.0][1], "margin": base_m},
    )
