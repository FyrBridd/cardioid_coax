"""The passive crossover topology: which parts exist for a set of options, the netlist they form in passive
and semi-active modes, the driver loads, and each way's terminal voltage per input volt; and the fit of the
part values at one pair of crossover points, snapped to standard values; the search over crossover points,
polarities and options, the part-tolerance check, and the series-woofer fallback.

Passive: a bass low-pass on `in`, and a coax high-pass from `in` to node `cx`, which feeds the coax LF
(low-pass, optional notch, L-pad) and the coax HF (high-pass, L-pad, optional series RLC shunt).
Semi-active: the coax network alone, with `cx` being `in` (the DSP supplies the split, the bass runs active).
"""
from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field, replace
from typing import Literal

import numpy as np
from scipy.optimize import least_squares

from speaker.crossover import MATCH_BANDS, match_gains
from speaker.network import Part, inductor_dcr, neighbours, snap, solve
from speaker.scores import Score, score_all
from speaker.spinorama import spinorama
from speaker.system import Ways, coax_section
from speaker.xo_fit import BASS_COAX_GRID, LF_HF_GRID, SEARCH_KEYS, pick

Mode = Literal["passive", "semi_active"]


@dataclass(frozen=True)
class Topology:
    bass_order: int = 2  # 2 or 3
    zobel: bool = False
    lf_order: int = 2  # 2 or 3
    lf_notch: bool = False
    lf_shelf: bool = False  # the coax LF contour: R_lc with L_lc across it, in series before the L-pad
    hf_shunt: bool = False
    bass_inverted: bool = False
    hf_inverted: bool = False
    bass_wiring: Literal["parallel", "series"] = "parallel"


_PREFIX_BOUNDS = {"L": (0.05e-3, 20e-3), "C": (0.5e-6, 400e-6), "R": (0.1, 50.0)}  # H, F, ohm


def value_names(topo: Topology, mode: Mode) -> list[str]:
    """Names of the fitted values, in a fixed order; semi-active has no bass and no coax high-pass."""
    names: list[str] = []
    if mode == "passive":
        names += ["L_b1", "C_b1"] + (["L_b2"] if topo.bass_order == 3 else [])
        names += ["R_bz", "C_bz"] if topo.zobel else []
        names += ["C_h1", "L_h1"]
    names += ["R_ls", "R_lp", "L_l1", "C_l1"] + (["L_l2"] if topo.lf_order == 3 else [])
    names += ["L_ln", "C_ln", "R_ln"] if topo.lf_notch else []
    names += ["R_lc", "L_lc"] if topo.lf_shelf else []
    names += ["C_t1", "L_t1", "C_t2", "R_ts", "R_tp"]
    names += ["R_tr", "L_tr", "C_tr"] if topo.hf_shunt else []
    return names


_ALL = value_names(Topology(bass_order=3, zobel=True, lf_order=3, lf_notch=True, lf_shelf=True, hf_shunt=True),
                   "passive")
BOUNDS: dict[str, tuple[float, float]] = {
    n: ((0.5, 100.0) if n in ("R_lp", "R_tp") else _PREFIX_BOUNDS[n[0]]) for n in _ALL}


def build(values: dict[str, float], topo: Topology, mode: Mode) -> list[Part]:
    """The netlist: unit-source node `in`, loads `bass`, `coax_lf` and `coax_hf` as "Z" parts."""
    parts: list[Part] = []

    def r(name, a, b):
        parts.append(Part(name, "R", a, b, values[name]))

    def c(name, a, b):
        parts.append(Part(name, "C", a, b, values[name]))

    def ind(name, a, b):
        v = values[name]
        parts.append(Part(name, "L", a, b, v, dcr=inductor_dcr(v, cored=name.startswith("L_b"))))

    if mode == "passive":
        ind("L_b1", "in", "b1")
        c("C_b1", "b1", "0")
        end = "b1"
        if topo.bass_order == 3:
            ind("L_b2", "b1", "b2")
            end = "b2"
        if topo.zobel:
            r("R_bz", end, "bz")
            c("C_bz", "bz", "0")
        parts.append(Part("bass", "Z", end, "0", 1.0))
        c("C_h1", "in", "cx")
        ind("L_h1", "cx", "0")
        cx = "cx"
    else:
        cx = "in"

    # Coax LF: low-pass, notch, contour, then the L-pad next to the driver. The contour (a resistor with an inductor
    # across it) sits just before the L-pad, whose input stays near the driver's |Z|: a steady load, so its shelf
    # depth and corner move little with the driver's impedance, and the low-pass ahead still sees it as a load.
    ind("L_l1", cx, "l1")
    c("C_l1", "l1", "0")
    end = "l1"
    if topo.lf_order == 3:
        ind("L_l2", "l1", "l2")
        end = "l2"
    if topo.lf_notch:
        ind("L_ln", end, "ln")
        c("C_ln", end, "ln")
        r("R_ln", end, "ln")
        end = "ln"
    if topo.lf_shelf:
        r("R_lc", end, "lc")
        ind("L_lc", end, "lc")
        end = "lc"
    r("R_ls", end, "lp")
    r("R_lp", "lp", "0")
    parts.append(Part("coax_lf", "Z", "lp", "0", 1.0))

    # Coax HF: high-pass, L-pad, optional series RLC across the driver.
    c("C_t1", cx, "t1")
    ind("L_t1", "t1", "0")
    c("C_t2", "t1", "t2")
    r("R_ts", "t2", "tp")
    r("R_tp", "tp", "0")
    if topo.hf_shunt:
        r("R_tr", "tp", "tr1")
        ind("L_tr", "tr1", "tr2")
        c("C_tr", "tr2", "0")
    parts.append(Part("coax_hf", "Z", "tp", "0", 1.0))
    return parts


def driver_loads(z: dict[str, np.ndarray], topo: Topology) -> dict[str, np.ndarray]:
    """Load impedances by netlist name from `speaker.vituixcad.impedances`; two woofers in series see 4x one pair."""
    bass = z["bass_pair"] * (4.0 if topo.bass_wiring == "series" else 1.0)
    return {"bass": bass, "coax_lf": z["coax_lf"], "coax_hf": z["coax_hf"]}


def bass_scale(topo: Topology) -> float:
    """Series woofers each see half the terminal voltage."""
    return 0.5 if topo.bass_wiring == "series" else 1.0


def transfers(parts: list[Part], f: np.ndarray, loads: dict[str, np.ndarray],
              topo: Topology) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Each way's terminal voltage per input volt (polarity applied), and the impedance the amplifier sees."""
    sol = solve(parts, f, loads)
    flip = {"bass": topo.bass_inverted, "coax_hf": topo.hf_inverted}
    out = {}
    for p in parts:
        if p.kind == "Z":
            v = sol.across(p)
            out[p.name] = -v if flip.get(p.name, False) else v
    return out, sol.z_in


Z_GRID = np.geomspace(20.0, 20000.0, 300)  # Hz, where the impedance floor is checked
FIT_LO = 100.0  # Hz, bottom of the fitted band
FIT_TOP = 10000.0  # Hz, top of the fitted band
ABOVE_WEIGHT = 0.25  # weight of the fit residual above the trusted limit
FIT_SMOOTH = 1 / 3  # octave
SHUNT_Q_MAX = 2.0  # the HF's RLC shunt, sqrt(L/C)/R
Z_PENALTY = 10.0  # residual per ohm under the floor, at each Z_GRID point
Q_PENALTY = 10.0  # residual per unit of shunt Q over SHUNT_Q_MAX
CROSS_WEIGHT = 3.0  # residual per dB between the two ways' 23° levels at their crossover point
BASS_UNDER_DB = 10.0  # dB, how far the bass must sit under the coax LF at 23° from 2 fb to BASS_UNDER_TOP
BASS_UNDER_TOP = 600.0  # Hz, top of the band the cardioid coax section must own (main spec §2 item 4)
LEAK_WEIGHT = 3.0  # residual per dB the bass sits less than BASS_UNDER_DB under the coax LF, at each point
LEVEL_WEIGHT = 3.0  # residual per dB the fitted level sits under the woofers' own level less LEVEL_SLACK_DB
LEVEL_SLACK_DB = 1.0  # dB, how far under the woofers' own level the fit may go free
NOTCH_HZ = 2700.0  # Hz, the coax LF cone breakup
SHELF_DB = 7.5  # dB, the coax LF contour's starting cut above its corner (the LF is ~9 dB down at 150 Hz)
SHELF_CORNER = 1.5  # the contour's starting corner, times fb
SHUNT_HZ = float(np.sqrt(2000.0 * 3300.0))  # Hz, centre of the coax HF's 2-3.3 kHz hump
STOP_DB = -12.0  # dB, the most a branch's |H| may reach in its far stopband, re |H| at the crossover point bounding it
STOP_RATIO = 4.0  # the far stopband starts this factor beyond the crossover point
# residual per dB of |H| over STOP_DB, at each Z_GRID point of a far stopband. On round 3's bass260-580 (whole search,
# 2026-10-30) weights 3 and 10 left the bass 0.7 and 1.8 dB over the rule near 4 fb; 15, 20, 30, 50 and 100 held it
# within 0.3 dB. 20 is one step above the lowest that held, so a small change in the driver data does not tip it over.
STOP_WEIGHT = 20.0
ELEC_GAIN_MAX_DB = 1.0  # dB, the most voltage gain a branch may give its driver on Z_GRID (only a resonance gains)
PEAK_TOL_DB = 1.0  # dB, how far the smoothed 23° sum may sit above the fitted level below FIT_LO (a roll-off is free)
# residual per dB over ELEC_GAIN_MAX_DB at each Z_GRID point of each branch, and per dB over PEAK_TOL_DB at each point
# below FIT_LO. On round 3's bass260-580 (whole search, 2026-10-30) no weight met the cap: the woofers look capacitive
# above their ~40 Hz resonance, and with the 3.5 ohm floor the bass low-pass gains at least ~+1.1 dB near 60 Hz.
# Weights 20-300 left it +2.0 to +3.6 dB (flat23 3.3-7.3 dB); 500-10000 gave +1.1 to +1.7 dB, the sum within 1 dB
# of its level below 100 Hz, but the bass 6-18 dB over its STOP_DB line (every weight checked, 50-10000). 1000 gave
# the best flat23 (1.42 dB; +1.4 dB, stopband 6.3 dB over). The outcome moves a lot with the weight: 500 gave flat23
# 6.06 dB, 2000 1.99 dB.
PEAK_WEIGHT = 1000.0
_EDGE = 1e-6  # keeps start values strictly inside the bounds (log scale)
_LSQ_TOL = 1e-4  # least_squares ftol and xtol: finer log-value steps are lost in the ~20 % snap to standard values


def stopbands(fb: float, fh: float, mode: Mode = "passive") -> list[tuple[str, float, np.ndarray]]:
    """Each branch's far stopbands as (branch, the crossover point bounding it, mask on Z_GRID): the bass from
    STOP_RATIO x fb up; the coax LF from STOP_RATIO x fh up and, in passive mode, from fb / STOP_RATIO down; the coax
    HF from fh / STOP_RATIO down. Semi-active has no bass branch, and the DSP owns the coax's lower edge."""
    out = []
    if mode == "passive":
        out.append(("bass", fb, Z_GRID >= STOP_RATIO * fb))
        out.append(("coax_lf", fb, Z_GRID <= fb / STOP_RATIO))
    out.append(("coax_lf", fh, Z_GRID >= STOP_RATIO * fh))
    out.append(("coax_hf", fh, Z_GRID <= fh / STOP_RATIO))
    return out


def stop_grid(fb: float, fh: float) -> np.ndarray:
    """The frequencies stop_excess_db reads transfers on: Z_GRID, then fb and fh."""
    return np.concatenate([Z_GRID, [fb, fh]])


def stop_excess_db(h: dict[str, np.ndarray], fb: float, fh: float, mode: Mode = "passive") -> np.ndarray:
    """dB by which each branch's |H| sits over STOP_DB re its |H| at the crossover point bounding the band, at every
    Z_GRID point of each of stopbands(fb, fh, mode) in turn (0 where it is under). h: transfers on stop_grid(fb, fh).
    This is the network alone, with the driver loads; the acoustic data play no part (the bass field, say, is zero
    above the trusted limit, where a resonance in its network would not show in the response)."""
    n = len(Z_GRID)
    ref = {fb: n, fh: n + 1}
    out = []
    for k, fc, m in stopbands(fb, fh, mode):
        a = np.abs(h[k])
        db = 20 * np.log10(np.maximum(a[:n][m], 1e-300) / max(float(a[ref[fc]]), 1e-300))
        out.append(np.maximum(0.0, db - STOP_DB))
    return np.concatenate(out)


def gain_excess_db(h: dict[str, np.ndarray]) -> np.ndarray:
    """dB by which each branch's |H| (h: transfers on Z_GRID) sits over ELEC_GAIN_MAX_DB, branch by branch in h's
    order (0 where under). A passive network gives a driver more than its input voltage only by resonating, so this
    holds whatever the acoustic data say (the fitted band starts at FIT_LO, and the bass is zero above tl)."""
    return np.concatenate([np.maximum(0.0, 20 * np.log10(np.maximum(np.abs(x), 1e-300)) - ELEC_GAIN_MAX_DB)
                           for x in h.values()])


def largest_gain_db(h: dict[str, np.ndarray]) -> tuple[float, str, int]:
    """The largest |H| in dB over every branch of h (transfers on any grid), which branch, and at which index."""
    db = {k: 20 * np.log10(np.maximum(np.abs(x), 1e-300)) for k, x in h.items()}
    k = max(db, key=lambda n: float(db[n].max()))
    j = int(np.argmax(db[k]))
    return float(db[k][j]), k, j


def _fit_band(f: np.ndarray, tl: float) -> np.ndarray:
    """The residual's point weights: 1 on [FIT_LO, tl], ABOVE_WEIGHT on (tl, FIT_TOP], 0 elsewhere."""
    return np.where((f >= FIT_LO) & (f <= tl), 1.0, np.where((f > tl) & (f <= FIT_TOP), ABOVE_WEIGHT, 0.0))


def _fit_level(p: np.ndarray, ps: np.ndarray, band: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """The residual's weights (band, 0 where the sum's power p or its smoothed power ps is exactly 0), the smoothed dB
    (0 where weightless) and their weighted mean: the level fit_passive_at fits to."""
    w = np.where((p > 0) & (ps > 0), band, 0.0)
    db = 10 * np.log10(np.where(w > 0, ps, 1.0))
    return w, db, float(np.sum(w * db) / max(np.sum(w), 1e-30))


def _low_rise(f: np.ndarray, p: np.ndarray, ps: np.ndarray, level: float) -> np.ndarray:
    """dB the smoothed sum (power ps; p unsmoothed) sits above level at each point of f under FIT_LO; -inf where the
    sum is exactly zero (skipped)."""
    lo = f < FIT_LO
    ok = (p[lo] > 0) & (ps[lo] > 0)
    return np.where(ok, 10 * np.log10(np.where(ok, ps[lo], 1.0)) - level, -np.inf)


@dataclass(frozen=True)
class PassiveFit:
    bass_coax: float  # Hz
    lf_hf: float  # Hz
    topo: Topology
    values: dict[str, float]  # snapped to standard values, in value_names order
    continuous: dict[str, float]  # the continuous optimum the values were snapped from
    z_floor: float  # ohm
    z_min: float  # ohm, min |z_in| on Z_GRID at the snapped values
    z_ok: bool  # z_min >= z_floor
    cost: float  # sum of squares of the weighted residual and the penalties, at the snapped values


def _mag(z_at: dict[str, np.ndarray], load: str, f0: float) -> float:
    """|Z| of a load (given on Z_GRID) at f0, interpolated on log f."""
    return float(np.interp(np.log(f0), np.log(Z_GRID), np.abs(z_at[load])))


def _l_pad(r: float, pad_db: float) -> tuple[float, float]:
    """Series and parallel resistors that attenuate by pad_db into a load r and still look like r."""
    a = 10 ** (-max(pad_db, 0.0) / 20)
    return r * (1 - a), (r * a / (1 - a) if a < 1 else np.inf)


def start_values(topo: Topology, fb: float, fh: float, z_at: dict[str, np.ndarray],
                 pad_db: dict[str, float] | None = None) -> dict[str, float]:
    """The fit's start (passive mode): textbook Butterworth filters at fb and fh, each on its load's |Z| at that
    point (z_at: driver_loads on Z_GRID), L-pads that attenuate the coax LF and HF by pad_db (dB, default 0)
    and still present that |Z|, the notch tuned to NOTCH_HZ and the shunt to SHUNT_HZ at Q 1. Every value
    is clipped into BOUNDS."""
    pad = pad_db or {}
    wb, wh = 2 * np.pi * fb, 2 * np.pi * fh
    rb, rcx = _mag(z_at, "bass", fb), _mag(z_at, "coax_lf", fb)
    rl, rt = _mag(z_at, "coax_lf", fh), _mag(z_at, "coax_hf", fh)
    v: dict[str, float] = {}
    if topo.bass_order == 3:
        v.update(L_b1=1.5 * rb / wb, C_b1=4 / 3 / (wb * rb), L_b2=0.5 * rb / wb)
    else:
        v.update(L_b1=np.sqrt(2) * rb / wb, C_b1=1 / (np.sqrt(2) * wb * rb))
    v.update(R_bz=rb, C_bz=1 / (wb * rb))  # Zobel: the bass load's |Z|, corner at fb
    v.update(C_h1=1 / (np.sqrt(2) * wb * rcx), L_h1=np.sqrt(2) * rcx / wb)
    if topo.lf_order == 3:
        v.update(L_l1=1.5 * rl / wh, C_l1=4 / 3 / (wh * rl), L_l2=0.5 * rl / wh)
    else:
        v.update(L_l1=np.sqrt(2) * rl / wh, C_l1=1 / (np.sqrt(2) * wh * rl))
    lf_pad = pad.get("coax_lf", 0.0)
    if topo.lf_shelf:  # SHELF_DB of the coax LF's pad moves into the contour, its corner at SHELF_CORNER x fb
        rc = _mag(z_at, "coax_lf", 2 * fb)
        v.update(R_lc=rc * (10 ** (SHELF_DB / 20) - 1))
        v.update(L_lc=v["R_lc"] / (2 * np.pi * SHELF_CORNER * fb))
        lf_pad = max(lf_pad - SHELF_DB, 0.0)
    v["R_ls"], v["R_lp"] = _l_pad(rl, lf_pad)
    wn, rn = 2 * np.pi * NOTCH_HZ, _mag(z_at, "coax_lf", NOTCH_HZ)
    v.update(L_ln=rn / wn, C_ln=1 / (wn * rn), R_ln=rn)  # tank reactance = the load at resonance, 6 dB deep
    v.update(C_t1=2 / 3 / (wh * rt), L_t1=0.75 * rt / wh, C_t2=2 / (wh * rt))
    v["R_ts"], v["R_tp"] = _l_pad(rt, pad.get("coax_hf", 0.0))
    ws, rs = 2 * np.pi * SHUNT_HZ, _mag(z_at, "coax_hf", SHUNT_HZ)
    v.update(R_tr=rs, L_tr=rs / ws, C_tr=1 / (ws * rs))
    out = {}
    for n in value_names(topo, "passive"):
        lo, hi = BOUNDS[n]
        out[n] = float(np.clip(v[n], lo * np.exp(_EDGE), hi * np.exp(-_EDGE)))
    return out


def passive_total(ways: Ways, fit: PassiveFit, z: Callable[[np.ndarray], dict[str, np.ndarray]],
                  mode: Mode = "passive") -> np.ndarray:
    """(nf, nd) complex sum of the ways through the network at fit.values; semi_active sums the coax only."""
    t = fit.topo
    h, _ = transfers(build(fit.values, t, mode), ways.f, driver_loads(z(ways.f), t), t)
    return _sum(ways, h, t)


def _sum(ways: Ways, h: dict[str, np.ndarray], topo: Topology) -> np.ndarray:
    """(nf, nd) sum of the ways in h, each through its terminal voltage h (on ways.f)."""
    scale = {"bass": bass_scale(topo)}
    return sum(scale.get(k, 1.0) * h[k][:, None] * ways.ways[k] for k in h)


def _smoother(f: np.ndarray, fraction: float) -> np.ndarray:
    """The matrix speaker.curves.smooth applies to powers on the grid f."""
    lf = np.log2(f)
    m = (np.abs(lf[:, None] - lf[None, :]) <= fraction / 2) * np.gradient(lf)
    return m / m.sum(axis=1)[:, None]


def _pads_db(ways: Ways, topo: Topology) -> dict[str, float]:
    """Attenuation (dB) that brings the coax LF and HF down to the quietest way at 23° (match_gains' bands;
    a way with no points in its band is taken as level with the coax LF); the bass has no pad."""
    g = match_gains(ways, {"bass": 0.0, "coax_hf": 0.0})
    level = {"bass": -g["bass"] + 20 * np.log10(bass_scale(topo)), "coax_lf": 0.0, "coax_hf": -g["coax_hf"]}
    floor = min(level.values())
    return {k: level[k] - floor for k in ("coax_lf", "coax_hf")}


def _woofer_level(f: np.ndarray, bass: np.ndarray, sm: np.ndarray) -> float | None:
    """Mean over MATCH_BANDS["bass"] of the smoothed (sm) dB of the unfiltered bass way at 23° (bass, complex Pa on
    f, bass_scale applied); None when the band holds no point where the way is non-zero."""
    ps = sm @ np.abs(bass) ** 2
    m = (f >= MATCH_BANDS["bass"][0]) & (f <= MATCH_BANDS["bass"][1]) & (ps > 0)
    return float(np.mean(10 * np.log10(ps[m]))) if m.any() else None


def woofer_level_db(ways: Ways, topo: Topology) -> float | None:
    """The woofers' own level at 23° (dB re 1 Pa peak, FIT_SMOOTH-smoothed, mean over MATCH_BANDS["bass"],
    wired as topo), the level fit_passive_at anchors to; None as _woofer_level."""
    f = np.asarray(ways.f, dtype=float)
    return _woofer_level(f, bass_scale(topo) * ways.ways["bass"][:, ways.obs.at23_idx], _smoother(f, FIT_SMOOTH))


def low_rise_db(ways: Ways, fit: PassiveFit, z: Callable[[np.ndarray], dict[str, np.ndarray]],
                tl: float) -> tuple[float, float] | None:
    """The largest rise (dB) of the FIT_SMOOTH-smoothed 23° sum (passive mode) above the level fit_passive_at fits,
    over the points of ways.f under FIT_LO, and where (Hz); None when no such point has a non-zero sum."""
    f = np.asarray(ways.f, dtype=float)
    p = np.abs(passive_total(ways, fit, z)[:, ways.obs.at23_idx]) ** 2
    ps = _smoother(f, FIT_SMOOTH) @ p
    rise = _low_rise(f, p, ps, _fit_level(p, ps, _fit_band(f, tl))[2])
    if not np.isfinite(rise).any():
        return None
    j = int(np.argmax(rise))
    return float(rise[j]), float(f[f < FIT_LO][j])


def fit_passive_at(ways: Ways, z: Callable[[np.ndarray], dict[str, np.ndarray]], tl: float, fb: float, fh: float,
                   topo: Topology, z_floor: float = 3.5) -> PassiveFit:
    """The passive network's values (passive mode, topo's options and polarity) at crossover points fb, fh.

    The residual is the FIT_SMOOTH-smoothed dB of the summed 23° response minus its weighted mean (the level is
    free), weighted 1 on [FIT_LO, tl] and ABOVE_WEIGHT on (tl, FIT_TOP], 0 where the sum is exactly 0. To it are
    appended: LEVEL_WEIGHT times the dB the level (that weighted mean) sits under woofer_level_db less
    LEVEL_SLACK_DB, which stops the fit buying flatness or the floor with level (a lossy pad, a high-impedance bass
    branch); LEAK_WEIGHT per dB (times the point's weight) that the bass sits less than BASS_UNDER_DB under the coax
    LF over [2 fb, BASS_UNDER_TOP] (each smoothed alike, at 23°; points where either way is exactly zero skipped),
    which keeps the non-cardioid woofers out of the band the cardioid coax section must own (flatness alone lets
    the bass fill the coax LF's low-end droop there); CROSS_WEIGHT times the dB between the bass and the coax LF
    (smoothed alike, at 23°) at fb, and
    between the coax LF and HF at fh, which holds the crossovers at the points (flatness alone does not: on
    flat ways the sections can stagger anywhere and still sum flat); Z_PENALTY per ohm of |z_in| under z_floor
    at each Z_GRID point; STOP_WEIGHT per dB that a branch's |H| (the network on the driver loads, electrical, not
    acoustic) sits over STOP_DB re its |H| at the crossover point, at each Z_GRID point of its far stopbands
    (stopbands, stop_excess_db), which keeps a branch from resonating where the acoustic data cannot see it (the
    bass is zero above tl, and the smoothing hides narrow peaks); PEAK_WEIGHT per dB that any branch's |H| sits over
    ELEC_GAIN_MAX_DB at each Z_GRID point (gain_excess_db: a passive network only gains by resonating), and PEAK_WEIGHT
    per dB that the smoothed 23° sum sits more than PEAK_TOL_DB over the level at each point of ways.f under FIT_LO
    (one-sided: a roll-off there costs nothing; points where the sum is exactly zero skipped), which together keep a
    resonance out of the band under FIT_LO that the residual does not see (the sealed woofers roll off there by
    design, so the band is not lowered); and, with the RLC shunt, Q_PENALTY per unit of its Q over SHUNT_Q_MAX.

    The log values are fitted inside BOUNDS by least squares from start_values (L-pads from the level
    difference), then snapped. Coordinate-descent passes then try each value's two standard neighbours (inside
    BOUNDS) and keep the better one if it lowers the cost without breaking the floor, until a pass changes
    nothing. Should the snap itself break the floor, a move that shrinks the shortfall is kept first, whatever
    it costs: the floor is a hard constraint."""
    names = value_names(topo, "passive")
    f = np.asarray(ways.f, dtype=float)
    nf, i = len(f), ways.obs.at23_idx
    zf, zg = driver_loads(z(f), topo), driver_loads(z(Z_GRID), topo)
    zs = driver_loads(z(np.array([fb, fh])), topo)
    grid = np.concatenate([f, stop_grid(fb, fh)])
    loads = {k: np.concatenate([zf[k], zg[k], zs[k]]) for k in zf}
    resp = {k: ways.ways[k][:, i] * (bass_scale(topo) if k == "bass" else 1.0)
            for k in ("bass", "coax_lf", "coax_hf")}
    band = _fit_band(f, tl)
    if not band.any():
        raise ValueError(f"no frequencies between {FIT_LO:g} and {max(tl, FIT_TOP):g} Hz to fit on")
    sm = _smoother(f, FIT_SMOOTH)
    target = _woofer_level(f, resp["bass"], sm)
    under = np.flatnonzero((f >= 2 * fb) & (f <= BASS_UNDER_TOP) & (np.abs(resp["bass"]) > 0)
                           & (np.abs(resp["coax_lf"]) > 0) & (band > 0))
    # the two grid points either side of each crossover point, and its log-f position between them
    cross = []
    for a, b, fc in (("bass", "coax_lf", fb), ("coax_lf", "coax_hf", fh)):
        if f[0] <= fc <= f[-1] and nf > 1:
            j = int(np.clip(np.searchsorted(f, fc), 1, nf - 1))
            cross.append((a, b, slice(j - 1, j + 1), np.log(fc / f[j - 1]) / np.log(f[j] / f[j - 1])))

    def evaluate(v: dict[str, float]) -> tuple[np.ndarray, float]:
        h, z_in = transfers(build(v, topo, "passive"), grid, loads, topo)
        each = {k: h[k][:nf] * resp[k] for k in h}
        p = np.abs(sum(each.values())) ** 2
        ps = sm @ p
        w, db, level = _fit_level(p, ps, band)
        r = w * (db - level)
        # no level term when the woofers have no level to anchor to (no points in their band, or exactly zero)
        lv = [] if target is None else [LEVEL_WEIGHT * max(0.0, target - LEVEL_SLACK_DB - level)]
        pb, pc = sm[under] @ np.abs(each["bass"]) ** 2, sm[under] @ np.abs(each["coax_lf"]) ** 2
        ok = (pb > 0) & (pc > 0)
        gap = np.where(ok, 10 * np.log10(np.where(ok, pb, 1.0) / np.where(ok, pc, 1.0)) + BASS_UNDER_DB, 0.0)
        bu = LEAK_WEIGHT * band[under] * np.maximum(0.0, gap)
        xr = np.zeros(2)
        for n, (a, b, js, t) in enumerate(cross):
            pa, pb = sm[js] @ np.abs(each[a]) ** 2, sm[js] @ np.abs(each[b]) ** 2
            if np.all(pa > 0) and np.all(pb > 0):
                d = 10 * np.log10(pa / pb)
                xr[n] = CROSS_WEIGHT * ((1 - t) * d[0] + t * d[1])
        zm = np.abs(z_in[nf:nf + len(Z_GRID)])
        pen = [Z_PENALTY * np.maximum(0.0, z_floor - zm),
               STOP_WEIGHT * stop_excess_db({k: x[nf:] for k, x in h.items()}, fb, fh),
               PEAK_WEIGHT * gain_excess_db({k: x[nf:nf + len(Z_GRID)] for k, x in h.items()}),
               PEAK_WEIGHT * np.maximum(0.0, _low_rise(f, p, ps, level) - PEAK_TOL_DB)]
        if topo.hf_shunt:
            pen.append([Q_PENALTY * max(0.0, np.sqrt(v["L_tr"] / v["C_tr"]) / v["R_tr"] - SHUNT_Q_MAX)])
        return np.concatenate([r, xr, lv, bu, *pen]), float(zm.min())

    lo = np.log([BOUNDS[n][0] for n in names])
    hi = np.log([BOUNDS[n][1] for n in names])
    start = start_values(topo, fb, fh, zg, _pads_db(ways, topo))
    x0 = np.log([start[n] for n in names])
    x = least_squares(lambda y: evaluate(dict(zip(names, np.exp(y))))[0], x0, bounds=(lo, hi),
                      ftol=_LSQ_TOL, xtol=_LSQ_TOL).x
    continuous = {n: float(v) for n, v in zip(names, np.exp(x))}

    def rank(v: dict[str, float]) -> tuple[float, float, float]:
        """(shortfall under the floor, cost, z_min): lower shortfall first, then lower cost."""
        r, zm = evaluate(v)
        return max(0.0, z_floor - zm), float(r @ r), zm

    vals = {n: snap(n[0], v) for n, v in continuous.items()}
    best = rank(vals)
    changed = True
    while changed:
        changed = False
        for n in names:
            moves = [(rank(t), t) for t in ({**vals, n: c} for c in neighbours(n[0], vals[n]))
                     if t[n] != vals[n] and BOUNDS[n][0] <= t[n] <= BOUNDS[n][1]]
            if moves:
                k, t = min(moves, key=lambda m: m[0][:2])
                if k[:2] < best[:2]:
                    vals, best, changed = t, k, True
    _, cost, z_min = best
    return PassiveFit(fb, fh, topo, vals, continuous, z_floor, z_min, z_min >= z_floor, cost)


OPTIONS = ("lf_shelf", "zobel", "bass_order", "lf_order", "lf_notch", "hf_shunt")  # in this order; an order "on" is 3
OPTION_TEXT = {"lf_shelf": "the coax LF contour (R_lc with L_lc across it)",
               "zobel": "the Zobel across the woofers", "bass_order": "a 3rd-order bass low-pass",
               "lf_order": "a 3rd-order coax LF low-pass", "lf_notch": f"the coax LF notch at {NOTCH_HZ:g} Hz",
               "hf_shunt": "the RLC shunt across the coax HF"}
OPTION_GAIN_DB = 0.1  # dB, the least trusted flat23 improvement that keeps an option
OPTION_SLACK_DB = 0.25  # dB, the most any other search score may get worse for it
_DB_EPS = 1e-9  # dB, so a change of exactly a threshold counts as reaching it
TOL_DRAWS = 200  # random draws of the part values
TOL_SPREAD = 0.05  # each L and C uniform within +-5 %
TOL_SEED = 0
TOL_Z_MIN = 3.3  # ohm, no draw may go lower
TOL_FLAT_DB = 0.5  # dB, the 95th-percentile flat23 may be at most this far above the nominal one
FLOOR_STEP = 0.1  # ohm, the floor is raised by this while the tolerance check fails
FLOOR_MAX = 3.8  # ohm, and no further


def _score(ways: Ways, total: np.ndarray, tl: float) -> list[Score]:
    spin = spinorama(total, ways.obs, f=ways.f)
    return score_all(ways.f, total, coax_section(ways), ways.obs, spin, trusted_limit=tl,
                     puck_blend=ways.meta.get("puck_blend"), coax_far=ways.meta.get("coax_far_h0_h180"))


def passive_scores(ways: Ways, fit: PassiveFit, z: Callable[[np.ndarray], dict[str, np.ndarray]], tl: float,
                   mode: Mode = "passive") -> list[Score]:
    """score_all on the passive total over all of ways, with the arguments xo_fit.fit_active passes; scores
    above tl come out low-confidence, and only the trusted ones enter the search."""
    return _score(ways, passive_total(ways, fit, z, mode), tl)


def _trusted(scores: list[Score]) -> dict[str, float]:
    """The trusted, finite SEARCH_KEYS values."""
    return {s.key: float(s.value) for s in scores
            if s.confidence == "trusted" and s.key in SEARCH_KEYS and math.isfinite(s.value)}


@dataclass(frozen=True)
class Tolerance:
    z_min: float  # ohm, min |z_in| on Z_GRID of the worst draw
    flat23_nominal: float  # dB, trusted, at the nominal values
    flat23_p95: float  # dB, 95th percentile of the draws
    whole_rear: tuple[float, float]  # dB, 5th and 95th percentile of the draws
    passed: bool  # z_min >= TOL_Z_MIN and flat23_p95 <= flat23_nominal + TOL_FLAT_DB
    # each draw's min |z_in| (ohm) and trusted flat23 (dB, nan where not trusted), for the report's plot
    draws_z_min: tuple[float, ...] = field(default=(), compare=False, repr=False)
    draws_flat23: tuple[float, ...] = field(default=(), compare=False, repr=False)


def tolerance(ways: Ways, z: Callable[[np.ndarray], dict[str, np.ndarray]], tl: float, fit: PassiveFit,
              n: int = 200, spread: float = 0.05, seed: int = 0) -> Tolerance:
    """n draws of fit's network (passive mode), every L and C scaled by its own uniform factor in
    [1 - spread, 1 + spread] (np.random.default_rng(seed)); inductor DCR stays at the nominal part's and
    resistors are exact. Each draw is evaluated as it stands, not refitted: min |z_in| on Z_GRID and the
    trusted flat23 and whole_rear (nan where not trusted, which fails the check)."""
    t = fit.topo
    parts = build(fit.values, t, "passive")
    f = np.asarray(ways.f, dtype=float)
    nf = len(f)
    grid = np.concatenate([f, Z_GRID])
    zf, zg = driver_loads(z(f), t), driver_loads(z(Z_GRID), t)
    loads = {k: np.concatenate([zf[k], zg[k]]) for k in zf}
    lc = [j for j, p in enumerate(parts) if p.kind in ("L", "C")]
    factors = np.random.default_rng(seed).uniform(1 - spread, 1 + spread, size=(n, len(lc)))

    def run(ps: list[Part]) -> tuple[float, float, float]:
        h, z_in = transfers(ps, grid, loads, t)
        v = _trusted(_score(ways, _sum(ways, {k: x[:nf] for k, x in h.items()}, t), tl))
        return float(np.abs(z_in[nf:]).min()), v.get("flat23", math.nan), v.get("whole_rear", math.nan)

    _, flat0, _ = run(parts)
    draws = []
    for row in factors:
        ps = list(parts)
        for j, a in zip(lc, row):
            ps[j] = replace(parts[j], value=parts[j].value * float(a))
        draws.append(run(ps))
    zs, flat, rear = np.array(draws).reshape(n, 3).T
    z_min, p95 = float(zs.min()), float(np.percentile(flat, 95))
    passed = z_min >= TOL_Z_MIN and math.isfinite(flat0) and math.isfinite(p95) and p95 <= flat0 + TOL_FLAT_DB
    return Tolerance(z_min, flat0, p95, (float(np.percentile(rear, 5)), float(np.percentile(rear, 95))),
                     bool(passed), tuple(float(x) for x in zs), tuple(float(x) for x in flat))


@dataclass(frozen=True)
class PassiveResult:
    fit: PassiveFit
    tolerance: Tolerance
    series: PassiveResult | None  # the woofers-in-series alternative, only when fit is not z_ok
    notes: list[str]  # one plain sentence per decision, in the order taken


def _option_verdict(old: dict[str, float], new: dict[str, float]) -> str:
    """Why an option's trusted SEARCH_KEYS scores new do not beat old ("" when they do): flat23 must improve by
    OPTION_GAIN_DB, and no other score may get worse by more than OPTION_SLACK_DB (whole_rear higher is better,
    every other key lower)."""
    gain = old["flat23"] - new["flat23"]
    if gain < OPTION_GAIN_DB - _DB_EPS:
        how = "no better" if gain <= 0 else f"less than the {OPTION_GAIN_DB:g} dB gain needed"
        return f"flat23 {old['flat23']:.2f} → {new['flat23']:.2f} dB, {how}"
    for k in SEARCH_KEYS:
        worse = old[k] - new[k] if k == "whole_rear" else new[k] - old[k]
        if k != "flat23" and worse > OPTION_SLACK_DB + _DB_EPS:
            return f"{k} worse by {worse:.2f} dB (at most {OPTION_SLACK_DB:g} dB allowed)"
    return ""


def _polarity(t: Topology) -> str:
    return (f"bass {'inverted' if t.bass_inverted else 'normal'}, "
            f"coax HF {'inverted' if t.hf_inverted else 'normal'}")


def _grid(ways: Ways, z: Callable[[np.ndarray], dict[str, np.ndarray]], tl: float, lf_hf: list[float],
          base: Topology, z_floor: float) -> tuple[PassiveFit, dict[str, float], int, bool]:
    """At every grid pair, the best of the four polarity pairs by (shortfall under the floor, cost), the key
    fit_passive_at's own descent uses; pick among those with every trusted SEARCH_KEYS score (grid order,
    bass_coax outer, lf_hf inner), only the z_ok ones when there are any. Returns the fit, its scores, the
    candidate count and whether any candidate was z_ok."""
    cands = []
    for fb in BASS_COAX_GRID:
        for fh in lf_hf:
            fit = min((fit_passive_at(ways, z, tl, fb, fh, replace(base, bass_inverted=bi, hf_inverted=hi), z_floor)
                       for bi in (False, True) for hi in (False, True)),
                      key=lambda x: (max(0.0, x.z_floor - x.z_min), x.cost))
            v = _trusted(passive_scores(ways, fit, z, tl))
            if len(v) == len(SEARCH_KEYS):
                cands.append(((fit, v), v))
    if not cands:
        raise ValueError(f"no passive crossover candidate has all {len(SEARCH_KEYS)} trusted search scores")
    ok = [c for c in cands if c[0][0].z_ok]
    fit, v = pick(ok or cands)
    return fit, v, len(cands), bool(ok)


def _tolerance_why(tol: Tolerance) -> str:
    why = []
    if not tol.z_min >= TOL_Z_MIN:
        why.append(f"worst draw {tol.z_min:.2f} Ω, under {TOL_Z_MIN:g} Ω")
    if not tol.flat23_p95 <= tol.flat23_nominal + TOL_FLAT_DB:
        why.append(f"flat23 {tol.flat23_p95:.2f} dB at the 95th percentile against {tol.flat23_nominal:.2f} dB "
                   f"nominal")
    return "; ".join(why) or "no reason recorded"


def _search(ways: Ways, z: Callable[[np.ndarray], dict[str, np.ndarray]], tl: float, lf_hf: list[float],
            base: Topology, z_floor: float) -> tuple[PassiveResult, bool]:
    """Steps 1-3 of search_passive for one bass wiring, and whether any grid candidate cleared the floor."""
    fit, v, count, any_ok = _grid(ways, z, tl, lf_hf, base, z_floor)
    notes = [f"Chose {fit.bass_coax:g} Hz bass to coax and {fit.lf_hf:g} Hz coax LF to HF from {count} grid pairs, "
             f"woofers in {base.bass_wiring} ({_polarity(fit.topo)}): flat23 {v['flat23']:.2f} dB, "
             f"whole_rear {v['whole_rear']:.2f} dB."]

    for opt in OPTIONS:
        t = replace(fit.topo, **{opt: 3 if opt.endswith("_order") else True})
        new = fit_passive_at(ways, z, tl, fit.bass_coax, fit.lf_hf, t, fit.z_floor)
        nv = _trusted(passive_scores(ways, new, z, tl))
        if len(nv) < len(SEARCH_KEYS):
            why = "not every search score is trusted with it"
        elif fit.z_ok and not new.z_ok:
            why = f"it takes the impedance to {new.z_min:.2f} Ω, under the {fit.z_floor:g} Ω floor"
        elif not fit.z_ok and not new.z_ok and new.z_min < fit.z_min:
            why = (f"it lowers the impedance from {fit.z_min:.2f} Ω to {new.z_min:.2f} Ω, further under the "
                   f"{fit.z_floor:g} Ω floor")
        elif new.z_ok and not fit.z_ok:  # the floor outranks the scores
            notes.append(f"Kept {OPTION_TEXT[opt]}: it lifts the impedance from {fit.z_min:.2f} Ω to "
                         f"{new.z_min:.2f} Ω, clearing the {fit.z_floor:g} Ω floor (flat23 {v['flat23']:.2f} → "
                         f"{nv['flat23']:.2f} dB).")
            fit, v = new, nv
            continue
        else:
            why = _option_verdict(v, nv)
        if why:
            notes.append(f"Left out {OPTION_TEXT[opt]}: {why}.")
        else:
            notes.append(f"Kept {OPTION_TEXT[opt]}: flat23 {v['flat23']:.2f} → {nv['flat23']:.2f} dB.")
            fit, v = new, nv

    tol = tolerance(ways, z, tl, fit, n=TOL_DRAWS, spread=TOL_SPREAD, seed=TOL_SEED)
    stop = ""  # why the floor stopped rising, for the note
    while not tol.passed:
        if not fit.z_ok:
            stop = ", which the network already misses, so it was not raised"
            break
        if fit.z_floor >= FLOOR_MAX - 1e-9:
            stop = ", the highest tried"
            break
        floor = round(fit.z_floor + FLOOR_STEP, 6)
        new = fit_passive_at(ways, z, tl, fit.bass_coax, fit.lf_hf, fit.topo, floor)
        if not new.z_ok:
            stop = (f"; a refit at {floor:g} Ω reaches only {new.z_min:.2f} Ω, so the {fit.z_floor:g} Ω "
                    f"network is kept")
            break
        notes.append(f"Raised the impedance floor to {floor:g} Ω: the tolerance check failed at "
                     f"{fit.z_floor:g} Ω ({_tolerance_why(tol)}).")
        fit = new
        tol = tolerance(ways, z, tl, fit, n=TOL_DRAWS, spread=TOL_SPREAD, seed=TOL_SEED)
    if tol.passed:
        notes.append(f"The tolerance check passes at the {fit.z_floor:g} Ω floor: worst draw {tol.z_min:.2f} Ω, "
                     f"flat23 {tol.flat23_nominal:.2f} dB nominal and {tol.flat23_p95:.2f} dB at the 95th "
                     f"percentile.")
    else:
        notes.append(f"The tolerance check still fails at the {fit.z_floor:g} Ω floor{stop}: worst draw "
                     f"{tol.z_min:.2f} Ω (limit {TOL_Z_MIN:g} Ω), flat23 "
                     f"{tol.flat23_p95 - tol.flat23_nominal:+.2f} dB at the 95th percentile on the nominal "
                     f"{tol.flat23_nominal:.2f} dB (limit +{TOL_FLAT_DB:g} dB).")
    return PassiveResult(fit, tol, None, notes), any_ok


def search_passive(ways: Ways, z: Callable[[np.ndarray], dict[str, np.ndarray]], tl: float,
                   lf_hf_min: float | None, z_floor: float = 3.5) -> PassiveResult:
    """The passive network: (1) at every grid pair (LF_HF_GRID floored at lf_hf_min) the polarity pair with the
    least shortfall under the floor, then the lowest cost, scored by passive_scores, the points chosen with
    xo_fit.pick (among the z_ok candidates when there are any); (2) at those points, each of OPTIONS in turn,
    kept when it lifts a network that misses the floor over it, or when _option_verdict allows and the floor
    still holds (and, for a network already under the floor, the option does not lower its z_min further); (3)
    while the tolerance check fails,
    the network meets its floor and the floor is under FLOOR_MAX, a refit at the same points and topology with
    the floor FLOOR_STEP higher, stopping (and keeping the last network that met its floor, with its tolerance
    result) at a refit that misses the raised floor; (4) if the final network's z_min is under z_floor (never
    for a missed raised floor), steps 1-3 again with the woofers in series, attached as `series` (should that search
    raise ValueError, the parallel network is returned with a note saying why). Each step leaves a sentence in
    notes."""
    lf_hf = [fh for fh in LF_HF_GRID if lf_hf_min is None or fh >= lf_hf_min]
    if not lf_hf:
        raise ValueError(f"no LF-HF grid point is at or above the driver's {lf_hf_min:g} Hz; "
                         f"the grid stops at {LF_HF_GRID[-1]:g} Hz")
    res, any_ok = _search(ways, z, tl, lf_hf, Topology(), z_floor)
    if res.fit.z_min >= z_floor:
        return res
    try:
        series, _ = _search(ways, z, tl, lf_hf, Topology(bass_wiring="series"), z_floor)
    except ValueError as e:  # keep the parallel network: it is still the answer for the woofers as they are
        return replace(res, notes=[*res.notes, f"The chosen parallel-woofer network does not clear {z_floor:g} Ω (its "
                                               f"lowest is {res.fit.z_min:.2f} Ω), and the woofers-in-series search "
                                               f"failed: {e}."])
    tail = "" if series.fit.z_min >= z_floor else f", though it does not clear {z_floor:g} Ω either"
    if any_ok:
        head = f"The chosen parallel-woofer network does not clear {z_floor:g} Ω"
    else:
        head = f"No parallel-woofer network in the search clears {z_floor:g} Ω"
    note = f"{head} (its lowest is {res.fit.z_min:.2f} Ω); the woofers-in-series alternative is attached{tail}."
    return replace(res, series=series, notes=[*res.notes, note])
