"""The passive crossover topology: which parts exist for a set of options, the netlist they form in passive
and semi-active modes, the driver loads, and each way's terminal voltage per input volt; and the fit of the
part values at one pair of crossover points, snapped to standard values.

Passive: a bass low-pass on `in`, and a coax high-pass from `in` to node `cx`, which feeds the coax LF
(low-pass, optional notch, L-pad) and the coax HF (high-pass, L-pad, optional series RLC shunt).
Semi-active: the coax network alone, with `cx` being `in` (the DSP supplies the split, the bass runs active).
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.optimize import least_squares

from speaker.crossover import match_gains
from speaker.network import Part, inductor_dcr, neighbours, snap, solve
from speaker.system import Ways

Mode = Literal["passive", "semi_active"]


@dataclass(frozen=True)
class Topology:
    bass_order: int = 2  # 2 or 3
    zobel: bool = False
    lf_order: int = 2  # 2 or 3
    lf_notch: bool = False
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
    names += ["C_t1", "L_t1", "C_t2", "R_ts", "R_tp"]
    names += ["R_tr", "L_tr", "C_tr"] if topo.hf_shunt else []
    return names


_ALL = value_names(Topology(bass_order=3, zobel=True, lf_order=3, lf_notch=True, hf_shunt=True), "passive")
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

    # Coax LF: low-pass, notch, then the L-pad next to the driver.
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
NOTCH_HZ = 2700.0  # Hz, the coax LF cone breakup
SHUNT_HZ = float(np.sqrt(2000.0 * 3300.0))  # Hz, centre of the coax HF's 2-3.3 kHz hump
_EDGE = 1e-6  # keeps start values strictly inside the bounds (log scale)
_LSQ_TOL = 1e-4  # least_squares ftol and xtol: finer log-value steps are lost in the ~20 % snap to standard values


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
    v["R_ls"], v["R_lp"] = _l_pad(rl, pad.get("coax_lf", 0.0))
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
    scale = {"bass": bass_scale(t)}
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


def fit_passive_at(ways: Ways, z: Callable[[np.ndarray], dict[str, np.ndarray]], tl: float, fb: float, fh: float,
                   topo: Topology, z_floor: float = 3.5) -> PassiveFit:
    """The passive network's values (passive mode, topo's options and polarity) at crossover points fb, fh.

    The residual is the FIT_SMOOTH-smoothed dB of the summed 23° response minus its weighted mean (the level is
    free), weighted 1 on [FIT_LO, tl] and ABOVE_WEIGHT on (tl, FIT_TOP], 0 where the sum is exactly 0. To it are
    appended: CROSS_WEIGHT times the dB between the bass and the coax LF (smoothed alike, at 23°) at fb, and
    between the coax LF and HF at fh, which holds the crossovers at the points (flatness alone does not: on
    flat ways the sections can stagger anywhere and still sum flat); Z_PENALTY per ohm of |z_in| under z_floor
    at each Z_GRID point; and, with the RLC shunt, Q_PENALTY per unit of its Q over SHUNT_Q_MAX.

    The log values are fitted inside BOUNDS by least squares from start_values (L-pads from the level
    difference), then snapped. Coordinate-descent passes then try each value's two standard neighbours (inside
    BOUNDS) and keep the better one if it lowers the cost without breaking the floor, until a pass changes
    nothing. Should the snap itself break the floor, a move that shrinks the shortfall is kept first, whatever
    it costs: the floor is a hard constraint."""
    names = value_names(topo, "passive")
    f = np.asarray(ways.f, dtype=float)
    nf, i = len(f), ways.obs.at23_idx
    zf, zg = driver_loads(z(f), topo), driver_loads(z(Z_GRID), topo)
    grid = np.concatenate([f, Z_GRID])
    loads = {k: np.concatenate([zf[k], zg[k]]) for k in zf}
    resp = {k: ways.ways[k][:, i] * (bass_scale(topo) if k == "bass" else 1.0)
            for k in ("bass", "coax_lf", "coax_hf")}
    band = np.where((f >= FIT_LO) & (f <= tl), 1.0, np.where((f > tl) & (f <= FIT_TOP), ABOVE_WEIGHT, 0.0))
    if not band.any():
        raise ValueError(f"no frequencies between {FIT_LO:g} and {max(tl, FIT_TOP):g} Hz to fit on")
    sm = _smoother(f, FIT_SMOOTH)
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
        w = np.where((p > 0) & (ps > 0), band, 0.0)
        db = 10 * np.log10(np.where(w > 0, ps, 1.0))
        r = w * (db - np.sum(w * db) / max(np.sum(w), 1e-30))
        xr = np.zeros(2)
        for n, (a, b, js, t) in enumerate(cross):
            pa, pb = sm[js] @ np.abs(each[a]) ** 2, sm[js] @ np.abs(each[b]) ** 2
            if np.all(pa > 0) and np.all(pb > 0):
                d = 10 * np.log10(pa / pb)
                xr[n] = CROSS_WEIGHT * ((1 - t) * d[0] + t * d[1])
        zm = np.abs(z_in[nf:])
        pen = [Z_PENALTY * np.maximum(0.0, z_floor - zm)]
        if topo.hf_shunt:
            pen.append([Q_PENALTY * max(0.0, np.sqrt(v["L_tr"] / v["C_tr"]) / v["R_tr"] - SHUNT_Q_MAX)])
        return np.concatenate([r, xr, *pen]), float(zm.min())

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
