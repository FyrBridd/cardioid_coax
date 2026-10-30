"""``speaker passive``: the passive crossover for a design, its semi-active variant and the active fit beside them,
with the parts list, schematics, a VituixCAD netlist and an HTML report for the builder.

The fields and ways come from speaker.directivity_report.design_ways (cached BEM bands), the driver impedances
from speaker.vituixcad.impedances. The passive network is speaker.passive.search_passive, the semi-active DSP stage
speaker.semi_active.fit_semi_active and the active comparison speaker.xo_fit.fit_active, all on the same ways and
trusted limit, and all scored with the same speaker.scores.score_all. Plots are in speaker.directivity_plots.
"""
from __future__ import annotations

import dataclasses
import datetime
import html
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from . import directivity_plots as plots
from . import passive as P
from .crossover import apply_active, way_transfer
from .directivity_report import _CSS, _commit, _h, _jsonable, _params_sha, _row, _spl, design_ways
from .evaluate import Design, _datasheet_fallback, ways_scores
from .network import Part, solve
from .passive import (
    BOUNDS, ELEC_GAIN_MAX_DB, FIT_LO, OPTION_TEXT, PEAK_TOL_DB, STOP_DB, TOL_FLAT_DB, TOL_SPREAD, TOL_Z_MIN, Z_GRID,
    PassiveFit, PassiveResult, Topology, bass_scale, build, driver_loads, largest_gain_db, low_rise_db, passive_scores,
    passive_total, search_passive, stop_grid, stopbands, transfers, woofer_level_db,
)
from .schematic import LINK_TEXT, draw_links, draw_network, fmt_value, link_states, links_text, terminal_text
from .scores import Score
from .semi_active import _transfer, fit_semi_active, semi_active_text, semi_total, semi_ways
from .spinorama import P_REF_PEAK, spinorama
from .system import Ways
from .vituixcad import impedances
from .xo_fit import fit_active, passive_warnings

Z_FLOOR = 3.5  # ohm, the lowest the SA550 may see (main spec §7.2)
RATINGS_W = (5, 10, 15, 20, 30, 50)  # standard resistor power ratings
BUY_NEAR = 0.85  # a resistor whose 1/3-voltage dissipation reaches this share of its rating: buy 2-2.5 times it
CAP_RATINGS_V = (50, 63, 100, 160, 250)  # standard capacitor voltage ratings
CAP_MARGIN = 1.5  # a capacitor's voltage rating is at least this times its peak voltage at full amplifier voltage
TOLERANCE_TEXT = "5 %"
ELECTROLYTIC_MIN = 30e-6  # F, bipolar electrolytics are allowed above this, outside the coax network
SHORTFALL_DB = 0.25  # dB, a score difference at least this large is reported
_BASS_NETWORK = ("L_b", "C_b", "R_bz", "R_bc")  # name prefixes of the bass low-pass, its Zobel and its compensator
# parts at the end of their range where they barely act: L-pad resistors (a series one small, a parallel one
# large) and the coax HF's second series capacitor (large: nearly a plain wire, the high-pass close to 2nd order)
_FADING = {("R_ls", "lower"), ("R_ts", "lower"), ("R_lp", "upper"), ("R_tp", "upper"), ("C_t2", "upper")}
# score keys where a higher value is better; for every other key lower is better
HIGHER_BETTER = ("cardioid_hard", "cardioid_hard_far", "cardioid_aim", "rear_hemisphere", "whole_rear")
MEANING = {
    "flat23": "the response at the listening seat (23°) is less even",
    "window_spread": "the curves from 0° to 30° follow the seat's curve less closely, so the sound changes more "
                     "as you move",
    "presence_lift": "there is more 2-5 kHz lift off-axis than at the seat, in the range you hear best",
    "di_smooth_bass": "the directivity changes less smoothly around the bass-to-coax crossover",
    "di_smooth_lfhf": "the directivity changes less smoothly around the coax LF-to-HF crossover",
    "rolloff_monotonic": "more bands get louder again further off-axis",
    "whole_rear": "the whole speaker sends more sound backwards at 150-300 Hz",
}
SENS_BAND = (300.0, 1000.0)  # Hz, where the system sensitivity is taken
BASS_BAND = (20.0, 500.0)  # Hz, where the bass-band minimum impedance is taken
_BRANCH_TEXT = {"bass": "bass", "coax_lf": "coax LF", "coax_hf": "coax HF"}
MODES = {"passive": "Passive", "semi_active": "Semi-active", "active": "Active"}


def _log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- parts

def _rating(w: float) -> int:
    """The first of RATINGS_W at or above w (W); above the series, the next multiple of its top value."""
    for r in RATINGS_W:
        if w <= r * (1 + 1e-9):
            return r
    return int(RATINGS_W[-1] * math.ceil(w / RATINGS_W[-1]))


def power_ratings(parts: list[Part], z: dict[str, np.ndarray], f: np.ndarray,
                  v_rms: float) -> dict[str, tuple[float, float, int]]:
    """Per resistor: the largest sine dissipation over f (W) at v_rms / 3 and at v_rms on the input, and the
    standard rating for the first (_rating). z: each driver load's impedance on f, by its part name."""
    sol = solve(parts, f, z)
    out = {}
    for p in parts:
        if p.kind == "R":
            per_v2 = float(np.max(np.abs(sol.across(p)) ** 2 / p.value))  # W per V² rms at the input
            third = per_v2 * (v_rms / 3) ** 2
            out[p.name] = (third, per_v2 * v_rms ** 2, _rating(third))
    return out


def peak_currents(parts: list[Part], z: dict[str, np.ndarray], f: np.ndarray, v_rms: float) -> dict[str, float]:
    """Per inductor: the peak current (A) of a sine of v_rms on the input (the peak, √2 times the rms current), the
    largest over f. z: each driver load's impedance on f, by its part name."""
    sol = solve(parts, f, z)
    w = 2 * np.pi * np.asarray(f, dtype=float)
    return {p.name: float(np.sqrt(2) * v_rms * np.max(np.abs(sol.across(p) / (p.dcr + 1j * w * p.value))))
            for p in parts if p.kind == "L"}


def capacitor_stress(parts: list[Part], z: dict[str, np.ndarray], f: np.ndarray,
                     v_rms: float) -> dict[str, tuple[float, float]]:
    """Per capacitor: the peak voltage across it (V, √2 times the rms) and its rms (ripple) current (A) with a sine of
    v_rms on the input, each the largest over f. z: each driver load's impedance on f, by its part name."""
    sol = solve(parts, f, z)
    w = 2 * np.pi * np.asarray(f, dtype=float)
    out = {}
    for p in parts:
        if p.kind == "C":
            v = np.abs(sol.across(p))
            out[p.name] = (float(np.sqrt(2) * v_rms * v.max()), float(v_rms * np.max(v * w * p.value)))
    return out


def cap_voltage_rating(peak_v: float) -> int:
    """The first of CAP_RATINGS_V at or above CAP_MARGIN times peak_v (V); above the series, whole multiples of its
    top value."""
    need = CAP_MARGIN * peak_v
    for r in CAP_RATINGS_V:
        if need <= r * (1 + 1e-9):
            return r
    return int(CAP_RATINGS_V[-1] * math.ceil(need / CAP_RATINGS_V[-1]))


def _in_bass_network(name: str) -> bool:
    return name.startswith(_BASS_NETWORK) or name == "bass"


def _num(v: float) -> str:
    return f"{v:.4g}"


def _value_unit(p: Part) -> tuple[str, str]:
    """The value in the parts list's units (ASCII, for spreadsheets): mH, uF, ohm."""
    if p.kind == "L":
        return _num(p.value * 1e3), "mH"
    if p.kind == "C":
        return _num(p.value * 1e6), "uF"
    return _num(p.value), "ohm"


def _type(p: Part, rating: int | None, third_w: float | None = None, peak_v: float | None = None) -> str:
    """The part's type with its rating: a resistor's power rating (and, when third_w, its dissipation at 1/3 of the
    amp's voltage, reaches BUY_NEAR of it, the 2-2.5 times larger part to buy), a capacitor's voltage rating from
    its peak voltage at full voltage (peak_v)."""
    if p.kind == "L":
        return "cored inductor" if p.name.startswith("L_b") else "air-core inductor"
    if p.kind == "C":
        kind = ("bipolar electrolytic capacitor" if _in_bass_network(p.name) and p.value > ELECTROLYTIC_MIN
                else "polypropylene film capacitor")
        return kind if peak_v is None else f"{kind} ({cap_voltage_rating(peak_v)} V)"
    if rating is not None and third_w is not None and third_w >= BUY_NEAR * rating * (1 - 1e-9):
        return f"resistor ({rating} W; buy {2 * rating:g}-{2.5 * rating:g} W)"
    return f"resistor ({rating} W)"


def _loads(fit: PassiveFit, d: Design, f: np.ndarray) -> dict[str, np.ndarray]:
    return driver_loads(impedances(d, f), fit.topo)


BOTH_MODES = ("passive", "semi_active")


def _mode_ratings(fit: PassiveFit, d: Design,
                  modes: tuple[str, ...] = BOTH_MODES) -> dict[str, dict[str, tuple[float, float, int]]]:
    """power_ratings in each mode (semi-active: the coax channel at the same amp voltage), on Z_GRID."""
    loads = _loads(fit, d, Z_GRID)
    v = d.p.drive.amp_v_rms
    return {m: power_ratings(build(fit.values, fit.topo, m), loads, Z_GRID, v) for m in modes}


def _mode_currents(fit: PassiveFit, d: Design, modes: tuple[str, ...] = BOTH_MODES) -> dict[str, dict[str, float]]:
    """peak_currents at the amp's full voltage in each mode, on Z_GRID."""
    loads = _loads(fit, d, Z_GRID)
    return {m: peak_currents(build(fit.values, fit.topo, m), loads, Z_GRID, d.p.drive.amp_v_rms) for m in modes}


def _mode_caps(fit: PassiveFit, d: Design, modes: tuple[str, ...] = BOTH_MODES) -> dict[str, dict]:
    """capacitor_stress at the amp's full voltage in each mode, on Z_GRID."""
    loads = _loads(fit, d, Z_GRID)
    return {m: capacitor_stress(build(fit.values, fit.topo, m), loads, Z_GRID, d.p.drive.amp_v_rms) for m in modes}


def _part_rows(fit: PassiveFit, d: Design, modes: tuple[str, ...] = BOTH_MODES) -> list[dict]:
    """One row per part of the passive netlist (build order, drivers left out), the resistor rating, an inductor's
    peak current and a capacitor's peak voltage and ripple current being the largest of the modes' (the series
    alternative is rated in passive mode only: semi-active uses the parallel network's coax boards)."""
    ratings = _mode_ratings(fit, d, modes)
    currents = _mode_currents(fit, d, modes)
    caps = _mode_caps(fit, d, modes)
    rows = []
    for p in build(fit.values, fit.topo, "passive"):
        if p.kind == "Z":
            continue
        rating = max((r[p.name][2] for r in ratings.values() if p.name in r), default=None)
        third = max((r[p.name][0] for r in ratings.values() if p.name in r), default=None)
        value, unit = _value_unit(p)
        peak = max((c[p.name] for c in currents.values() if p.name in c), default=None)
        peak_v = max((c[p.name][0] for c in caps.values() if p.name in c), default=None)
        ripple = max((c[p.name][1] for c in caps.values() if p.name in c), default=None)
        rows.append({"part": p, "value": value, "unit": unit, "type": _type(p, rating, third, peak_v),
                     "rating": rating, "modes": {m: r.get(p.name) for m, r in ratings.items()}, "peak_a": peak,
                     "peak_v": peak_v, "ripple_a": ripple})
    return rows


def bound_notes(fit: PassiveFit) -> list[str]:
    """A sentence for each value at (or snapped just outside) a BOUNDS limit of the search."""
    out = []
    for n, v in fit.values.items():
        lo, hi = BOUNDS[n]
        c = fit.continuous.get(n, v)
        if v >= hi * (1 - 1e-3) or c >= hi * 0.99:
            end, limit = "upper", hi
        elif v <= lo * (1 + 1e-3) or c <= lo * 1.01:
            end, limit = "lower", lo
        else:
            continue
        s = f"{n} is at the {fmt_value(n[0], limit)} limit of the search (its {end} end)"
        if n in ("L_b1", "L_b2") and end == "upper":
            s += (": the network may be buying the impedance floor with bass loss; compare flat23 with the "
                  "active fit's.")
        elif (n, end) in _FADING:
            s += ": at this value the part barely acts (the fit wanted it further out of the way than the search "
            s += "allows), so this costs little."
        else:
            s += ": the fit wanted to go further than the search allows, so the response near this part is a "
            s += "compromise."
        out.append(s)
    return out


def _z_min(parts: list[Part], loads: dict[str, np.ndarray], f: np.ndarray,
           band: tuple[float, float] | None = None) -> tuple[float, float]:
    """min |Z_in| (ohm) over f (only inside band, Hz, when given) and where (Hz)."""
    z = np.abs(solve(parts, f, loads).z_in)
    m = np.ones(len(f), bool) if band is None else (f >= band[0]) & (f <= band[1])
    i = int(np.flatnonzero(m)[np.argmin(z[m])])
    return float(z[i]), float(f[i])


def _peak_text(p: Part, peak: float | None) -> str:
    """An inductor's peak current at full voltage, and the heat in its DCR at that current (sine, W)."""
    if peak is None:
        return "-"
    return (f"{peak:.2f} A" + (" (cored: saturation current at least this)" if p.name.startswith("L_b") else "")
            + f"; DCR heat {(peak / np.sqrt(2)) ** 2 * p.dcr:.1f} W at full V")


def _cap_text(r: dict) -> str:
    """A capacitor's peak voltage and ripple current at full voltage."""
    return "-" if r.get("peak_v") is None else f"{r['peak_v']:.1f} V peak, {r['ripple_a']:.2f} A rms"


def stopband_lines(fit: PassiveFit, d: Design) -> list[str]:
    """One sentence per far stopband of passive.stopbands: the worst |H| there re the branch's |H| at its crossover
    point, against STOP_DB (passive mode, nominal values, the driver loads)."""
    t, fb, fh = fit.topo, fit.bass_coax, fit.lf_hf
    g = stop_grid(fb, fh)
    h, _ = transfers(build(fit.values, t, "passive"), g, _loads(fit, d, g), t)
    n, out = len(Z_GRID), []
    for k, fc, m in stopbands(fb, fh):
        a = 20 * np.log10(np.abs(h[k]))
        rel = a[:n][m] - a[n if fc == fb else n + 1]
        j = int(np.argmax(rel))
        where = f"from {Z_GRID[m][0]:.0f} Hz up" if Z_GRID[m][0] > fc else f"up to {Z_GRID[m][-1]:.0f} Hz"
        out.append(f"electrical stopband, {_BRANCH_TEXT[k]} {where}: worst {rel[j]:+.1f} dB at {Z_GRID[m][j]:.0f} Hz "
                   f"re its level at {fc:g} Hz (limit {STOP_DB:g} dB){'' if rel[j] <= STOP_DB + 0.05 else ' OVER'}")
    return out


def peak_line(fit: PassiveFit, d: Design, tl: float, ways: Ways | None = None) -> str:
    """One sentence: the largest voltage gain any branch gives its driver (passive mode, nominal values, the driver
    loads, Z_GRID) against ELEC_GAIN_MAX_DB, and, with ways, the largest rise of the smoothed 23° sum above the fitted
    level below FIT_LO against PEAK_TOL_DB."""
    t = fit.topo
    h, _ = transfers(build(fit.values, t, "passive"), Z_GRID, _loads(fit, d, Z_GRID), t)
    g, k, j = largest_gain_db(h)
    out = (f"largest branch voltage gain: {g:+.1f} dB ({_BRANCH_TEXT[k]}) at {Z_GRID[j]:.0f} Hz (limit "
           f"+{ELEC_GAIN_MAX_DB:g} dB){'' if g <= ELEC_GAIN_MAX_DB + 0.05 else ' OVER'}")
    if ways is None:
        return out
    rise = low_rise_db(ways, fit, lambda x: impedances(d, x), tl)
    if rise is None:
        return out + f"; no response data below {FIT_LO:g} Hz"
    r, at = rise
    return (out + f"; largest rise of the 23° sum above the fitted level below {FIT_LO:g} Hz: {r:+.1f} dB at "
                  f"{at:.0f} Hz (limit +{PEAK_TOL_DB:g} dB){'' if r <= PEAK_TOL_DB + 0.05 else ' OVER'}")


def _draws(tol) -> int:
    return len(tol.draws_z_min) or P.TOL_DRAWS


def _wiring(t: Topology) -> str:
    return f"woofers in {t.bass_wiring}"


def _polarity_text(t: Topology) -> str:
    def one(name: str, inv: bool) -> str:
        return (f"{name} inverted (its + terminal to the network's ground side, its - terminal to the network "
                f"output)" if inv else f"{name} normal")
    return f"{one('bass', t.bass_inverted)}; {one('coax HF', t.hf_inverted)}; coax LF normal"


def _options_text(t: Topology) -> str:
    on = [OPTION_TEXT[k] for k, v in (("bass_comp", t.bass_comp), ("lf_shelf", t.lf_shelf), ("zobel", t.zobel),
                                      ("bass_order", t.bass_order == 3), ("lf_order", t.lf_order == 3),
                                      ("lf_notch", t.lf_notch), ("hf_shunt", t.hf_shunt)) if v]
    return ", ".join(on) if on else ("none (2nd-order bass and coax LF, no compensator, contour, Zobel, notch or "
                                     "shunt)")


def _table(head: list[str], rows: list[list[str]]) -> list[str]:
    w = [max(len(r[i]) for r in [head, *rows]) for i in range(len(head))]
    return ["  ".join(c.ljust(w[i]) for i, c in enumerate(r)).rstrip() for r in [head, *rows]]


_DRIVER_TEXT = {"bass": "the woofers", "coax_lf": "coax LF driver", "coax_hf": "coax HF driver"}


def _driver_nodes(p: Part, t: Topology) -> str:
    """'+ node to - node' of a driver; an inverted one has its + on the ground side."""
    inv = {"bass": t.bass_inverted, "coax_hf": t.hf_inverted}.get(p.name, False)
    a, b = (p.b, p.a) if inv else (p.a, p.b)
    return f"+ {a}, - {b}" + (" (inverted)" if inv else "")


def _netlist_table(fit: PassiveFit, d: Design, mode: str, rows: list[dict]) -> list[str]:
    """The mode's netlist in build order, drivers included (value "driver", their + and - nodes)."""
    by_name = {r["part"].name: r for r in rows}
    out = []
    for p in build(fit.values, fit.topo, mode):
        if p.kind == "Z":
            out.append([p.name, "driver", _DRIVER_TEXT[p.name], "-", "-", "-", "-", _driver_nodes(p, fit.topo), "-",
                        "-"])
            continue
        r = by_name[p.name]
        m = r["modes"][mode] if p.kind == "R" else None
        a = "cx (COAX)" if mode == "semi_active" and p.a == "in" else p.a  # semi-active: the amp is on COAX at cx
        out.append([p.name, fmt_value(p.kind, p.value), r["type"], TOLERANCE_TEXT,
                    f"{p.dcr:.2f} Ω" if p.kind == "L" else "-",
                    f"{m[0]:.2f} W" if m else "-", f"{m[1]:.1f} W" if m else "-", f"{a}-{p.b}",
                    _cap_text(r), _peak_text(p, r["peak_a"])])
    return _table(["part", "value", "type", "tolerance", "DCR", "power at 1/3 V", "power at full V", "nodes",
                   "capacitor at full V", "peak current at full V"], out)


def _terminals(m: dict) -> str:
    """One mode's terminals of link_states as a sentence: where each + and - lands, and what it feeds."""
    return "; ".join(f"{terminal_text(x)}, feeding {x['feeds']}" for x in m["terminals"])


def _links_lines(t: Topology) -> list[str]:
    out = []
    for mode, m in link_states(t).items():
        out.append(f"{mode} ({m['pairs']}): {_terminals(m)}. Links: {links_text(m)}. Polarity: {m['polarity']}.")
    out += ["('either': that link is not in the signal path in this mode.)", ""]
    out += [f"{k}: {v}" for k, v in LINK_TEXT.items()]
    return out


def _network_block(res: PassiveResult, d: Design, title: str, *, semi: bool = True,
                   sens: str | None = None, peak: str | None = None) -> list[str]:
    """One network's section of network.txt; semi False (the series alternative) replaces the semi-active netlist
    with a pointer to the parallel network's boards. sens: the sensitivity sentence, when known; peak: peak_line."""
    fit, tol, t = res.fit, res.tolerance, res.fit.topo
    v = d.p.drive.amp_v_rms
    loads = _loads(fit, d, Z_GRID)
    z, zf = _z_min(build(fit.values, t, "passive"), loads, Z_GRID)
    zbb, zbf = _z_min(build(fit.values, t, "passive"), loads, Z_GRID, BASS_BAND)
    zs, zsf = _z_min(build(fit.values, t, "semi_active"), loads, Z_GRID)
    zb = np.abs(loads["bass"])
    rows = _part_rows(fit, d, BOTH_MODES if semi else ("passive",))
    margin = z - Z_FLOOR
    out = [
        f"== {title}: {_wiring(t)} ==",
        f"crossover points: bass to coax {fit.bass_coax:g} Hz, coax LF to HF {fit.lf_hf:g} Hz",
        f"polarity: {_polarity_text(t)}",
        f"options kept: {_options_text(t)}",
        f"impedance floor used in the fit: {fit.z_floor:g} Ω (the spec minimum is {Z_FLOOR:g} Ω)",
        f"minimum impedance: {z:.2f} Ω at {zf:.0f} Hz ({'+' if margin >= 0 else ''}{margin:.2f} Ω against "
        f"{Z_FLOOR:g} Ω) | flat23 {tol.flat23_nominal:.2f} dB (trusted, nominal values)",
        f"minimum impedance in the bass band ({BASS_BAND[0]:g}-{BASS_BAND[1]:g} Hz): {zbb:.2f} Ω at {zbf:.0f} Hz",
        *stopband_lines(fit, d),
        *([peak] if peak else []),
        f"tolerance check ({_draws(tol)} draws, ±{TOL_SPREAD * 100:g} % on every L and C): "
        f"{'passed' if tol.passed else 'FAILED'}; worst draw {tol.z_min:.2f} Ω (limit {TOL_Z_MIN:g} Ω); flat23 "
        f"{tol.flat23_nominal:.2f} dB nominal, {tol.flat23_p95:.2f} dB at the 95th percentile (limit "
        f"+{TOL_FLAT_DB:g} dB); whole_rear {tol.whole_rear[0]:.2f}-{tol.whole_rear[1]:.2f} dB (5th-95th percentile)",
        *([sens] if sens else []),
    ]
    out += [f"WARNING: {s}" for s in bound_notes(fit)]
    out += ["", f"-- Passive mode (one amp channel, {v:g} V rms, on IN) --", *_netlist_table(fit, d, "passive", rows)]
    if semi:
        semi_m = link_states(t)["Semi-active"]
        out += ["", f"-- Semi-active mode ({_terminals(semi_m)}; links {links_text(semi_m)}; {v:g} V rms assumed on "
                    f"COAX) --",
                f"coax channel minimum impedance {zs:.2f} Ω at {zsf:.0f} Hz; bass channel (bare {_wiring(t)}) "
                f"{float(zb.min()):.2f} Ω at {float(Z_GRID[np.argmin(zb)]):.0f} Hz",
                *_netlist_table(fit, d, "semi_active", rows)]
    else:
        out += ["", "-- Semi-active mode --",
                "Semi-active mode uses the parallel network's coax boards; semi_active.txt is fitted to those."]
    full_m = link_states(t)["Fully active"]
    out += ["", f"-- Fully active mode: no network parts in circuit ({_terminals(full_m)}; links "
                f"{links_text(full_m)}) --",
            "", "-- Links to fit for each mode --", *_links_lines(t),
            "", "-- How the search got here --", *(f"- {n}" for n in res.notes), ""]
    return out


def _header(d: Design, tl: float, what: str) -> list[str]:
    v = d.p.drive.amp_v_rms
    return [
        f"# {what} (SIMULATED design, not measured)",
        f"# design file {d.path.name}, sha256 {_params_sha(d)[:8]}; trusted limit {tl:.0f} Hz",
        "# one network per speaker; parts.csv lists the parts for the pair",
        "# values are standard values: capacitors E12, resistors E24, inductors in 0.1 mH steps below 1 mH and "
        "0.5 mH steps from 1 mH",
        f"# tolerance: {TOLERANCE_TEXT} on every part (the tolerance check draws ±{TOL_SPREAD * 100:g} % on every "
        "L and C)",
        "# DCR: the coil resistance the model assumes (air core 0.5 Ω × √(L / 1 mH), cored 0.14 Ω × √(L / 1 mH)); "
        "a coil with a lower DCR is fine, a higher one changes the response",
        "# capacitors: polypropylene film in the coax network; bipolar electrolytic only in the bass network above "
        f"{ELECTROLYTIC_MIN * 1e6:g} µF; voltage rating (in the type) the first of "
        f"{', '.join(str(r) for r in CAP_RATINGS_V)} V at or above {CAP_MARGIN:g} × the peak voltage at full amp "
        "voltage (sine, 20 Hz-20 kHz, the larger of the modes); the ripple current is the rms current at full voltage",
        f"# power: the largest sine dissipation 20 Hz-20 kHz with 1/3 of the amp's {v:g} V rms on the input (about a "
        "tenth of full power, roughly what music averages) and with the full voltage; a resistor's rating is the "
        f"first of {', '.join(str(r) for r in RATINGS_W)} W at or above the 1/3-voltage figure; for the parallel "
        "network the larger of the passive and semi-active modes, for the series alternative passive mode only "
        f"(semi-active uses the parallel network's coax boards); where the 1/3-voltage figure reaches {BUY_NEAR:.0%} of "
        "the rating, buy 2-2.5 times the rating: a sustained low bass note is close to a sine; inductors: the heat "
        "in the DCR at the peak current's frequency, at full voltage",
        "# nodes: in = the IN terminal (amp +), 0 = ground (amp -); see vituixcad_netlist.txt for every node",
        "",
    ]


def network_text(result: PassiveResult, d: Design, tl: float, *, ways: Ways | None = None,
                 active_total: np.ndarray | None = None) -> str:
    """The network sheet: per network (the chosen one, and the woofers-in-series alternative when attached), the
    points, polarity, options, the floor used, the minimum impedance beside flat23, the stopbands, the largest branch
    gain (and with ways the rise below FIT_LO), the tolerance result, the sensitivity (with ways), any value at a
    search limit, the netlist of each mode with type, tolerance, DCR and power, the links per mode and the search
    notes. active_total (with ways): the active crossover's total, whose
    sensitivity is given for comparison."""
    out = _header(d, tl, "Passive crossover network")
    sens = {}
    if ways is not None:
        def z(x: np.ndarray) -> dict[str, np.ndarray]:
            return impedances(d, x)
        for r in [result, *([result.series] if result.series is not None else [])]:
            sens[r.fit.topo.bass_wiring] = sensitivity_text(ways, passive_total(ways, r.fit, z), r.fit.topo)
        if active_total is not None:
            out += [f"# active crossover, for comparison: {sensitivity_text(ways, active_total, Topology())}", ""]
    out += _network_block(result, d, "Chosen network" if result.series is None else
                          f"Network with the woofers in parallel (does not clear {Z_FLOOR:g} Ω)",
                          sens=sens.get(result.fit.topo.bass_wiring), peak=peak_line(result.fit, d, tl, ways))
    if result.series is not None:
        out += _network_block(result.series, d, "Alternative: the woofers in series", semi=False,
                              sens=sens.get(result.series.fit.topo.bass_wiring),
                              peak=peak_line(result.series.fit, d, tl, ways))
    return "\n".join(out)


def parts_csv(result: PassiveResult, d: Design) -> str:
    """The parts for a pair of speakers: part, value, unit, type (with a capacitor's voltage rating and a
    resistor's power rating), tolerance, DCR (ohm, inductors), peak current (A, inductors: the sine peak at the amp's
    full voltage, the larger of the modes'; a cored coil's saturation current must be at least this), peak voltage
    (V) and ripple current (A rms) of a capacitor at full voltage (the larger of the modes'), rating (W, resistors),
    quantity, and the network (parallel or series woofers) it belongs to."""
    lines = ["part,value,unit,type,tolerance,dcr_ohm,peak_a,peak_v,ripple_a,rating_w,qty_pair,network"]
    for res, modes in [(result, BOTH_MODES), *([(result.series, ("passive",))] if result.series is not None else [])]:
        for r in _part_rows(res.fit, d, modes):
            p = r["part"]
            lines.append(",".join([p.name, r["value"], r["unit"], r["type"], TOLERANCE_TEXT,
                                   f"{p.dcr:.3f}" if p.kind == "L" else "",
                                   f"{r['peak_a']:.2f}" if r["peak_a"] is not None else "",
                                   f"{r['peak_v']:.1f}" if r["peak_v"] is not None else "",
                                   f"{r['ripple_a']:.2f}" if r["ripple_a"] is not None else "", str(r["rating"] or ""), "2",
                                   res.fit.topo.bass_wiring]))
    return "\n".join(lines) + "\n"


_NODES = [
    "# nodes: in = the IN terminal (amp +); 0 = ground (amp -, and every driver's return through its two-pole link "
    "K3, K4 or K5)",
    "#   b1, b2 = bass low-pass (the woofers hang on the last); bz = between the Zobel's R and C",
    "#   bc1, bc2 = inside the impedance-peak compensator across the woofers (R_bc, L_bc, C_bc in series)",
    "#   cx = after the coax high-pass (C_h1 from in, L_h1 to ground); the coax LF and HF networks start here",
    "#   l1, l2 = coax LF low-pass; ln = after the notch (L_ln, C_ln, R_ln in parallel); lp = coax LF driver +",
    "#   t1, t2 = coax HF high-pass; tp = coax HF driver +; tr1, tr2 = inside the RLC shunt across the HF",
    "# drivers (kind Z): listed + node then - node; an inverted driver is listed with its terminals swapped",
    "# inductors: Rs = the series resistance (DCR) the model assumes; enter it in VituixCAD's coil",
    "# semi-active mode: the same parts with C_h1 and L_h1 out of circuit and the coax amp channel on node cx",
]


def vituixcad_netlist(result: PassiveResult) -> str:
    """The passive netlist, one line per part in build order, `name kind value node_a node_b` (and Rs for coils),
    per network (the woofers-in-series alternative too, when attached), under a header that explains the nodes."""
    out = ["# Passive crossover netlist for typing into VituixCAD (SIMULATED design, not measured)", *_NODES]
    for res in [result, *([result.series] if result.series is not None else [])]:
        t = res.fit.topo
        inv = {"bass": t.bass_inverted, "coax_hf": t.hf_inverted}
        out += ["", f"# network: {_wiring(t)} (bass to coax {res.fit.bass_coax:g} Hz, coax LF to HF "
                    f"{res.fit.lf_hf:g} Hz)"]
        for p in build(res.fit.values, t, "passive"):
            if p.kind == "Z":
                a, b = (p.b, p.a) if inv.get(p.name) else (p.a, p.b)
                out.append(f"{p.name} Z driver {a} {b}")
                continue
            value, unit = _value_unit(p)
            line = f"{p.name} {p.kind} {value}{unit} {p.a} {p.b}"
            out.append(line + (f" Rs={p.dcr:.3f}ohm" if p.kind == "L" else ""))
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------- report

def sensitivity(ways: Ways, total: np.ndarray) -> float:
    """Mean dB SPL at 23° (2 m, 2.83 V on the input) of total over SENS_BAND."""
    m = (ways.f >= SENS_BAND[0]) & (ways.f <= SENS_BAND[1])
    return float(np.mean(_spl(total[m, ways.obs.at23_idx])))


def sensitivity_text(ways: Ways, total: np.ndarray, topo: Topology) -> str:
    """One sentence: the sensitivity of total and how far it sits below the woofers' own level (wired as topo),
    the level the passive fit anchors to (before its LEVEL_SLACK_DB)."""
    s = sensitivity(ways, total)
    out = f"Sensitivity {s:.1f} dB at 2 m for 2.83 V (23°, {SENS_BAND[0]:g}-{SENS_BAND[1]:g} Hz)"
    w = woofer_level_db(ways, topo)
    if w is None:
        return out + "."
    gap = round(float(w - 20 * np.log10(P_REF_PEAK) - s), 1) + 0.0  # + 0.0: no "-0.0"
    where = f"{gap:.1f} dB below" if gap >= 0 else f"{-gap:.1f} dB above"
    return out + f", {where} what the woofers alone give ({_wiring(topo)})."



def _trusted(scores: list[Score]) -> dict[str, Score]:
    return {s.key: s for s in scores if s.confidence == "trusted" and math.isfinite(s.value)}


def shortfalls(worse: list[Score], ref: list[Score], name: str, ref_name: str) -> list[str]:
    """Plain sentences for each trusted score where `worse` is behind `ref` by SHORTFALL_DB or more (one count
    for rolloff_monotonic)."""
    a, b = _trusted(worse), _trusted(ref)
    out = []
    for k, s in a.items():
        if k not in b:
            continue
        gap = (b[k].value - s.value) if k in HIGHER_BETTER else (s.value - b[k].value)
        if gap < (1 if s.unit == "count" else SHORTFALL_DB) - 1e-9:
            continue
        better = "higher" if k in HIGHER_BETTER else "lower"
        u = "" if s.unit == "count" else f" {s.unit}"
        out.append(f"{s.label} ({k}): {name} {s.value:.2f}{u} against {ref_name} {b[k].value:.2f}{u}, "
                   f"{gap:.2f}{u} worse ({better} is better)"
                   + (f": {MEANING[k]}" if k in MEANING else "")
                   + (" (an active crossover can delay the bass to line it up with the cardioid coax; a passive "
                      "network cannot)" if k == "whole_rear" and name == "passive" else "") + ".")
    return out


def _way_curves(f: np.ndarray, ways: Ways, terms: dict[str, np.ndarray], total: np.ndarray) -> dict[str, np.ndarray]:
    """Each way at 23° through its transfer, in dB SPL (NaN where the way is exactly zero), and the sum."""
    i = ways.obs.at23_idx
    out = {}
    for k, h in terms.items():
        p = h * ways.ways[k][:, i]
        out[k] = _spl(np.where(np.abs(p) > 0, p, np.nan))
    out["sum"] = _spl(total[:, i])
    return out


def _score_table(cols: dict[str, list[Score]]) -> str:
    tr = {m: _trusted(s) for m, s in cols.items()}
    keys = list(dict.fromkeys(k for s in tr.values() for k in s))
    e = html.escape
    rows = []
    for k in keys:
        first = next(t[k] for t in tr.values() if k in t)
        cells = [f"<td><code>{e(k)}</code><br>{e(first.label)}"
                 f"{' (higher is better)' if k in HIGHER_BETTER else ''}</td>", f"<td>{e(first.unit)}</td>"]
        for m in cols:
            s = tr[m].get(k)
            cells.append(f'<td class="num {e(s.status)}">{s.value:.2f}</td>' if s else "<td>-</td>")
        p, a = tr["passive"].get(k), tr["active"].get(k)
        behind = (a.value - p.value if k in HIGHER_BETTER else p.value - a.value) if p and a else None
        cells.append(f'<td class="num">{behind:+.2f}</td>' if behind is not None else "<td>-</td>")
        rows.append(_row(cells))
    head = "".join(f"<th>{e(MODES.get(m, m))}</th>" for m in cols)
    return (f"<table><tr><th>Score</th><th>Unit</th>{head}<th>Passive behind active</th></tr>" + "".join(rows)
            + "</table>")


def _parts_html(fit: PassiveFit, d: Design, modes: tuple[str, ...] = BOTH_MODES) -> str:
    e = html.escape
    rows = []
    for r in _part_rows(fit, d, modes):
        p = r["part"]
        m = [x for x in r["modes"].values() if x]
        power = (f"{max(x[0] for x in m):.2f} W / {max(x[1] for x in m):.1f} W" if m else "-")
        rows.append(_row([f"<td>{e(p.name)}</td>", f'<td class="num">{e(fmt_value(p.kind, p.value))}</td>',
                          f"<td>{e(r['type'])}</td>", f"<td>{TOLERANCE_TEXT}</td>",
                          f'<td class="num">{p.dcr:.2f} Ω</td>' if p.kind == "L" else "<td>-</td>",
                          f'<td class="num">{power}</td>', f"<td>{e(p.a)} – {e(p.b)}</td>",
                          f"<td>{e(_cap_text(r))}</td>", f"<td>{e(_peak_text(p, r['peak_a']))}</td>"]))
    return ("<table><tr><th>Part</th><th>Value</th><th>Type</th><th>Tolerance</th><th>DCR</th>"
            "<th>Power at 1/3 V / full V</th><th>Nodes</th><th>Capacitor at full V</th><th>Peak current at full V</th>"
            "</tr>" + "".join(rows)
            + "</table>")


def _links_html(t: Topology) -> str:
    e = html.escape
    ks = list(LINK_TEXT)
    rows = [_row([f"<td><b>{e(mode)}</b><br>{e(m['pairs'])}</td>",
                  "<td>" + "<br>".join(f"{e(terminal_text(x))}: {e(x['feeds'])}" for x in m["terminals"])
                  + f"<br><i>{e(m['polarity'][0].upper() + m['polarity'][1:])}.</i></td>",
                  *(f"<td>{e(m['K'][k])}</td>" for k in ks)]) for mode, m in link_states(t).items()]
    return ("<table><tr><th>Mode</th><th>Terminals</th>" + "".join(f"<th>{k}</th>" for k in ks) + "</tr>"
            + "".join(rows) + "</table><ul>" + "".join(f"<li><b>{k}</b>: {e(v)}</li>" for k, v in LINK_TEXT.items())
            + "</ul>")


def _summary(res: PassiveResult, d: Design, label: str, sens: str | None = None, peak: str | None = None) -> str:
    """One network's key facts as an HTML list: points, polarity, options, z_min beside flat23, tolerance, and the
    sensitivity sentence and peak_line when given."""
    e = html.escape
    fit, tol, t = res.fit, res.tolerance, res.fit.topo
    z, zf = _z_min(build(fit.values, t, "passive"), _loads(fit, d, Z_GRID), Z_GRID)
    zbb, zbf = _z_min(build(fit.values, t, "passive"), _loads(fit, d, Z_GRID), Z_GRID, BASS_BAND)
    ok = z >= Z_FLOOR
    items = [
        f"<b>{e(label)}</b>: {e(_wiring(t))}; bass to coax {fit.bass_coax:g} Hz, coax LF to HF {fit.lf_hf:g} Hz; "
        f"{e(_polarity_text(t))}.",
        f"Options kept: {e(_options_text(t))}.",
        f'Minimum impedance <span class="{"pass" if ok else "fail"}">{z:.2f} Ω at {zf:.0f} Hz</span> '
        f"({z - Z_FLOOR:+.2f} Ω against the {Z_FLOOR:g} Ω minimum; the fit used a {fit.z_floor:g} Ω floor), "
        f"and <b>flat23 {tol.flat23_nominal:.2f} dB</b> (read the two together: a network can clear the floor by "
        "making the bass branch high-impedance, at a cost in flatness). In the bass band "
        f"({BASS_BAND[0]:g}-{BASS_BAND[1]:g} Hz) the minimum is {zbb:.2f} Ω at {zbf:.0f} Hz.",
        f'Part tolerance: <span class="{"pass" if tol.passed else "fail"}">'
        f"{'passed' if tol.passed else 'failed'}</span> (worst draw {tol.z_min:.2f} Ω; flat23 "
        f"{tol.flat23_p95:.2f} dB at the 95th percentile).",
        *([e(peak[0].upper() + peak[1:] + ".")] if peak else []),
        *([e(sens)] if sens else []),
    ]
    items += [f'<span class="warn">{e(s)}</span>' for s in bound_notes(fit)]
    return "<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>"


def _tolerance_html(res: PassiveResult) -> str:
    tol = res.tolerance
    rows = [("Result", f'<span class="{"pass" if tol.passed else "fail"}">{"passed" if tol.passed else "failed"}'
                       "</span>"),
            ("Worst draw's minimum impedance", f"{tol.z_min:.2f} Ω (no draw may go under {TOL_Z_MIN:g} Ω)"),
            ("flat23, nominal values", f"{tol.flat23_nominal:.2f} dB"),
            ("flat23, 95th percentile", f"{tol.flat23_p95:.2f} dB (at most {TOL_FLAT_DB:g} dB above nominal)"),
            ("whole_rear, 5th-95th percentile", f"{tol.whole_rear[0]:.2f} – {tol.whole_rear[1]:.2f} dB"),
            ("Floor the fit used", f"{res.fit.z_floor:g} Ω")]
    return "<table class='meta'>" + "".join(_row([f"<td>{k}</td>", f"<td>{v}</td>"]) for k, v in rows) + "</table>"


def _semi_html(sa, topo: Topology) -> str:
    rows = []
    e = html.escape
    for ch in ("bass", "coax"):
        wf = getattr(sa, ch)
        eq = "<br>".join(f"{p.f0:.0f} Hz, {p.gain_db:+.1f} dB, Q {p.q:.2f}" for p in wf.peq) or "none"
        rows.append(_row([f"<td>{ch}</td>", f'<td class="num">{wf.gain_db:+.2f} dB</td>',
                          f'<td class="num">{wf.delay_s * 1e3:.2f} ms</td>', f'<td class="num">{wf.polarity:+d}</td>',
                          f"<td>{eq}</td>"]))
    m = link_states(topo)["Semi-active"]
    return (f"<p>The DSP splits bass and coax with an LR4 at {sa.bass_coax:g} Hz; the coax channel feeds the same "
            f"passive coax boards ({e(_terminals(m))}; links {e(links_text(m))}). "
            f"{e(m['polarity'][0].upper() + m['polarity'][1:])}. Delays are relative.</p>"
            "<table><tr><th>Channel</th><th>Gain</th><th>Delay</th><th>Polarity</th><th>PEQ (f0, gain, Q)</th></tr>"
            + "".join(rows) + '</table><p>Setup sheet with biquads: <a href="semi_active.txt">semi_active.txt</a>.</p>')


def run_passive(params_path: Path | None, out_dir: Path, *, trusted_limit: tuple[float, str] | None) -> Path:
    """Design the passive crossover for the design at params_path and write out_dir/report.html with network.txt,
    parts.csv, semi_active.txt, vituixcad_netlist.txt, scores.json, the schematics and plots/. Returns the report.

    The fields come from design_ways (out_dir's mesh_check.json sets the trusted limit when trusted_limit is
    None; the command line always passes it). The passive search takes a minute or two on a real design."""
    out_dir = Path(out_dir)
    t0 = time.perf_counter()
    dw = design_ways(params_path, trusted_limit=trusted_limit, out_dir=out_dir)
    d, f, ways, tl, tl_text = dw.d, dw.f, dw.ways, dw.tl, dw.tl_text

    def z(x: np.ndarray) -> dict[str, np.ndarray]:
        return impedances(d, x)

    _log(f"passive: searching the network ({len(f)} frequencies, trusted to {tl:.0f} Hz; a minute or two)")
    t1 = time.perf_counter()
    res = search_passive(ways, z, tl, d.coax.hf_xo_min, z_floor=Z_FLOOR)
    _log(f"passive: search done, {time.perf_counter() - t1:.0f} s"
         + (" (with the woofers-in-series alternative)" if res.series is not None else ""))
    t1 = time.perf_counter()
    sa = fit_semi_active(ways, z, tl, res.fit)
    _log(f"semi-active: DSP stage at {sa.bass_coax:g} Hz, {time.perf_counter() - t1:.0f} s")
    t1 = time.perf_counter()
    xo = fit_active(ways, tl, fallback_db=_datasheet_fallback(), lf_hf_min=d.coax.hf_xo_min)
    _log(f"active: {xo.bass_coax:g} / {xo.lf_hf:g} Hz, {time.perf_counter() - t1:.0f} s")

    def scores_of(total: np.ndarray) -> list[Score]:
        return ways_scores(f, ways, total, spinorama(total, ways.obs, f=f), tl)

    fit, t = res.fit, res.fit.topo
    sw = semi_ways(ways, z, fit)
    totals = {"passive": passive_total(ways, fit, z), "semi_active": semi_total(sw, sa),
              "active": apply_active(ways, xo)}
    scores = {"passive": passive_scores(ways, fit, z, tl), "semi_active": scores_of(totals["semi_active"]),
              "active": scores_of(totals["active"])}
    series_scores = None if res.series is None else passive_scores(ways, res.series.fit, z, tl)
    sens = {"passive": sensitivity_text(ways, totals["passive"], t),
            "active": sensitivity_text(ways, totals["active"], Topology())}
    if res.series is not None:
        sens["series"] = sensitivity_text(ways, passive_total(ways, res.series.fit, z), res.series.fit.topo)
    peaks = {"passive": peak_line(fit, d, tl, ways),
             **({"series": peak_line(res.series.fit, d, tl, ways)} if res.series is not None else {})}

    _log("passive: writing outputs")
    out_dir.mkdir(parents=True, exist_ok=True)
    sha = _params_sha(d)
    (out_dir / "network.txt").write_text(network_text(res, d, tl, ways=ways, active_total=totals["active"]))
    (out_dir / "parts.csv").write_text(parts_csv(res, d))
    (out_dir / "semi_active.txt").write_text(semi_active_text(sa, tl, design=d.path.name, sha=sha[:8], topo=t))
    (out_dir / "vituixcad_netlist.txt").write_text(vituixcad_netlist(res))
    (out_dir / "scores.json").write_text(json.dumps(
        _jsonable({m: [dataclasses.asdict(s) for s in v] for m, v in scores.items()}), indent=2, allow_nan=False))
    schem = {}
    for r, suffix in [(res, ""), *([(res.series, "_series")] if res.series is not None else [])]:
        parts = build(r.fit.values, r.fit.topo, "passive")
        bass = [p for p in parts if _in_bass_network(p.name)]
        coax = [p for p in parts if not _in_bass_network(p.name)]
        w, rt = _wiring(r.fit.topo), r.fit.topo
        inv = {n for n, i in (("bass", rt.bass_inverted), ("coax_hf", rt.hf_inverted)) if i}
        draw_network(bass, f"Bass low-pass ({w})", out_dir / f"schematic_bass{suffix}.svg", inverted=inv)
        draw_network(coax, "Coax: high-pass to cx, then the LF and HF networks",
                     out_dir / f"schematic_coax{suffix}.svg", inverted=inv)
        schem[suffix] = (f"schematic_bass{suffix}.svg", f"schematic_coax{suffix}.svg")
    draw_links(t, out_dir / "links.svg")

    # plots
    pd, obs, blend = out_dir / "plots", ways.obs, ways.meta.get("puck_blend")
    h, _ = transfers(build(fit.values, t, "passive"), f, driver_loads(z(f), t), t)
    hs, _ = transfers(build(fit.values, t, "semi_active"), f, driver_loads(z(f), t), t)
    terms = {
        "passive": {"bass": bass_scale(t) * h["bass"], "coax_lf": h["coax_lf"], "coax_hf": h["coax_hf"]},
        "semi_active": {"bass": bass_scale(t) * _transfer(f, sa, "bass"),
                        "coax_lf": _transfer(f, sa, "coax") * hs["coax_lf"],
                        "coax_hf": _transfer(f, sa, "coax") * hs["coax_hf"]},
        "active": {k: way_transfer(f, xo, k) for k in ("bass", "coax_lf", "coax_hf")},
    }
    for m in MODES:
        plots.crossover_23(pd / f"crossover_23_{m}.png", f, _way_curves(f, ways, terms[m], totals[m]), tl, tl_text,
                           blend=blend, title=f"{MODES[m]}: each way at 23° after its crossover, and the sum")
    lv = _spl(totals["passive"])
    i0, i180 = _h(obs, 0), _h(obs, -180)
    for plane, ang, idx in (("H", obs.h_angles, obs.h_idx), ("V", obs.v_angles, obs.v_idx)):
        plots.polar_map(pd / f"polar_{plane.lower()}_passive.png", f, ang, lv[:, idx] - lv[:, [i0]],
                        "Horizontal" if plane == "H" else "Vertical", tl, tl_text, blend=blend)
    fr = f <= min(f[-1], 2000.0)
    rej = {MODES[m]: (_spl(x[:, i0]) - _spl(x[:, i180]))[fr] for m, x in totals.items()}
    fcs = {MODES["passive"]: fit.bass_coax, MODES["semi_active"]: sa.bass_coax, MODES["active"]: xo.bass_coax}
    if res.series is not None:
        st = passive_total(ways, res.series.fit, z)
        rej["Passive, woofers in series"] = (_spl(st[:, i0]) - _spl(st[:, i180]))[fr]
    plots.whole_rear_compare(pd / "whole_rear_compare.png", f[fr], rej, fcs)
    loads = driver_loads(z(Z_GRID), t)
    zc = {f"Passive, whole speaker ({_wiring(t)})": np.abs(solve(build(fit.values, t, "passive"), Z_GRID, loads).z_in),
          "Semi-active, coax channel": np.abs(solve(build(fit.values, t, "semi_active"), Z_GRID, loads).z_in),
          "Semi-active, bass channel (bare woofers)": np.abs(loads["bass"])}
    if res.series is not None:
        ts = res.series.fit.topo
        zc["Passive, whole speaker (woofers in series)"] = np.abs(
            solve(build(res.series.fit.values, ts, "passive"), Z_GRID, driver_loads(z(Z_GRID), ts)).z_in)
    plots.impedance(pd / "impedance.png", Z_GRID, zc, Z_FLOOR, f"{Z_FLOOR:g} Ω minimum for the SA550",
                    dashed=("Semi-active, bass channel (bare woofers)",))
    g = stop_grid(fit.bass_coax, fit.lf_hf)
    he, _ = transfers(build(fit.values, t, "passive"), g, driver_loads(z(g), t), t)
    n = len(Z_GRID)
    names = {"bass": f"bass ({_wiring(t)})", "coax_lf": "coax LF", "coax_hf": "coax HF"}
    db = {k: 20 * np.log10(np.abs(x)) for k, x in he.items()}
    limits = [(names[k], float(Z_GRID[m][0]), float(Z_GRID[m][-1]), float(db[k][n if fc == fit.bass_coax else n + 1])
               + STOP_DB) for k, fc, m in stopbands(fit.bass_coax, fit.lf_hf)]
    plots.electrical(pd / "electrical.png", Z_GRID, {names[k]: v[:n] for k, v in db.items()}, limits)
    tol = res.tolerance
    plots.tolerance_spread(pd / "tolerance.png", np.array(tol.draws_flat23), tol.flat23_nominal, tol.flat23_p95,
                           TOL_FLAT_DB,
                           np.array(tol.draws_z_min), TOL_Z_MIN)

    zb = np.abs(loads["bass"])
    page = _html(d, res=res, sa=sa, xo=xo, scores=scores, series_scores=series_scores, sens=sens, peaks=peaks,
                 tl=tl, tl_text=tl_text,
                 schem=schem, bass_min=(float(zb.min()), float(Z_GRID[np.argmin(zb)])),
                 wall=time.perf_counter() - t0)
    path = out_dir / "report.html"
    path.write_text(page)
    _log(f"passive: done in {time.perf_counter() - t0:.0f} s")
    return path


def _html(d: Design, *, res: PassiveResult, sa, xo, scores: dict[str, list[Score]], sens: dict[str, str],
          peaks: dict[str, str], series_scores: list[Score] | None, tl: float, tl_text: str,
          schem: dict[str, tuple[str, str]], bass_min: tuple[float, float], wall: float) -> str:
    e = html.escape
    fit, t = res.fit, res.fit.topo

    def fig(src: str, caption: str) -> str:
        return f'<img src="{e(src)}" alt="{e(Path(src).stem)}"><div class="cap">{caption}</div>'

    head = [("Params", f"{e(str(d.path))}<br><code>sha256 {_params_sha(d)}</code>"),
            ("Code commit", e(_commit())),
            ("Date", e(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))),
            ("Trusted 3D limit", e(tl_text)),
            ("Run time", f"{wall:.0f} s")]
    short = shortfalls(scores["passive"], scores["active"], "passive", "active")
    short_semi = shortfalls(scores["semi_active"], scores["active"], "semi-active", "active")
    warns = passive_warnings(xo, tl)

    def ul(items: list[str], empty: str) -> str:
        return "<ul>" + "".join(f"<li>{e(i)}</li>" for i in items) + "</ul>" if items else f"<p>{e(empty)}</p>"

    bz, bf = bass_min
    semi_warn = (f"In semi-active mode the bass amplifier channel drives the bare woofer pair ({_wiring(t)}) with no "
                 f"network in between. Its impedance falls to {bz:.2f} Ω at {bf:.0f} Hz"
                 + (f", under the {Z_FLOOR:g} Ω the SA550 is limited to" if bz < Z_FLOOR else "")
                 + ". Whatever amplifier drives that channel must handle this load. This is a warning, not a pass or "
                   "fail: that amplifier is not chosen yet.")
    series = ""
    if res.series is not None:
        sr = res.series
        series = (
            "<h2>Series fallback: the woofers in series</h2>"
            f"<p>The woofers-in-parallel network does not keep the speaker at or above {Z_FLOOR:g} Ω, so the best "
            "network with the two woofers wired in series is given too. In series the pair is four times the "
            "impedance but 6 dB quieter, so the coax is padded down further and the speaker plays about 6 dB less "
            "loud. You choose.</p>"
            + _summary(sr, d, "Series alternative", sens.get("series"), peaks.get("series"))
            + "<p>The semi-active results above use the parallel network's coax boards; the series alternative's "
              "coax boards differ (more padding).</p>"
            + (_score_table({"passive": series_scores, "active": scores["active"]}) if series_scores else "")
            + fig(schem["_series"][0], "Bass low-pass, woofers in series.")
            + fig(schem["_series"][1], "Coax network of the series alternative.")
            + _parts_html(sr.fit, d, ("passive",))
            + "<h3>How the search got here</h3>" + ul(sr.notes, "no notes"))
    parts = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        f"<title>Passive crossover</title><style>{_CSS}</style></head><body>",
        "<h1>Passive crossover report</h1>",
        "<p>A simulated design, not a measured one: the driver data are datasheet and model curves. The passive "
        "network is fitted to the same target as the active crossover and scored the same way, so the three modes "
        "can be compared directly.</p>",
        "<table class='meta'>" + "".join(_row([f"<td>{k}</td>", f"<td>{v}</td>"]) for k, v in head) + "</table>",
        "<h2>Summary</h2>",
        _summary(res, d, "Passive network", sens["passive"], peaks["passive"]),
        f"<p><b>Semi-active</b>: DSP split at {sa.bass_coax:g} Hz ahead of the same coax boards. "
        f"<b>Active</b> (for comparison): {xo.bass_coax:g} Hz and {xo.lf_hf:g} Hz. {e(sens['active'])}</p>",
        (f'<p class="warn">The woofers-in-parallel network does not clear {Z_FLOOR:g} Ω; see the series fallback '
         "below.</p>" if res.series is not None else ""),
        "<h2>Scores</h2>",
        "<p>Trusted scores only (up to the trusted limit); every score, low-confidence ones included, is in "
        "scores.json. Lower is better unless marked. 'Passive behind active': how much worse passive is "
        "(positive = passive worse, negative = passive better), whichever direction is better for the score.</p>",
        _score_table(scores),
        "<h2>Where passive falls short of active</h2>",
        ul(short, f"No trusted score is {SHORTFALL_DB:g} dB (or one band, for the count) or more behind active."),
        "<h3>Semi-active against active</h3>",
        ul(short_semi, f"No trusted semi-active score is {SHORTFALL_DB:g} dB (or one band, for the count) or more "
                       "behind active."),
        "<h3>What the active crossover does that a passive network cannot copy</h3>",
        ul(warns, "Nothing: the active fit asks for nothing a passive network cannot supply."),
        "<h2>Impedance</h2>",
        fig("plots/impedance.png", f"What the amplifier sees, modelled at nominal part values. The dotted line is "
                                   f"the {Z_FLOOR:g} Ω minimum. Dashed: the bare woofers, which the bass channel "
                                   "drives in semi-active mode."),
        f'<p class="warn"><b>Semi-active bass channel:</b> {e(semi_warn)}</p>',
        "<h2>Electrical transfer</h2>",
        fig("plots/electrical.png", f"Each branch's terminal voltage per input volt through the passive network, on "
                                    f"the modelled driver loads; the dashed lines are the far-stopband limits the fit "
                                    f"holds ({STOP_DB:g} dB re the branch's level at its crossover point, from 4 times "
                                    "past it), which the acoustic data cannot see above the trusted limit."),
        "<ul>" + "".join(f"<li>{e(x)}</li>" for x in stopband_lines(fit, d)) + "</ul>",
        "<h2>23° response</h2>",
        fig("plots/crossover_23_passive.png", "Passive: each way at the listening angle (23°) through the network, "
                                              "and the sum, smoothed to 1/6 octave."),
        fig("plots/crossover_23_semi_active.png", "Semi-active: the bass through the DSP, the coax LF and HF through "
                                                  "the DSP and the coax boards."),
        fig("plots/crossover_23_active.png", "Active, for comparison."),
        "<h2>Polar maps (passive)</h2>",
        fig("plots/polar_h_passive.png", "Horizontal: level re on-axis at 2 m. The dashed line is the trusted "
                                         "3D limit."),
        fig("plots/polar_v_passive.png", "Vertical (positive angles up), normalised to on-axis."),
        "<h2>Rear rejection</h2>",
        fig("plots/whole_rear_compare.png", "Info only: the whole speaker, front minus back, for each mode. A passive "
                                            "network cannot delay the bass, so the rear rejection near the "
                                            "bass-to-coax crossover changes."),
        "<h2>Tolerance</h2>",
        f"<p>{_draws(res.tolerance)} random draws with every inductor and capacitor off by up to "
        f"±{TOL_SPREAD * 100:g} % (resistors exact, coil DCR fixed).</p>",
        _tolerance_html(res),
        fig("plots/tolerance.png", "Left: flat23 over the draws. Right: each draw's minimum impedance."),
        "<h2>Schematics and links</h2>",
        fig(schem[""][0], "Bass low-pass, from IN to the woofers."),
        fig(schem[""][1], "Coax: the high-pass from IN to node cx, then the coax LF and HF networks."),
        fig("links.svg", "Terminals and removable links K1-K5 for each mode."),
        _links_html(t),
        "<h2>Parts (one speaker; order two of each)</h2>",
        _parts_html(fit, d),
        "<p>Power: the largest sine dissipation with 1/3 and with all of the amp's "
        f"{d.p.drive.amp_v_rms:g} V rms on the input, the larger of passive and semi-active mode. "
        'Files: <a href="network.txt">network.txt</a> (netlist per mode, links, notes), '
        '<a href="parts.csv">parts.csv</a> (the pair), '
        '<a href="vituixcad_netlist.txt">vituixcad_netlist.txt</a> (to type into VituixCAD).</p>',
        "<h2>Semi-active DSP</h2>",
        _semi_html(sa, t),
        series,
        "<h2>How the search got here</h2>",
        ul(res.notes, "no notes"),
        "</body></html>",
    ]
    return "\n".join(parts)
