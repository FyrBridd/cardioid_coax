"""``speaker directivity``: the 3D directivity run, its HTML report, and the mesh-convergence check.

The run solves the cabinet BEM (and, unless fast, the coax puck) one BEM band per call, so each
band is cached on its own and an interrupted run resumes where it stopped. Ways, crossover,
spinorama and scores then follow speaker.system, speaker.crossover, speaker.spinorama and
speaker.scores. Plots are in speaker.directivity_plots.
"""
from __future__ import annotations

import dataclasses
import datetime
import hashlib
import html
import json
import math
import subprocess
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import directivity_plots as plots
from .bem.mesh import cabinet_mesh
from .bem.run import BAND_EDGES, PUCK_F_RANGE, SourceField, cabinet_field, freq_grid, puck_field
from .bem.solver import BEM_MEMORY_BUDGET, estimate_bytes
from .cardioid import front_back, simulate
from .constants import C
from .crossover import MATCH_BANDS, FixedCrossover, apply_crossover
from .curves import interp_db, load_curve, smooth
from .drivers import Driver, load_driver
from .layout import Layout, cardioid_params, compute_layout
from .params import DEFAULT_PATH, Params, load_params
from .report import FB_AIM_BAND, FB_BAND
from .scores import LW_RING, Score, score_all
from .spinorama import P_REF_PEAK, Spin, spinorama
from .system import COAX_CURVES, V_DRIVE, Ways, coax_section, vent_scaled_section, way_responses

CACHE_DIR = Path("out/cache/bem")
CAB_H_MAX = 0.05  # m, cabinet mesh size cap of the full run and the mesh check
FULL_GRID = (50.0, 20000.0, 12)  # f_lo, f_hi, points per octave
FAST_GRID = (50.0, 1000.0, 6)
FAST_H_MAX = 0.08  # m
COARSE = 6.0  # elements per wavelength of the run the mesh check verifies
MAX_DEV_DB = 0.5  # mesh check: allowed |dB| change on refinement
UNVERIFIED_LIMIT = 3000.0  # Hz, trusted limit without a mesh check
VENT_SPREAD = (0.7, 1.3)  # vent resistance factors for the felt's unknown resistance
LUMPED_DIFF_DB = 3.0  # report both cardioid models where they differ by more
_REPO = Path(__file__).resolve().parents[2]


@dataclass
class _Design:
    path: Path
    text: str
    p: Params
    coax: Driver
    bass: Driver
    layout: Layout


def _design(params_path: Path | None) -> _Design:
    path = DEFAULT_PATH if params_path is None else Path(params_path)
    p = load_params(path)
    coax, bass = load_driver(p.coax), load_driver(p.bass)
    return _Design(path, path.read_text(), p, coax, bass, compute_layout(p, coax, bass))


def _concat(parts: list[SourceField]) -> SourceField:
    """Join per-band fields along frequency; meta["bands"] lists every band, array metas are joined."""
    first = parts[0]
    meta = {k: np.concatenate([x.meta[k] for x in parts]) for k in first.meta if k != "bands"}
    meta["bands"] = [b for x in parts for b in x.meta["bands"]]
    return SourceField(np.concatenate([x.f for x in parts]), first.names, first.obs,
                       np.concatenate([x.p for x in parts]), meta)


def _banded(solve: Callable[[np.ndarray], SourceField], f: np.ndarray, label: str) -> SourceField:
    """solve() once per BEM band of f, one progress line per band on stderr."""
    band = np.searchsorted(BAND_EDGES, f)
    parts = []
    for b in np.unique(band):
        t0 = time.perf_counter()
        parts.append(solve(f[band == b]))
        m = parts[-1].meta["bands"][0]
        print(f"{label}: band <= {BAND_EDGES[b]:g} Hz, {m['triangles']} elements, "
              f"{time.perf_counter() - t0:.1f} s", file=sys.stderr, flush=True)
    return _concat(parts)


def _datasheet_fallback() -> dict[str, float]:
    """coax_hf gain from the datasheet on-axis curves (LF over its matching band minus HF over
    its own), used only when the grid stops below the HF matching band."""
    def mean_db(way: str, band: tuple[float, float]) -> float:
        return float(np.mean(interp_db(load_curve(f"{COAX_CURVES}_{way}_spl"), np.geomspace(*band, 64))))
    return {"coax_hf": mean_db("lf", MATCH_BANDS["coax_lf"]) - mean_db("hf", MATCH_BANDS["coax_hf"])}


def _evaluate(d: _Design, cab: SourceField, puck: SourceField | None, f: np.ndarray, tl: float,
              xo: FixedCrossover) -> tuple[Ways, np.ndarray, dict[str, float], Spin]:
    ways = way_responses(d.p, d.coax, d.bass, d.layout, cab, puck, f, trusted_limit=tl)
    total, gains = apply_crossover(ways, xo, fallback_db=_datasheet_fallback())
    return ways, total, gains, spinorama(total, ways.obs, f=f)


def _spl(p: np.ndarray) -> np.ndarray:
    return 20 * np.log10(np.abs(p) / P_REF_PEAK)


def _h(obs, angle: float) -> int:
    return int(obs.h_idx[np.isclose(obs.h_angles, angle)][0])


def _jsonable(x):
    """Plain JSON types: tuples to lists, numpy scalars and arrays to Python, non-finite floats to None."""
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.generic):
        x = x.item()
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


# ---------------------------------------------------------------- mesh check

def _fine_dofs(p: Params, coax: Driver, bass: Driver, h: float) -> int:
    """P1 node count of the cabinet mesh at size h."""
    return len(np.unique(cabinet_mesh(p, coax, bass, h=h).triangles))


def mesh_check(params_path: Path | None, out_dir: Path, *, f_max: float = 4000.0, fine: float = 8.0) -> Path:
    """Compare the cabinet run at COARSE and fine elements per wavelength up to f_max.

    Per 1/3-octave band: the max |dB| difference over the grid points and over the at23, LW, ER,
    SP and H 0-30 curves after the crossover. The trusted limit is the upper edge of the highest
    band that, with every band below it, differs by <= MAX_DEV_DB, capped at the top solved
    frequency. A fine mesh over BEM_MEMORY_BUDGET stops the check at that band. If no band
    converges, the JSON is still written and ValueError is raised.
    """
    d = _design(params_path)
    f = freq_grid(FULL_GRID[0], f_max, FULL_GRID[2])
    band = np.searchsorted(BAND_EDGES, f)
    if band.max() >= len(BAND_EDGES):
        raise ValueError(f"f_max {f_max:g} Hz is above the last band edge, {BAND_EDGES[-1]:g} Hz")
    kw = dict(params_text=d.text, h_max=CAB_H_MAX, cache_dir=CACHE_DIR)
    coarse, refined, elements, capped = [], [], {}, False
    for b in np.unique(band):
        edge, fb = BAND_EDGES[b], f[band == b]
        h_c, h_f = (min(CAB_H_MAX, C / (edge * pw)) for pw in (COARSE, fine))
        same = h_f == h_c
        if not same:
            need = estimate_bytes(_fine_dofs(d.p, d.coax, d.bass, h_f))
            if need > BEM_MEMORY_BUDGET:
                capped = True
                print(f"mesh check: stopped at band <= {edge:g} Hz, its fine mesh needs about "
                      f"{need / 1e9:.1f} GB (budget {BEM_MEMORY_BUDGET / 1e9:.1f} GB)", file=sys.stderr)
                break
        t0 = time.perf_counter()
        c = cabinet_field(d.p, d.coax, d.bass, fb, per_wavelength=COARSE, **kw)
        r = c if same else cabinet_field(d.p, d.coax, d.bass, fb, per_wavelength=fine, **kw)
        coarse.append(c)
        refined.append(r)
        elements[f"{edge:g}"] = {"same_mesh": same, "coarse": c.meta["bands"][0], "fine": r.meta["bands"][0]}
        print(f"mesh check: band <= {edge:g} Hz, {c.meta['bands'][0]['triangles']} / "
              f"{r.meta['bands'][0]['triangles']} elements, {time.perf_counter() - t0:.1f} s",
              file=sys.stderr, flush=True)
    if not coarse:
        raise ValueError("mesh check: the first band's fine mesh is over the memory budget")

    fc, fr = _concat(coarse), _concat(refined)
    curves = []
    for fld in (fc, fr):
        _, total, _, spin = _evaluate(d, fld, None, fld.f, float(fld.f[-1]), FixedCrossover())
        lv = _spl(total)
        curves.append(np.stack([spin.at23, spin.lw, spin.er, spin.sp,
                                *(lv[:, _h(fld.obs, a)] for a in LW_RING)]))
    dev = np.abs(curves[0] - curves[1]).max(axis=0)
    k = np.floor(3 * np.log2(fc.f / 1000.0) + 0.5).astype(int)  # nominal 1/3-octave band index
    band_dev, limit, ok = {}, 0.0, True
    for q in np.unique(k):
        v = float(dev[k == q].max())
        band_dev[f"{1000.0 * 2.0 ** (q / 3):.0f}"] = v
        ok = ok and v <= MAX_DEV_DB
        if ok:
            limit = 1000.0 * 2.0 ** ((q + 0.5) / 3)
    out = {
        "trusted_limit_hz": min(limit, float(fc.f[-1])),  # never past the top solved frequency
        "band_max_dev_db": band_dev,
        "per_wavelength": [COARSE, float(fine)],
        "elements": elements,
        "capped_by_memory": capped,
        "f_max_hz": float(f_max),
        "h_max_m": CAB_H_MAX,
        "max_dev_db": MAX_DEV_DB,
    }
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "mesh_check.json"
    path.write_text(json.dumps(_jsonable(out), indent=2, allow_nan=False))
    if limit == 0.0:
        first = next(iter(band_dev.items()))
        raise ValueError(f"mesh check: no 1/3-octave band converged; the lowest ({first[0]} Hz) already differs "
                         f"by {first[1]:.2f} dB (> {MAX_DEV_DB:g} dB) between the meshes. See {path}")
    return path


# ---------------------------------------------------------------- report run

def _commit() -> str:
    def git(*a: str) -> str:
        return subprocess.run(["git", *a], cwd=_REPO, check=True, capture_output=True, text=True).stdout.strip()
    try:
        head = git("rev-parse", "--short", "HEAD")
        return head + (" (with uncommitted changes)" if git("status", "--porcelain", "--", "src/") else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _trusted_value(scores: list[Score], key: str) -> float:
    return next(s.value for s in scores if s.key == key and s.confidence == "trusted")


def _trusted_limit(mc: dict | None, *, cab_top: float, f_top: float, puck: bool) -> tuple[float, str]:
    """The trusted limit and its header text, from the mesh check or UNVERIFIED_LIMIT.

    When the grid goes past the cabinet field (cab_top), a limit above it is capped there. A limit
    below what the run needs (2 x the puck's lowest frequency for the cabinet-to-puck blend, or the
    top of a puck-less grid) raises ValueError.
    """
    tl = float(mc["trusted_limit_hz"]) if mc else UNVERIFIED_LIMIT
    note = "" if mc else " (unverified)"
    if f_top > cab_top and tl > cab_top:
        note += f" (from {tl:.0f} Hz, capped at the top of the cabinet field)"
        tl = cab_top
    need = 2 * PUCK_F_RANGE[0] if puck else f_top
    if tl < need:
        why = "the cabinet-to-puck blend starts at half of it" if puck else "a run without the puck stops at it"
        raise ValueError(f"trusted 3D limit {tl:.0f} Hz is below {need:.0f} Hz, which this run needs ({why})")
    return tl, f"{tl:.0f} Hz{note}"


def run_directivity(params_path: Path | None, out_dir: Path, *, fast: bool = False,
                    f_max_3d: float = 4000.0, per_wavelength: float = 6.0) -> Path:
    """Simulate, score and plot the design; write out_dir/report.html, scores.json and plots/."""
    d = _design(params_path)
    out_dir = Path(out_dir)
    mc_path = out_dir / "mesh_check.json"
    mc = json.loads(mc_path.read_text()) if mc_path.exists() else None
    f = freq_grid(*(FAST_GRID if fast else FULL_GRID))
    f_cab = f if fast else f[f <= f_max_3d]
    tl, tl_text = _trusted_limit(mc, cab_top=float(f_cab[-1]), f_top=float(f[-1]), puck=not fast)
    kw = dict(params_text=d.text, per_wavelength=per_wavelength, cache_dir=CACHE_DIR)
    t0 = time.perf_counter()
    if fast:
        cab = _banded(lambda x: cabinet_field(d.p, d.coax, d.bass, x, h_max=FAST_H_MAX, **kw), f_cab, "cabinet")
        puck = None
    else:
        cab = _banded(lambda x: cabinet_field(d.p, d.coax, d.bass, x, h_max=CAB_H_MAX, **kw), f_cab, "cabinet")
        fp = f[(f >= PUCK_F_RANGE[0]) & (f <= PUCK_F_RANGE[1])]
        puck = _banded(lambda x: puck_field(d.p, d.coax, x, **kw), fp, "puck")
    xo = FixedCrossover()
    ways, total, gains, spin = _evaluate(d, cab, puck, f, tl, xo)
    obs, sec = ways.obs, coax_section(ways)
    scores = score_all(f, total, sec, obs, spin, trusted_limit=tl, xo=xo)

    # Cardioid: the vent resistance spread (same fields), and the lumped model for comparison.
    cp = cardioid_params(d.p, d.layout, d.coax)
    spread = {s: vent_scaled_section(ways, cab, d.coax, cp, s) for s in VENT_SPREAD}
    spread_scores = {1.0: scores}
    spread_scores.update({s: score_all(f, total, x, obs, spin, trusted_limit=tl, xo=xo) for s, x in spread.items()})
    vent = {k: {s: _trusted_value(sc, k) for s, sc in spread_scores.items()} for k in ("cardioid_hard", "cardioid_aim")}
    i0, i180 = _h(obs, 0), _h(obs, -180)
    rej = lambda x: _spl(x[:, i0]) - _spl(x[:, i180])  # noqa: E731
    fr = f[f <= min(f[-1], 1000.0)]
    lumped = front_back(simulate(d.coax, cp, fr, V_DRIVE), cp)
    bem = rej(sec)[: len(fr)]
    m = (fr >= FB_BAND[0]) & (fr <= FB_AIM_BAND[1])
    max_diff = float(np.abs(bem - lumped)[m].max())

    plot_dir = out_dir / "plots"
    lv = _spl(total)
    for plane, ang, idx in (("H", obs.h_angles, obs.h_idx), ("V", obs.v_angles, obs.v_idx)):
        plots.polar_map(plot_dir / f"polar_{plane.lower()}.png", f, ang, lv[:, idx] - lv[:, [i0]],
                        "Horizontal" if plane == "H" else "Vertical", tl, tl_text)
    plots.spinorama(plot_dir / "spinorama.png", spin, f, tl, tl_text)
    ref = smooth(f, lv[:, obs.at23_idx], 1 / 6)
    plots.window_family(plot_dir / "window_family.png", f,
                        {a: smooth(f, lv[:, _h(obs, a)], 1 / 6) - ref for a in LW_RING}, tl, tl_text)
    title = "Coax section rear rejection (cone + vents, before the crossover)"
    if max_diff > LUMPED_DIFF_DB:
        title += f"\nlumped vs BEM differ by > {LUMPED_DIFF_DB:g} dB (up to {max_diff:.1f} dB in 150-500 Hz)"
    lo, hi = (rej(spread[s])[: len(fr)] for s in VENT_SPREAD)
    plots.rear_rejection(plot_dir / "rear_rejection.png", fr, bem, lumped, (np.minimum(lo, hi), np.maximum(lo, hi)),
                         f"x{VENT_SPREAD[0]:g}-x{VENT_SPREAD[1]:g}", title)
    plots.whole_rear(plot_dir / "whole_rear.png", f[f <= min(f[-1], 2000.0)],
                     rej(total)[f <= min(f[-1], 2000.0)], xo.bass_coax)

    (out_dir / "scores.json").write_text(
        json.dumps(_jsonable([dataclasses.asdict(s) for s in scores]), indent=2, allow_nan=False))
    page = _html(d, fast=fast, mc=mc, tl_text=tl_text, scores=scores, gains=gains, xo=xo, ways=ways,
                 vent=vent, cp=cp, max_diff=max_diff, fields={"cabinet": cab, "puck": puck},
                 wall=time.perf_counter() - t0)
    path = out_dir / "report.html"
    path.write_text(page)
    return path


# ---------------------------------------------------------------- HTML

_CSS = """
body { font-family: system-ui, sans-serif; max-width: 1100px; margin: 1.5em auto; padding: 0 16px;
       color: #222; background: #fff; line-height: 1.45; }
h1 { margin-bottom: 0.2em; } h2 { margin-top: 1.6em; border-bottom: 1px solid #ccc; }
table { border-collapse: collapse; font-size: 0.9em; margin: 0.5em 0; }
th, td { border: 1px solid #ccc; padding: 3px 7px; text-align: left; vertical-align: top; }
td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
td.formula { font-size: 0.85em; max-width: 420px; }
tr.low td { color: #999; }
.pass { color: #1a7f37; } .fail { color: #c62828; font-weight: bold; } .warn { color: #b26a00; }
tr.low .pass, tr.low .fail, tr.low .warn { color: #999; }
img { max-width: 100%; display: block; margin: 0.6em 0 0.2em; }
.cap { font-size: 0.9em; color: #555; margin-bottom: 1em; }
.meta td:first-child { font-weight: 600; }
"""


def _row(cells: list[str], cls: str = "") -> str:
    return f"<tr{f' class={cls!r}' if cls else ''}>" + "".join(cells) + "</tr>"


def _html(d: _Design, *, fast, mc, tl_text, scores, gains, xo, ways, vent, cp, max_diff, fields, wall) -> str:
    e = html.escape
    sha = hashlib.sha256(d.path.read_bytes()).hexdigest()
    head = [("Params", f"{e(str(d.path))}<br><code>sha256 {sha}</code>"),
            ("Code commit", e(_commit())),
            ("Date", e(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))),
            ("Design", "Round 0 — current design (baseline)"),
            ("Trusted 3D limit", e(tl_text)),
            ("Mode", "fast (coarse mesh, 50 Hz-1 kHz, no puck)" if fast else "full")]

    def fig(name: str, caption: str) -> str:
        return f'<img src="plots/{name}.png" alt="{e(name)}"><div class="cap">{caption}</div>'

    srows = []
    for s in scores:
        srows.append(_row([f"<td><code>{e(s.key)}</code><br>{e(s.label)}</td>", f'<td class="num">{s.value:.2f}</td>',
                           f"<td>{e(s.unit)}</td>", f'<td class="{e(s.status)}">{e(s.status)}</td>',
                           f'<td class="num">{s.band[0]:.0f}-{s.band[1]:.0f} Hz</td>', f"<td>{e(s.confidence)}</td>",
                           f'<td class="formula">{e(s.formula)}</td>'], "low" if s.confidence == "low" else ""))

    r = cp.vent_resistance
    if mc:
        how = (f"from mesh_check.json: {mc['per_wavelength'][0]:g} vs {mc['per_wavelength'][1]:g} elements per "
               f"wavelength agree within {mc.get('max_dev_db', MAX_DEV_DB):g} dB up to {tl_text}")
        if mc.get("capped_by_memory"):
            how += "; the check stopped early because the finer mesh would not fit in memory"
    else:
        how = f"{tl_text}: no mesh check in this output directory (run with --mesh-check)"
    mismatch = ways.meta["lf_join_mismatch_db"]
    notes = [
        f"<b>Trusted 3D limit</b> {how}. Scores are split at the limit; rows above it are grey "
        "(low confidence) and do not rank designs.",
        f"<b>Vent resistance</b> (the felt as built is unknown, ±30%; R = {r:.0f} Pa s/m): "
        + "; ".join(f"{k} (trusted) is " + ", ".join(f"{v[s]:.1f} dB at R × {s:g}" for s in (VENT_SPREAD[0], 1.0, VENT_SPREAD[1]))
                    for k, v in vent.items())
        + ". Only the lumped vent-to-cone volume-velocity ratio is taken at the scaled resistance, at every "
        "frequency the vents radiate (up to the trusted limit); the coax LF strength and the 3D vent fields "
        "are those of the nominal run, since the fields do not depend on the resistance.",
        f"<b>Cardioid models</b>: the 3D coax-section rear rejection and the lumped model differ by up to "
        f"{max_diff:.1f} dB over 150-500 Hz; both are plotted. The 3D model is the authority for the cardioid "
        "gate (spec §4); the lumped model ignores cabinet diffraction.",
        f"<b>Provisional cone geometry</b>: the coax cone depth ({d.coax.cone_depth * 1000:.0f} mm), HF exit "
        f"({d.coax.hf_exit_d * 1000:.1f} mm) and cutout ({d.coax.cutout_d * 1000:.0f} mm) are estimates until the "
        "driver is measured with a profile gauge. Above the trusted limit the coax's horn behaviour rests on them.",
        f"<b>LF join</b>: the lumped and datasheet coax LF volume velocities differ by {mismatch:+.1f} dB at "
        "400 Hz; the coax LF strength is blended between them across 300-500 Hz.",
    ]
    if fast:
        notes.append(f"<b>Fast mode</b>: coarse mesh (h_max {FAST_H_MAX * 1000:.0f} mm), 1/6-octave grid to 1 kHz, "
                     f"no puck. The coax HF gain ({gains['coax_hf']:+.1f} dB) comes from the datasheet on-axis curves, "
                     "because the grid stops below its 2-8 kHz matching band.")

    run_rows = []
    for name, fld in fields.items():
        for b in (fld.meta["bands"] if fld is not None else []):
            run_rows.append(_row([f"<td>{name}</td>", f'<td class="num">{b["edge"]:g} Hz</td>',
                                  f'<td class="num">{b["freqs"]}</td>', f'<td class="num">{b["h"] * 1000:.1f} mm</td>',
                                  f'<td class="num">{b["nodes"]}</td>', f'<td class="num">{b["triangles"]}</td>',
                                  f'<td class="num">{b["mesh_s"] + b["solve_s"]:.0f} s</td>']))

    parts = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        f"<title>Directivity round 0</title><style>{_CSS}</style></head><body>",
        "<h1>Directivity report — Round 0</h1>",
        "<table class='meta'>" + "".join(_row([f"<td>{k}</td>", f"<td>{v}</td>"]) for k, v in head) + "</table>",
        "<h2>Curves</h2>",
        fig("polar_h", "Horizontal polar map: level re on-axis at 2 m. The dashed line is the trusted 3D limit."),
        fig("polar_v", "Vertical polar map (positive angles up), normalised to on-axis."),
        fig("spinorama", "Spinorama at 2 m for 2.83 V into the coax LF way after the crossover. "
                         "Listening window: H 0/10/20/30° plus az 23° at el ±10°. DIs on the right axis."),
        fig("window_family", "H 0-30° in 5° steps minus the 23° listener curve, 1/6-octave smoothed; "
                             "±1 dB band and 2-5 kHz shaded."),
        fig("rear_rejection", "Coax section alone (cone plus vents, no bass, before the crossover filter): "
                              "front minus back. Shaded: 3D result with the vent strength taken at vent resistance ×0.7 to ×1.3."),
        fig("whole_rear", "Info only: the whole speaker after the crossover, front minus back."),
        "<h2>Scores</h2>",
        "<table><tr><th>Score</th><th>Value</th><th>Unit</th><th>Status</th><th>Band</th><th>Confidence</th>"
        "<th>Formula</th></tr>" + "".join(srows) + "</table>",
        "<h2>Confidence</h2><ul>" + "".join(f"<li>{n}</li>" for n in notes) + "</ul>",
        "<h2>Crossover</h2>",
        f"<p>Fixed and provisional: acoustic LR4 at {xo.bass_coax:g} Hz (bass to coax) and {xo.lf_hf:g} Hz "
        "(coax LF to HF). No EQ, no delays. Each way's gain matches its level at 23° to the coax LF way.</p>",
        "<table><tr><th>Way</th><th>Gain</th><th>Matching band</th></tr>"
        + "".join(_row([f"<td>{k}</td>", f'<td class="num">{gains[k]:+.2f} dB</td>',
                        f'<td class="num">{MATCH_BANDS[k][0]:g}-{MATCH_BANDS[k][1]:g} Hz</td>'])
                  for k in ("bass", "coax_lf", "coax_hf")) + "</table>",
        f"<h2>Run</h2><p>BEM bands (cached results keep their original times). Wall time {wall:.0f} s.</p>",
        "<table><tr><th>Field</th><th>Band edge</th><th>Freqs</th><th>h</th><th>Nodes</th><th>Elements</th>"
        "<th>Time</th></tr>" + "".join(run_rows) + "</table>",
        "</body></html>",
    ]
    return "\n".join(parts)
