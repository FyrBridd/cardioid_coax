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
from pathlib import Path

import numpy as np

from . import directivity_plots as plots
from .bem import run as bem_run
from .bem.mesh import cabinet_mesh
from .bem.run import (BAND_EDGES, CAB_H_COAX, CAB_H_VENT, PUCK_F_RANGE, SourceField, cabinet_field, puck_field,
                      recess_ratio)
from .bem.solver import BEM_MEMORY_BUDGET, estimate_bytes
from .cardioid import front_back, simulate
from .constants import C
from .crossover import MATCH_BANDS, FixedCrossover
from .curves import smooth
from .drivers import Driver
from .evaluate import CACHE_DIR, Design, banded, concat, grid, load_design, responses, ways_scores
from .layout import cardioid_params
from .params import Params
from .report import FB_AIM_BAND, FB_BAND
from .scores import LW_RING, REAR_H, ROLL_ANGLES, Score, score_all
from .spinorama import P_REF_PEAK
from .system import BLEND_MIN_OCT, PUCK_HOLD_F, V_DRIVE, coax_section, recess_db, vent_scaled_section

CAB_H_MAX = 0.05  # m, cabinet mesh size cap of the full run and the mesh check
FULL_GRID = (50.0, 20000.0, 12)  # f_lo, f_hi, points per octave
FAST_GRID = (50.0, 1000.0, 6)
FAST_H_MAX = 0.08  # m
COARSE = 6.0  # elements per wavelength of the run the mesh check verifies
FINE = 8.0  # elements per wavelength of the mesh check's finer run; all mesh sizes scale by COARSE / FINE
REJ_BAND = (150.0, 500.0)  # Hz, where the mesh check compares the coax-section rear rejection
MAX_DEV_DB = 0.5  # mesh check: allowed |dB| change on refinement
UNVERIFIED_LIMIT = 3000.0  # Hz, trusted limit without a mesh check
VENT_SPREAD = (0.7, 1.3)  # vent resistance factors for the felt's unknown resistance
LUMPED_DIFF_DB = 3.0  # report both cardioid models where they differ by more
_REPO = Path(__file__).resolve().parents[2]


def _params_sha(d: Design) -> str:
    return hashlib.sha256(d.path.read_bytes()).hexdigest()


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

def _fine_dofs(p: Params, coax: Driver, bass: Driver, h: float, h_coax: float, h_vent: float) -> int:
    """P1 node count of the cabinet mesh with these sizes."""
    return len(np.unique(cabinet_mesh(p, coax, bass, h=h, h_coax=h_coax, h_vent=h_vent).triangles))


def _mesh_sizes(edge: float, per_wavelength: float) -> tuple[float, float, float]:
    """(h, h_coax, h_vent) of the band below `edge`: COARSE gives the full run's mesh; a denser
    per_wavelength scales h_max, h_coax and h_vent by COARSE / per_wavelength as well."""
    s = COARSE / per_wavelength
    return min(CAB_H_MAX * s, C / (edge * per_wavelength)), CAB_H_COAX * s, CAB_H_VENT * s


_COMPARED = ("the 23° curve, listening window, early reflections and sound power, and H 0-30° (5° steps), after "
             "the crossover, over the whole grid; H ±40-90° (10° steps, both sides) after the crossover, over the "
             "whole grid; and the coax-section rear rejection L(0°) - L(θ) for θ = ±120-180° (5° steps, before "
             "the crossover) over 150-500 Hz")


def _check_curves(d: Design, fld: SourceField) -> np.ndarray:
    """(n_curves, nf) dB curves the mesh check compares; NaN where a curve is not compared."""
    ways, total, _, spin = responses(d, fld, None, fld.f, float(fld.f[-1]), FixedCrossover())
    obs, lv = fld.obs, _spl(total)
    sec = _spl(coax_section(ways))
    out_rej = (fld.f < REJ_BAND[0]) | (fld.f > REJ_BAND[1])
    rej = [sec[:, _h(obs, 0)] - sec[:, _h(obs, a)] for a in REAR_H]
    for r in rej:
        r[out_rej] = np.nan
    return np.stack([spin.at23, spin.lw, spin.er, spin.sp, *(lv[:, _h(obs, a)] for a in LW_RING),
                     *(lv[:, _h(obs, sg * a)] for a in ROLL_ANGLES for sg in (1, -1)), *rej])


def mesh_check(params_path: Path | None, out_dir: Path, *, f_max: float = 4000.0, fine: float = FINE) -> Path:
    """Compare the cabinet run at COARSE and fine elements per wavelength up to f_max.

    The coarse run is the full run's mesh; the fine run scales h_max, h_coax and h_vent by
    COARSE / fine too (_mesh_sizes). Per 1/3-octave band: the max |dB| difference over the grid
    points and over the curves in _COMPARED. The trusted limit is the upper edge of the highest
    band that, with every band below it, differs by <= MAX_DEV_DB, capped at the top solved
    frequency. A fine mesh over BEM_MEMORY_BUDGET stops the check at that band. If no band
    converges, the JSON is still written and ValueError is raised.
    """
    d = load_design(params_path)
    g = grid(*FULL_GRID)
    f = g[g <= f_max]
    band = np.searchsorted(BAND_EDGES, f)
    if band.max() >= len(BAND_EDGES):
        raise ValueError(f"f_max {f_max:g} Hz is above the last band edge, {BAND_EDGES[-1]:g} Hz")
    kw = dict(cache_dir=CACHE_DIR)
    coarse, refined, elements, capped = [], [], {}, False
    for b in np.unique(band):
        edge, fb = BAND_EDGES[b], f[band == b]
        h_f, hc_f, hv_f = _mesh_sizes(edge, fine)
        need = estimate_bytes(_fine_dofs(d.p, d.coax, d.bass, h_f, hc_f, hv_f))
        if need > BEM_MEMORY_BUDGET:
            capped = True
            print(f"mesh check: stopped at band <= {edge:g} Hz, its fine mesh needs about "
                  f"{need / 1e9:.1f} GB (budget {BEM_MEMORY_BUDGET / 1e9:.1f} GB)", file=sys.stderr)
            break
        t0 = time.perf_counter()
        c = cabinet_field(d.p, d.coax, d.bass, fb, per_wavelength=COARSE, h_max=CAB_H_MAX, **kw)
        r = cabinet_field(d.p, d.coax, d.bass, fb, per_wavelength=fine, h_max=CAB_H_MAX * COARSE / fine,
                          h_coax=hc_f, h_vent=hv_f, **kw)
        coarse.append(c)
        refined.append(r)
        elements[f"{edge:g}"] = {"coarse": c.meta["bands"][0], "fine": r.meta["bands"][0]}
        print(f"mesh check: band <= {edge:g} Hz, {c.meta['bands'][0]['triangles']} / "
              f"{r.meta['bands'][0]['triangles']} elements, {time.perf_counter() - t0:.1f} s",
              file=sys.stderr, flush=True)
    if not coarse:
        raise ValueError("mesh check: the first band's fine mesh is over the memory budget")

    fc, fr = concat(coarse), concat(refined)
    dev = np.nanmax(np.abs(_check_curves(d, fc) - _check_curves(d, fr)), axis=0)
    k = np.floor(3 * np.log2(fc.f / 1000.0) + 0.5).astype(int)  # nominal 1/3-octave band index
    band_dev, limit, ok = {}, 0.0, True
    for q in np.unique(k):
        v = float(dev[k == q].max())
        band_dev[f"{1000.0 * 2.0 ** (q / 3):.0f}"] = v
        ok = ok and v <= MAX_DEV_DB
        if ok:
            limit = 1000.0 * 2.0 ** ((q + 0.5) / 3)
    s = COARSE / fine
    out = {
        "trusted_limit_hz": min(limit, float(fc.f[-1])),  # never past the top solved frequency
        "band_max_dev_db": band_dev,
        "per_wavelength": [COARSE, float(fine)],
        "h_max_m": CAB_H_MAX,
        "mesh_sizes_m": {"coarse": {"h_max": CAB_H_MAX, "h_coax": CAB_H_COAX, "h_vent": CAB_H_VENT},
                         "fine": {"h_max": CAB_H_MAX * s, "h_coax": CAB_H_COAX * s, "h_vent": CAB_H_VENT * s}},
        "compared": _COMPARED,
        "elements": elements,
        "capped_by_memory": capped,
        "f_max_hz": float(f_max),
        "max_dev_db": MAX_DEV_DB,
        "params_sha256": _params_sha(d),
        "code_version": bem_run.code_version(),
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


def _usable_mesh_check(mc: dict | None, d: Design, *, h_max: float, per_wavelength: float) -> tuple[dict | None, str]:
    """mc if it was made for this run (params, code version, mesh densities, h_max), else (None, why not)."""
    if mc is None:
        return None, "no mesh check in this output directory (run with --mesh-check)"
    want = {"params_sha256": ("params", _params_sha(d)), "code_version": ("code version", bem_run.code_version()),
            "per_wavelength": ("per_wavelength", [float(per_wavelength), FINE]), "h_max_m": ("h_max", h_max)}
    bad = [f"{what} {mc.get(k)!r} vs {v!r}" if k in ("per_wavelength", "h_max_m") else what
           for k, (what, v) in want.items() if mc.get(k) != v]
    if bad:
        return None, "mesh_check.json is not for this run (" + "; ".join(bad) + "), so it is ignored"
    return mc, ""


def _trusted_limit(mc: dict | None, *, cab_top: float, f_top: float, puck_lo: float | None,
                   given: tuple[float, str] | None = None) -> tuple[float, str]:
    """The trusted limit and its header text, from `given` (Hz, note), the mesh check or UNVERIFIED_LIMIT.

    given replaces the mesh check's limit and text, as "<Hz> Hz (<note>)", and is capped as below.
    With the puck (puck_lo, its lowest frequency, given), the cabinet-to-puck blend runs above the
    limit, so a limit less than BLEND_MIN_OCT below the top of the cabinet field (cab_top) is capped
    there. Without it, a grid past the cabinet field caps the limit at cab_top. A limit below what
    the run needs (puck_lo, where the blend starts, or the top of a puck-less grid) raises ValueError.
    """
    if given is not None:
        tl, note = float(given[0]), f" ({given[1]})"
    else:
        tl = float(mc["trusted_limit_hz"]) if mc else UNVERIFIED_LIMIT
        note = "" if mc else " (unverified)"
    if puck_lo is not None:
        cap = cab_top / 2 ** BLEND_MIN_OCT
        if tl > cap:
            note += (f" (from {tl:.0f} Hz, capped {BLEND_MIN_OCT:.2f} octave below the top of the "
                     "cabinet field to leave room for the cabinet-to-puck blend)")
            tl = cap
    elif f_top > cab_top and tl > cab_top:
        note += f" (from {tl:.0f} Hz, capped at the top of the cabinet field)"
        tl = cab_top
    need = puck_lo if puck_lo is not None else f_top
    if tl < need:
        why = ("the cabinet-to-puck blend starts at it" if puck_lo is not None
               else "a run without the puck stops at it")
        raise ValueError(f"trusted 3D limit {tl:.0f} Hz is below {need:.0f} Hz, which this run needs ({why})")
    return tl, f"{tl:.0f} Hz{note}"


def run_directivity(params_path: Path | None, out_dir: Path, *, fast: bool = False,
                    f_max_3d: float = 4000.0, per_wavelength: float = 6.0,
                    trusted_limit: tuple[float, str] | None = None) -> Path:
    """Simulate, score and plot the design; write out_dir/report.html, scores.json, gains.json (the way gains in dB) and plots/.

    trusted_limit (Hz, note) replaces the mesh check's limit and its header text; it is still capped
    to leave room for the cabinet-to-puck blend (_trusted_limit).
    """
    d = load_design(params_path)
    out_dir = Path(out_dir)
    mc_path = out_dir / "mesh_check.json"
    h_max = FAST_H_MAX if fast else CAB_H_MAX
    if trusted_limit is None:
        mc, mc_why = _usable_mesh_check(json.loads(mc_path.read_text()) if mc_path.exists() else None, d,
                                        h_max=h_max, per_wavelength=per_wavelength)
    else:
        mc, mc_why = None, ""
    f = grid(*(FAST_GRID if fast else FULL_GRID))
    f_cab = f if fast else f[f <= f_max_3d]
    fp = None
    if not fast:  # the puck's grid points, ending on the hold frequency when the grid goes past it
        fp = f[(f >= PUCK_F_RANGE[0]) & (f <= PUCK_F_RANGE[1])]
        if f[-1] > PUCK_HOLD_F:
            fp = np.union1d(fp, [PUCK_HOLD_F])
    tl, tl_text = _trusted_limit(mc, cab_top=float(f_cab[-1]), f_top=float(f[-1]),
                                 puck_lo=None if fp is None else float(fp[0]), given=trusted_limit)
    kw = dict(per_wavelength=per_wavelength, cache_dir=CACHE_DIR)
    t0 = time.perf_counter()
    cab = banded(lambda x: cabinet_field(d.p, d.coax, d.bass, x, h_max=h_max, **kw), f_cab, "cabinet")
    puck = None if fp is None else banded(lambda x: puck_field(d.p, d.coax, x, **kw), fp, "puck")
    rec = None
    if not fast:
        t1 = time.perf_counter()
        rec = recess_ratio(d.coax, cache_dir=CACHE_DIR)
        print(f"recess ratio: {len(rec[0])} frequencies, {time.perf_counter() - t1:.1f} s", file=sys.stderr, flush=True)
    xo = FixedCrossover()
    ways, total, gains, spin = responses(d, cab, puck, f, tl, xo, lf_recess=rec)
    obs, sec = ways.obs, coax_section(ways)
    blend, far = ways.meta["puck_blend"], ways.meta["coax_far_h0_h180"]
    scores = ways_scores(f, ways, total, spin, tl, xo)

    # Cardioid: the vent resistance spread (same fields), and the lumped model for comparison.
    cp = cardioid_params(d.p, d.layout, d.coax)
    spread = {s: vent_scaled_section(ways, cab, d.coax, cp, s) for s in VENT_SPREAD}
    spread_scores = {1.0: scores}
    spread_scores.update({s: score_all(f, total, x, obs, spin, trusted_limit=tl, xo=xo, puck_blend=blend)
                          for s, x in spread.items()})
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
                        "Horizontal" if plane == "H" else "Vertical", tl, tl_text, blend=blend)
    plots.spinorama(plot_dir / "spinorama.png", spin, f, tl, tl_text, blend=blend)
    ref = smooth(f, lv[:, obs.at23_idx], 1 / 6)
    plots.window_family(plot_dir / "window_family.png", f,
                        {a: smooth(f, lv[:, _h(obs, a)], 1 / 6) - ref for a in LW_RING}, tl, tl_text, blend=blend)
    title = "Coax section rear rejection (cone + vents, before the crossover)"
    if max_diff > LUMPED_DIFF_DB:
        title += f"\nlumped vs BEM differ by > {LUMPED_DIFF_DB:g} dB (up to {max_diff:.1f} dB in 150-500 Hz)"
    lo, hi = (rej(spread[s])[: len(fr)] for s in VENT_SPREAD)
    far_rej = None if far is None else (_spl(far[:, 0]) - _spl(far[:, 1]))[: len(fr)]
    plots.rear_rejection(plot_dir / "rear_rejection.png", fr, bem, lumped, (np.minimum(lo, hi), np.maximum(lo, hi)),
                         f"x{VENT_SPREAD[0]:g}-x{VENT_SPREAD[1]:g}", title, far=far_rej)
    plots.whole_rear(plot_dir / "whole_rear.png", f[f <= min(f[-1], 2000.0)],
                     rej(total)[f <= min(f[-1], 2000.0)], xo.bass_coax)
    if rec is not None:
        fs = np.geomspace(rec[0][0], rec[0][-1], 200)
        plots.lf_reference(plot_dir / "lf_reference.png", rec[0], 20 * np.log10(np.abs(rec[1])), fs,
                           recess_db(rec, fs))

    (out_dir / "scores.json").write_text(
        json.dumps(_jsonable([dataclasses.asdict(s) for s in scores]), indent=2, allow_nan=False))
    (out_dir / "gains.json").write_text(json.dumps(_jsonable({k: float(v) for k, v in gains.items()}), indent=2,
                                                   allow_nan=False))
    page = _html(d, fast=fast, mc=mc, mc_why=mc_why, tl_text=tl_text, scores=scores, gains=gains, xo=xo, ways=ways,
                 vent=vent, cp=cp, max_diff=max_diff, fields={"cabinet": cab, "puck": puck}, lf_recess=rec,
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


def _html(d: Design, *, fast, mc, mc_why, tl_text, scores, gains, xo, ways, vent, cp, max_diff, fields, wall,
          lf_recess=None) -> str:
    e = html.escape
    sha = _params_sha(d)
    head = [("Params", f"{e(str(d.path))}<br><code>sha256 {sha}</code>"),
            ("Code commit", e(_commit())),
            ("Date", e(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))),
            ("Design", "Round 0 — current design (baseline)"),
            ("Trusted 3D limit", e(tl_text) + (f"<br>{e(mc_why)}" if mc_why else "")),
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
        pw = mc["per_wavelength"]
        how = (f"from mesh_check.json: {pw[0]:g} vs {pw[1]:g} elements per wavelength (the finer mesh also has "
               f"h_max, h_coax and h_vent scaled by {pw[0]:g}/{pw[1]:g}) agree within "
               f"{mc.get('max_dev_db', MAX_DEV_DB):g} dB up to {tl_text}. Compared: {mc.get('compared', '')}")
        if mc.get("capped_by_memory"):
            how += "; the check stopped early because the finer mesh would not fit in memory"
    else:
        how = f"{tl_text}: {mc_why}" if mc_why else tl_text
    mismatch = ways.meta["lf_join_mismatch_db"]
    rec_lo, rec_hi = (float(lf_recess[0][0]), float(lf_recess[0][-1])) if lf_recess is not None else (0.0, 0.0)
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
        "<b>Datasheet de-embedding</b>: above the LF join the coax LF volume velocity is the datasheet on-axis "
        f"response divided by that of an infinite-baffle piston of the diaphragm diameter "
        f"({d.coax.diaphragm_d * 1000:.0f} mm) at 1 m"
        + (", times the recess ratio R. The datasheet already contains the cone's own recess, which the 3D "
           "field adds again; dividing by the piston alone would double-count it. R is the BEM on-axis "
           "pressure at 1 m of the cone in its own frame recess, front-mounted on a 1 m diameter disc baffle, "
           "over that of a flat piston of the diaphragm on the same disc, both per unit volume velocity, so most "
           f"of the disc's edge diffraction cancels. R runs {rec_lo:.0f}-{rec_hi:.0f} Hz, is 1/3-octave smoothed, and holds "
           "its end values outside that range (plots/lf_reference.png). <i>Assumption:</i> the datasheet was "
           f"measured front-mounted on a large baffle, through a frame opening of the cutout "
           f"({d.coax.cutout_d * 1000:.0f} mm) with the surround {d.coax.flange_t * 1000:.0f} mm (the flange) "
           "behind the mounting face."
           if ways.meta["lf_reference"] == "recess" else
           ": flat piston (fast run). The datasheet already contains the cone's own recess, which the 3D field "
           "adds again, so this may double-count it.")
        + " The coax HF volume velocity uses "
        + ("the puck's on-axis response at 1 m as its reference, 1/3-octave smoothed as an excess over the "
           "HF exit's infinite-baffle piston; above the puck's top the reference continues as that piston, "
           "shifted to meet the puck."
           if ways.meta["hf_reference"] == "puck" else
           f"an infinite-baffle piston of the HF exit ({d.coax.hf_exit_d * 1000:.1f} mm) at 1 m as its reference "
           "(no puck in this run)."),
        "<b>Coax displacement and filler volumes</b> are estimates (spec §7.1.1/§8): the coax is a basket "
        "cylinder until the driver is measured. The chamber volume, and with it the lumped cardioid model and "
        "the vent strengths, rest on them.",
    ]
    if ways.meta["puck_blend"] is not None:
        lo, hi = ways.meta["puck_blend"]
        checked = float(mc["trusted_limit_hz"]) if mc else 0.0
        verified = ("The mesh check verified the cabinet fields across the blend" if checked >= hi else
                    f"The cabinet fields in the blend above {max(checked, lo):.0f} Hz are not mesh-checked")
        notes.append(f"<b>Cabinet-to-puck blend</b>: from {lo:.0f} Hz to {hi:.0f} Hz, above the trusted limit, the "
                     "coax LF and HF fields cross-fade from the cabinet BEM to the puck BEM (the coax on a round "
                     "cylinder, with no cabinet facets or edges), log-linear in frequency, and the vents fade "
                     f"out; above it the puck alone. The trusted range is the cabinet BEM alone. {verified}. "
                     "The range is hatched on the polar maps, spinorama and window family. Each score's "
                     "detail in scores.json gives puck_fraction (the share of its band's grid points with any "
                     "puck weight) and puck_mean_weight.")
    far = [s for s in scores if s.key == "cardioid_hard_far"]
    if far:
        near = [s for s in scores if s.key == "cardioid_hard" and s.confidence == far[0].confidence]
        notes.append("<b>Far field</b>: The gate is scored at 2 m (spec §4); far field is shown for reference. "
                     f"At 50 m the coax section's 150-300 Hz rear rejection is {far[0].value:.1f} dB"
                     + (f", against {near[0].value:.1f} dB at 2 m." if near else ".")
                     + " The rear rejection plot shows both.")
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
        fig("polar_h", "Horizontal polar map: level re on-axis at 2 m. The dashed line is the trusted 3D limit."
                + (" Hatched: the cabinet-to-puck blend." if ways.meta["puck_blend"] is not None else "")),
        fig("polar_v", "Vertical polar map (positive angles up), normalised to on-axis."),
        fig("spinorama", "Spinorama at 2 m for 2.83 V into the coax LF way after the crossover. "
                         "Listening window: H 0/10/20/30° plus az 23° at el ±10°. DIs on the right axis."),
        fig("window_family", "H 0-30° in 5° steps minus the 23° listener curve, 1/6-octave smoothed; "
                             "±1 dB band and 2-5 kHz shaded."),
        fig("rear_rejection", "Coax section alone (cone plus vents, no bass, before the crossover filter): "
                              "front minus back. Shaded: 3D result with the vent strength taken at vent resistance ×0.7 to ×1.3."),
        fig("whole_rear", "Info only: the whole speaker after the crossover, front minus back."),
        (fig("lf_reference", "Coax LF de-embedding: the recess ratio R (cone in its frame recess over a flat "
                             "piston, both on a 1 m diameter disc baffle, on axis at 1 m) and the smoothed curve "
                             "applied to the LF reference.") if lf_recess is not None else ""),
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
