"""One evaluation of a design: fields, vent tune, ways, crossover, spinorama and scores.

evaluate() runs a params file at an EvalSettings preset (a quick vent screen or a wider screen), with
the vents tuned or as they stand. Every score it returns is trusted: the grid stops at the trusted
limit, and the puck (when the preset has one) only sets the coax HF reference, with no blend.
A preset with fit_crossover scores the design through the active crossover fitted to it
(speaker.xo_fit.fit_active); the others use the fixed LR4 crossover.
The design helpers and the shared response pipeline here are used by speaker.directivity_report too.
"""
from __future__ import annotations

import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .bem.run import BAND_EDGES, PUCK_F_RANGE, SourceField, cabinet_field, freq_grid, puck_field, recess_ratio
from .crossover import MATCH_BANDS, ActiveCrossover, FixedCrossover, apply_active, apply_crossover
from .curves import interp_db, load_curve
from .drivers import Driver, load_driver
from .layout import Layout, compute_layout
from .params import DEFAULT_PATH, Params, load_params
from .scores import Score, score_all
from .spinorama import Spin, spinorama
from .system import COAX_CURVES, Ways, coax_section, way_responses
from .vent_tune import VentTune, tune_vents, with_vents
from .xo_fit import fit_active

CACHE_DIR = Path("out/cache/bem")
SCORE_EDGES = (150.0, 300.0, 500.0, 2000.0, 5000.0)  # Hz, score band edges added to the output grids
SCORING_VERSION = 2  # 1: the fixed crossover in every preset; 2: screens score through the fitted crossover


@dataclass
class Design:
    path: Path
    text: str
    p: Params
    coax: Driver
    bass: Driver
    layout: Layout


def grid(f_lo: float, f_hi: float, ppo: int) -> np.ndarray:
    """freq_grid(f_lo, f_hi, ppo) with the score band edges inside its range added."""
    g = freq_grid(f_lo, f_hi, ppo)
    return np.union1d(g, [e for e in SCORE_EDGES if g[0] <= e <= g[-1]])


def load_design(path: Path | None) -> Design:
    """The params at path (None: params/default.toml), the drivers they name and their layout."""
    path = DEFAULT_PATH if path is None else Path(path)
    p = load_params(path)
    coax, bass = load_driver(p.coax), load_driver(p.bass)
    return Design(path, path.read_text(), p, coax, bass, compute_layout(p, coax, bass))


def concat(parts: list[SourceField]) -> SourceField:
    """Join per-band fields along frequency; meta["bands"] lists every band, array metas are joined."""
    first = parts[0]
    meta = {k: np.concatenate([x.meta[k] for x in parts]) for k in first.meta if k != "bands"}
    meta["bands"] = [b for x in parts for b in x.meta["bands"]]
    return SourceField(np.concatenate([x.f for x in parts]), first.names, first.obs,
                       np.concatenate([x.p for x in parts]), meta)


def banded(solve: Callable[[np.ndarray], SourceField], f: np.ndarray, label: str) -> SourceField:
    """solve() once per BEM band of f, one progress line per band on stderr."""
    band = np.searchsorted(BAND_EDGES, f)
    parts = []
    for b in np.unique(band):
        t0 = time.perf_counter()
        parts.append(solve(f[band == b]))
        m = parts[-1].meta["bands"][0]
        print(f"{label}: band <= {BAND_EDGES[b]:g} Hz, {m['triangles']} elements, "
              f"{time.perf_counter() - t0:.1f} s", file=sys.stderr, flush=True)
    return concat(parts)


def _datasheet_fallback() -> dict[str, float]:
    """coax_hf gain from the datasheet on-axis curves (LF over its matching band minus HF over
    its own), used only when the grid stops below the HF matching band."""
    def mean_db(way: str, band: tuple[float, float]) -> float:
        return float(np.mean(interp_db(load_curve(f"{COAX_CURVES}_{way}_spl"), np.geomspace(*band, 64))))
    return {"coax_hf": mean_db("lf", MATCH_BANDS["coax_lf"]) - mean_db("hf", MATCH_BANDS["coax_hf"])}


def responses(d: Design, cab: SourceField, puck: SourceField | None, f: np.ndarray, tl: float,
              xo: FixedCrossover | None, lf_recess: tuple[np.ndarray, np.ndarray] | None = None,
              ) -> tuple[Ways, np.ndarray, dict[str, float], Spin, FixedCrossover | ActiveCrossover]:
    """Ways, the crossed-over total, the way gains, the spinorama of d on the grid f and the crossover used.

    xo None: the active crossover fitted to the ways (fit_active), whose gains are returned.
    """
    ways = way_responses(d.p, d.coax, d.bass, d.layout, cab, puck, f, trusted_limit=tl, lf_recess=lf_recess)
    if xo is None:
        used: FixedCrossover | ActiveCrossover = fit_active(ways, tl, fallback_db=_datasheet_fallback())
        total, gains = apply_active(ways, used), used.gains_db
    else:
        used = xo
        total, gains = apply_crossover(ways, xo, fallback_db=_datasheet_fallback())
    return ways, total, gains, spinorama(total, ways.obs, f=f), used


def ways_scores(f: np.ndarray, ways: Ways, total: np.ndarray, spin: Spin, tl: float,
                xo: FixedCrossover | ActiveCrossover) -> list[Score]:
    """score_all on the nominal coax section, with the run's puck blend and far field."""
    return score_all(f, total, coax_section(ways), ways.obs, spin, trusted_limit=tl, xo=xo,
                     puck_blend=ways.meta["puck_blend"], coax_far=ways.meta["coax_far_h0_h180"])


@dataclass(frozen=True)
class EvalSettings:
    name: str
    f_lo: float  # Hz
    f_hi: float | None  # Hz; None: the trusted limit
    ppo: int  # grid points per octave
    per_wavelength: float  # BEM elements per wavelength
    h_max: float  # m, cabinet mesh size cap
    puck_top: float | None  # Hz; a puck from 500 Hz to here, used only for the HF reference
    recess: bool  # de-embed the coax LF datasheet with the recess ratio
    fit_crossover: bool = False  # score through the fitted active crossover, not the fixed one


VENTS = EvalSettings("vents", 120.0, 600.0, 12, 6.0, 0.05, None, True, False)
SCREEN = EvalSettings("screen", 120.0, None, 12, 5.0, 0.05, 4000.0, True, True)
SETTINGS = {"vents": VENTS, "screen": SCREEN}


@dataclass
class Evaluation:
    design: Design  # with the tuned vents when tuned (path and text stay the params file's)
    settings: EvalSettings
    f: np.ndarray
    ways: Ways
    total: np.ndarray
    spin: Spin
    scores: list[Score]
    gains: dict
    xo: FixedCrossover | ActiveCrossover  # the crossover used
    tune: VentTune
    tuned: bool
    trusted_limit: float
    bands: list[dict]  # every BEM band solved or read from the cache, with its "field"
    wall_s: float


def evaluate(path: Path, settings: EvalSettings, *, trusted_limit: float, tune: bool,
             cache_dir: Path = CACHE_DIR) -> Evaluation:
    """Score the design at path on settings' grid, which stops at the trusted limit.

    The vents are always tuned (tune_vents, on the cabinet field and recess ratio). With tune, the
    design is rebuilt with the tuned resistance and net volume; without it, the params' own values
    are scored and the tune is kept for its base margin. Every score is trusted.
    """
    t0 = time.perf_counter()
    d = load_design(path)
    tl = float(trusted_limit)
    f = grid(settings.f_lo, tl if settings.f_hi is None else min(settings.f_hi, tl), settings.ppo)
    kw = dict(per_wavelength=settings.per_wavelength, cache_dir=cache_dir)
    cab = banded(lambda x: cabinet_field(d.p, d.coax, d.bass, x, h_max=settings.h_max, **kw), f, "cabinet")
    puck = None
    if settings.puck_top is not None:
        fp = grid(PUCK_F_RANGE[0], settings.puck_top, settings.ppo)
        puck = banded(lambda x: puck_field(d.p, d.coax, x, **kw), fp, "puck")
    rec = recess_ratio(d.coax, cache_dir=cache_dir) if settings.recess else None

    vt = tune_vents(d.p, d.coax, d.bass, cab, lf_recess=rec)
    if tune:
        q = with_vents(d.p, vt.resistance, vt.net_volume)
        d = Design(d.path, d.text, q, d.coax, d.bass, compute_layout(q, d.coax, d.bass))

    ways, total, gains, spin, xo = responses(d, cab, puck, f, tl, None if settings.fit_crossover else FixedCrossover(),
                                             lf_recess=rec)
    scores = ways_scores(f, ways, total, spin, tl, xo)
    bands = [{"field": name, **b} for name, fld in (("cabinet", cab), ("puck", puck)) if fld is not None
             for b in fld.meta["bands"]]
    return Evaluation(d, settings, f, ways, total, spin, scores, gains, xo, vt, tune, tl, bands,
                      time.perf_counter() - t0)
