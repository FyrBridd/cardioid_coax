"""Spec §4 scores: directivity smoothness and cardioid rejection, graded per confidence part.

Inputs are complex PEAK pressures (Pa); every level is SPL = 20 log10(|p| / √2 / 20 µPa).
Each band-based score is computed over band ∩ [f[0], trusted_limit] (confidence "trusted") and
band ∩ [trusted_limit, f[-1]] (confidence "low"); a part narrower than 1/3 octave is skipped.
Every score's detail carries puck_fraction (the share of the part's grid points where the puck
field takes any part in the coax fields) and puck_mean_weight (the mean puck weight over them).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Literal

import numpy as np

from speaker.bem.run import Observation
from speaker.crossover import FixedCrossover
from speaker.curves import smooth
from speaker.spinorama import P_REF_PEAK, Spin

MIN_OCT = 1 / 3
FULL = (150.0, 16000.0)
LW_RING = (0, 5, 10, 15, 20, 25, 30)  # H angles of the presence-lift LW curves
ROLL_ANGLES = (40, 50, 60, 70, 80, 90)
REAR_H = (*range(120, 180, 5), *range(-120, -185, -5))  # ring angles; 180 is the -180 entry


@dataclass
class Score:
    key: str
    label: str
    value: float
    unit: str
    status: Literal["pass", "fail", "warn", "info"]
    band: tuple[float, float]
    confidence: Literal["trusted", "low"]
    formula: str
    detail: dict = field(default_factory=dict)


def _spl(p: np.ndarray) -> np.ndarray:
    return 20 * np.log10(np.abs(p) / P_REF_PEAK)


def _third_bands(f: np.ndarray, db: np.ndarray, m: np.ndarray):
    """Nominal 1/3-octave bands (centre 1000·2^(k/3), edges ×2^(±1/6)): centres and power-mean dB
    of the masked grid points in each, empty bands skipped."""
    k = np.floor(3 * np.log2(f[m] / 1000.0) + 0.5).astype(int)
    keys = np.unique(k)
    lev = [10 * np.log10(np.mean(10 ** (db[m][k == q] / 10))) for q in keys]
    return 1000.0 * 2.0 ** (keys / 3), np.array(lev)


def _runs(f: np.ndarray, ok: np.ndarray, m: np.ndarray) -> list[tuple[float, float]]:
    """Contiguous (lo, hi) frequency ranges of the masked grid points where ok holds."""
    out, start, prev = [], None, None
    for i in np.flatnonzero(m):
        if ok[i]:
            start = f[i] if start is None else start
            prev = f[i]
        elif start is not None:
            out.append((float(start), float(prev))); start = None
    if start is not None:
        out.append((float(start), float(prev)))
    return out


class _Ctx:
    def __init__(self, f, total, coax_sec, obs: Observation, spin: Spin, coax_far=None):
        self.f, self.obs, self.spin = f, obs, spin
        self.total = _spl(total)
        self.coax = _spl(coax_sec)
        self.far = None if coax_far is None else _spl(coax_far)  # (nf, 2): H 0°, H 180° at 50 m
        self.h = {int(a): int(i) for a, i in zip(obs.h_angles, obs.h_idx)}
        self._sm: dict[int, np.ndarray] = {}

    def sm(self, idx: int) -> np.ndarray:  # 1/6-octave smoothed level of one direction
        if idx not in self._sm:
            self._sm[idx] = smooth(self.f, self.total[:, idx], 1 / 6)
        return self._sm[idx]

    def lw_curves(self) -> list[np.ndarray]:
        return [self.sm(self.h[a]) for a in LW_RING] + [self.sm(i) for i in self.obs.lw_extra_idx]


def _flat23(c: _Ctx, m):
    cen, lev = _third_bands(c.f, c.spin.at23, m)
    if not len(lev):
        return None
    dev = lev - lev.mean()
    return float(np.abs(dev).max()), "info", {"rms": float(np.sqrt(np.mean(dev**2))), "bands": cen.tolist()}


def _window_spread(c: _Ctx, m):
    ref = c.sm(c.obs.at23_idx)
    dev = np.stack([c.sm(c.h[a]) - ref for a in LW_RING])
    worst = float(np.abs(dev[:, m]).max())
    ok = np.abs(dev).max(axis=0) <= 1.0
    return worst, "pass" if worst <= 1.0 else "warn", {"bands_within_1db": _runs(c.f, ok, m)}


def _presence_lift(c: _Ctx, m):
    # the c23 reference is the 2-5 kHz mean restricted to this part; the 1-2 kHz mean spans the whole grid
    pm = m & (c.f >= 2000) & (c.f <= 5000)
    if not pm.any():
        return None
    ref23 = c.sm(c.obs.at23_idx)[pm].mean()
    lo = (c.f >= 1000) & (c.f <= 2000)
    worst = -np.inf
    for cur in c.lw_curves():
        terms = cur - ref23
        if lo.any():
            terms = np.maximum(terms, cur - cur[lo].mean())
        worst = max(worst, float(terms[pm].max()))
    return worst, "warn" if worst > 0.5 else "pass", {}


def _di_smooth(fc: float):
    def run(c: _Ctx, m):
        m = m & (c.f >= fc / 2) & (c.f <= 2 * fc)
        if m.sum() < 3:
            return None
        x, res = np.log10(c.f[m]), {}
        for name, di in (("spdi", c.spin.spdi), ("erdi", c.spin.erdi)):
            fit = np.polyval(np.polyfit(x, di[m], 1), x)
            res[name] = float(np.abs(di[m] - fit).max())
        return max(res.values()), "info", res
    return run


def _rolloff(c: _Ctx, m):
    hits: dict[float, list[str]] = {}
    for sign, side in ((1, "+"), (-1, "-")):
        levels = {a: _third_bands(c.f, c.total[:, c.h[sign * a]], m) for a in ROLL_ANGLES}
        cen = levels[ROLL_ANGLES[0]][0]
        mat = np.stack([levels[a][1] for a in ROLL_ANGLES])  # (angles, bands); same bands for all
        for j in np.flatnonzero((np.diff(mat, axis=0) > 0.5).any(axis=0)):
            hits.setdefault(float(cen[j]), []).append(side)
    bands = [{"centre_hz": k, "sides": v} for k, v in sorted(hits.items())]
    return float(len(bands)), "pass" if not bands else "warn", {"bands": bands}


def _rejection(curve: str, status: Callable[[float], str], angles=(-180,)):  # 180 deg is the -180 ring entry
    def run(c: _Ctx, m):
        lv = getattr(c, curve)
        d = np.stack([lv[:, c.h[0]] - lv[:, c.h[a]] for a in angles])[:, m]  # (angles, part)
        k, j = np.unravel_index(int(d.argmin()), d.shape)
        return float(d[k, j]), status(float(d[k, j])), {"worst_angle": angles[k], "worst_hz": float(c.f[m][j])}
    return run


def _far_rejection(c: _Ctx, m):
    if c.far is None:
        return None
    d = (c.far[:, 0] - c.far[:, 1])[m]
    if np.isnan(d).all():
        return None
    j = int(np.nanargmin(d))
    return float(d[j]), "info", {"worst_hz": float(c.f[m][j])}


def _puck_weight(f: np.ndarray, blend: tuple[float, float] | None) -> np.ndarray:
    """Weight of the puck field in the coax fields: 0 below the blend, log ramp across it, 1 above."""
    if blend is None:
        return np.zeros_like(f)
    lo, hi = blend
    return np.clip(np.log(f / lo) / np.log(hi / lo), 0.0, 1.0)


# key, label, formula text, band (or a function of xo), computation
_SPECS = [
    ("flat23", "23° H flatness", FULL, _flat23,
     "The 23° H response (the exact listener direction) in 1/3-octave band levels (power mean of the grid points "
     "in each band) over 150 Hz-16 kHz, minus the mean of those band levels. Value = max |dev|; detail rms holds the RMS."),
    ("window_spread", "Window spread vs 23°", FULL, _window_spread,
     "For H 0, 5, ..., 30, the 1/6-octave-smoothed level minus the 1/6-octave-smoothed 23° curve. Value = worst |dev| "
     "over 150 Hz-16 kHz; pass if <= 1 dB, else warn. detail bands_within_1db lists the contiguous (lo, hi) ranges "
     "where every angle is within +/-1 dB."),
    ("presence_lift", "Presence lift 2-5 kHz", (2000.0, 5000.0), _presence_lift,
     "Over the LW curves (H 0-30 step 5 plus the two az 23 el +/-10 directions), each 1/6-octave smoothed: value = max "
     "over f in 2-5 kHz of max(c(f) - mean of c23 over the grid points in 2-5 kHz within this confidence part, c(f) - mean over 1-2 kHz of c). Warn if > 0.5 dB, else pass."),
    ("di_smooth_bass", "DI smoothness at bass-coax crossover", None, None,
     "For SPDI and ERDI over [fc/2, 2fc] with fc the bass-coax crossover: the max |residual| from a least-squares "
     "line in log f. Value = the worse of the two. Info."),
    ("di_smooth_lfhf", "DI smoothness at LF-HF crossover", None, None,
     "For SPDI and ERDI over [fc/2, 2fc] with fc the LF-HF crossover: the max |residual| from a least-squares "
     "line in log f. Value = the worse of the two. Info."),
    ("rolloff_monotonic", "Off-axis rolloff monotonic", FULL, _rolloff,
     "For each 1/3-octave band (power mean of the grid points) over 150 Hz-16 kHz, on the H +40...+90 side and the "
     "H -40...-90 side separately: the band counts once if on either side any step to the next larger |angle| "
     "(40, 50, ..., 90) rises by more than 0.5 dB. Value = count of such bands; pass if 0, else warn. detail bands "
     "lists the band centres with the side(s)."),
    ("cardioid_hard", "Cardioid rejection 150-300 Hz", (150.0, 300.0),
     _rejection("coax", lambda v: "pass" if v >= 10 else "fail"),
     "The coax section alone: min over 150-300 Hz of L(0 deg) - L(180 deg). Pass if >= 10 dB, else fail."),
    ("cardioid_hard_far", "Cardioid rejection 150-300 Hz at 50 m", (150.0, 300.0), _far_rejection,
     "The coax section alone at 50 m (far field): min over 150-300 Hz of L(0 deg) - L(180 deg), the plain level "
     "ratio of the two points (the 1/r asymmetry is negligible at 50 m). Info; the gate is the 2 m score."),
    ("cardioid_aim", "Cardioid rejection 300-500 Hz", (300.0, 500.0),
     _rejection("coax", lambda v: "pass" if v >= 6 else "warn"),
     "The coax section alone: min over 300-500 Hz of L(0 deg) - L(180 deg). Pass if >= 6 dB, else warn."),
    ("rear_hemisphere", "Rear hemisphere rejection", (150.0, 500.0),
     _rejection("coax", lambda v: "info", REAR_H),
     "The coax section alone: min over 150-500 Hz and over H +/-120...180 (5 deg steps) of L(0 deg) - L(theta). Info."),
    ("whole_rear", "Whole-system rear rejection", (150.0, 300.0),
     _rejection("total", lambda v: "info"),
     "The total after the crossover: min over 150-300 Hz of L(0 deg) - L(180 deg). Info."),
]


def score_all(f, total, coax_sec, obs: Observation, spin: Spin, *, trusted_limit: float,
              xo: FixedCrossover, puck_blend: tuple[float, float] | None = None,
              coax_far: np.ndarray | None = None) -> list[Score]:
    """Every §4 score, one per key per confidence part wide enough (>= 1/3 octave).

    puck_blend: (lo, hi) Hz of the cabinet -> puck crossfade (Ways.meta["puck_blend"]).
    coax_far: (nf, 2) coax section at 50 m on H 0° and 180° (Ways.meta["coax_far_h0_h180"]);
    without it cardioid_hard_far is skipped.
    """
    f = np.asarray(f, dtype=float)
    c = _Ctx(f, total, coax_sec, obs, spin, coax_far)
    w = _puck_weight(f, puck_blend)
    out: list[Score] = []
    for key, label, band, run, formula in _SPECS:
        if key.startswith("di_smooth"):
            fc = xo.bass_coax if key == "di_smooth_bass" else xo.lf_hf
            band, run = (fc / 2, 2 * fc), _di_smooth(fc)
        lo, hi = max(band[0], f[0]), min(band[1], f[-1])
        parts = (("trusted", lo, min(hi, trusted_limit)), ("low", max(lo, trusted_limit), hi))
        for conf, plo, phi in parts:
            if phi < plo * 2**MIN_OCT * (1 - 1e-9):
                continue
            m = (f >= plo) & (f <= phi)
            r = run(c, m)
            if r is not None:
                puck = {"puck_fraction": float(np.mean(w[m] > 0)), "puck_mean_weight": float(np.mean(w[m]))}
                out.append(Score(key, label, r[0], "count" if key == "rolloff_monotonic" else "dB",
                                 r[1], (float(plo), float(phi)), conf, formula, {**r[2], **puck}))
    return out
