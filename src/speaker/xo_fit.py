"""Active crossover fit: per-way parametric EQ towards a flat level, delay and polarity alignment,
and the search over crossover points that minimises the directivity key."""
from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
from scipy.optimize import least_squares

from speaker.crossover import (
    REF_WAY, ActiveCrossover, Peq, WayFilter, apply_active, match_gains, peq, ref_db, way_transfer,
)
from speaker.curves import smooth
from speaker.ranking import directivity_key
from speaker.scores import score_all
from speaker.spinorama import spinorama
from speaker.system import Ways, coax_section

N_PEQ = {"bass": 2, "coax_lf": 4, "coax_hf": 5}  # filters per way
GAIN_RANGE = (-15.0, 12.0)  # dB
Q_RANGE = (0.5, 6.0)
EQ_SMOOTH = 1 / 3  # octave, smoothing applied before fitting
EQ_TOL_DB = 0.5  # dB, the fit stops adding filters once the smoothed response is this close to ref
MIN_GAIN = 0.5  # dB, filters weaker than this are dropped
MIN_BAND_OCT = 1 / 3  # narrower bands get no EQ
EQ_TOP = 16000.0  # Hz
BASS_EQ_LO = 80.0  # Hz
_EDGE = 1e-6  # keeps start values strictly inside the bounds

BASS_COAX_GRID = (150.0, 170.0, 190.0, 210.0, 230.0, 250.0)  # Hz
LF_HF_GRID = (1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0)  # Hz
BASS_DELAYS = np.arange(-3.0, 3.0001, 0.05) * 1e-3  # s
HF_DELAYS = np.arange(-1.0, 1.0001, 0.01) * 1e-3  # s
ALIGN_SMOOTH = 1 / 6  # octave
SEARCH_KEYS = ("window_spread", "presence_lift", "di_smooth_bass", "di_smooth_lfhf", "flat23")


def eq_band(way: str, fb: float, fh: float) -> tuple[float, float]:
    """Band (Hz) a way is equalised over; never below the way's own high-pass point."""
    if way == "bass":
        return (BASS_EQ_LO, 2 * fb)
    if way == "coax_lf":
        return (fb, 2 * fh)
    return (fh, EQ_TOP)


def _peq_db(f: np.ndarray, f0: float, gain_db: float, q: float) -> np.ndarray:
    return 20 * np.log10(np.abs(peq(f, Peq(f0, gain_db, q))))


def fit_peq(
    f: np.ndarray, db: np.ndarray, ref: float, band: tuple[float, float], n: int
) -> tuple[Peq, ...]:
    """The fewest peaking filters (at most n) that bring smooth(f, db, EQ_SMOOTH) within EQ_TOL_DB of ref
    over band (clipped to f), or the best n if none do. Filters are added one at a time at the worst
    residual, with all refitted jointly after each.

    Returns () when the clipped band (its edges, not its points) spans less than MIN_BAND_OCT.
    """
    f = np.asarray(f, dtype=float)
    lo, hi = max(band[0], f[0]), min(band[1], f[-1])
    if hi <= lo or np.log2(hi / lo) < MIN_BAND_OCT:
        return ()
    m = (f >= lo) & (f <= hi)
    fb = f[m]
    r = smooth(f, db, EQ_SMOOTH)[m] - ref
    if m.sum() < 2 or n < 1:
        return ()

    llo, lhi = np.log(lo), np.log(hi)
    lb1 = np.array([llo, GAIN_RANGE[0], np.log(Q_RANGE[0])])
    ub1 = np.array([lhi, GAIN_RANGE[1], np.log(Q_RANGE[1])])

    def total(x: np.ndarray) -> np.ndarray:
        out = r.copy()
        for i in range(len(x) // 3):
            out += _peq_db(fb, float(np.exp(x[3 * i])), x[3 * i + 1], float(np.exp(x[3 * i + 2])))
        return out

    x = np.empty(0)
    for k in range(1, n + 1):
        rem = total(x)
        j = int(np.argmax(np.abs(rem)))
        new = np.array([np.log(fb[j]), np.clip(-rem[j], *GAIN_RANGE), np.log(2.0)])
        pad = _EDGE * (ub1 - lb1)
        new = np.clip(new, lb1 + pad, ub1 - pad)
        x0 = np.concatenate([x, new])
        x = least_squares(total, x0, bounds=(np.tile(lb1, k), np.tile(ub1, k))).x
        if np.max(np.abs(total(x))) <= EQ_TOL_DB:
            break

    eq = [Peq(float(np.exp(x[3 * i])), float(x[3 * i + 1]), float(np.exp(x[3 * i + 2])))
          for i in range(len(x) // 3)]
    return tuple(sorted((p for p in eq if abs(p.gain_db) >= MIN_GAIN), key=lambda p: p.f0))


def _smooth_cols(f: np.ndarray, db: np.ndarray, fraction: float) -> np.ndarray:
    """smooth() applied to each column of db (nf, k)."""
    lf = np.log2(f)
    inside = np.abs(lf[:, None] - lf[None, :]) <= fraction / 2
    m = inside * np.gradient(lf)
    return 10.0 * np.log10(m @ 10.0 ** (db / 10.0) / m.sum(axis=1)[:, None])


def _align(f: np.ndarray, h: np.ndarray, h_ref: np.ndarray, ref: float, fc: float, tl: float,
           delays: np.ndarray) -> tuple[float, int]:
    """Delay and polarity of h (23° complex, filtered) that minimise the RMS of the smoothed sum
    with h_ref against ref over [fc/2, 2 fc] ∩ [f[0], tl]. Ties keep the first (polarity +1 first)."""
    m = (f >= max(fc / 2, f[0])) & (f <= min(2 * fc, tl))
    if not m.any():
        raise ValueError(f"no frequencies to align over around {fc:g} Hz")
    shift = np.exp(-2j * np.pi * f[:, None] * delays[None, :])  # (nf, nd)
    best = (math.inf, 0.0, 1)
    for pol in (1, -1):
        s = pol * h[:, None] * shift + h_ref[:, None]
        db = 20 * np.log10(np.maximum(np.abs(s), 1e-30))
        rms = np.sqrt(np.mean((_smooth_cols(f, db, ALIGN_SMOOTH)[m] - ref) ** 2, axis=0))
        j = int(np.argmin(rms))
        if rms[j] < best[0]:
            best = (float(rms[j]), float(delays[j]), pol)
    return best[1], best[2]


def fit_at(w: Ways, tl: float, fb: float, fh: float,
           fallback_db: dict[str, float] | None = None) -> ActiveCrossover:
    """The active crossover at fixed points fb, fh: matched gains, per-way EQ towards the coax LF
    band level at 23°, then the bass and the coax HF delay and polarity aligned to the coax LF way.
    Delays are shifted so the smallest is 0."""
    f, i = w.f, w.obs.at23_idx
    gains = match_gains(w, fallback_db)
    ref = ref_db(w)
    ways = {}
    for k in ("bass", "coax_lf", "coax_hf"):
        top = min(f[-1], tl) if k == "bass" else f[-1]
        m = f <= top
        db = 20 * np.log10(np.abs(w.ways[k][m, i])) + gains[k]
        ways[k] = WayFilter(gains[k], 0.0, 1, fit_peq(f[m], db, ref, eq_band(k, fb, fh), N_PEQ[k]))
    xo = ActiveCrossover(fb, fh, ways, float(f[-1]))
    h = {k: way_transfer(f, xo, k) * w.ways[k][:, i] for k in ways}
    d = {REF_WAY: 0.0}
    pol = {REF_WAY: 1}
    d["bass"], pol["bass"] = _align(f, h["bass"], h[REF_WAY], ref, fb, tl, BASS_DELAYS)
    d["coax_hf"], pol["coax_hf"] = _align(f, h["coax_hf"], h[REF_WAY], ref, fh, tl, HF_DELAYS)
    d0 = min(d.values())
    return replace(xo, ways={k: replace(v, delay_s=d[k] - d0, polarity=pol[k]) for k, v in ways.items()})


def fit_active(w: Ways, tl: float, fallback_db: dict[str, float] | None = None) -> ActiveCrossover:
    """fit_at on every grid point pair; the candidate with the smallest directivity_key over its
    trusted SEARCH_KEYS scores wins. Ties keep grid order (bass_coax outer, lf_hf inner)."""
    best, best_key = None, None
    for fb in BASS_COAX_GRID:
        for fh in LF_HF_GRID:
            xo = fit_at(w, tl, fb, fh, fallback_db)
            total = apply_active(w, xo)
            spin = spinorama(total, w.obs, f=w.f)
            scores = score_all(w.f, total, coax_section(w), w.obs, spin, trusted_limit=tl, xo=xo,
                               puck_blend=w.meta.get("puck_blend"), coax_far=w.meta.get("coax_far_h0_h180"))
            v = {s.key: float(s.value) for s in scores
                 if s.confidence == "trusted" and s.key in SEARCH_KEYS and math.isfinite(s.value)}
            if len(v) < len(SEARCH_KEYS):
                continue
            key = directivity_key(v)
            if best_key is None or key < best_key:
                best, best_key = xo, key
    if best is None:
        raise ValueError("no crossover candidate has all five trusted directivity scores")
    return best
