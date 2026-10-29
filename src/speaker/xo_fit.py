"""Active crossover fit: per-way parametric EQ towards a flat level, delay and polarity alignment,
and the search over crossover points: the best whole-system rear rejection, then directivity order."""
from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
from scipy.optimize import least_squares

from speaker.crossover import (
    REF_WAY, ActiveCrossover, Peq, WayFilter, apply_active, biquad_deviation, biquads, deviation_bands, match_gains,
    peq, ref_db, way_transfer,
)
from speaker.curves import smooth
from speaker.ranking import directivity_order
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
TRUSTED_WEIGHTS = (1.0, 4.0, 16.0, 64.0)  # weights of the trusted part in the above-limit EQ fit, tried in turn

BASS_COAX_GRID = (150.0, 170.0, 190.0, 210.0, 230.0, 250.0)  # Hz
LF_HF_GRID = (1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0)  # Hz
BASS_DELAYS = np.arange(-3.0, 3.0001, 0.05) * 1e-3  # s
HF_DELAYS = np.arange(-1.0, 1.0001, 0.01) * 1e-3  # s
ALIGN_SMOOTH = 1 / 6  # octave
SEARCH_KEYS = ("window_spread", "presence_lift", "di_smooth_bass", "di_smooth_lfhf", "flat23", "whole_rear")
WHOLE_REAR_TOL_DB = 1.0  # dB, candidates further below the best whole_rear are out of the search
EXPORT_FS = (48000.0, 96000.0)  # Hz, sample rates of the biquad export
WAYS = ("bass", "coax_lf", "coax_hf")
PASSIVE_BOOST_DB = 6.0  # dB, the most summed EQ boost a passive network can supply
HF_RESONANCE_HZ = 750.0  # Hz, the coax HF resonance (main spec §6.1)
HF_RESONANCE_OCT_HZ = 1500.0  # Hz, one octave above it: a lower lf_hf point sits too close
PASSIVE_DELAY_S = 0.1e-3  # s, the most coax HF-to-LF delay a passive network can add
BIQUAD_TOL_DB = 0.25  # dB, a 48 kHz cascade missing a way by more (trusted band or 2-5 kHz) recommends 96 kHz
LIMIT_TOL_DB = 0.01  # dB, a PEQ gain this close to a GAIN_RANGE end is at the fitter's limit
AT_LIMIT = "at the fitter's limit"


def eq_band(way: str, fb: float, fh: float) -> tuple[float, float]:
    """Band (Hz) a way is equalised over; never below the way's own high-pass point."""
    if way == "bass":
        return (BASS_EQ_LO, 2 * fb)
    if way == "coax_lf":
        return (fb, 2 * fh)
    return (fh, EQ_TOP)


def _peq_db(f: np.ndarray, f0: float, gain_db: float, q: float) -> np.ndarray:
    return 20 * np.log10(np.abs(peq(f, Peq(f0, gain_db, q))))


def _fit(fb: np.ndarray, r: np.ndarray, f0_lo: float, f0_hi: float, n: int, weight: np.ndarray,
         done) -> tuple[Peq, ...]:
    """Up to n peaking filters with f0 in [f0_lo, f0_hi] added one at a time at the worst residual r (dB on
    fb) inside that range, all refitted jointly after each to minimise weight * (r + their dB); stops once
    done(residual) holds. Filters weaker than MIN_GAIN are dropped; sorted by f0."""
    llo, lhi = np.log(f0_lo), np.log(f0_hi)
    lb1 = np.array([llo, GAIN_RANGE[0], np.log(Q_RANGE[0])])
    ub1 = np.array([lhi, GAIN_RANGE[1], np.log(Q_RANGE[1])])
    start = (fb >= f0_lo) & (fb <= f0_hi)

    def total(x: np.ndarray) -> np.ndarray:
        out = r.copy()
        for i in range(len(x) // 3):
            out += _peq_db(fb, float(np.exp(x[3 * i])), x[3 * i + 1], float(np.exp(x[3 * i + 2])))
        return out

    x = np.empty(0)
    for k in range(1, n + 1):
        rem = total(x)
        j = int(np.argmax(np.where(start, np.abs(rem), -1.0)))
        new = np.array([np.log(fb[j]), np.clip(-rem[j], *GAIN_RANGE), np.log(2.0)])
        pad = _EDGE * (ub1 - lb1)
        new = np.clip(new, lb1 + pad, ub1 - pad)
        x0 = np.concatenate([x, new])
        x = least_squares(lambda y: weight * total(y), x0, bounds=(np.tile(lb1, k), np.tile(ub1, k))).x
        if done(total(x)):
            break

    eq = [Peq(float(np.exp(x[3 * i])), float(x[3 * i + 1]), float(np.exp(x[3 * i + 2])))
          for i in range(len(x) // 3)]
    return tuple(sorted((p for p in eq if abs(p.gain_db) >= MIN_GAIN), key=lambda p: p.f0))


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
    if m.sum() < 2 or n < 1:
        return ()
    r = smooth(f, db, EQ_SMOOTH)[m] - ref
    return _fit(f[m], r, lo, hi, n, np.ones(int(m.sum())), lambda t: np.max(np.abs(t)) <= EQ_TOL_DB)


def fit_peq_above(
    f: np.ndarray, db: np.ndarray, ref: float, band: tuple[float, float], tl: float, n: int
) -> tuple[Peq, ...]:
    """At most n peaking filters with f0 in (tl, band top] that bring smooth(f, db, EQ_SMOOTH) towards ref
    over the part of band (clipped to f) above the trusted limit tl.

    The trusted part of the band enters the residual too, as the change the new filters make there (its own
    EQ is already fitted), weighted equally with the part above tl at first and by each of TRUSTED_WEIGHTS
    in turn until the filters move every point of f at or below tl by less than EQ_TOL_DB; () if none does,
    or when the part above tl spans less than MIN_BAND_OCT.
    """
    f = np.asarray(f, dtype=float)
    lo, hi = max(band[0], f[0]), min(band[1], f[-1])
    f0_lo = max(lo, tl * (1 + 1e-9))
    if n < 1 or hi <= f0_lo or np.log2(hi / f0_lo) < MIN_BAND_OCT:
        return ()
    m = (f >= lo) & (f <= hi)
    fb = f[m]
    up = fb > tl
    if up.sum() < 2:
        return ()
    r = np.where(up, smooth(f, db, EQ_SMOOTH)[m] - ref, 0.0)
    ft = f[f <= tl]
    for wt in TRUSTED_WEIGHTS:
        eq = _fit(fb, r, f0_lo, hi, n, np.where(up, 1.0, wt), lambda t: np.max(np.abs(t[up])) <= EQ_TOL_DB)
        if not eq:
            return ()
        moved = sum(_peq_db(ft, p.f0, p.gain_db, p.q) for p in eq)
        if np.max(np.abs(moved), initial=0.0) < EQ_TOL_DB:
            return eq
    return ()


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


def _clip(w: Ways, tl: float) -> Ways:
    """w on its frequencies at or below tl: every way, and every meta array (or dict of them) on the grid."""
    m = w.f <= tl
    if m.all():
        return w
    n = len(w.f)

    def cut(v):
        if isinstance(v, np.ndarray) and v.ndim and v.shape[0] == n:
            return v[m]
        return {k: cut(x) for k, x in v.items()} if isinstance(v, dict) else v
    return Ways(w.f[m], w.obs, {k: v[m] for k, v in w.ways.items()}, cut(w.meta))


def fit_at(w: Ways, tl: float, fb: float, fh: float,
           fallback_db: dict[str, float] | None = None) -> ActiveCrossover:
    """The active crossover at fixed points fb, fh, from the data at or below the trusted limit tl only:
    matched gains (each match band clipped to tl), per-way EQ towards the coax LF band level at 23°, then
    the bass and the coax HF delay and polarity aligned to the coax LF way. Delays are shifted so the
    smallest is 0. Nothing above tl enters it (spec §4); fit_above_limit adds the EQ above it."""
    w = _clip(w, tl)
    f, i = w.f, w.obs.at23_idx
    gains = match_gains(w, fallback_db)
    ref = ref_db(w)
    ways = {}
    for k in ("bass", "coax_lf", "coax_hf"):
        db = 20 * np.log10(np.abs(w.ways[k][:, i])) + gains[k]
        ways[k] = WayFilter(gains[k], 0.0, 1, fit_peq(f, db, ref, eq_band(k, fb, fh), N_PEQ[k]))
    xo = ActiveCrossover(fb, fh, ways, float(f[-1]))
    h = {k: way_transfer(f, xo, k) * w.ways[k][:, i] for k in ways}
    d = {REF_WAY: 0.0}
    pol = {REF_WAY: 1}
    d["bass"], pol["bass"] = _align(f, h["bass"], h[REF_WAY], ref, fb, tl, BASS_DELAYS)
    d["coax_hf"], pol["coax_hf"] = _align(f, h["coax_hf"], h[REF_WAY], ref, fh, tl, HF_DELAYS)
    d0 = min(d.values())
    return replace(xo, ways={k: replace(v, delay_s=d[k] - d0, polarity=pol[k]) for k, v in ways.items()})


def fit_above_limit(w: Ways, tl: float, xo: ActiveCrossover) -> ActiveCrossover:
    """xo with low-confidence PEQs added above the trusted limit tl: for each coax way whose eq_band reaches
    above tl, at most its remaining N_PEQ filters with f0 in (tl, band top] (fit_peq_above), fitted on the
    data above tl at the way's gain and trusted EQ. They move the response at or below tl by less than
    EQ_TOL_DB. Points, gains, delays, polarities and the trusted EQ are kept; fit_top becomes f[-1]."""
    f, i = w.f, w.obs.at23_idx
    if f[-1] <= tl:
        return xo
    ref = ref_db(_clip(w, tl))
    ways = {}
    for k, wf in xo.ways.items():
        add: tuple[Peq, ...] = ()
        n = N_PEQ[k] - len(wf.peq)
        band = eq_band(k, xo.bass_coax, xo.lf_hf)
        if k != "bass" and n > 0 and band[1] > tl:
            db = 20 * np.log10(np.abs(w.ways[k][:, i])) + wf.gain_db
            db = db + sum((_peq_db(f, p.f0, p.gain_db, p.q) for p in wf.peq), np.zeros(len(f)))
            add = fit_peq_above(f, db, ref, band, tl, n)
        ways[k] = replace(wf, peq=tuple(sorted(wf.peq + add, key=lambda p: p.f0)))
    return replace(xo, ways=ways, fit_top=float(f[-1]))


def pick(cands: list[tuple[object, dict[str, float]]]) -> object:
    """The first of (candidate, SEARCH_KEYS values) pairs whose whole_rear is within WHOLE_REAR_TOL_DB of
    the best, in directivity order. Ties keep the input order."""
    best = max(v["whole_rear"] for _, v in cands)
    keep = [c for c in cands if c[1]["whole_rear"] >= best - WHOLE_REAR_TOL_DB]
    return directivity_order(keep, lambda c: c[1])[0][0]


def fit_active(w: Ways, tl: float, fallback_db: dict[str, float] | None = None,
               lf_hf_min: float | None = None) -> ActiveCrossover:
    """fit_at on every grid point pair (lf_hf no lower than lf_hf_min, the driver's recommended HF
    crossover), each scored on the data at or below tl; pick chooses among those with every trusted
    SEARCH_KEYS score (ties keep grid order, bass_coax outer, lf_hf inner). The winner then gets its EQ
    above tl (fit_above_limit)."""
    lf_hf = [fh for fh in LF_HF_GRID if lf_hf_min is None or fh >= lf_hf_min]
    if not lf_hf:
        raise ValueError(f"no LF-HF grid point is at or above the driver's {lf_hf_min:g} Hz; "
                         f"the grid stops at {LF_HF_GRID[-1]:g} Hz")
    wt = _clip(w, tl)
    cands = []
    for fb in BASS_COAX_GRID:
        for fh in lf_hf:
            xo = fit_at(wt, tl, fb, fh, fallback_db)
            total = apply_active(wt, xo)
            spin = spinorama(total, wt.obs, f=wt.f)
            scores = score_all(wt.f, total, coax_section(wt), wt.obs, spin, trusted_limit=tl,
                               puck_blend=wt.meta.get("puck_blend"), coax_far=wt.meta.get("coax_far_h0_h180"))
            v = {s.key: float(s.value) for s in scores
                 if s.confidence == "trusted" and s.key in SEARCH_KEYS and math.isfinite(s.value)}
            if len(v) == len(SEARCH_KEYS):
                cands.append((xo, v))
    if not cands:
        raise ValueError(f"no crossover candidate has all {len(SEARCH_KEYS)} trusted search scores")
    return fit_above_limit(w, tl, pick(cands))


def at_limit(p: Peq) -> bool:
    """Whether p's gain sits at an end of GAIN_RANGE (within LIMIT_TOL_DB)."""
    return min(abs(p.gain_db - g) for g in GAIN_RANGE) <= LIMIT_TOL_DB


def max_boost(xo: ActiveCrossover) -> tuple[float, str, float]:
    """The largest gain plus summed PEQ (dB) of any way over 20 Hz-20 kHz, its way and frequency (Hz)."""
    g = np.geomspace(20.0, 20000.0, 2000)
    best = (-math.inf, "", 0.0)
    for way in WAYS:
        wf = xo.ways[way]
        db = wf.gain_db + sum((20 * np.log10(np.abs(peq(g, p))) for p in wf.peq), np.zeros(len(g)))
        k = int(np.argmax(db))
        if db[k] > best[0]:
            best = (float(db[k]), way, float(g[k]))
    return best


def biquad_advice(xo: ActiveCrossover, tl: float) -> tuple[dict[float, dict[str, dict[str, float]]], str]:
    """biquad_deviation at every EXPORT_FS, and a recommendation ("" for none): use 96 kHz when the 48 kHz
    cascade misses any way by more than BIQUAD_TOL_DB in the trusted band or in 2-5 kHz."""
    dev = {fs: biquad_deviation(xo, fs, tl) for fs in EXPORT_FS}
    worst = max(((v, way, band) for way, d in dev[48000.0].items() for band, v in d.items()
                 if band in ("trusted", "2-5 kHz") and math.isfinite(v)), default=(0.0, "", ""))
    if worst[0] <= BIQUAD_TOL_DB:
        return dev, ""
    return dev, (f"use 96 kHz: at 48 kHz the biquads miss the analog {worst[1]} way by up to {worst[0]:.2f} dB "
                 f"({worst[2]}), more than {BIQUAD_TOL_DB:g} dB")


def _num(v: float) -> str:
    return f"{v:.2f}" if math.isfinite(v) else "-"


def crossover_text(xo: ActiveCrossover, tl: float, *, design: str, sha: str) -> str:
    """The crossover as a DSP setup sheet: '#' comment lines, per way one tab-separated row each of
    gain, delay_ms, polarity and every PEQ, then per EXPORT_FS sample rate the biquad sections and how far
    they miss the analog ways. design and sha name the design file and the first characters of its sha256."""
    boost, bway, bf = max_boost(xo)
    dev, advice = biquad_advice(xo, tl)
    out = [
        "# Active crossover (SIMULATED design, not measured)",
        f"# design file {design}, sha256 {sha}; trusted limit {tl:.0f} Hz",
        f"# points: bass to coax {xo.bass_coax:g} Hz, coax LF to HF {xo.lf_hf:g} Hz",
        "# acoustic target: LR4 at each point; the electrical filters are the LR4 sections plus per-way PEQ",
        "# delays are relative, the smallest is 0",
        "# gains are relative (coax LF 0 dB): normalise them in the DSP so the largest boost does not clip; "
        f"the largest gain plus PEQ boost of any way is {boost:+.2f} dB ({bway} at {bf:.0f} Hz)",
        f"# PEQ above {tl:.0f} Hz rests on low-confidence data (marked low); a PEQ marked {AT_LIMIT!r} has its "
        f"gain at {GAIN_RANGE[0]:g} or {GAIN_RANGE[1]:+g} dB, where the fit stopped",
        *([f"# recommendation: {advice}"] if advice else []),
        "# way\titem\tvalue: gain dB, delay ms, polarity +1/-1, peq f0 Hz, gain dB, Q",
    ]
    for way in WAYS:
        wf = xo.ways[way]
        out += [f"{way}\tgain\t{wf.gain_db:+.2f}", f"{way}\tdelay_ms\t{wf.delay_s * 1e3:.3f}",
                f"{way}\tpolarity\t{wf.polarity:+d}"]
        out += ["\t".join([way, "peq", f"{p.f0:.1f}", f"{p.gain_db:+.2f}", f"{p.q:.3f}",
                            *(["low"] if p.f0 > tl else []), *([AT_LIMIT] if at_limit(p) else [])])
                for p in wf.peq]
    for fs in EXPORT_FS:
        bands = deviation_bands(fs, tl)
        out += ["", f"# fs = {fs:g}",
                "# way\tlabel\tb0\tb1\tb2\ta1\ta2: H(z) = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2), "
                "in cascade per way; gain, delay and polarity are not in them",
                "# some DSP tools (miniDSP-style) expect -a1 and -a2: negate those two columns for them",
                f"# max |dB| these sections (with gain and polarity, no delay) miss the analog way by, where the "
                f"analog way is within 20 dB of its peak ('-': none there):",
                "# way\t" + "\t".join(f"{n} ({lo:.0f}-{hi:.0f} Hz)" for n, (lo, hi) in bands.items())]
        out += ["# " + "\t".join([way, *(_num(dev[fs][way][n]) for n in bands)]) for way in WAYS]
        for way, secs in biquads(xo, fs).items():
            out += ["\t".join([way, label, *(f"{c:.10g}" for c in coeffs)]) for label, coeffs in secs]
    return "\n".join(out) + "\n"


def peak_boost(xo: ActiveCrossover, way: str) -> tuple[float, float]:
    """The largest summed PEQ boost (dB, 0 for none) of a way in its eq_band, and the frequency (Hz) of it."""
    lo, hi = eq_band(way, xo.bass_coax, xo.lf_hf)
    if not xo.ways[way].peq or hi <= lo:
        return 0.0, lo
    g = np.geomspace(lo, hi, 400)
    db = sum(20 * np.log10(np.abs(peq(g, p))) for p in xo.ways[way].peq)
    k = int(np.argmax(db))
    return max(float(db[k]), 0.0), float(g[k])


def passive_warnings(xo: ActiveCrossover, tl: float) -> list[str]:
    """Where a passive network could not follow this active crossover, one sentence each: a way whose summed
    EQ boosts more than PASSIVE_BOOST_DB in its eq_band, lf_hf within an octave of the HF resonance, a coax
    HF-LF delay over PASSIVE_DELAY_S, and any PEQ above the trusted limit tl."""
    out = []
    for way in WAYS:
        boost, at = peak_boost(xo, way)
        if boost > PASSIVE_BOOST_DB:
            out.append(f"{way} needs {boost:.1f} dB of summed EQ boost at {at:.0f} Hz, more than a passive "
                       f"network can supply ({PASSIVE_BOOST_DB:g} dB).")
    if xo.lf_hf < HF_RESONANCE_OCT_HZ:
        out.append(f"The coax LF to HF point, {xo.lf_hf:g} Hz, is within an octave of the {HF_RESONANCE_HZ:g} Hz "
                   f"HF resonance (the limit is {HF_RESONANCE_OCT_HZ:g} Hz).")
    dt = abs(xo.ways["coax_hf"].delay_s - xo.ways["coax_lf"].delay_s)
    if dt > PASSIVE_DELAY_S:
        out.append(f"The coax HF and LF delays differ by {dt * 1e3:.2f} ms, more than a passive "
                   f"network can add ({PASSIVE_DELAY_S * 1e3:.1f} ms).")
    for way in WAYS:
        for p in xo.ways[way].peq:
            if p.f0 > tl:
                out.append(f"{way} has a {p.gain_db:+.1f} dB filter at {p.f0:.0f} Hz, above the trusted limit "
                           f"{tl:.0f} Hz, so it rests on low-confidence data.")
    return out
