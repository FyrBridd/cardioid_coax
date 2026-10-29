"""Per-way EQ fit: the band each way is equalised over, and a parametric-EQ fit towards a flat level."""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from speaker.crossover import Peq, peq
from speaker.curves import smooth

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
