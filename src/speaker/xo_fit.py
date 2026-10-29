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
    """Up to n peaking filters flattening smooth(f, db, EQ_SMOOTH) towards ref over band (clipped to f).

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
    x0, lb, ub = [], [], []
    rem = r.copy()
    for _ in range(n):
        f0 = fb[np.argmax(np.abs(rem))]
        gain = float(np.clip(-rem[np.argmax(np.abs(rem))], *GAIN_RANGE))
        rem = rem + _peq_db(fb, f0, gain, 2.0)
        x0 += [np.log(f0), gain, np.log(2.0)]
        lb += [llo, GAIN_RANGE[0], np.log(Q_RANGE[0])]
        ub += [lhi, GAIN_RANGE[1], np.log(Q_RANGE[1])]
    lb, ub = np.array(lb), np.array(ub)
    pad = _EDGE * (ub - lb)
    x0 = np.clip(x0, lb + pad, ub - pad)

    def resid(x: np.ndarray) -> np.ndarray:
        out = r.copy()
        for i in range(n):
            out += _peq_db(fb, float(np.exp(x[3 * i])), x[3 * i + 1], float(np.exp(x[3 * i + 2])))
        return out

    x = least_squares(resid, x0, bounds=(lb, ub)).x
    eq = [Peq(float(np.exp(x[3 * i])), float(x[3 * i + 1]), float(np.exp(x[3 * i + 2]))) for i in range(n)]
    return tuple(sorted((p for p in eq if abs(p.gain_db) >= MIN_GAIN), key=lambda p: p.f0))
