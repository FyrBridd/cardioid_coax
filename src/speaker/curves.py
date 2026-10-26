"""Traced datasheet curves (SPL at 1 m, 2.83 V; impedance) and spectral helpers."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

CURVES_DIR = Path(__file__).resolve().parents[2] / "data" / "curves"

_EXT_OCT = 4  # octaves the end slopes are extended before the level is held
_FFT_N = 1 << 21  # minimum-phase FFT length
_SLOPE_OCT = 1 / 3  # octaves at each end used to fit the end slope


@dataclass(frozen=True)
class Curve:
    f: np.ndarray  # Hz
    value: np.ndarray  # dB SPL or ohm, per `unit`
    unit: str
    source: str


def load_curve(stem: str) -> Curve:
    """Read data/curves/<stem>.csv: `#` header lines (source, unit), then f,value rows."""
    path = CURVES_DIR / f"{stem}.csv"
    meta: dict[str, str] = {}
    rows: list[tuple[float, float]] = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            key, _, val = line[1:].partition(":")
            meta[key.strip()] = val.strip()
        elif line and line[0].isdigit():
            f, v = line.split(",")
            rows.append((float(f), float(v)))
    if "unit" not in meta or "source" not in meta:
        raise ValueError(f"{path}: missing '# unit:' or '# source:' header")
    arr = np.array(rows)
    return Curve(f=arr[:, 0], value=arr[:, 1], unit=meta["unit"], source=meta["source"])


def interp_db(curve: Curve, f: np.ndarray) -> np.ndarray:
    """Curve value at f, interpolated in log f; held at the end values outside the trace."""
    return np.interp(np.log(np.asarray(f, dtype=float)), np.log(curve.f), curve.value)


def smooth(f: np.ndarray, db: np.ndarray, fraction: float) -> np.ndarray:
    """Power average over +/- fraction/2 octave around each point of the grid f."""
    f = np.asarray(f, dtype=float)
    lf = np.log2(f)
    w = np.gradient(lf)  # each point's share of the log-f axis
    inside = np.abs(lf[:, None] - lf[None, :]) <= fraction / 2
    p = 10.0 ** (np.asarray(db, dtype=float) / 10.0)
    return 10.0 * np.log10((inside * w) @ p / (inside @ w))


def _end_slope(lf: np.ndarray, db: np.ndarray, at_top: bool) -> float:
    """dB per octave fitted over the last _SLOPE_OCT octaves at one end."""
    m = lf >= lf[-1] - _SLOPE_OCT if at_top else lf <= lf[0] + _SLOPE_OCT
    if m.sum() < 2:
        m = slice(-2, None) if at_top else slice(0, 2)
    return float(np.polyfit(lf[m], db[m], 1)[0])


def minimum_phase(f: np.ndarray, db: np.ndarray) -> np.ndarray:
    """Minimum phase (rad, e^{+jwt} convention) of the magnitude db(f).

    ln|H| is put on a uniform FFT grid whose Nyquist frequency lies _EXT_OCT octaves
    above the data. Past the data it follows the end slopes for _EXT_OCT octaves and is
    held beyond (below the data, down to DC). Its Hilbert transform is taken via the
    folded real cepstrum and the phase interpolated back to f.
    """
    f = np.asarray(f, dtype=float)
    db = np.asarray(db, dtype=float)
    lf = np.log2(f)
    lo_slope, hi_slope = _end_slope(lf, db, False), _end_slope(lf, db, True)

    n = _FFT_N
    fs = 2.0 * 2.0**_EXT_OCT * f[-1]
    fg = np.arange(n // 2 + 1) * fs / n
    lg = np.log2(np.maximum(fg, fs / n))  # DC takes the first bin's level
    g = np.interp(lg, lf, db)
    below, above = lg < lf[0], lg > lf[-1]
    g[below] = db[0] + lo_slope * (np.maximum(lg[below], lf[0] - _EXT_OCT) - lf[0])
    g[above] = db[-1] + hi_slope * (np.minimum(lg[above], lf[-1] + _EXT_OCT) - lf[-1])

    ln_mag = g * np.log(10.0) / 20.0
    c = np.fft.irfft(ln_mag, n)  # real cepstrum
    c[1 : n // 2] *= 2.0  # fold onto positive quefrency: the minimum-phase cepstrum
    c[n // 2 + 1 :] = 0.0
    phase = np.fft.rfft(c).imag
    return np.interp(f, fg, phase)
