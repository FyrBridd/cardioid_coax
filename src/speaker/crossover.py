"""Fixed provisional crossover: ideal acoustic LR4 at two frequencies, one gain per way.

e^{+jωt} convention, so s = j f/fc. No EQ and no delays.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from speaker.system import Ways

# Level-matching bands (Hz): the reference way and the band each way is matched over.
REF_WAY = "coax_lf"
MATCH_BANDS = {"coax_lf": (300.0, 1000.0), "coax_hf": (2000.0, 8000.0), "bass": (80.0, 150.0)}


def lr4(f: np.ndarray, fc: float, kind: Literal["low", "high"]) -> np.ndarray:
    """Linkwitz-Riley 4th order (squared 2nd-order Butterworth), complex."""
    s = 1j * np.asarray(f, dtype=float) / fc
    d = s * s + np.sqrt(2.0) * s + 1.0
    h = 1.0 / d if kind == "low" else (s * s) / d
    return h * h


@dataclass(frozen=True)
class FixedCrossover:
    bass_coax: float = 170.0  # Hz
    lf_hf: float = 1400.0  # Hz


def _band_db(p: np.ndarray, f: np.ndarray, band: tuple[float, float]) -> float:
    m = (f >= band[0]) & (f <= band[1])
    if not m.any():
        raise ValueError(f"no frequencies in level-matching band {band}")
    return float(np.mean(20 * np.log10(np.abs(p[m]))))


def _has_points(f: np.ndarray, band: tuple[float, float]) -> bool:
    return bool(((f >= band[0]) & (f <= band[1])).any())


def apply_crossover(
    w: Ways, xo: FixedCrossover, *, fallback_db: dict[str, float] | None = None
) -> tuple[np.ndarray, dict[str, float]]:
    """Sum the filtered, level-matched ways: (nf, nd) complex total and per-way gains in dB.

    Gains compare the unfiltered ways at the listener direction (dB mean over each band)
    against the coax LF way. fallback_db gives the gain of a way whose matching band holds
    no grid points; without one, such a way raises ValueError.
    """
    f, i = w.f, w.obs.at23_idx
    fb = fallback_db or {}
    ref = _band_db(w.ways[REF_WAY][:, i], f, MATCH_BANDS[REF_WAY])
    gains = {
        k: fb[k] if k in fb and not _has_points(f, band) else ref - _band_db(w.ways[k][:, i], f, band)
        for k, band in MATCH_BANDS.items()
    }
    gains[REF_WAY] = 0.0
    g = {k: 10 ** (v / 20) for k, v in gains.items()}
    filt = {
        "bass": lr4(f, xo.bass_coax, "low"),
        "coax_lf": lr4(f, xo.bass_coax, "high") * lr4(f, xo.lf_hf, "low"),
        "coax_hf": lr4(f, xo.lf_hf, "high"),
    }
    total = sum(g[k] * filt[k][:, None] * w.ways[k] for k in filt)
    return total, gains
