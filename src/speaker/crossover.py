"""Crossovers: the fixed provisional one (ideal acoustic LR4 at two frequencies, one gain per way, no EQ
and no delays) and the active one (LR4 sections, per-way gain, delay, polarity and parametric EQ, fitted by
speaker.xo_fit), with its digital biquad sections and how far they miss the analog filters.

e^{+jωt} convention, so s = j f/fc.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
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


def ref_db(w: Ways) -> float:
    """Level of the reference way at the listener direction: dB mean over its matching band."""
    return _band_db(w.ways[REF_WAY][:, w.obs.at23_idx], w.f, MATCH_BANDS[REF_WAY])


def match_gains(w: Ways, fallback_db: dict[str, float] | None = None) -> dict[str, float]:
    """Per-way gains in dB that match each unfiltered way at the listener direction (dB mean over
    its band) to the coax LF way. fallback_db gives the gain of a way whose matching band holds no
    grid points; without one, such a way raises ValueError.
    """
    f, i = w.f, w.obs.at23_idx
    fb = fallback_db or {}
    ref = ref_db(w)
    gains = {
        k: fb[k] if k in fb and not _has_points(f, band) else ref - _band_db(w.ways[k][:, i], f, band)
        for k, band in MATCH_BANDS.items()
    }
    gains[REF_WAY] = 0.0
    return gains


def apply_crossover(
    w: Ways, xo: FixedCrossover, *, fallback_db: dict[str, float] | None = None
) -> tuple[np.ndarray, dict[str, float]]:
    """Sum the filtered, level-matched ways: (nf, nd) complex total and per-way gains in dB.

    Gains as match_gains gives them.
    """
    f = w.f
    gains = match_gains(w, fallback_db)
    g = {k: 10 ** (v / 20) for k, v in gains.items()}
    filt = {
        "bass": lr4(f, xo.bass_coax, "low"),
        "coax_lf": lr4(f, xo.bass_coax, "high") * lr4(f, xo.lf_hf, "low"),
        "coax_hf": lr4(f, xo.lf_hf, "high"),
    }
    total = sum(g[k] * filt[k][:, None] * w.ways[k] for k in filt)
    return total, gains


@dataclass(frozen=True)
class Peq:
    f0: float  # Hz
    gain_db: float
    q: float


def peq(f: np.ndarray, p: Peq) -> np.ndarray:
    """Analog RBJ peaking filter, complex: (s^2 + s A/Q + 1) / (s^2 + s/(A Q) + 1), s = j f/f0."""
    s = 1j * np.asarray(f, dtype=float) / p.f0
    a = 10 ** (p.gain_db / 40)
    return (s * s + s * a / p.q + 1.0) / (s * s + s / (a * p.q) + 1.0)


@dataclass(frozen=True)
class WayFilter:
    gain_db: float
    delay_s: float
    polarity: int  # +1 or -1
    peq: tuple[Peq, ...]


@dataclass(frozen=True)
class ActiveCrossover:
    bass_coax: float  # Hz
    lf_hf: float  # Hz
    ways: dict[str, WayFilter]  # "bass", "coax_lf", "coax_hf"
    fit_top: float  # Hz, top of the data the EQ was fitted on

    @property
    def gains_db(self) -> dict[str, float]:
        return {k: v.gain_db for k, v in self.ways.items()}


def _edges(xo: ActiveCrossover, way: str) -> list[tuple[str, float]]:
    """The LR4 edges of a way as (kind, fc)."""
    return {
        "bass": [("low", xo.bass_coax)],
        "coax_lf": [("high", xo.bass_coax), ("low", xo.lf_hf)],
        "coax_hf": [("high", xo.lf_hf)],
    }[way]


def way_transfer(f: np.ndarray, xo: ActiveCrossover, way: str) -> np.ndarray:
    """Complex transfer of one way: gain, polarity, LR4 edges, PEQs and delay."""
    wf = xo.ways[way]
    f = np.asarray(f, dtype=float)
    h = wf.polarity * 10 ** (wf.gain_db / 20) * np.ones(len(f), complex)
    for kind, fc in _edges(xo, way):
        h = h * lr4(f, fc, kind)
    for p in wf.peq:
        h = h * peq(f, p)
    return h * np.exp(-2j * np.pi * f * wf.delay_s)


def apply_active(w: Ways, xo: ActiveCrossover) -> np.ndarray:
    """Sum the filtered ways: (nf, nd) complex total."""
    return sum(way_transfer(w.f, xo, k)[:, None] * w.ways[k] for k in ("bass", "coax_lf", "coax_hf"))


Biquad = tuple[float, float, float, float, float]  # b0, b1, b2, a1, a2 with a0 = 1


def _rbj_butterworth(kind: Literal["low", "high"], fc: float, fs: float) -> Biquad:
    w0 = 2 * np.pi * fc / fs
    c, al = np.cos(w0), np.sin(w0) / (2 * np.sqrt(0.5))
    b0 = (1 - c) / 2 if kind == "low" else (1 + c) / 2
    b1 = 1 - c if kind == "low" else -(1 + c)
    a0 = 1 + al
    return (b0 / a0, b1 / a0, b0 / a0, -2 * c / a0, (1 - al) / a0)


def _rbj_peaking(p: Peq, fs: float) -> Biquad:
    w0 = 2 * np.pi * p.f0 / fs
    c, al, a = np.cos(w0), np.sin(w0) / (2 * p.q), 10 ** (p.gain_db / 40)
    a0 = 1 + al / a
    return ((1 + al * a) / a0, -2 * c / a0, (1 - al * a) / a0, -2 * c / a0, (1 - al / a) / a0)


def biquads(xo: ActiveCrossover, fs: float) -> dict[str, list[tuple[str, Biquad]]]:
    """Per way, the labelled digital biquad sections (LR4 edges as two Butterworth sections, then PEQs)."""
    out: dict[str, list[tuple[str, Biquad]]] = {}
    for way in ("bass", "coax_lf", "coax_hf"):
        secs: list[tuple[str, Biquad]] = []
        for kind, fc in _edges(xo, way):
            c = _rbj_butterworth(kind, fc, fs)
            secs += [(f"lr4_{'lp' if kind == 'low' else 'hp'} {fc:g} Hz ({i}/2)", c) for i in (1, 2)]
        secs += [(f"peq {p.f0:g} Hz {p.gain_db:g} dB Q{p.q:g}", _rbj_peaking(p, fs)) for p in xo.ways[way].peq]
        out[way] = secs
    return out


def biquad_response(c: tuple[float, ...], f: np.ndarray, fs: float) -> np.ndarray:
    """H(z) of one biquad (b0, b1, b2, a1, a2) at z = exp(j 2 pi f / fs)."""
    zi = np.exp(-2j * np.pi * np.asarray(f, dtype=float) / fs)
    return (c[0] + c[1] * zi + c[2] * zi**2) / (1 + c[3] * zi + c[4] * zi**2)


DEVIATION_FLOOR_DB = -20.0  # dB re the analog way's peak: quieter frequencies are not compared
DEVIATION_TOP = 20000.0  # Hz, and at most 0.45 fs


def deviation_bands(fs: float, tl: float) -> dict[str, tuple[float, float]]:
    """The bands (Hz) biquad_deviation reports, by name."""
    return {"trusted": (20.0, tl), "2-5 kHz": (2000.0, 5000.0),
            "above 5 kHz": (5000.0, min(DEVIATION_TOP, 0.45 * fs))}


def biquad_deviation(xo: ActiveCrossover, fs: float, tl: float) -> dict[str, dict[str, float]]:
    """Per way and deviation_bands band, the largest |dB| difference between the way's digital cascade at fs
    (its biquads, gain and polarity) and way_transfer, both without the delay, at the frequencies where the
    analog way is within DEVIATION_FLOOR_DB of its peak over 20 Hz to the top band's end. NaN where no
    frequency of a band qualifies."""
    bands = deviation_bands(fs, tl)
    f = np.geomspace(20.0, bands["above 5 kHz"][1], 2000)
    out: dict[str, dict[str, float]] = {}
    for way, secs in biquads(xo, fs).items():
        wf = xo.ways[way]
        digital = wf.polarity * 10 ** (wf.gain_db / 20) * np.prod([biquad_response(c, f, fs) for _, c in secs], axis=0)
        analog = way_transfer(f, replace(xo, ways={**xo.ways, way: replace(wf, delay_s=0.0)}), way)
        a = np.abs(analog)
        ok = a >= 10 ** (DEVIATION_FLOOR_DB / 20) * a.max()
        err = np.abs(20 * np.log10(np.abs(digital[ok]) / a[ok]))
        fo = f[ok]
        out[way] = {}
        for name, (lo, hi) in bands.items():
            m = (fo >= lo) & (fo <= hi)
            out[way][name] = float(err[m].max()) if m.any() else float("nan")
    return out
