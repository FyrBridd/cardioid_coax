"""Semi-active mode: a DSP in front of the passive coax boards. The coax LF and HF run through the passive
coax network (the same boards as passive mode, with the coax high-pass taken out by links K1 and K2) from one amp
channel;
the bass runs bare from another. The DSP does the bass-to-coax crossover: an LR4 split at one point, the bass
gain, delay, polarity and up to N_PEQ["bass"] PEQs, and the coax gain, delay and up to N_PEQ_COAX PEQs (the
passive network shapes the coax too, but cannot fill the cardioid chamber's low-end droop).

e^{+jωt} convention, as speaker.crossover.
"""
from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from speaker.crossover import MATCH_BANDS, WayFilter, lr4, peq, rbj_butterworth, rbj_peaking
from speaker.passive import PassiveFit, Topology, bass_scale, build, driver_loads, transfers
from speaker.schematic import link_states, links_text, terminal_text
from speaker.scores import score_all
from speaker.spinorama import spinorama
from speaker.system import Ways, coax_section
from speaker.xo_fit import (
    AT_LIMIT, BASS_COAX_GRID, BASS_DELAYS, EXPORT_FS, N_PEQ, SEARCH_KEYS, align, at_limit, eq_band, fit_peq, pick,
)

CHANNELS = ("bass", "coax")
N_PEQ_COAX = N_PEQ["coax_lf"]  # coax channel PEQs, as many as the active fit gives the coax LF


@dataclass(frozen=True)
class SemiActive:
    bass_coax: float  # Hz
    bass: WayFilter  # gain, delay, polarity, up to N_PEQ["bass"] PEQs
    coax: WayFilter  # gain, delay, polarity +1, up to N_PEQ_COAX PEQs


def semi_ways(ways: Ways, z: Callable[[np.ndarray], dict[str, np.ndarray]], fit: PassiveFit) -> Ways:
    """The two DSP channels' acoustic ways: "coax", the LF and HF through the semi_active network (polarity
    included) summed, and "bass", the woofers with bass_scale applied. ways.meta is kept."""
    t = fit.topo
    h, _ = transfers(build(fit.values, t, "semi_active"), ways.f, driver_loads(z(ways.f), t), t)
    coax = h["coax_lf"][:, None] * ways.ways["coax_lf"] + h["coax_hf"][:, None] * ways.ways["coax_hf"]
    return Ways(ways.f, ways.obs, {"bass": bass_scale(t) * ways.ways["bass"], "coax": coax}, ways.meta)


def _transfer(f: np.ndarray, sa: SemiActive, channel: str) -> np.ndarray:
    """Complex transfer of one DSP channel: gain, polarity, its LR4 edge, PEQs and delay."""
    wf = getattr(sa, channel)
    f = np.asarray(f, dtype=float)
    h = wf.polarity * 10 ** (wf.gain_db / 20) * lr4(f, sa.bass_coax, "low" if channel == "bass" else "high")
    for p in wf.peq:
        h = h * peq(f, p)
    return h * np.exp(-2j * np.pi * f * wf.delay_s)


def semi_total(sw: Ways, sa: SemiActive) -> np.ndarray:
    """(nf, nd) complex sum of the two filtered channels."""
    return sum(_transfer(sw.f, sa, k)[:, None] * sw.ways[k] for k in CHANNELS)


def _level(p: np.ndarray, f: np.ndarray, band: tuple[float, float]) -> float | None:
    """dB mean of |p| over band, None when no point of f is in it."""
    m = (f >= band[0]) & (f <= band[1])
    return float(np.mean(20 * np.log10(np.abs(p[m])))) if m.any() else None


def _fit_at(sw: Ways, tl: float, fb: float, lf_hf: float, fallback_db: float | None) -> SemiActive:
    """The DSP stage at bass-to-coax point fb from the data at or below tl only: coax at 0 dB, the bass gain
    matching it over MATCH_BANDS["bass"] (fallback_db when that band holds no point), the coax PEQs (on the coax
    way before its LR4 high-pass, over eq_band("coax_lf", fb, lf_hf) up to tl) and the bass PEQs towards the coax
    level over MATCH_BANDS["coax_lf"], then the bass delay and polarity aligned to the coax. Delays are shifted so
    the smaller is 0."""
    m = sw.f <= tl
    f, i = sw.f[m], sw.obs.at23_idx
    b, c = sw.ways["bass"][m, i], sw.ways["coax"][m, i]
    ref = _level(c, f, MATCH_BANDS["coax_lf"])
    if ref is None:
        raise ValueError(f"no frequencies in level-matching band {MATCH_BANDS['coax_lf']}")
    lb = _level(b, f, MATCH_BANDS["bass"])
    gain = fallback_db if lb is None and fallback_db is not None else ref - (lb if lb is not None else math.nan)
    if not math.isfinite(gain):
        raise ValueError(f"no frequencies in level-matching band {MATCH_BANDS['bass']}")
    # as xo_fit.fit_at: the way before its LR4 high-pass (gain 0 dB), which _transfer then applies on top
    ceq = fit_peq(f, 20 * np.log10(np.abs(c)), ref, eq_band("coax_lf", fb, lf_hf), N_PEQ_COAX)
    eq = fit_peq(f, 20 * np.log10(np.abs(b)) + gain, ref, eq_band("bass", fb, fb), N_PEQ["bass"])
    sa = SemiActive(fb, WayFilter(gain, 0.0, 1, eq), WayFilter(0.0, 0.0, 1, ceq))
    d, pol = align(f, _transfer(f, sa, "bass") * b, _transfer(f, sa, "coax") * c, ref, fb, tl, BASS_DELAYS)
    return SemiActive(fb, WayFilter(gain, max(d, 0.0), pol, eq), WayFilter(0.0, max(-d, 0.0), 1, ceq))


def fit_semi_active(ways: Ways, z: Callable[[np.ndarray], dict[str, np.ndarray]], tl: float, fit: PassiveFit,
                    fallback_db: float | None = None) -> SemiActive:
    """_fit_at on every BASS_COAX_GRID point (fitted on the data at or below tl), each scored by score_all on
    semi_total over all of the semi ways, with the coax section and meta of the original ways (so the cardioid
    scores are those of the other modes); pick chooses among those with every trusted SEARCH_KEYS score (grid
    order). fallback_db is the bass gain when its match band holds no grid points."""
    sw = semi_ways(ways, z, fit)
    cands = []
    for fb in BASS_COAX_GRID:
        sa = _fit_at(sw, tl, fb, fit.lf_hf, fallback_db)
        total = semi_total(sw, sa)
        spin = spinorama(total, ways.obs, f=ways.f)
        scores = score_all(ways.f, total, coax_section(ways), ways.obs, spin, trusted_limit=tl,
                           puck_blend=ways.meta.get("puck_blend"), coax_far=ways.meta.get("coax_far_h0_h180"))
        v = {s.key: float(s.value) for s in scores
             if s.confidence == "trusted" and s.key in SEARCH_KEYS and math.isfinite(s.value)}
        if len(v) == len(SEARCH_KEYS):
            cands.append((sa, v))
    if not cands:
        raise ValueError(f"no semi-active candidate has all {len(SEARCH_KEYS)} trusted search scores")
    return pick(cands)


def _sections(sa: SemiActive, fs: float) -> dict[str, list]:
    """Per channel, the labelled digital biquad sections (the LR4 edge as two Butterworth sections, then PEQs)."""
    out = {}
    for ch, kind in (("bass", "low"), ("coax", "high")):
        c = rbj_butterworth(kind, sa.bass_coax, fs)
        out[ch] = [(f"lr4_{'lp' if kind == 'low' else 'hp'} {sa.bass_coax:g} Hz ({i}/2)", c) for i in (1, 2)]
        out[ch] += [(f"peq {p.f0:g} Hz {p.gain_db:g} dB Q{p.q:g}", rbj_peaking(p, fs)) for p in getattr(sa, ch).peq]
    return out


def semi_active_text(sa: SemiActive, tl: float, *, design: str, sha: str, topo: Topology = Topology()) -> str:
    """The DSP stage as a setup sheet: '#' comment lines (the terminals, links and polarity reference of
    schematic.link_states(topo)["Semi-active"]), per channel one tab-separated row each of gain, delay_ms, polarity
    and every PEQ, then per EXPORT_FS sample rate the biquad sections. design and sha name the design file and the
    first characters of its sha256; topo is the passive network's (its coax boards' wiring)."""
    m = link_states(topo)["Semi-active"]
    out = [
        "# Semi-active DSP stage (SIMULATED design, not measured)",
        f"# design file {design}, sha256 {sha}; trusted limit {tl:.0f} Hz",
        f"# point: bass to coax {sa.bass_coax:g} Hz",
        f"# two amp channels on {m['pairs']}, links {links_text(m)}:",
        *(f"#   {terminal_text(t)}, feeding {t['feeds']}" for t in m["terminals"]),
        f"# {m['polarity']}.",
        "# acoustic target: LR4 at the point; the electrical filters are the LR4 sections plus each channel's PEQs",
        "# delays are relative, the smaller is 0",
        "# gains are relative (coax 0 dB): normalise them in the DSP so the largest boost does not clip",
        f"# a PEQ marked {AT_LIMIT!r} has its gain at the fitter's end of range, where the fit stopped",
        "# channel\titem\tvalue: gain dB, delay ms, polarity +1/-1, peq f0 Hz, gain dB, Q",
    ]
    for ch in CHANNELS:
        wf = getattr(sa, ch)
        out += [f"{ch}\tgain\t{wf.gain_db:+.2f}", f"{ch}\tdelay_ms\t{wf.delay_s * 1e3:.3f}",
                f"{ch}\tpolarity\t{wf.polarity:+d}"]
        out += ["\t".join([ch, "peq", f"{p.f0:.1f}", f"{p.gain_db:+.2f}", f"{p.q:.3f}",
                            *([AT_LIMIT] if at_limit(p) else [])]) for p in wf.peq]
    for fs in EXPORT_FS:
        out += ["", f"# fs = {fs:g}",
                "# channel\tlabel\tb0\tb1\tb2\ta1\ta2: H(z) = (b0 + b1 z^-1 + b2 z^-2) / (1 + a1 z^-1 + a2 z^-2), "
                "in cascade per channel; gain, delay and polarity are not in them",
                "# some DSP tools (miniDSP-style) expect -a1 and -a2: negate those two columns for them"]
        for ch, secs in _sections(sa, fs).items():
            out += ["\t".join([ch, label, *(f"{c:.10g}" for c in coeffs)]) for label, coeffs in secs]
    return "\n".join(out) + "\n"
