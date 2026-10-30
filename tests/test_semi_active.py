"""Semi-active mode: the DSP bass-to-coax stage ahead of the passive coax boards."""
from __future__ import annotations

import numpy as np
import pytest
from fakes import SAMPLE, omni_ways

from speaker.crossover import Peq, WayFilter
from speaker.passive import PassiveFit, Topology, value_names
from speaker.semi_active import SemiActive, fit_semi_active, semi_active_text, semi_total, semi_ways

F = np.geomspace(50.0, 5000.0, 120)


def _resistive():
    return lambda f: {k: np.full(len(f), v, complex) for k, v in
                      (("coax_lf", 6.0), ("coax_hf", 6.0), ("bass_pair", 3.3))}


def _fit(topo=Topology()):
    return PassiveFit(200.0, 1600.0, topo, {n: SAMPLE[n] for n in value_names(topo, "passive")}, {},
                      3.5, 4.0, True, 0.0)


def test_semi_ways_are_the_coax_through_the_network_and_the_scaled_bass():
    w = omni_ways(F)
    w.meta["marker"] = 1
    sw = semi_ways(w, _resistive(), _fit(Topology(bass_wiring="series")))
    assert set(sw.ways) == {"bass", "coax"} and sw.meta is w.meta
    np.testing.assert_allclose(sw.ways["bass"], 0.5 * w.ways["bass"])
    assert np.all(np.isfinite(sw.ways["coax"]))


def test_semi_total_sums_the_two_filtered_channels():
    w = omni_ways(F)
    sw = semi_ways(w, _resistive(), _fit())
    sa = SemiActive(200.0, WayFilter(0.0, 0.0, 1, ()), WayFilter(0.0, 0.0, 1, ()))
    from speaker.crossover import lr4
    want = lr4(F, 200.0, "low")[:, None] * sw.ways["bass"] + lr4(F, 200.0, "high")[:, None] * sw.ways["coax"]
    np.testing.assert_allclose(semi_total(sw, sa), want)


def test_semi_active_recovers_a_known_bass_delay(monkeypatch):
    import speaker.semi_active as S
    monkeypatch.setattr(S, "BASS_COAX_GRID", (200.0,))
    tl = 5000.0
    base = fit_semi_active(omni_ways(F), _resistive(), tl, _fit())
    w = omni_ways(F)
    w.ways["bass"] = w.ways["bass"] * np.exp(-2j * np.pi * F * 0.8e-3)[:, None]
    sa = fit_semi_active(w, _resistive(), tl, _fit())
    shift = lambda s: s.coax.delay_s - s.bass.delay_s
    assert shift(sa) - shift(base) == pytest.approx(0.8e-3, abs=0.05e-3)
    assert min(sa.bass.delay_s, sa.coax.delay_s) == 0.0
    assert sa.bass.polarity == base.bass.polarity == 1 and sa.coax.polarity == 1 and sa.coax.peq == ()


def test_bass_gain_matches_the_coax_level_through_the_network():
    w = omni_ways(F)
    w.ways["bass"] = w.ways["bass"] * 10 ** (-6 / 20)
    sw = semi_ways(w, _resistive(), _fit())
    m = (F >= 300) & (F <= 1000)
    i = w.obs.at23_idx
    coax_db = np.mean(20 * np.log10(np.abs(sw.ways["coax"][m, i])))
    low = (F >= 80) & (F <= 150)
    bass_db = np.mean(20 * np.log10(np.abs(sw.ways["bass"][low, i])))
    sa = fit_semi_active(w, _resistive(), 5000.0, _fit())
    assert sa.coax.gain_db == 0.0 and sa.coax.peq == ()
    assert sa.bass.gain_db == pytest.approx(coax_db - bass_db, abs=1e-9) and len(sa.bass.peq) <= 2


def test_semi_active_text_has_two_channels_and_biquads_at_both_rates():
    sa = SemiActive(200.0, WayFilter(1.5, 0.25e-3, -1, (Peq(100.0, 2.0, 1.5),)), WayFilter(0.0, 0.0, 1, ()))
    txt = semi_active_text(sa, 2851.0, design="round-3.toml", sha="abcd1234")
    assert "bass\tdelay_ms" in txt and "coax\tgain" in txt and "# fs = 48000" in txt and "# fs = 96000" in txt
    assert "high-pass links" in txt and "bass\tpeq\t100.0" in txt and "round-3.toml" in txt
    assert txt.count("lr4_lp 200 Hz") == 4 and txt.count("lr4_hp 200 Hz") == 4
