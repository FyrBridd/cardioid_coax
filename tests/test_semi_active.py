"""Semi-active mode: the DSP bass-to-coax stage ahead of the passive coax boards."""
from __future__ import annotations

import numpy as np
import pytest
from fakes import SAMPLE, omni_ways

from dataclasses import replace

from speaker.crossover import Peq, WayFilter
from speaker.curves import smooth
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
    assert sa.bass.polarity == base.bass.polarity == 1 and sa.coax.polarity == 1


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
    assert sa.coax.gain_db == 0.0
    assert sa.bass.gain_db == pytest.approx(coax_db - bass_db, abs=1e-9) and len(sa.bass.peq) <= 2


def test_semi_active_text_has_two_channels_and_biquads_at_both_rates():
    sa = SemiActive(200.0, WayFilter(1.5, 0.25e-3, -1, (Peq(100.0, 2.0, 1.5),)), WayFilter(0.0, 0.0, 1, ()))
    txt = semi_active_text(sa, 2851.0, design="round-3.toml", sha="abcd1234")
    assert "bass\tdelay_ms" in txt and "coax\tgain" in txt and "# fs = 48000" in txt and "# fs = 96000" in txt
    assert "high-pass links" in txt and "bass\tpeq\t100.0" in txt and "round-3.toml" in txt
    assert txt.count("lr4_lp 200 Hz") == 4 and txt.count("lr4_hp 200 Hz") == 4


def test_semi_active_coax_channel_gets_peqs_for_a_low_end_droop(monkeypatch):
    import speaker.semi_active as S
    monkeypatch.setattr(S, "BASS_COAX_GRID", (150.0,))
    w = omni_ways(F)
    # first-order low shelf: -8 dB well below 300 Hz, 0 dB well above
    g = 10 ** (-8 / 20)
    shelf = (g + 1j * F / 300.0) / (1 + 1j * F / 300.0)
    w.ways["coax_lf"] = w.ways["coax_lf"] * shelf[:, None]
    sa = fit_semi_active(w, _resistive(), 5000.0, _fit())
    assert len(sa.coax.peq) > 0
    bare = replace(sa, coax=replace(sa.coax, peq=()))
    sw = semi_ways(w, _resistive(), _fit())
    m = (F >= sa.bass_coax) & (F <= 1000.0)

    def spread(x):
        db = smooth(F, 20 * np.log10(np.abs(semi_total(sw, x)[:, w.obs.at23_idx])), 1 / 3)[m]
        return db.max() - db.min()

    assert spread(sa) < spread(bare)


def test_semi_active_text_lists_coax_peqs():
    sa = SemiActive(200.0, WayFilter(1.5, 0.25e-3, -1, ()), WayFilter(0.0, 0.0, 1, (Peq(150.0, 6.0, 1.2),)))
    txt = semi_active_text(sa, 2851.0, design="round-3.toml", sha="abcd1234")
    assert "coax\tpeq\t150.0\t+6.00\t1.200" in txt
    assert txt.count("coax\tpeq 150 Hz 6 dB Q1.2") == 2


def test_semi_active_coax_peqs_lie_in_the_active_fits_band_below_tl(monkeypatch):
    import speaker.semi_active as S
    from speaker.xo_fit import eq_band
    monkeypatch.setattr(S, "BASS_COAX_GRID", (150.0,))
    w = omni_ways(F)
    g = 10 ** (-8 / 20)
    w.ways["coax_lf"] = w.ways["coax_lf"] * ((g + 1j * F / 300.0) / (1 + 1j * F / 300.0))[:, None]
    fit, tl = _fit(), 5000.0
    sa = fit_semi_active(w, _resistive(), tl, fit)
    lo, hi = eq_band("coax_lf", sa.bass_coax, fit.lf_hf)
    assert sa.coax.peq and all(lo <= p.f0 <= hi and p.f0 < tl for p in sa.coax.peq)


def test_semi_active_coax_peqs_are_fitted_before_the_lr4_high_pass():
    import speaker.semi_active as S
    from speaker.system import Ways
    w = omni_ways(F)
    flat = Ways(F, w.obs, {"bass": w.ways["bass"], "coax": w.ways["coax_lf"]}, {})
    sa = S._fit_at(flat, 5000.0, 150.0, 1600.0, None)
    # a flat coax way needs no EQ: the LR4 high-pass alone must not provoke a boost near fb
    assert sa.coax.peq == ()
