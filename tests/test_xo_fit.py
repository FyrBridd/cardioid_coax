import numpy as np
import pytest

from speaker.bem.run import observation_set
from speaker.crossover import Peq, peq
from speaker.system import Ways
from speaker.xo_fit import BASS_COAX_GRID, GAIN_RANGE, LF_HF_GRID, Q_RANGE, eq_band, fit_active, fit_at, fit_peq


def test_fit_removes_a_single_bump():
    f = np.geomspace(100, 10000, 160)
    bump = Peq(1500.0, 5.0, 2.0)
    db = 80 + 20 * np.log10(np.abs(peq(f, bump)))
    eq = fit_peq(f, db, 80.0, (500.0, 5000.0), 3)
    after = db + sum(20 * np.log10(np.abs(peq(f, p))) for p in eq)
    m = (f >= 500) & (f <= 5000)
    assert np.max(np.abs(after[m] - 80)) < 0.5
    assert any(abs(p.f0 - 1500) < 150 and p.gain_db == pytest.approx(-5.0, abs=0.7) for p in eq)


def test_boost_and_cut_stay_in_range():
    f = np.geomspace(100, 10000, 160)
    db = np.where((f > 1000) & (f < 2000), 50.0, 80.0)  # a 30 dB hole
    for p in fit_peq(f, db, 80.0, (500.0, 5000.0), 4):
        assert GAIN_RANGE[0] <= p.gain_db <= GAIN_RANGE[1] and Q_RANGE[0] <= p.q <= Q_RANGE[1]


def test_short_band_gets_no_eq():
    f = np.geomspace(100, 2851, 60)
    assert fit_peq(f, np.full(len(f), 70.0), 80.0, (2600.0, 16000.0), 5) == ()


def test_eq_bands_never_reach_below_a_high_pass():
    assert eq_band("coax_lf", 170, 1400) == (170, 2800)
    assert eq_band("coax_hf", 170, 1400) == (1400, 16000.0)
    assert eq_band("bass", 170, 1400) == (80.0, 340)


def test_a_single_bump_needs_one_filter():
    f = np.geomspace(100, 10000, 160)
    db = 80 + 20 * np.log10(np.abs(peq(f, Peq(1500.0, 5.0, 2.0))))
    eq = fit_peq(f, db, 80.0, (500.0, 5000.0), 3)
    assert len(eq) == 1
    assert abs(eq[0].f0 - 1500) < 150 and eq[0].gain_db == pytest.approx(-5.0, abs=0.7)


def omni_ways():
    o = observation_set(2.0)
    f = np.geomspace(50, 16000, 200)
    return Ways(f, o, {k: np.ones((len(f), len(o.dirs)), complex) for k in ("bass", "coax_lf", "coax_hf")}, {})


def test_alignment_recovers_a_known_delay_and_inverted_polarity():
    # the LR4 sections have phase of their own, so pin the delay added on top of the flat fit
    xo0 = fit_at(omni_ways(), 16000.0, 170.0, 1400.0)
    d0 = xo0.ways["coax_lf"].delay_s - xo0.ways["coax_hf"].delay_s
    w = omni_ways()
    w.ways["coax_hf"] = -w.ways["coax_hf"] * np.exp(-2j * np.pi * w.f * 0.3e-3)[:, None]
    xo1 = fit_at(w, 16000.0, 170.0, 1400.0)
    d1 = xo1.ways["coax_lf"].delay_s - xo1.ways["coax_hf"].delay_s
    assert xo0.ways["coax_hf"].polarity == 1 and xo1.ways["coax_hf"].polarity == -1
    assert d1 - d0 == pytest.approx(0.3e-3, abs=0.02e-3)
    assert min(x.delay_s for x in xo1.ways.values()) == 0.0


def test_flat_omni_ways_need_no_eq():
    xo = fit_at(omni_ways(), 16000.0, 170.0, 1400.0)
    assert all(x.peq == () for x in xo.ways.values())


def test_search_stays_in_the_spec_ranges_and_returns_a_grid_point():
    xo = fit_active(omni_ways(), 2851.0)
    assert xo.bass_coax in BASS_COAX_GRID and xo.lf_hf in LF_HF_GRID


def test_directivity_key_is_public():
    from speaker.ranking import directivity_key
    v = dict(window_spread=1.0, presence_lift=1.0, di_smooth_bass=0.5, di_smooth_lfhf=0.5, flat23=2.0)
    assert directivity_key(v) == (8, 1.0, 2.0)


def test_passive_warnings():
    from speaker.crossover import ActiveCrossover, WayFilter
    from speaker.xo_fit import passive_warnings
    eq = (Peq(3000.0, 8.0, 1.0),)
    xo = ActiveCrossover(170.0, 1200.0, {"bass": WayFilter(0, 0, 1, ()), "coax_lf": WayFilter(0, 0.4e-3, 1, ()),
                                         "coax_hf": WayFilter(-5, 0.0, 1, eq)}, 20000.0)
    w = " ".join(passive_warnings(xo, 2851.0))
    assert "coax_hf" in w and "8.0 dB" in w      # boost over 6 dB
    assert "750 Hz" in w                          # 1.2 kHz is within an octave of the HF resonance
    assert "0.40 ms" in w                         # LF-HF delay a passive network can't add
    assert "3000 Hz" in w                         # PEQ above the trusted limit
    flat = {k: WayFilter(0.0, 0.0, 1, ()) for k in ("bass", "coax_lf", "coax_hf")}
    assert passive_warnings(ActiveCrossover(170.0, 1600.0, flat, 20000.0), 2851.0) == []
