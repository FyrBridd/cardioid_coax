import numpy as np
import pytest

from speaker.crossover import Peq, peq
from speaker.xo_fit import GAIN_RANGE, Q_RANGE, eq_band, fit_peq


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
