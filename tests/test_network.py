import numpy as np
import pytest

from speaker.network import Part, inductor_dcr, neighbours, snap, solve


def test_rc_low_pass_matches_the_analytic_transfer():
    f = np.geomspace(20, 20000, 200)
    parts = [Part("R1", "R", "in", "o", 8.0), Part("C1", "C", "o", "0", 10e-6)]
    s = solve(parts, f, {})
    expect = 1 / (1 + 2j * np.pi * f * 8.0 * 10e-6)
    assert np.allclose(s.v["o"], expect, rtol=1e-9)
    assert np.allclose(s.z_in, 8.0 + 1 / (2j * np.pi * f * 10e-6), rtol=1e-9)


def test_inductor_dcr_is_in_series_with_the_inductance():
    f = np.array([1.0, 1000.0])
    s = solve([Part("L1", "L", "in", "o", 1e-3, dcr=0.5), Part("D", "Z", "o", "0", 0.0)], f,
              {"D": np.full(2, 8.0 + 0j)})
    assert np.allclose(s.z_in, 8.5 + 2j * np.pi * f * 1e-3)
    assert np.allclose(s.across(Part("D", "Z", "o", "0", 0.0)), 8.0 / s.z_in)


def test_snap_never_returns_zero_and_lands_on_the_series():
    assert snap("C", 9.4e-6) == pytest.approx(10e-6)
    assert snap("C", 1e-9) > 0
    assert snap("R", 3.05) == pytest.approx(3.0)
    assert snap("L", 0.43e-3) == pytest.approx(0.4e-3)
    assert snap("L", 2.3e-3) == pytest.approx(2.5e-3)
    assert snap("L", 1e-7) == pytest.approx(0.1e-3)
    lo, hi = neighbours("C", 10e-6)
    assert lo == pytest.approx(8.2e-6) and hi == pytest.approx(12e-6)


def test_inductor_dcr_model():
    assert inductor_dcr(1e-3, cored=False) == pytest.approx(0.5)
    assert inductor_dcr(4e-3, cored=True) == pytest.approx(0.28)


def test_unknown_load_names_itself():
    with pytest.raises(ValueError, match="Dx"):
        solve([Part("Dx", "Z", "in", "0", 0.0)], np.array([100.0]), {})
