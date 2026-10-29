import math, numpy as np, pytest
from speaker.bass import sealed_impedance
from speaker.drivers import load_driver
from speaker.cardioid import CardioidParams, electrical_impedance, simulate, front_back, tune_resistance, polar

D = load_driver("bms_10c262")
F = np.geomspace(100, 1000, 200)
VENTS = ((-0.30, 0.16, 0.0), (-0.30, -0.16, 0.0))

def P(r):
    return CardioidParams(0.010, 0.012, 0.018, r, VENTS)

def test_sealed_is_omni_and_matches_box():
    res = simulate(D, P(math.inf), F, 2.83)
    assert np.allclose(front_back(res, P(math.inf)), 0.0, atol=1e-6)
    assert np.all(res.u_vent == 0)

def test_open_back_is_symmetric_front_back():
    res = simulate(D, P(0.0), F, 2.83)
    assert np.all(np.isfinite(res.u_vent))
    assert np.allclose(front_back(res, P(0.0)), 0.0, atol=0.5)

def test_tuned_vents_reach_10db_150_to_300():
    r = tune_resistance(D, P(100.0))
    res = simulate(D, P(r), F, 2.83)
    fb = front_back(res, P(r))
    band = (F >= 150) & (F <= 300)
    assert fb[band].min() >= 10.0
    assert 30 < r < 400

def test_excursion_tiny_at_1w():
    res = simulate(D, P(100.0), np.array([150.0]), 2.83)
    assert res.excursion[0] < 0.5e-3

def test_polar_shape():
    res = simulate(D, P(100.0), F, 2.83)
    assert polar(res, P(100.0), np.arange(0, 181, 10)).shape == (200, 19)

def test_bad_frequency_raises():
    with pytest.raises(ValueError):
        simulate(D, P(100.0), np.array([0.0]), 2.83)


def test_closed_chamber_matches_a_sealed_box():
    f = np.geomspace(20, 2000, 100)
    p = CardioidParams(0.016, 0.0168, 0.018, math.inf, ((0.0, 0.1, 0.0),), 12)
    assert np.allclose(electrical_impedance(D, p, f), sealed_impedance(D, 0.016, f), rtol=1e-9)
