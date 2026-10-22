import numpy as np, pytest
from speaker.params import load_params
from speaker.room import boundary_response, axial_modes, sbir_null

R = load_params().room

def test_front_wall_null_matches_quarter_wave():
    d = 0.6
    f = np.linspace(80, 260, 2000)
    h = boundary_response((d, 2.3, 1.0), (d + 3.0, 2.3, 1.0), f, R, surfaces=["front"])
    assert f[np.argmin(np.abs(h))] == pytest.approx(sbir_null(d), rel=0.03)
    assert sbir_null(0.6) == pytest.approx(142.9, abs=0.1)
    assert sbir_null(0.4) == pytest.approx(214.4, abs=0.1)

def test_floor_bounce_null():               # ear and coax at 0.95 m, 3 m apart -> 311 Hz
    f = np.linspace(200, 450, 2000)
    h = boundary_response((1.0, 2.3, 0.95), (4.0, 2.3, 0.95), f, R, surfaces=["floor"])
    assert f[np.argmin(np.abs(h))] == pytest.approx(311.2, rel=0.02)

def test_axial_modes():
    m = axial_modes(R)
    assert m["x"] == pytest.approx([42.9, 85.8, 128.6], abs=0.1)
    assert m["y"][:2] == pytest.approx([37.3, 74.6], abs=0.1)
    assert m["z"][0] == pytest.approx(85.8, abs=0.1)

def test_no_surfaces_is_direct_only():
    f = np.array([100.0, 1000.0])
    h = boundary_response((0.5, 1, 1), (3.5, 1, 1), f, R, surfaces=[])
    assert np.abs(h) == pytest.approx([1/3.0, 1/3.0])

def test_nonpositive_frequency_raises():
    with pytest.raises(ValueError):
        boundary_response((1, 1, 1), (2, 1, 1), np.array([0.0]), R)
