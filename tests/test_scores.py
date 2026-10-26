import numpy as np, pytest
from speaker.bem.run import observation_set
from speaker.crossover import FixedCrossover
from speaker.scores import score_all
from speaker.spinorama import spinorama

O = observation_set(2.0); F = np.geomspace(50, 16000, 300)


def _run(pattern, coax_pattern=None, limit=3000.0):
    p = np.stack([pattern(O.dirs, f) for f in F]).astype(complex)
    c = p if coax_pattern is None else np.stack([coax_pattern(O.dirs, f) for f in F]).astype(complex)
    return {s.key + "/" + s.confidence: s for s in score_all(F, p, c, O, spinorama(p, O), trusted_limit=limit, xo=FixedCrossover())}


def _az(d):
    return np.degrees(np.arctan2(d[:, 0], -d[:, 1]))


def test_omni_passes_window_and_fails_cardioid():
    s = _run(lambda d, f: np.ones(len(d)))
    assert s["window_spread/trusted"].value == pytest.approx(0, abs=1e-6) and s["window_spread/trusted"].status == "pass"
    assert s["cardioid_hard/trusted"].status == "fail" and s["cardioid_hard/trusted"].value == pytest.approx(0, abs=1e-6)
    assert s["window_spread/low"].band == (3000.0, 16000.0)


def test_ideal_cardioid_passes_hard():
    card = lambda d, f: 0.5 * (1 - d[:, 1]) + 1e-3          # forward is -y
    s = _run(lambda d, f: np.ones(len(d)), coax_pattern=card)
    assert s["cardioid_hard/trusted"].status == "pass" and s["cardioid_hard/trusted"].value > 40


def test_presence_lift_flags_on_axis_hot_2_to_5k():
    def pat(d, f):
        az = np.degrees(np.arctan2(d[:, 0], -d[:, 1]))
        hot = (2000 <= f <= 5000) * (np.abs(az) < 7) * 10 ** (2 / 20)   # 0 deg 2 dB hot in 2-5 kHz
        return np.where(hot > 0, hot, 1.0)
    s = _run(pat, limit=8000.0)
    assert s["presence_lift/trusted"].status == "warn" and s["presence_lift/trusted"].value == pytest.approx(2.0, abs=0.3)


def test_rolloff_counts_rising_bands():
    def pat(d, f):
        az = np.abs(np.degrees(np.arctan2(d[:, 0], -d[:, 1])))
        return np.where((az >= 60) & (az <= 70) & (f > 1000), 2.0, 1.0)
    assert _run(pat)["rolloff_monotonic/trusted"].value >= 1


def test_flat23_and_window_spread_use_the_exact_23_direction():
    # only the exact 23 deg direction is 3 dB hot above 1 kHz; the 20 and 25 deg ring entries stay flat
    pat = lambda d, f: np.where((f > 1000) & (np.abs(_az(d) - 23) < 0.5), 10 ** (3 / 20), 1.0)
    s = _run(pat, limit=20000.0)
    assert s["window_spread/trusted"].value == pytest.approx(3.0, abs=0.05)
    assert s["window_spread/trusted"].status == "warn"
    assert s["flat23/trusted"].value > 1.0
    assert "window_spread/low" not in s  # the low part is empty above the grid


def test_parts_narrower_than_a_third_octave_are_skipped():
    s = _run(lambda d, f: np.ones(len(d)), limit=15000.0)
    assert "window_spread/low" not in s and "window_spread/trusted" in s
    assert "cardioid_hard/low" not in s


def test_rear_hemisphere_reports_worst_rear_angle():
    card = lambda d, f: 0.5 * (1 - d[:, 1]) + 1e-3
    s = _run(lambda d, f: np.ones(len(d)), coax_pattern=card)
    r = s["rear_hemisphere/trusted"]
    assert r.status == "info" and 0 < r.value < s["cardioid_hard/trusted"].value
    assert s["whole_rear/trusted"].value == pytest.approx(0, abs=1e-6)


def test_window_spread_includes_the_0_degree_curve():
    pat = lambda d, f: np.where(np.abs(_az(d)) < 0.5, 10 ** (3 / 20), 1.0)  # only 0 deg is 3 dB hot
    s = _run(pat, limit=20000.0)["window_spread/trusted"]
    assert s.value == pytest.approx(3.0, abs=0.05) and s.status == "warn"


def test_presence_lift_trusted_ignores_data_above_the_limit():
    dip = lambda d, f: np.where((f > 3300) & (np.abs(_az(d) - 23) < 0.5), 10 ** (-3 / 20), 1.0)
    flat = _run(lambda d, f: np.ones(len(d)), limit=3000.0)["presence_lift/trusted"]
    dipped = _run(dip, limit=3000.0)["presence_lift/trusted"]
    assert dipped.value == pytest.approx(flat.value, abs=1e-6) and dipped.status == flat.status == "pass"
