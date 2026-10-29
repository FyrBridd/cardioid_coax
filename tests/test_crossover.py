import numpy as np
import pytest

from speaker.bem.run import observation_set
from speaker.crossover import (
    ActiveCrossover,
    FixedCrossover,
    Peq,
    WayFilter,
    apply_active,
    apply_crossover,
    biquad_response,
    biquads,
    lr4,
    peq,
    way_transfer,
)
from speaker.system import Ways


def test_lr4_sums_to_allpass():
    f = np.geomspace(20, 20000, 300)
    s = lr4(f, 1400, "low") + lr4(f, 1400, "high")
    assert np.allclose(np.abs(s), 1.0, atol=1e-9)
    assert 20 * np.log10(abs(lr4(np.array([1400.0]), 1400, "low")[0])) == pytest.approx(-6.02, abs=0.01)


def _ways(fn):
    o = observation_set(2.0)
    f = np.geomspace(50, 16000, 200)
    d = o.dirs
    w = {k: np.broadcast_to(fn(d)[None, :], (len(f), len(d))).astype(complex) for k in ("bass", "coax_lf", "coax_hf")}
    return Ways(f, o, w, {})


def test_level_match_and_omni_spin():
    from speaker.spinorama import spinorama

    w = _ways(lambda d: np.ones(len(d)))
    w.ways["coax_hf"] *= 10.0
    total, gains = apply_crossover(w, FixedCrossover())
    assert gains["coax_hf"] == pytest.approx(-20.0, abs=0.01) and gains["coax_lf"] == 0.0
    s = spinorama(total, w.obs)
    assert np.allclose(s.lw, s.on_axis, atol=1e-6) and np.allclose(s.sp, s.on_axis, atol=1e-6)
    assert np.allclose(s.spdi, 0.0, atol=1e-6)


def test_level_match_uses_exact_23_direction():
    w = _ways(lambda d: np.ones(len(d)))
    o = w.obs
    hf = np.ones(len(o.dirs))
    hf[o.at23_idx] = 10.0  # 20 deg and 25 deg stay at 1
    w.ways["coax_hf"] = np.broadcast_to(hf, w.ways["coax_hf"].shape).astype(complex)
    _, gains = apply_crossover(w, FixedCrossover())
    assert gains["coax_hf"] == pytest.approx(-20.0, abs=1e-6)


def test_each_way_gets_its_filter():
    o = observation_set(2.0)
    f = np.geomspace(20, 20000, 200)
    w = Ways(f, o, {k: np.ones((len(f), len(o.dirs)), complex) for k in ("bass", "coax_lf", "coax_hf")}, {})
    total, gains = apply_crossover(w, FixedCrossover())
    assert gains == {"coax_lf": 0.0, "coax_hf": 0.0, "bass": 0.0}
    t = total[:, 0]
    lo, hi = f < 40, f > 8000  # one way dominates, the others are 40+ dB down
    assert np.allclose(t[lo], lr4(f[lo], 170, "low"), atol=5e-3)
    assert np.allclose(t[hi], lr4(f[hi], 1400, "high"), atol=5e-3)
    fm = f[(f > 400) & (f < 700)]  # all three ways contribute, each through its own filter
    ref = lr4(fm, 170, "low") + lr4(fm, 170, "high") * lr4(fm, 1400, "low") + lr4(fm, 1400, "high")
    assert np.allclose(t[(f > 400) & (f < 700)], ref, atol=1e-12)


def test_fallback_gain_only_for_a_way_with_no_points_in_its_band():
    o = observation_set(2.0)
    f = np.geomspace(50, 1000, 30)  # nothing in the coax_hf 2-8 kHz matching band
    w = Ways(f, o, {k: np.ones((len(f), len(o.dirs)), complex) for k in ("bass", "coax_lf", "coax_hf")}, {})
    w.ways["bass"] *= 2.0
    with pytest.raises(ValueError, match="no frequencies"):
        apply_crossover(w, FixedCrossover())
    _, gains = apply_crossover(w, FixedCrossover(), fallback_db={"coax_hf": -7.0, "bass": 3.0})
    assert gains["coax_hf"] == -7.0
    assert gains["bass"] == pytest.approx(-6.02, abs=0.01)  # its band has points: matched, fallback unused


def test_peq_has_its_gain_at_f0_and_none_far_away():
    f = np.array([20.0, 2000.0, 20000.0])
    h = peq(f, Peq(2000.0, -6.0, 2.0))
    assert 20 * np.log10(abs(h[1])) == pytest.approx(-6.0, abs=1e-9)
    assert abs(20 * np.log10(abs(h[0]))) < 0.05 and abs(20 * np.log10(abs(h[2]))) < 0.05


def _flat_xo(**kw):
    w = {k: WayFilter(0.0, 0.0, 1, ()) for k in ("bass", "coax_lf", "coax_hf")}
    w.update(kw)
    return ActiveCrossover(170.0, 1400.0, w, 20000.0)


def test_active_without_eq_or_delay_matches_the_fixed_crossover():
    o = observation_set(2.0)
    f = np.geomspace(20, 20000, 200)
    w = Ways(f, o, {k: np.ones((len(f), len(o.dirs)), complex) for k in ("bass", "coax_lf", "coax_hf")}, {})
    fixed, _ = apply_crossover(w, FixedCrossover())
    assert np.allclose(apply_active(w, _flat_xo()), fixed, atol=1e-12)


def test_delay_and_polarity():
    f = np.array([1000.0])
    xo = _flat_xo(coax_hf=WayFilter(-3.0, 1e-4, -1, ()))
    h = way_transfer(f, xo, "coax_hf") / lr4(f, 1400.0, "high")
    assert abs(h[0]) == pytest.approx(10 ** (-3 / 20))
    assert np.angle(h[0]) == pytest.approx(np.angle(-np.exp(-2j * np.pi * 1000.0 * 1e-4)))


@pytest.mark.parametrize("fs, cerr, dberr", [(48000.0, 0.05, 0.2), (96000.0, 0.015, 0.05)])
def test_biquads_match_the_analog_filters_below_5k(fs, cerr, dberr):
    xo = _flat_xo(coax_hf=WayFilter(0.0, 0.0, 1, (Peq(2000.0, 6.0, 2.0),)))
    f = np.geomspace(100, 5000, 50)
    for way in ("bass", "coax_lf", "coax_hf"):
        digital = np.prod([biquad_response(c, f, fs) for _, c in biquads(xo, fs)[way]], axis=0)
        analog = way_transfer(f, xo, way)
        assert np.max(np.abs(digital - analog)) / np.max(np.abs(analog)) < cerr
        m = np.abs(analog) >= 0.1 * np.max(np.abs(analog))
        assert np.max(np.abs(20 * np.log10(np.abs(digital[m] / analog[m])))) < dberr
