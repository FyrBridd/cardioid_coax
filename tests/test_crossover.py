import numpy as np
import pytest

from speaker.bem.run import observation_set
from speaker.crossover import FixedCrossover, apply_crossover, lr4
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
