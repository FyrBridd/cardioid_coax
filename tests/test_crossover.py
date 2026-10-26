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
