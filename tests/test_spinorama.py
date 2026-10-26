import numpy as np
import pytest

from speaker.bem.run import observation_set
from speaker.spinorama import spinorama


def test_dipole_power_di():
    # |cos| pattern along the forward axis: DI = 10 log10(3) = 4.77 dB on axis.
    o = observation_set(2.0)
    p = np.abs(-o.dirs[:, 1])
    s = spinorama(np.broadcast_to(p, (5, len(p))).astype(complex), o)
    assert (s.on_axis - s.sp)[0] == pytest.approx(4.77, abs=0.2)


def test_er_parts_and_listener_use_their_directions():
    # Unit pressure on the rear ring entries only: rear part is 0 dB-ish, others silent.
    o = observation_set(2.0)
    p = np.full(len(o.dirs), 1e-9)
    rear = [o.h_idx[np.isclose(o.h_angles, a)][0] for a in (90, -90, -180)]
    p[rear] = 1.0
    p[o.at23_idx] = 2.0
    s = spinorama(np.broadcast_to(p, (2, len(p))).astype(complex), o)
    spl1 = 20 * np.log10(1 / np.sqrt(2) / 20e-6)
    assert s.er_parts["rear"][0] == pytest.approx(spl1, abs=1e-6)
    assert s.er_parts["front"][0] < spl1 - 100
    assert s.at23[0] == pytest.approx(20 * np.log10(2 / np.sqrt(2) / 20e-6), abs=1e-6)
    assert set(s.er_parts) == {"floor", "ceiling", "front", "side", "rear"}


def _ring_idx(o, angs, ring="h"):
    a, ix = (o.h_angles, o.h_idx) if ring == "h" else (o.v_angles, o.v_idx)
    return [int(ix[np.flatnonzero(np.isclose(a, x))[0]]) for x in angs]


def _spin(o, p):
    return spinorama(np.broadcast_to(p, (2, len(p))).astype(complex), o)


def _spl(p2):
    return 10 * np.log10(p2 / (20e-6 * np.sqrt(2)) ** 2)


def test_lw_uses_only_its_six_directions():
    o = observation_set(2.0)
    base = _spin(o, np.ones(len(o.dirs))).lw[0]
    for angs, changes in (([-20], False), ([-10, -30, 40, -180], False), ([20], True), ([0], True)):
        p = np.ones(len(o.dirs))
        p[_ring_idx(o, angs)] = 10.0
        assert (abs(_spin(o, p).lw[0] - base) > 1e-6) == changes, angs
    p = np.ones(len(o.dirs))
    p[list(o.lw_extra_idx)] = 10.0  # both az 23 el +-10 directions are in the window
    assert _spin(o, p).lw[0] == pytest.approx(_spl((4 + 2 * 100) / 6), abs=1e-6)


def test_er_is_equal_weight_power_mean_of_five_parts():
    o = observation_set(2.0)
    p = np.ones(len(o.dirs))
    p[_ring_idx(o, (-20, -30, -40), "v")] = 3.0  # floor only
    s = _spin(o, p)
    assert s.er_parts["floor"][0] == pytest.approx(_spl(9.0), abs=1e-6)
    assert s.er[0] == pytest.approx(_spl((9 + 4) / 5), abs=1e-6)


def test_sp_uses_sphere_weights():
    o = observation_set(2.0)
    p = np.ones(len(o.dirs))
    k = int(np.argmax(o.sphere_w))
    p[o.sphere_idx[k]] = 2.0
    assert _spin(o, p).sp[0] == pytest.approx(_spl(1 + 3 * o.sphere_w[k]), abs=1e-6)


def test_spin_carries_f():
    o = observation_set(2.0)
    f = np.array([100.0, 200.0])
    s = spinorama(np.ones((2, len(o.dirs)), complex), o, f=f)
    assert np.array_equal(s.f, f)
