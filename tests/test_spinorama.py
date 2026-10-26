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
