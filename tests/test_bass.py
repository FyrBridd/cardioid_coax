import numpy as np, pytest
from speaker.drivers import load_driver
from speaker.bass import sealed_alignment, sealed_spl, max_spl, sealed_excursion, sealed_impedance
from speaker.curves import minimum_phase

D = load_driver("dayton_rss265ho_44")

@pytest.mark.parametrize("vb,fc,qtc,f3", [(0.025, 42.7, 0.70, 43.4), (0.030, 40.6, 0.66, 43.6)])
def test_alignment_matches_spec(vb, fc, qtc, f3):
    a = sealed_alignment(D, vb)
    assert a.fc == pytest.approx(fc, abs=0.3)
    assert a.qtc == pytest.approx(qtc, abs=0.01)
    assert a.f3 == pytest.approx(f3, abs=0.5)

def test_inductor_dcr_raises_qtc():
    assert sealed_alignment(D, 0.025, series_r=0.4).qtc > sealed_alignment(D, 0.025).qtc

def test_passband_sensitivity():          # datasheet 84.5 dB @ 2.83 V series-wired
    spl = sealed_spl(D, 0.030, np.array([150.0]), 2.83)
    assert spl[0] == pytest.approx(84.7, abs=1.0)

def test_xmax_limited_spl_40hz():          # half-space 1 m, single driver, hand calc 106.0 dB
    spl = max_spl(D, 0.030, np.array([40.0]), None)
    assert spl[0] == pytest.approx(106.0, abs=0.3)

def test_two_drivers_add_6db():
    f = np.array([60.0])
    assert (sealed_spl(D, 0.03, f, 2.83, n_drivers=2) - sealed_spl(D, 0.03, f, 2.83))[0] == pytest.approx(6.02, abs=0.05)

def test_excursion_below_xmax_at_amp_limit_80hz():
    assert sealed_excursion(D, 0.03, np.array([80.0]), 25.0)[0] < D.xmax

def test_zero_frequency_rejected():
    with pytest.raises(ValueError):
        sealed_spl(D, 0.03, np.array([0.0, 50.0]), 2.83)

def test_series_r_lowers_passband_spl():
    f = np.array([150.0])
    drop = sealed_spl(D, 0.03, f, 2.83) - sealed_spl(D, 0.03, f, 2.83, series_r=0.8)
    assert drop[0] == pytest.approx(20 * np.log10((D.re + 0.8) / D.re), abs=0.05)  # ~1.0 dB

def test_series_r_changes_excursion():
    f = np.array([20.0, 40.0, 80.0, 150.0])
    x0 = sealed_excursion(D, 0.03, f, 25.0)
    x1 = sealed_excursion(D, 0.03, f, 25.0, series_r=0.8)
    assert np.all(x1 < x0)
    # Well above Fc the cone is mass-controlled and current-limited: x ~ 1 / (Re + series_r).
    assert x0[3] / x1[3] == pytest.approx((D.re + 0.8) / D.re, rel=0.01)

def test_max_spl_is_min_of_voltage_and_xmax_limits():
    f = np.geomspace(20, 300, 50)
    x_lim = max_spl(D, 0.03, f, None)
    v_lim = sealed_spl(D, 0.03, f, 40.0, series_r=0.8)
    out = max_spl(D, 0.03, f, 40.0, series_r=0.8)
    assert np.allclose(out, np.minimum(x_lim, v_lim))
    assert np.any(x_lim < v_lim) and np.any(v_lim < x_lim)  # both branches active over 20-300 Hz


def test_sealed_impedance_peaks_at_fc():
    f = np.geomspace(10, 1000, 2000)
    z = sealed_impedance(D, 0.028, f)
    assert f[np.argmax(np.abs(z))] == pytest.approx(sealed_alignment(D, 0.028).fc, rel=0.02)
    assert D.re < abs(z[0]) < 1.5 * D.re  # well below fc the motional part is small


def test_minimum_phase_reproduces_the_lumped_phase():
    # validates the method used on the datasheet traces: impedance is minimum phase
    f = np.geomspace(10, 20000, 48 * 11)
    z = sealed_impedance(D, 0.028, f)
    ph = minimum_phase(f, 20 * np.log10(np.abs(z)))
    m = (f >= 30) & (f <= 300)
    assert np.max(np.abs(np.degrees(ph[m] - np.angle(z[m])))) < 5.0
