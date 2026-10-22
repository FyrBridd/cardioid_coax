import numpy as np, pytest
from speaker.drivers import load_driver
from speaker.bass import sealed_alignment, sealed_spl, max_spl, sealed_excursion

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
