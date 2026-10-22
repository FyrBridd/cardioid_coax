import pytest
from speaker.drivers import load_driver

def test_coax_derived():
    d = load_driver("bms_10c262")
    assert d.qts == pytest.approx(0.367, abs=0.005)
    assert d.cms == pytest.approx(1.758e-4, rel=0.01)
    assert d.mms == pytest.approx(0.0341, rel=0.01)
    assert d.bl == pytest.approx(14.33, rel=0.01)

def test_dayton_series_derived():
    d = load_driver("dayton_rss265ho_44")
    assert d.qts == pytest.approx(0.450, abs=0.005)
    assert d.bl == pytest.approx(19.31, rel=0.01)   # datasheet 2*9.39=18.78, within 3%
    assert d.mms == pytest.approx(0.168, rel=0.01)  # datasheet 161 g
    assert d.re == pytest.approx(6.4)

def test_unknown_driver():
    with pytest.raises(FileNotFoundError):
        load_driver("nope")


def test_driver_names():
    assert load_driver("dayton_rss265ho_44").name == "Dayton RSS265HO-44"
    assert load_driver("bms_10c262").name == "BMS 10C262"
