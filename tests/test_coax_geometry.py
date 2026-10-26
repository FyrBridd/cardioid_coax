import pytest

from speaker.coax_geometry import cavity_profile, flare_radii
from speaker.drivers import load_driver
from speaker.params import load_params

P = load_params()
CX = load_driver("bms_10c262")


def test_new_driver_fields():
    assert CX.diaphragm_d == pytest.approx(0.202)  # datasheet "Effect. Diaphragm Diameter"
    assert CX.cone_depth == pytest.approx(0.065)  # PROVISIONAL estimate
    assert CX.cone_inner_d == pytest.approx(0.0636)  # voice coil diameter
    assert CX.hf_exit_d == pytest.approx(0.0254)  # PROVISIONAL
    assert load_driver("dayton_rss265ho_44").diaphragm_d == pytest.approx(0.2118, abs=5e-4)  # from Sd 352.3 cm2


def test_flare_radii_match_cad_constants():
    r_throat, r_mouth, r_out = flare_radii(P, CX)
    assert r_throat == pytest.approx(0.115) and r_out == pytest.approx(0.1449)
    assert r_mouth == pytest.approx(0.1419)


def test_cavity_profile_is_continuous_and_ordered():
    segs = cavity_profile(P, CX)
    assert [n for n, _ in segs] == ["flare", "bore", "surround", "cone", "pole", "hf_exit"]
    pts = [pt for _, s in segs for pt in s]
    for (n0, a), (n1, b) in zip(segs, segs[1:]):
        assert a[-1] == pytest.approx(b[0])
    assert pts[0] == pytest.approx((0.1419, 0.0)) and pts[-1] == pytest.approx((0.0, 0.036 + 0.065))


def test_cavity_profile_requires_cone_fields():
    with pytest.raises(ValueError, match="cone"):
        cavity_profile(P, load_driver("dayton_rss265ho_44"))
