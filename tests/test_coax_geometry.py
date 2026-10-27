import dataclasses

import pytest

from speaker.coax_geometry import cavity_profile, flare_radii, recess_profile
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


def test_recess_profile_is_the_frame_opening_then_the_shared_cone():
    segs = recess_profile(CX)
    assert [n for n, _ in segs] == ["bore", "surround", "cone", "pole", "hf_exit"]
    for (_, a), (_, b) in zip(segs, segs[1:]):
        assert a[-1] == pytest.approx(b[0])
    prof = dict(segs)
    assert sum(prof["bore"], ()) == pytest.approx((CX.cutout_d / 2, 0.0, CX.cutout_d / 2, CX.flange_t))
    assert prof["cone"][0] == pytest.approx((CX.diaphragm_d / 2, CX.flange_t))
    assert prof["hf_exit"][-1] == pytest.approx((0.0, CX.flange_t + CX.cone_depth))
    # The cone, pole and HF exit are the cavity's, moved forward from baffle_t to flange_t.
    cav = dict(cavity_profile(P, CX))
    dy = P.cabinet.baffle_t - CX.flange_t
    for n in ("cone", "pole", "hf_exit"):
        assert [v for r, y in prof[n] for v in (r, y + dy)] == pytest.approx(sum(cav[n], ()))


def test_recess_profile_requires_cone_fields():
    with pytest.raises(ValueError, match="cone"):
        recess_profile(load_driver("dayton_rss265ho_44"))


def _with(**kw):
    c = dataclasses.replace(P.cabinet, **kw)
    return dataclasses.replace(P, cabinet=dataclasses.replace(c, baffle_t=c.outer_baffle_t + c.panel_t))


@pytest.mark.parametrize("profile", ["ellipse", "cone"])
def test_flare_runs_over_the_outer_baffle_and_the_bore_over_the_inner(profile):
    p = _with(outer_baffle_t=0.030, flare_profile=profile)
    prof = dict(cavity_profile(p, CX))
    r_throat, r_mouth, _ = flare_radii(p, CX)
    assert prof["flare"][0] == pytest.approx((r_mouth, 0.0))
    assert prof["flare"][-1] == pytest.approx((r_throat, 0.030))
    assert sum(prof["bore"], ()) == pytest.approx((r_throat, 0.030, r_throat, 0.048))
    assert prof["surround"][0] == pytest.approx((r_throat, 0.048))


def test_cone_flare_is_one_straight_segment():
    p = _with(flare_profile="cone")
    r_throat, r_mouth, _ = flare_radii(p, CX)
    flare = dict(cavity_profile(p, CX))["flare"]
    assert sum(flare, ()) == pytest.approx((r_mouth, 0.0, r_throat, P.cabinet.outer_baffle_t))


def test_ellipse_flare_is_tangent_to_the_front_and_the_bore():
    flare = dict(cavity_profile(_with(outer_baffle_t=0.030), CX))["flare"]
    assert len(flare) == 17
    (r0, y0), (r1, y1) = flare[:2]
    assert abs(y1 - y0) < 0.1 * abs(r1 - r0)  # leaves the front face nearly radially
    (r0, y0), (r1, y1) = flare[-2:]
    assert abs(r1 - r0) < 0.1 * abs(y1 - y0)  # meets the bore nearly axially


def test_unknown_flare_profile_is_rejected():
    with pytest.raises(ValueError, match="flare_profile"):
        cavity_profile(_with(flare_profile="horn"), CX)
