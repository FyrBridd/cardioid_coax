import dataclasses
import itertools

import pytest
from speaker.params import load_params
from speaker.drivers import load_driver
from speaker.layout import compute_layout
from speaker.cad.panels import build_panels, panel_solid, assembly, flare_insert

P = load_params(); CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")
L = compute_layout(P, CX, B)

def test_panel_names():
    names = {pn.name for pn in build_panels(L, P, CX, B)}
    assert names == {"baffle_outer", "baffle_inner", "side_left", "side_right",
                     "back", "top", "bottom", "divider"}

def test_assembly_bbox_mm():
    bb = assembly(L, P, CX, B).toCompound().BoundingBox()
    assert (bb.xlen, bb.ylen, bb.zlen) == pytest.approx((320, 350, L.height * 1000), abs=0.5)

def test_side_has_six_vent_slots():
    side = next(pn for pn in build_panels(L, P, CX, B) if pn.name == "side_left")
    slots = [f for f in side.features if f.shape == "slot"]
    assert len(slots) == 6 and all(f.corner_r_mm >= 3 for f in slots)

def test_inner_corner_radii_at_least_3mm():
    for pn in build_panels(L, P, CX, B):
        for f in pn.features:
            if f.shape in ("rect", "slot"):
                assert f.corner_r_mm >= 3.0

def test_flare_insert_dims():
    bb = flare_insert(L, P, CX).val().BoundingBox()
    assert max(bb.xlen, bb.ylen) == pytest.approx(289.8, abs=0.2)
    # R7: insert occupies the outer baffle layer only (panel_t long).
    assert bb.zlen == pytest.approx(18.0, abs=0.2)

def _assert_no_intersections(l, p):
    # Every placed part of the assembly (8 panels + flare insert), pairwise.
    solids = {ch.name: ch.obj.val().moved(ch.loc) for ch in assembly(l, p, CX, B).children}
    assert len(solids) == 9
    for (na, a), (nb, b) in itertools.combinations(solids.items(), 2):
        vol = a.intersect(b).Volume()
        assert vol < 1.0, f"{na} and {nb} intersect by {vol:.1f} mm^3"

def test_panels_do_not_intersect():
    _assert_no_intersections(L, P)

def test_divider_uses_divider_t():
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, divider_t=0.025))
    l = compute_layout(p, CX, B)
    div = next(pn for pn in build_panels(l, p, CX, B) if pn.name == "divider")
    assert div.thickness_mm == pytest.approx(25.0)
    bb = panel_solid(div).val().moved(div.placement).BoundingBox()
    assert (bb.zmin, bb.zmax) == pytest.approx((l.divider_z * 1000, l.divider_z * 1000 + 25.0), abs=1e-6)
    _assert_no_intersections(l, p)

def test_baffle_must_be_two_panel_layers():
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, baffle_t=0.040))
    l = compute_layout(p, CX, B)
    with pytest.raises(ValueError, match=r"baffle_t 40\.0 mm.*panel_t 18\.0 mm"):
        build_panels(l, p, CX, B)
