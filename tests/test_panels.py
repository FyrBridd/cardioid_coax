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

def test_panels_do_not_intersect():
    # Every placed part of the assembly (8 panels + flare insert), pairwise.
    solids = {ch.name: ch.obj.val().moved(ch.loc) for ch in assembly(L, P, CX, B).children}
    assert len(solids) == 9
    for (na, a), (nb, b) in itertools.combinations(solids.items(), 2):
        vol = a.intersect(b).Volume()
        assert vol < 1.0, f"{na} and {nb} intersect by {vol:.1f} mm^3"
