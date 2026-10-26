import cadquery as cq
import pytest
from speaker.cad.fillers import filler_solids
from speaker.drivers import load_driver
from speaker.layout import compute_layout
from speaker.params import load_params

P = load_params(); CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")
L = compute_layout(P, CX, B)

def test_filler_volumes_match_layout():
    s = filler_solids(L)
    assert set(s) == {n for n, _ in L.fillers.corners} | {"filler_rear_block"}
    corners = sum(s[n].val().Volume() for n, _ in L.fillers.corners)
    assert corners == pytest.approx(L.fillers.corner_volume * 1e9, rel=1e-4)
    assert s["filler_rear_block"].val().Volume() == pytest.approx(L.fillers.block_volume * 1e9, rel=1e-6)

def test_corner_pieces_clear_coax_sweep():
    # The coax slides in from the back past the corner pieces; the rear block
    # comes out with the back panel, so only the corners must clear the sweep.
    cx, cz = (v * 1000 for v in L.coax_centre)
    sweep = cq.Solid.makeCylinder(CX.outer_d * 500, 296, cq.Vector(cx, 36, cz), cq.Vector(0, 1, 0))
    s = filler_solids(L)
    for name, _ in L.fillers.corners:
        assert s[name].val().intersect(sweep).Volume() < 1.0, name

def test_rear_block_clears_installed_coax():
    coax_rear_mm = (P.cabinet.baffle_t + CX.depth - CX.flange_t) * 1000  # 172 mm
    assert filler_solids(L)["filler_rear_block"].val().BoundingBox().ymin >= coax_rear_mm + 10 - 1e-6

def test_corner_pieces_stop_before_vents_and_block_clears_sides():
    s = filler_solids(L)
    for n, _ in L.fillers.corners:
        bb = s[n].val().BoundingBox()
        assert bb.ymin == pytest.approx(36, abs=1e-6) and bb.ymax == pytest.approx(205, abs=1e-6)
    bb = s["filler_rear_block"].val().BoundingBox()
    assert bb.xmin >= 18 + 25 - 1e-6 and bb.xmax <= 302 - 25 + 1e-6

def test_default_parts_fit_250mm_printer():
    for name, wp in filler_solids(L).items():
        bb = wp.val().BoundingBox()
        assert max(bb.xlen, bb.ylen, bb.zlen) < 250, name
