import dataclasses
import math

import pytest
from speaker.drivers import load_driver
from speaker.fillers import VOLUME_TOL, corner_profile, polygon_area
from speaker.layout import compute_layout
from speaker.params import load_params

P = load_params(); CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")

def _cab(**kw):
    return dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, **kw))

def _vents(**kw):
    return dataclasses.replace(P, vents=dataclasses.replace(P.vents, **kw))

def test_polygon_area_unit_square():
    assert polygon_area([(0, 0), (1, 0), (1, 1), (0, 1)]) == pytest.approx(1.0)
    assert polygon_area([(0, 1), (1, 1), (1, 0), (0, 0)]) == pytest.approx(1.0)

def test_corner_profile_stays_outside_coax_sweep_and_inside_walls():
    pts = corner_profile(0.160, 0.950, 0.302, 1.100, 0.133, 0.040)
    for x, z in pts:
        assert math.hypot(x - 0.160, z - 0.950) >= 0.133 - 1e-9
        assert 0.160 + 0.040 - 1e-9 <= x <= 0.302 + 1e-9
        assert 0.950 + 0.040 - 1e-9 <= z <= 1.100 + 1e-9

def test_default_fillers():
    l = compute_layout(P, CX, B)
    f = l.fillers
    assert [n for n, _ in f.corners] == ["filler_corner_top_left", "filler_corner_top_right",
                                         "filler_corner_bottom_left", "filler_corner_bottom_right"]
    assert f.corner_y == pytest.approx((0.036, 0.205))
    assert f.corner_volume == pytest.approx(3.85e-3, abs=0.05e-3)
    assert f.block_min[0] == pytest.approx(0.043) and f.block_max[0] == pytest.approx(0.277)
    assert f.block_min[1] == pytest.approx(0.205) and f.block_max[1] == pytest.approx(0.332)
    assert f.block_max[2] - f.block_min[2] == pytest.approx(0.167, abs=0.002)
    assert l.coax_chamber_empty == pytest.approx(0.01883, abs=2e-4)
    assert l.coax_chamber_as_built == pytest.approx(0.010, abs=VOLUME_TOL)
    assert l.filler_volume == pytest.approx(f.total_volume)

def test_large_coax_frame_rejected_for_corner_fillers():
    with pytest.raises(ValueError, match="corner filler"):
        compute_layout(P, dataclasses.replace(CX, outer_d=0.284), B)  # a = 139.4 > 137 mm

def test_vents_too_far_forward_rejected():
    with pytest.raises(ValueError, match="vent"):
        compute_layout(_vents(centre_from_front=0.100), CX, B)

def test_unreachable_small_target_rejected():
    with pytest.raises(ValueError, match="rear block"):
        compute_layout(_cab(coax_chamber_net_volume=0.004), CX, B)

def test_corners_alone_overshoot_rejected():
    with pytest.raises(ValueError, match="corner fillers alone"):
        compute_layout(_cab(coax_chamber_net_volume=0.016), CX, B)

@pytest.mark.parametrize("ear", [0.85, 1.05])
def test_fillers_follow_ear_height(ear):
    l = compute_layout(_cab(coax_centre_z=ear), CX, B)
    assert l.coax_chamber_as_built == pytest.approx(0.010, abs=VOLUME_TOL)
    for _, pts in l.fillers.corners:
        assert all(math.hypot(x - 0.160, z - ear) >= CX.outer_d / 2 for x, z in pts)
