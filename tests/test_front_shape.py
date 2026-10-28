import dataclasses
import math
import re
from pathlib import Path

import pytest

from speaker.drivers import load_driver
from speaker.fillers import polygon_area
from speaker.front_shape import CHAMFER_MAX, FAMILIES, FRONT_LAND, cut_depth, flat_width, section
from speaker.layout import compute_layout
from speaker.params import load_params

CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")
P = load_params()
W, D = P.cabinet.width, P.cabinet.depth
DEFAULT_TOML = Path("params/default.toml").read_text()

# A facet design that passes every rule: the mouth is narrowed so the 30 mm facets leave it room.
FACETS = dict(front_edge='"facets"', facet_angle_deg=20, facet_width_mm=30, facet_base_z_mm=700,
              facet_ramp_mm=50, flare_mouth_d_mm=250)


def toml(**keys) -> str:
    text = DEFAULT_TOML
    for k, v in keys.items():
        text, n = re.subn(rf"^{k}\s*=.*$", f"{k} = {v}", text, flags=re.M)
        assert n == 1, k
    return text


def layout(tmp_path, **keys):
    f = tmp_path / "p.toml"
    f.write_text(toml(**keys))
    return compute_layout(load_params(f), CX, B)


def cab(**kw):
    return dataclasses.replace(P.cabinet, **kw)


def test_constants():
    assert FRONT_LAND == 0.005 and CHAMFER_MAX == 0.045
    assert FAMILIES == ("square", "roundover", "chamfer", "facets")


def test_valid_facets_pass(tmp_path):
    layout(tmp_path, **FACETS)


@pytest.mark.parametrize("keys", [dict(front_edge='"roundover"', edge_size_mm=20, flare_mouth_d_mm=250),
                                  dict(front_edge='"chamfer"', edge_size_mm=8)])
def test_valid_edges_pass(tmp_path, keys):
    layout(tmp_path, **keys)


def test_a_family(tmp_path):
    with pytest.raises(ValueError, match=r"front_edge 'bevel' is not one of"):
        layout(tmp_path, front_edge='"bevel"')


@pytest.mark.parametrize("keys, msg", [
    (dict(edge_size_mm=5), r"edge_size 5\.0 mm must be 0 for square"),
    (dict(front_edge='"facets"', edge_size_mm=5, **{k: v for k, v in FACETS.items() if k != "front_edge"}),
     r"edge_size 5\.0 mm must be 0 for facets"),
    (dict(front_edge='"roundover"'), r"edge_size 0\.0 mm must be > 0 for roundover"),
    (dict(front_edge='"chamfer"', edge_size_mm=50), r"chamfer edge_size 50\.0 mm exceeds CHAMFER_MAX 45\.0 mm"),
])
def test_b_edge_size(tmp_path, keys, msg):
    with pytest.raises(ValueError, match=msg):
        layout(tmp_path, **keys)


@pytest.mark.parametrize("keys, msg", [
    (dict(FACETS, facet_angle_deg=50), r"facet_angle 50\.0° must be > 0 and ≤ 45°"),
    (dict(FACETS, facet_angle_deg=0), r"facet_angle 0\.0° must be > 0 and ≤ 45°"),
    (dict(FACETS, facet_width_mm=0), r"facet_width 0\.0 mm must be > 0"),
    (dict(FACETS, facet_ramp_mm=-1), r"facet_ramp -1\.0 mm must be ≥ 0"),
    (dict(front_edge='"chamfer"', edge_size_mm=10, facet_angle_deg=10),
     r"chamfer.*facet values must all be 0.*facet_angle 10\.0°"),
])
def test_c_facet_values(tmp_path, keys, msg):
    with pytest.raises(ValueError, match=msg):
        layout(tmp_path, **keys)


@pytest.mark.parametrize("keys", [dict(front_edge='"roundover"', edge_size_mm=100),
                                  dict(FACETS, facet_angle_deg=45, facet_width_mm=100)])
def test_d_treatment_reaches_interior(tmp_path, keys):
    with pytest.raises(ValueError, match=r"reaches the interior.*front corner"):
        layout(tmp_path, **keys)


@pytest.mark.parametrize("keys", [dict(front_edge='"chamfer"', edge_size_mm=20), dict(coax_offset_x_mm=11)])
def test_e_flat_width_at_coax(tmp_path, keys):
    with pytest.raises(ValueError, match=r"flat front width .* at the coax centre height .* below"):
        layout(tmp_path, **keys)


def test_e_offset_within_land_passes(tmp_path):
    assert layout(tmp_path, coax_offset_x_mm=8).coax_centre[0] == pytest.approx(W / 2 + 0.008)


def test_f_flat_width_at_bass(tmp_path):
    with pytest.raises(ValueError, match=r"flat front width 240\.0 mm at bass centre height 200\.0 mm is below .* 245\.0 mm"):
        layout(tmp_path, front_edge='"roundover"', edge_size_mm=40, flare_mouth_d_mm=200)


def test_g_coax_frame_clears_sides(tmp_path):
    with pytest.raises(ValueError, match=r"coax_offset_x 15\.0 mm.*side panels.*142\.0 mm"):
        layout(tmp_path, coax_offset_x_mm=15, flare_mouth_d_mm=250)


def test_g_negative_offset(tmp_path):
    with pytest.raises(ValueError, match=r"coax_offset_x -15\.0 mm.*side panels"):
        layout(tmp_path, coax_offset_x_mm=-15, flare_mouth_d_mm=250)


@pytest.mark.parametrize("keys, msg", [
    (dict(FACETS, facet_base_z_mm=600), r"facet_base_z 600\.0 mm is below .*upper bass cutout.* 642\.5 mm"),
    (dict(FACETS, facet_ramp_mm=150), r"facet_base_z \+ facet_ramp 850\.0 mm is above .* 825\.0 mm"),
])
def test_h_facet_heights(tmp_path, keys, msg):
    with pytest.raises(ValueError, match=msg):
        layout(tmp_path, **keys)


def test_i_vents_behind_facet_cut(tmp_path):
    with pytest.raises(ValueError, match=r"vent slots start at y=40\.0 mm.*facet cut.*45\.0 mm"):
        layout(tmp_path, **dict(FACETS, facet_angle_deg=45, facet_width_mm=40, flare_mouth_d_mm=200),
               centre_from_front_mm=90)


@pytest.mark.parametrize("keys, msg", [
    (dict(flare_profile='"horn"'), r"flare_profile 'horn' is not one of"),
    (dict(outer_baffle_t_mm=12), r"outer_baffle_t 12\.0 mm must be between panel_t 18\.0 mm and 3 x panel_t 54\.0 mm"),
    (dict(outer_baffle_t_mm=60), r"outer_baffle_t 60\.0 mm must be between"),
])
def test_j_flare_and_baffle(tmp_path, keys, msg):
    with pytest.raises(ValueError, match=msg):
        layout(tmp_path, **keys)


def test_k_bass_net_target(tmp_path):
    with pytest.raises(ValueError, match=r"bass net volume 45\.7 L is below the 50 L target"):
        layout(tmp_path, coax_centre_z_mm=850)


def test_baffle_t_is_outer_plus_panel(tmp_path):
    l = layout(tmp_path, outer_baffle_t_mm=30)
    assert l.interior_d == pytest.approx(D - 0.048 - 0.018)


def test_cut_depth_and_flat_width():
    assert cut_depth(cab()) == 0 and flat_width(cab(), 0.5) == pytest.approx(W)
    assert cut_depth(cab(front_edge="chamfer", edge_size=0.03)) == pytest.approx(0.03)
    assert flat_width(cab(front_edge="roundover", edge_size=0.02), 0.1) == pytest.approx(W - 0.04)
    fa = cab(front_edge="facets", facet_angle=math.radians(20), facet_width=0.04, facet_base_z=0.7, facet_ramp=0.1)
    assert cut_depth(fa) == pytest.approx(0.04 * math.tan(math.radians(20)))
    assert flat_width(fa, 0.69) == pytest.approx(W)
    assert flat_width(fa, 0.75) == pytest.approx(W - 0.04)   # half way up the ramp
    assert flat_width(fa, 0.9) == pytest.approx(W - 0.08)
    step = dataclasses.replace(fa, facet_ramp=0.0)
    assert flat_width(step, 0.699) == pytest.approx(W) and flat_width(step, 0.7) == pytest.approx(W - 0.08)


def test_section_areas():
    sq = cab(front_edge="square")
    assert polygon_area(section(sq, 0.5)) == pytest.approx(W * D)
    ch = cab(front_edge="chamfer", edge_size=0.03)
    assert polygon_area(section(ch, 0.5)) == pytest.approx(W * D - 0.03**2)          # two s²/2 corners
    ro = cab(front_edge="roundover", edge_size=0.03)
    assert polygon_area(section(ro, 0.5), ) == pytest.approx(W * D - 2 * (1 - math.pi / 4) * 0.03**2, rel=1e-3)
    fa = cab(front_edge="facets", facet_angle=math.radians(20), facet_width=0.04, facet_base_z=0.7, facet_ramp=0.1)
    d = 0.04 * math.tan(math.radians(20))
    assert polygon_area(section(fa, 0.69)) == pytest.approx(W * D)
    assert polygon_area(section(fa, 0.9)) == pytest.approx(W * D - 0.04 * d)        # two fw·d/2 triangles


@pytest.mark.parametrize("c", [cab(), cab(front_edge="chamfer", edge_size=0.03),
                               cab(front_edge="roundover", edge_size=0.03),
                               cab(front_edge="facets", facet_angle=0.3, facet_width=0.04, facet_base_z=0.7)])
def test_section_is_counter_clockwise_and_simple(c):
    pts = section(c, 0.9)
    n = len(pts)
    signed = sum(pts[i][0] * pts[(i + 1) % n][1] - pts[(i + 1) % n][0] * pts[i][1] for i in range(n)) / 2
    assert signed > 0
    assert len(set(pts)) == n                                       # no repeated vertices
    assert all(-1e-12 <= x <= W + 1e-12 and -1e-12 <= y <= D + 1e-12 for x, y in pts)
