import dataclasses
import math

import numpy as np
import pytest

from fakes import default_design
from speaker.bem.mesh import (DISC_R, PUCK_EXTRA_DEPTH, box_piston_mesh, cabinet_mesh, disc_baffle_mesh, edge_check,
                              mesh_volume, puck_mesh)
from speaker.coax_geometry import cavity_profile
from speaker.drivers import load_driver
from speaker.layout import compute_layout
from speaker.params import load_params

P = load_params()
CX = load_driver("bms_10c262")
B = load_driver("dayton_rss265ho_44")


def _area(m, name):
    i = m.tag_names.index(name) + 1
    v = m.vertices[m.triangles[m.tags == i]]
    return 0.5 * np.linalg.norm(np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=1).sum()


def test_cabinet_mesh_closed_and_tagged():
    m = cabinet_mesh(P, CX, B, h=0.05)
    edge_check(m)
    assert m.tag_names == ("coax_lf", "coax_hf", "vent_left", "vent_right", "bass_lower", "bass_upper")
    assert _area(m, "bass_lower") == pytest.approx(math.pi * 0.2118**2 / 4, rel=0.03)
    assert _area(m, "coax_hf") == pytest.approx(math.pi * 0.0254**2 / 4, rel=0.05)
    assert _area(m, "vent_left") == pytest.approx(0.100 * (11 * 0.022 + 0.010), rel=0.03)
    slant = math.hypot(0.065, (0.202 - 0.0636) / 2)
    assert _area(m, "coax_lf") == pytest.approx(math.pi * (0.101 + 0.0318) * slant, rel=0.03)


def test_vents_on_correct_sides_and_coax_region_refined():
    m = cabinet_mesh(P, CX, B, h=0.05)
    c = m.vertices[m.triangles].mean(axis=1)
    assert np.allclose(c[m.tags == 3][:, 0], 0.0, atol=1e-6)  # vent_left at x = 0
    assert np.allclose(c[m.tags == 4][:, 0], 0.320, atol=1e-6)  # vent_right at x = W
    e = np.linalg.norm(
        m.vertices[m.triangles[m.tags == 1][:, 0]] - m.vertices[m.triangles[m.tags == 1][:, 1]], axis=1
    )
    assert e.max() < 0.010 * 1.5


def test_puck_and_test_box_closed():
    edge_check(puck_mesh(P, CX, h=0.02, h_coax=0.01))
    edge_check(box_piston_mesh((0.2, 0.15, 0.3), 0.015, h=0.02))


def test_edge_check_rejects_open_mesh():
    from speaker.bem.solver import SurfaceMesh

    m = SurfaceMesh(np.eye(3), np.array([[0, 1, 2]]), np.zeros(1, int), ())
    with pytest.raises(ValueError, match="edge"):
        edge_check(m)


def test_edge_check_rejects_inward_and_inconsistent_normals():
    m = box_piston_mesh((0.2, 0.15, 0.3), 0.015, h=0.05)
    inward = dataclasses.replace(m, triangles=m.triangles[:, ::-1])
    with pytest.raises(ValueError, match="inward"):
        edge_check(inward)
    mixed = m.triangles.copy()
    mixed[0] = mixed[0, ::-1]
    with pytest.raises(ValueError, match="inconsistent"):
        edge_check(dataclasses.replace(m, triangles=mixed))


@pytest.mark.parametrize("recess", [True, False])
def test_disc_baffle_mesh_closed_and_tagged(recess):
    m = disc_baffle_mesh(load_driver("bms_10c262"), h=0.08, h_coax=0.02, recess=recess)
    assert m.tag_names == ("coax_lf",)
    assert np.count_nonzero(m.tags == 1) > 0  # edge_check already ran inside the mesher
    if recess:   # the cone, set back by the frame's flange
        slant = math.hypot(CX.cone_depth, (CX.diaphragm_d - CX.cone_inner_d) / 2)
        assert _area(m, "coax_lf") == pytest.approx(math.pi * (CX.diaphragm_d + CX.cone_inner_d) / 2 * slant, rel=0.03)
        y = m.vertices[m.triangles[m.tags == 1]][..., 1]
        assert y.min() == pytest.approx(CX.flange_t, abs=1e-6)
        assert y.max() == pytest.approx(CX.flange_t + CX.cone_depth, abs=1e-6)
    else:   # a flat piston of the diaphragm on the front face
        assert _area(m, "coax_lf") == pytest.approx(math.pi * CX.diaphragm_d**2 / 4, rel=0.03)
        assert np.allclose(m.vertices[m.triangles[m.tags == 1]][..., 1], 0.0, atol=1e-9)
    r = np.hypot(m.vertices[:, 0], m.vertices[:, 2])
    assert r.max() == pytest.approx(DISC_R, rel=1e-6)
    assert m.vertices[:, 1].max() == pytest.approx(CX.flange_t + CX.cone_depth + PUCK_EXTRA_DEPTH, abs=1e-6)


# --- the front-edge families, the coax offset, the outer baffle and the flare profile -------------

def _with(p, **kw):
    """p with cabinet fields replaced; baffle_t follows outer_baffle_t."""
    c = dataclasses.replace(p.cabinet, **kw)
    c = dataclasses.replace(c, baffle_t=c.outer_baffle_t + c.panel_t)
    return dataclasses.replace(p, cabinet=c)


BASE = dict(width=0.330, flare_mouth_d=0.240)  # room on the front for every family's rules
SHAPES = {
    "roundover": dict(front_edge="roundover", edge_size=0.020),
    "chamfer": dict(front_edge="chamfer", edge_size=0.030),
    "facets_step": dict(front_edge="facets", facet_angle=math.radians(30), facet_width=0.035, facet_base_z=0.700),
    "facets_ramp": dict(front_edge="facets", facet_angle=math.radians(30), facet_width=0.035, facet_base_z=0.700,
                        facet_ramp=0.100),
}
SQUARE = dict(front_edge="square", edge_size=0.0, facet_angle=0.0, facet_width=0.0, facet_base_z=0.0, facet_ramp=0.0)


def shaped_params(family):
    return _with(P, **BASE, **SHAPES[family])


def square_like(p):
    return _with(p, **SQUARE)


def fd(c):
    return c.facet_width * math.tan(c.facet_angle)


CUTS = {  # analytic volume removed by both sides together
    "roundover": lambda c, H: 2 * (1 - math.pi / 4) * c.edge_size**2 * H,
    "chamfer": lambda c, H: c.edge_size**2 * H,
    "facets_step": lambda c, H: c.facet_width * fd(c) * (H - c.facet_base_z),
    "facets_ramp": lambda c, H: c.facet_width * fd(c) * (H - c.facet_base_z - c.facet_ramp)
    + c.facet_width * fd(c) * c.facet_ramp / 3,
}


def _centre(m, name):
    """Area-weighted mean of the named source's triangle centroids."""
    v = m.vertices[m.triangles[m.tags == m.tag_names.index(name) + 1]]
    a = 0.5 * np.linalg.norm(np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=1)
    return (v.mean(axis=1) * a[:, None]).sum(axis=0) / a.sum()


@pytest.mark.parametrize("family", list(CUTS))
def test_front_edge_cut_volumes(family):
    p_x = shaped_params(family)
    p_sq = square_like(p_x)
    v_sq = mesh_volume(cabinet_mesh(p_sq, CX, B, h=0.04, h_coax=0.02))
    v_x = mesh_volume(cabinet_mesh(p_x, CX, B, h=0.04, h_coax=0.02))
    H = compute_layout(p_x, CX, B).height
    assert v_sq - v_x == pytest.approx(CUTS[family](p_x.cabinet, H), rel=0.03)


def test_coax_offset_moves_the_cone():
    p = _with(P, coax_offset_x=0.008)
    m = cabinet_mesh(p, CX, B, h=0.04, h_coax=0.02)
    assert _centre(m, "coax_lf")[0] == pytest.approx(P.cabinet.width / 2 + 0.008, abs=1e-3)
    assert _centre(m, "coax_hf")[0] == pytest.approx(P.cabinet.width / 2 + 0.008, abs=1e-3)


def test_thicker_outer_baffle_deepens_the_cavity():
    y = [_centre(cabinet_mesh(_with(P, outer_baffle_t=t), CX, B, h=0.04, h_coax=0.02), "coax_lf")[1]
         for t in (0.018, 0.030)]
    assert y[1] - y[0] == pytest.approx(0.012, abs=1e-3)


def _revolved(poly):
    """Volume swept by a closed (r, y) polygon about the axis: 2π ∫∫ r dA (Pappus)."""
    total = 0.0
    for (r0, y0), (r1, y1) in zip(poly, poly[1:] + poly[:1]):
        total += (r0 * y1 - r1 * y0) * (r0 + r1) / 6  # ∫∫ r dA over the edge's triangle with the origin
    return 2 * math.pi * abs(total)


def test_flare_profile_changes_cavity_volume_as_revolved():
    # h_coax 5 mm: at 20 mm the ellipse's ~50 mm profile gets about three chords and the mesh keeps
    # only 87 % of the difference (97 % at the run's 10 mm).
    p_el = _with(P, outer_baffle_t=0.036)
    p_co = _with(p_el, flare_profile="cone")
    v = {n: mesh_volume(cabinet_mesh(q, CX, B, h=0.04, h_coax=0.005)) for n, q in (("el", p_el), ("co", p_co))}
    t = p_el.cabinet.outer_baffle_t
    swept = {n: _revolved([(0.0, 0.0), *dict(cavity_profile(q, CX))["flare"], (0.0, t)])
             for n, q in (("el", p_el), ("co", p_co))}
    assert v["el"] - v["co"] == pytest.approx(swept["co"] - swept["el"], rel=0.03)


@pytest.mark.parametrize("profile", ["ellipse", "cone"])
@pytest.mark.parametrize("outer", [0.018, 0.030])
def test_puck_mesh_builds_each_flare(profile, outer):
    m = puck_mesh(_with(P, flare_profile=profile, outer_baffle_t=outer), CX, h=0.04, h_coax=0.02)
    edge_check(m)
    assert m.vertices[:, 1].max() == pytest.approx(outer + P.cabinet.panel_t + CX.cone_depth + PUCK_EXTRA_DEPTH,
                                                   abs=1e-6)


DEFAULT_TRIANGLES = 4046  # cabinet_mesh(default params, h=0.04, h_coax=0.02), captured before Task 4
DEFAULT_VOLUME = 0.12264353980870406  # m^3, the same mesh


def test_default_cabinet_mesh_unchanged():
    d = default_design()
    m = cabinet_mesh(d.p, d.coax, d.bass, h=0.04, h_coax=0.02)
    assert len(m.triangles) == DEFAULT_TRIANGLES
    assert mesh_volume(m) == pytest.approx(DEFAULT_VOLUME, rel=1e-9)
