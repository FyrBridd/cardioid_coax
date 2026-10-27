import dataclasses
import math

import numpy as np
import pytest

from speaker.bem.mesh import (DISC_R, PUCK_EXTRA_DEPTH, box_piston_mesh, cabinet_mesh, disc_baffle_mesh, edge_check,
                              puck_mesh)
from speaker.drivers import load_driver
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
