"""Closed, tagged surface meshes for the BEM solver, built with gmsh (OpenCASCADE). SI units (m).

Layout frame: origin bottom-left-front, x across the width, y backwards (front face
at y=0), z up. Triangles are oriented outward (positive signed volume). Tag 0 is
rigid; tag i >= 1 is source tag_names[i-1].
"""
from __future__ import annotations

import math
from collections import deque
from collections.abc import Iterator
from contextlib import contextmanager

import gmsh
import numpy as np

from speaker.bem.solver import SurfaceMesh
from speaker.coax_geometry import cavity_profile
from speaker.drivers import Driver
from speaker.layout import compute_layout
from speaker.params import Params

CAVITY_REFINE = 0.020  # m, h_coax applies within this distance of the coax cavity
PUCK_EXTRA_DEPTH = 0.030  # m, puck depth behind the cone's inner edge
PUCK_FILLET = 0.020  # m, front rim radius of the puck
FRONT_EXT = 0.010  # m, cavity extends this far in front of y=0 so the cut is not coplanar
MIN_CIRCLE_NODES = 24  # per circle: an inscribed polygon then loses about 1 % of a disc's area
_MATCH_TOL = 1e-5  # m on centre of mass; relative on area

CABINET_TAGS = ("coax_lf", "coax_hf", "vent_left", "vent_right", "bass_lower", "bass_upper")
PUCK_TAGS = ("coax_lf", "coax_hf")

Signature = tuple[tuple[float, float, float], float]  # (centre of mass, area)


def cabinet_mesh(
    p: Params, coax: Driver, bass: Driver, *, h: float, h_coax: float = 0.010, h_vent: float = 0.020
) -> SurfaceMesh:
    """Cabinet exterior less the coax cavity, with the coax, vent and bass sources imprinted."""
    lay = compute_layout(p, coax, bass)
    w, d = p.cabinet.width, p.cabinet.depth
    v = p.vents
    zs = [z for _, z, _, _ in lay.vent_slots]
    y0, y1 = v.centre_from_front - v.slot_w / 2, v.centre_from_front + v.slot_w / 2
    z0, z1 = min(zs) - v.slot_h / 2, max(zs) + v.slot_h / 2
    r_bass = bass.diaphragm_d / 2
    with _session():
        occ = gmsh.model.occ
        body = occ.addBox(0, 0, 0, w, d, lay.height)
        cavity, cavity_box, coax_sigs = _cavity(p, coax, lay.coax_centre)
        occ.cut([(3, body)], [(3, cavity)])
        vents, vent_sigs = [], []
        for x in (0.0, w):
            vents.append(_quad([(x, y0, z0), (x, y1, z0), (x, y1, z1), (x, y0, z1)]))
            vent_sigs.append(((x, (y0 + y1) / 2, (z0 + z1) / 2), (y1 - y0) * (z1 - z0)))
        discs, disc_sigs = [], []
        for x, z in sorted(lay.bass_centres, key=lambda c: c[1]):
            discs.append(occ.addDisk(x, 0, z, r_bass, r_bass, zAxis=[0, 1, 0], xAxis=[1, 0, 0]))
            disc_sigs.append(((x, 0.0, z), math.pi * r_bass**2))
        occ.fragment(occ.getEntities(3), [(2, s) for s in vents + discs])
        occ.synchronize()
        fields = [_cavity_field(cavity_box, h, h_coax)]
        for x in (0.0, w):
            fields.append(_box_field((x, y0, z0), (x, y1, z1), min(h, h_vent), h))
        return _mesh_and_read(h, fields, CABINET_TAGS, coax_sigs + vent_sigs + disc_sigs)


def puck_mesh(p: Params, coax: Driver, *, h: float, h_coax: float) -> SurfaceMesh:
    """A filleted cylinder of the cabinet's width, coaxial with the coax, less the coax cavity."""
    cx, cz = p.cabinet.width / 2, p.cabinet.coax_centre_z
    r = p.cabinet.width / 2
    depth = p.cabinet.baffle_t + coax.cone_depth + PUCK_EXTRA_DEPTH
    with _session():
        occ = gmsh.model.occ
        body = occ.addCylinder(cx, 0, cz, 0, depth, 0, r)
        rim = occ.getEntitiesInBoundingBox(cx - 2 * r, -1e-6, cz - 2 * r, cx + 2 * r, 1e-6, cz + 2 * r, dim=1)
        body = occ.fillet([body], [c for _, c in rim], [PUCK_FILLET])[0][1]
        cavity, cavity_box, coax_sigs = _cavity(p, coax, (cx, cz))
        occ.cut([(3, body)], [(3, cavity)])
        occ.synchronize()
        return _mesh_and_read(h, [_cavity_field(cavity_box, h, h_coax)], PUCK_TAGS, coax_sigs)


def box_piston_mesh(size: tuple[float, float, float], piston_r: float, *, h: float) -> SurfaceMesh:
    """Box [0, sx] x [0, sy] x [0, sz] with a flush piston centred on its front face (y=0)."""
    sx, sy, sz = size
    with _session():
        occ = gmsh.model.occ
        body = occ.addBox(0, 0, 0, sx, sy, sz)
        disc = occ.addDisk(sx / 2, 0, sz / 2, piston_r, piston_r, zAxis=[0, 1, 0], xAxis=[1, 0, 0])
        occ.fragment([(3, body)], [(2, disc)])
        occ.synchronize()
        sig = ((sx / 2, 0.0, sz / 2), math.pi * piston_r**2)
        return _mesh_and_read(h, [], ("piston",), [sig])


def edge_check(m: SurfaceMesh) -> None:
    """Raise ValueError unless the mesh is closed, consistently oriented and outward-facing."""
    t = np.asarray(m.triangles)
    directed = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    _, counts = np.unique(np.sort(directed, axis=1), axis=0, return_counts=True)
    if np.any(counts != 2):
        raise ValueError(
            f"mesh is not closed: {np.count_nonzero(counts != 2)} of {len(counts)} edges "
            f"are not shared by exactly two triangles"
        )
    if len(np.unique(directed, axis=0)) != len(directed):
        raise ValueError("triangle orientation is inconsistent: an edge is traversed twice in the same direction")
    vol = _signed_volume(m.vertices, t)
    if vol <= 0:
        raise ValueError(f"signed volume {vol:.3g} m^3 is not positive: normals point inward")


@contextmanager
def _session() -> Iterator[None]:
    gmsh.initialize(readConfigFiles=False, interruptible=False)
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        yield
    finally:
        gmsh.finalize()


def _cavity(
    p: Params, coax: Driver, centre: tuple[float, float]
) -> tuple[int, tuple[float, ...], list[Signature]]:
    """Revolve the cavity profile about the coax axis (+y through centre=(x, z)).

    Returns the volume tag, its bounding box, and the (centre, area) signatures of
    the cone (coax_lf) and the HF exit (coax_hf) surfaces.
    """
    occ = gmsh.model.occ
    cx, cz = centre
    prof = dict(cavity_profile(p, coax))
    segs = list(prof.values())
    r_mouth, y_end = segs[0][0][0], segs[-1][-1][1]

    def pt(r: float, y: float) -> int:
        return occ.addPoint(cx + r, y, cz)

    start = pt(r_mouth, -FRONT_EXT)
    tags = [pt(*segs[0][0])]
    curves = [occ.addLine(start, tags[0])]
    for seg in segs:
        inner = [pt(*q) for q in seg[1:-1]]
        tags.append(pt(*seg[-1]))
        ends = [tags[-2], *inner, tags[-1]]
        curves.append(occ.addSpline(ends) if inner else occ.addLine(*ends))
    axis_front = pt(0.0, -FRONT_EXT)
    curves += [occ.addLine(tags[-1], axis_front), occ.addLine(axis_front, start)]
    face = occ.addPlaneSurface([occ.addCurveLoop(curves)])
    vol = next(t for d, t in occ.revolve([(2, face)], cx, 0, cz, 0, 1, 0, 2 * math.pi) if d == 3)
    box = (cx - r_mouth, -FRONT_EXT, cz - r_mouth, cx + r_mouth, y_end, cz + r_mouth)
    sigs = [_revolved_signature(prof[n], cx, cz) for n in ("cone", "hf_exit")]
    return vol, box, sigs


def _revolved_signature(seg: tuple[tuple[float, float], ...], cx: float, cz: float) -> Signature:
    """Centre of mass and area of a straight (r, y) segment revolved about the axis."""
    (r1, y1), (r2, y2) = seg
    dr, dy = r2 - r1, y2 - y1
    area = math.pi * (r1 + r2) * math.hypot(dr, dy)
    y_c = (r1 * y1 + (r1 * dy + y1 * dr) / 2 + dr * dy / 3) / ((r1 + r2) / 2)
    return (cx, y_c, cz), area


def _quad(corners: list[tuple[float, float, float]]) -> int:
    occ = gmsh.model.occ
    pts = [occ.addPoint(*c) for c in corners]
    lines = [occ.addLine(pts[i], pts[(i + 1) % 4]) for i in range(4)]
    return occ.addPlaneSurface([occ.addCurveLoop(lines)])


def _cavity_field(box: tuple[float, ...], h: float, h_coax: float) -> int:
    """min(h, h_coax) within CAVITY_REFINE of the cavity surfaces, grading to h over a further h."""
    e = 1e-3  # OCC bounding boxes of curved surfaces are slightly loose
    lo, hi = np.array(box[:3]) - e, np.array(box[3:]) + e
    surfs = [t for _, t in gmsh.model.getEntitiesInBoundingBox(*lo, *hi, dim=2)]
    f = gmsh.model.mesh.field
    dist = f.add("Distance")
    f.setNumbers(dist, "SurfacesList", surfs)
    f.setNumber(dist, "Sampling", 100)
    th = f.add("Threshold")
    f.setNumber(th, "InField", dist)
    f.setNumber(th, "SizeMin", min(h, h_coax))
    f.setNumber(th, "SizeMax", h)
    f.setNumber(th, "DistMin", CAVITY_REFINE)
    f.setNumber(th, "DistMax", CAVITY_REFINE + h)
    return th


def _box_field(lo: tuple[float, float, float], hi: tuple[float, float, float], v_in: float, v_out: float) -> int:
    f = gmsh.model.mesh.field
    b = f.add("Box")
    e = 1e-6
    for k, (a, z) in zip("XYZ", zip(lo, hi)):
        f.setNumber(b, f"{k}Min", a - e)
        f.setNumber(b, f"{k}Max", z + e)
    f.setNumber(b, "VIn", v_in)
    f.setNumber(b, "VOut", v_out)
    return b


def _mesh_and_read(h: float, fields: list[int], names: tuple[str, ...], sigs: list[Signature]) -> SurfaceMesh:
    """Mesh the single volume's boundary and read it back with tags from the source signatures."""
    f = gmsh.model.mesh.field
    if fields:
        bg = f.add("Min")
        f.setNumbers(bg, "FieldsList", fields)
        f.setAsBackgroundMesh(bg)
    for name, value in (
        ("Mesh.Algorithm", 6),
        ("Mesh.MeshSizeMax", h),
        ("Mesh.MeshSizeFromPoints", 0),
        ("Mesh.MeshSizeFromCurvature", 0),
        ("Mesh.MeshSizeExtendFromBoundary", 0),
        ("Mesh.MinimumCircleNodes", MIN_CIRCLE_NODES),
    ):
        gmsh.option.setNumber(name, value)
    (vol,) = gmsh.model.getEntities(3)
    surfs = [t for _, t in gmsh.model.getBoundary([vol], combined=False, oriented=False)]
    tag_of = dict.fromkeys(surfs, 0)
    for i, sig in enumerate(sigs, start=1):
        tag_of[_find_surface(surfs, sig, names[i - 1])] = i
    gmsh.model.mesh.generate(2)

    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    row = np.empty(int(node_tags.max()) + 1, int)
    row[node_tags.astype(int)] = np.arange(len(node_tags))
    tris, tags = [], []
    for s, tag in tag_of.items():
        _, nodes = gmsh.model.mesh.getElementsByType(2, s)
        tris.append(row[nodes.astype(int)].reshape(-1, 3))
        tags.append(np.full(len(tris[-1]), tag))
    tri = np.concatenate(tris)
    used, tri = np.unique(tri, return_inverse=True)
    verts = coords.reshape(-1, 3)[used]
    m = SurfaceMesh(verts, _orient_outward(verts, tri.reshape(-1, 3)), np.concatenate(tags), names)
    edge_check(m)
    return m


def _find_surface(surfs: list[int], sig: Signature, what: str) -> int:
    com, area = sig
    for s in surfs:
        if np.allclose(gmsh.model.occ.getCenterOfMass(2, s), com, atol=_MATCH_TOL) and math.isclose(
            gmsh.model.occ.getMass(2, s), area, rel_tol=_MATCH_TOL
        ):
            return s
    raise RuntimeError(f"no surface for {what!r} with centre of mass {com} and area {area:.4g} m^2")


def _orient_outward(v: np.ndarray, tri: np.ndarray) -> np.ndarray:
    """Make orientation consistent over shared edges, then outward per connected component."""
    tri = tri.copy()
    by_edge: dict[tuple[int, int], list[int]] = {}
    for i, (a, b, c) in enumerate(tri):
        for e in ((a, b), (b, c), (c, a)):
            by_edge.setdefault((min(e), max(e)), []).append(i)
    seen = np.zeros(len(tri), bool)
    for seed in range(len(tri)):
        if seen[seed]:
            continue
        seen[seed] = True
        comp, queue = [seed], deque([seed])
        while queue:
            i = queue.popleft()
            a, b, c = tri[i]
            for x, y in ((a, b), (b, c), (c, a)):
                for j in by_edge[(min(x, y), max(x, y))]:
                    if seen[j]:
                        continue
                    if _has_directed_edge(tri[j], x, y):
                        tri[j] = tri[j, ::-1]
                    seen[j] = True
                    comp.append(j)
                    queue.append(j)
        if _signed_volume(v, tri[comp]) < 0:
            tri[comp] = tri[comp, ::-1]
    return tri


def _has_directed_edge(t: np.ndarray, x: int, y: int) -> bool:
    return (t[0] == x and t[1] == y) or (t[1] == x and t[2] == y) or (t[2] == x and t[0] == y)


def _signed_volume(v: np.ndarray, tri: np.ndarray) -> float:
    a, b, c = v[tri[:, 0]], v[tri[:, 1]], v[tri[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6)
