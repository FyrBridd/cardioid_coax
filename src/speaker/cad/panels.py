"""Parametric cabinet panels, the coax flare insert and the full assembly.

Units: everything here is millimetres; SI inputs are converted with ``_mm``.

Global frame (mm) is the Layout frame: origin at the bottom-left-front corner,
x across the width, y backwards (front face at y=0), z up.

Panel-local frame: each panel is a W x H x t slab occupying
[0, W] x [0, H] x [0, t]. Feature coordinates are 2D (u, v) in that frame,
origin at the panel's lower-left corner, *as seen looking at the local +Z face*
(the face at local z = t). That is the face the panel is machined from: pockets
are cut from it to ``depth_mm``. ``Panel.placement`` maps the local frame onto
the global frame (proper rotation + translation, no mirroring).

Joinery: baffle_outer is full width and height; the sides are full height and
run from behind the outer baffle to the rear; everything else fits between the
sides.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import cadquery as cq

from ..drivers import Driver
from ..layout import Layout
from ..params import Params

BOLT_HOLE_D = 6.0  # mm
POCKET_CLEARANCE = 2.0  # mm added to driver outer_d for the flange pocket
VENT_CORNER_R = 3.0  # mm, tool radius for vent slots
TERMINAL_W, TERMINAL_H, TERMINAL_CORNER_R = 100.0, 160.0, 3.0  # mm
FLARE_FIT_CLEARANCE = 0.2  # mm on diameter, push fit
FLARE_LIP = 3.0  # mm radial wall at the insert's front face
COAX_BOLTS, BASS_BOLTS = 8, 6


def _mm(x: float) -> float:
    return x * 1000.0


@dataclass(frozen=True)
class Feature:
    """One cut, in panel-local mm (see module docstring).

    circle: ``w_mm == h_mm`` is the diameter; ``corner_r_mm`` is unused (0).
    rect / slot: rounded rectangle ``w_mm`` x ``h_mm`` centred on (cx, cy)
    with corner radius ``corner_r_mm``; "slot" marks vent slots.
    """

    kind: Literal["through", "pocket"]
    depth_mm: float | None  # None for through-cuts
    shape: Literal["circle", "rect", "slot"]
    cx_mm: float
    cy_mm: float
    w_mm: float
    h_mm: float
    corner_r_mm: float = 0.0


@dataclass(frozen=True)
class Panel:
    name: str
    outline_w_mm: float
    outline_h_mm: float
    thickness_mm: float
    features: tuple[Feature, ...]
    placement: cq.Location


def _placement(origin: tuple[float, float, float], u: tuple, v: tuple) -> cq.Location:
    """Location mapping local x->u, y->v, z->u x v, with local (0,0,0) at origin."""
    normal = cq.Vector(*u).cross(cq.Vector(*v))
    return cq.Location(cq.Plane(origin=origin, xDir=u, normal=normal.toTuple()))


def _circle(kind: str, cx: float, cy: float, d: float, depth: float | None = None) -> Feature:
    return Feature(kind, depth, "circle", cx, cy, d, d)


def _driver_cuts(cx: float, cy: float, d: Driver, n_bolts: int) -> list[Feature]:
    """Rear-mount cuts on the inner baffle: bore, flange pocket, bolt circle."""
    feats = [
        _circle("through", cx, cy, _mm(d.cutout_d)),
        _circle("pocket", cx, cy, _mm(d.outer_d) + POCKET_CLEARANCE, _mm(d.flange_t)),
    ]
    r = _mm(d.outer_d + d.cutout_d) / 4  # bolt-circle diameter is (outer + cutout)/2
    for i in range(n_bolts):
        a = 2 * math.pi * (i + 0.5) / n_bolts
        feats.append(_circle("through", cx + r * math.cos(a), cy + r * math.sin(a), BOLT_HOLE_D))
    return feats


def build_panels(l: Layout, p: Params, coax: Driver, bass: Driver) -> list[Panel]:
    c = p.cabinet
    W, D, H = _mm(c.width), _mm(c.depth), _mm(l.height)
    t, bt = _mm(c.panel_t), _mm(c.baffle_t)
    iw = W - 2 * t  # between the sides
    ih = H - 2 * t  # between top and bottom
    cx_g, cz_g = (_mm(v) for v in l.coax_centre)
    bass_g = [(_mm(x), _mm(z)) for x, z in l.bass_centres]
    X, Y, Z = (1, 0, 0), (0, 1, 0), (0, 0, 1)
    NX = (-1, 0, 0)

    # baffle_outer: viewed from the front (local +Z = -y); local u = x, v = z.
    outer_feats = [_circle("through", cx_g, cz_g, _mm(c.flare_mouth_d))]
    outer_feats += [_circle("through", x, z, _mm(bass.cutout_d)) for x, z in bass_g]

    # baffle_inner: machined from its rear face (local +Z = +y), so local
    # u = -x (u = 302 - x), v = z - t.
    def inner_uv(x: float, z: float) -> tuple[float, float]:
        return W - t - x, z - t

    inner_feats = _driver_cuts(*inner_uv(cx_g, cz_g), coax, COAX_BOLTS)
    for x, z in bass_g:
        inner_feats += _driver_cuts(*inner_uv(x, z), bass, BASS_BOLTS)

    # sides start behind the outer baffle layer: u = y - t, v = z, local +Z = +x.
    side_feats = tuple(
        Feature("through", None, "slot", _mm(y) - t, _mm(z), _mm(w), _mm(h), VENT_CORNER_R)
        for y, z, w, h in l.vent_slots
    )

    # back: viewed from outside (local +Z = +y); u = 302 - x, v = z - t.
    bass_mid_z = (t + _mm(l.divider_z)) / 2  # mid-height of the bass section
    back_feats = (
        Feature("through", None, "rect", iw / 2, bass_mid_z - t,
                TERMINAL_W, TERMINAL_H, TERMINAL_CORNER_R),
    )

    return [
        Panel("baffle_outer", W, H, t, tuple(outer_feats), _placement((0, t, 0), X, Z)),
        Panel("baffle_inner", iw, ih, t, tuple(inner_feats), _placement((W - t, t, t), NX, Z)),
        Panel("side_left", D - t, H, t, side_feats, _placement((0, t, 0), Y, Z)),
        Panel("side_right", D - t, H, t, side_feats, _placement((W - t, t, 0), Y, Z)),
        Panel("top", iw, D - t, t, (), _placement((t, t, H - t), X, Y)),
        Panel("bottom", iw, D - t, t, (), _placement((t, t, 0), X, Y)),
        Panel("back", iw, ih, t, back_feats, _placement((W - t, D - t, t), NX, Z)),
        Panel("divider", iw, D - bt - t, t, (), _placement((t, bt, _mm(l.divider_z)), X, Y)),
    ]


def _cutter(f: Feature, thickness: float) -> cq.Solid:
    if f.kind == "through":
        z0, h = -1.0, thickness + 2.0
    else:
        z0, h = thickness - f.depth_mm, f.depth_mm + 1.0
    wp = cq.Workplane("XY", origin=(f.cx_mm, f.cy_mm, z0))
    if f.shape == "circle":
        return wp.circle(f.w_mm / 2).extrude(h).val()
    sk = cq.Sketch().rect(f.w_mm, f.h_mm).vertices().fillet(f.corner_r_mm)
    return wp.placeSketch(sk).extrude(h).val()


def panel_solid(panel: Panel) -> cq.Workplane:
    """The panel in its local frame (unplaced), with every feature cut."""
    slab = cq.Workplane("XY").box(
        panel.outline_w_mm, panel.outline_h_mm, panel.thickness_mm, centered=False
    )
    if not panel.features:
        return slab
    cutters = cq.Compound.makeCompound([_cutter(f, panel.thickness_mm) for f in panel.features])
    return slab.cut(cutters)


def flare_insert(l: Layout, p: Params, coax: Driver) -> cq.Workplane:
    """3D-printed flare for the outer baffle layer, axis along local +z.

    Rear face (throat, Ø coax.cutout_d) at z=0, front face (mouth) at
    z=panel_t. The bore is a quarter-ellipse, tangent to the axis at the throat
    and to the front face at the mouth, ending FLARE_LIP inside the outer wall.
    """
    length = _mm(p.cabinet.panel_t)
    r_out = (_mm(p.cabinet.flare_mouth_d) - FLARE_FIT_CLEARANCE) / 2
    r_throat = _mm(coax.cutout_d) / 2
    r_mouth = r_out - FLARE_LIP
    if r_mouth <= r_throat:
        raise ValueError(
            f"flare mouth Ø{2 * r_mouth:.1f} mm is not larger than the throat Ø{2 * r_throat:.1f} mm"
        )
    # Profile in the XZ half-plane (x = radius, z = axial), revolved about z.
    arc = cq.Edge.makeEllipse(
        r_mouth - r_throat, length, pnt=(r_mouth, 0, 0), dir=(0, -1, 0), xdir=(1, 0, 0),
        angle1=90, angle2=180,
    )  # from (r_mouth, 0, length) to (r_throat, 0, 0)
    pts = [(r_throat, 0, 0), (r_out, 0, 0), (r_out, 0, length), (r_mouth, 0, length)]
    lines = [cq.Edge.makeLine(cq.Vector(*a), cq.Vector(*b)) for a, b in zip(pts, pts[1:])]
    face = cq.Face.makeFromWires(cq.Wire.assembleEdges([*lines, arc]))
    solid = cq.Solid.revolve(face, 360, cq.Vector(0, 0, 0), cq.Vector(0, 0, 1))
    return cq.Workplane("XY").add(solid)


def assembly(l: Layout, p: Params, coax: Driver, bass: Driver) -> cq.Assembly:
    asm = cq.Assembly(name="cabinet")
    for pn in build_panels(l, p, coax, bass):
        asm.add(panel_solid(pn), name=pn.name, loc=pn.placement)
    # Insert sits in the outer baffle hole: local +z (front) -> global -y,
    # rear face on the inner baffle's front face (y = panel_t).
    cx, cz = (_mm(v) for v in l.coax_centre)
    insert_loc = _placement((cx, _mm(p.cabinet.panel_t), cz), (1, 0, 0), (0, 0, 1))
    asm.add(flare_insert(l, p, coax), name="flare_insert", loc=insert_loc)
    return asm
