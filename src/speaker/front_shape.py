"""The front-edge shape families and the rules that keep them buildable (spec §5).

Frame as in layout: x across the width, y backwards from the front face (y = 0), z up. The
treatment cuts the two vertical front edges, mirrored about x = width / 2; the back edges and
the top and bottom edges stay square. All SI.
"""
from __future__ import annotations

import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .drivers import Driver
    from .layout import Layout
    from .params import Cabinet, Params

FRONT_LAND = 0.005  # m of flat front kept beside every cutout
CHAMFER_MAX = 0.045  # m, the largest 45° chamfer leg built from edge strips
FAMILIES = ("square", "roundover", "chamfer", "facets")
FLARE_PROFILES = ("ellipse", "cone")
BASS_NET_TARGET = 0.050  # m^3, 25 L per driver (spec §5)
_EPS = 1e-9


def _mm(x: float) -> str:
    return f"{x * 1000:.1f} mm"


def _deg(a: float) -> str:
    return f"{math.degrees(a):.1f}°"


def cut_depth(c: Cabinet) -> float:
    """How far back from the front face the treatment reaches on the side face."""
    if c.front_edge in ("roundover", "chamfer"):
        return c.edge_size
    if c.front_edge == "facets":
        return c.facet_width * math.tan(c.facet_angle)
    return 0.0


def _facet_fraction(c: Cabinet, z: float) -> float:
    """0 below the facet base, 1 at and above base + ramp, linear across the ramp."""
    if z < c.facet_base_z:
        return 0.0
    if c.facet_ramp <= 0 or z >= c.facet_base_z + c.facet_ramp:
        return 1.0
    return (z - c.facet_base_z) / c.facet_ramp


def _legs(c: Cabinet, z: float) -> tuple[float, float]:
    """(across the front, back along the side) of one corner's cut at height z."""
    if c.front_edge in ("roundover", "chamfer"):
        return c.edge_size, c.edge_size
    if c.front_edge == "facets":
        u = _facet_fraction(c, z)
        return u * c.facet_width, u * cut_depth(c)
    return 0.0, 0.0


def flat_width(c: Cabinet, z: float) -> float:
    """Width of the flat front face at height z."""
    return c.width - 2 * _legs(c, z)[0]


def section(c: Cabinet, z: float, n_arc: int = 8) -> list[tuple[float, float]]:
    """The outer horizontal cross-section at height z, a counter-clockwise (x, y) polygon."""
    W, D = c.width, c.depth
    a, b = _legs(c, z)
    if a <= 0 or b <= 0:
        return [(0.0, 0.0), (W, 0.0), (W, D), (0.0, D)]
    if c.front_edge == "roundover":
        r = a
        # Right corner, centre (W - r, r): from (W - r, 0) round to (W, r).
        right = [(W - r + r * math.cos(t), r + r * math.sin(t))
                 for t in (-math.pi / 2 + (math.pi / 2) * k / n_arc for k in range(n_arc + 1))]
        # Left corner, centre (r, r): from (0, r) round to (r, 0).
        left = [(r + r * math.cos(t), r + r * math.sin(t))
                for t in (math.pi + (math.pi / 2) * k / n_arc for k in range(n_arc + 1))]
    else:  # chamfer and facets: a straight cut from the front face to the side face
        right = [(W - a, 0.0), (W, b)]
        left = [(0.0, b), (a, 0.0)]
    # Front edge left to right, the right corner, the back, the left corner back to the start.
    return [left[-1], *right, (W, D), (0.0, D), *left[:-1]]


def _inside(pts: list[tuple[float, float]], x: float, y: float) -> bool:
    """Point in a convex counter-clockwise polygon, boundary included."""
    n = len(pts)
    for i in range(n):
        (x0, y0), (x1, y1) = pts[i], pts[(i + 1) % n]
        if (x1 - x0) * (y - y0) - (y1 - y0) * (x - x0) < -1e-12:
            return False
    return True


def validate_front(p: Params, coax: Driver, bass: Driver, l: Layout) -> None:
    """Raise ValueError if the front shape, coax offset, outer baffle or flare cannot be built."""
    c, v = p.cabinet, p.vents
    fe = c.front_edge

    # (a) family
    if fe not in FAMILIES:
        raise ValueError(f"front_edge {fe!r} is not one of {FAMILIES}")

    # (b) edge size
    if fe in ("square", "facets") and c.edge_size != 0:
        raise ValueError(f"edge_size {_mm(c.edge_size)} must be 0 for {fe} (it sizes a roundover or chamfer)")
    if fe in ("roundover", "chamfer") and c.edge_size <= 0:
        raise ValueError(f"edge_size {_mm(c.edge_size)} must be > 0 for {fe}")
    if fe == "chamfer" and c.edge_size > CHAMFER_MAX + _EPS:
        raise ValueError(f"chamfer edge_size {_mm(c.edge_size)} exceeds CHAMFER_MAX {_mm(CHAMFER_MAX)}")

    # (c) facet values
    if fe == "facets":
        if not 0 < c.facet_angle <= math.radians(45) + _EPS:
            raise ValueError(f"facet_angle {_deg(c.facet_angle)} must be > 0 and ≤ 45°")
        if c.facet_width <= 0:
            raise ValueError(f"facet_width {_mm(c.facet_width)} must be > 0")
        if c.facet_ramp < 0:
            raise ValueError(f"facet_ramp {_mm(c.facet_ramp)} must be ≥ 0")
    elif (c.facet_angle, c.facet_width, c.facet_base_z, c.facet_ramp) != (0, 0, 0, 0):
        raise ValueError(
            f"front_edge {fe!r} takes no facets: the facet values must all be 0 (got facet_angle "
            f"{_deg(c.facet_angle)}, facet_width {_mm(c.facet_width)}, facet_base_z {_mm(c.facet_base_z)}, "
            f"facet_ramp {_mm(c.facet_ramp)})"
        )

    # (d) the treatment stays outside the interior, at its deepest height
    if fe != "square":
        z_full = c.facet_base_z + c.facet_ramp if fe == "facets" else 0.0
        pts = section(c, z_full)
        for x in (c.panel_t, c.width - c.panel_t):
            if not _inside(pts, x, c.baffle_t):
                raise ValueError(
                    f"the {fe} treatment reaches the interior: its front corner at x={_mm(x)}, "
                    f"y={_mm(c.baffle_t)} (panel_t, baffle_t) is outside the outer section "
                    f"(cut {_mm(_legs(c, z_full)[0])} across the front, {_mm(cut_depth(c))} back along the side)"
                )

    # (e) flat front around the flare mouth, with the offset
    off = abs(c.coax_offset_x)
    fw = flat_width(c, c.coax_centre_z)
    need = c.flare_mouth_d + 2 * off + 2 * FRONT_LAND
    if fw < need - _EPS:
        raise ValueError(
            f"flat front width {_mm(fw)} at the coax centre height {_mm(c.coax_centre_z)} is below "
            f"flare_mouth_d {_mm(c.flare_mouth_d)} + 2 x |coax_offset_x| {_mm(off)} + 2 x FRONT_LAND "
            f"{_mm(FRONT_LAND)} = {_mm(need)}"
        )

    # (f) flat front around each bass cutout
    need = bass.cutout_d + 2 * FRONT_LAND
    for _, z in l.bass_centres:
        fw = flat_width(c, z)
        if fw < need - _EPS:
            raise ValueError(
                f"flat front width {_mm(fw)} at bass centre height {_mm(z)} is below bass cutout_d "
                f"{_mm(bass.cutout_d)} + 2 x FRONT_LAND {_mm(FRONT_LAND)} = {_mm(need)}"
            )

    # (g) the coax frame clears the side panels
    if off + coax.outer_d / 2 > l.interior_w / 2 + _EPS:
        raise ValueError(
            f"coax_offset_x {_mm(c.coax_offset_x)} puts the coax frame into the side panels: "
            f"|offset| + coax outer_d / 2 = {_mm(off + coax.outer_d / 2)} exceeds interior_w / 2 "
            f"{_mm(l.interior_w / 2)}"
        )

    if fe == "facets":
        # (h) the facets start above the upper bass cutout and are full before the flare mouth
        cutout_top = max(z for _, z in l.bass_centres) + bass.cutout_d / 2
        lo = cutout_top + FRONT_LAND
        if c.facet_base_z < lo - _EPS:
            raise ValueError(
                f"facet_base_z {_mm(c.facet_base_z)} is below the top of the upper bass cutout "
                f"({_mm(cutout_top)}) + FRONT_LAND {_mm(FRONT_LAND)} = {_mm(lo)}: the bass section stays flat"
            )
        hi = c.coax_centre_z - c.flare_mouth_d / 2
        top = c.facet_base_z + c.facet_ramp
        if top > hi + _EPS:
            raise ValueError(
                f"facet_base_z + facet_ramp {_mm(top)} is above the bottom of the flare mouth, "
                f"coax_centre_z - flare_mouth_d / 2 = {_mm(hi)}"
            )
        # (i) the vent slots start behind the facet cut
        start = v.centre_from_front - v.slot_w / 2
        need = cut_depth(c) + FRONT_LAND
        if start < need - _EPS:
            raise ValueError(
                f"vent slots start at y={_mm(start)}, inside the facet cut: they must start at or behind "
                f"cut_depth {_mm(cut_depth(c))} + FRONT_LAND {_mm(FRONT_LAND)} = {_mm(need)}"
            )

    # (j) flare profile and outer baffle laminate
    if c.flare_profile not in FLARE_PROFILES:
        raise ValueError(f"flare_profile {c.flare_profile!r} is not one of {FLARE_PROFILES}")
    if not c.panel_t - _EPS <= c.outer_baffle_t <= 3 * c.panel_t + _EPS:
        raise ValueError(
            f"outer_baffle_t {_mm(c.outer_baffle_t)} must be between panel_t {_mm(c.panel_t)} and "
            f"3 x panel_t {_mm(3 * c.panel_t)} (a laminate of standard sheets)"
        )

    # (k) bass volume
    if l.bass_net < BASS_NET_TARGET:
        raise ValueError(
            f"bass net volume {l.bass_net * 1000:.1f} L is below the "
            f"{BASS_NET_TARGET * 1000:.0f} L target (25 L per driver)"
        )
