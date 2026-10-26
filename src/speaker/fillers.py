"""Printed fillers that bring the coax chamber down to its target net volume.

Four corner pieces sit outside the coax frame's sweep (so the coax still goes in
from the back) and stop short of the vent slots; one rear block behind the
magnet, on the back panel, makes up the rest. The corner profiles are polygons
so the CAD extrudes exactly the shape whose volume is reported here.

Frame as in layout: origin at bottom-left-front corner; x across the width,
y backwards, z up. All SI.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from .drivers import Driver
from .params import Params

CORNER_CLEARANCE = 0.003   # m, radial, beyond coax.outer_d / 2
CORNER_TRIM = 0.040        # m, pieces keep |dx| >= and |dz| >= this from the coax centre lines
CORNER_MIN_T = 0.005       # m, thinnest allowed wall between arc end and chamber wall
VENT_GAP = 0.015           # m, corner-piece rear face to first vent slot edge
MAGNET_GAP = 0.010         # m, coax rear to rear-block front face
SIDE_CHANNEL = 0.025       # m, rear block to each side panel
ARC_SEGMENTS = 64          # straight segments per corner arc
VOLUME_TOL = 0.05e-3       # m^3


@dataclass(frozen=True)
class Fillers:
    # (name, closed polygon of (x, z) vertices, global frame, first vertex not repeated)
    corners: tuple[tuple[str, tuple[tuple[float, float], ...]], ...]
    corner_y: tuple[float, float]  # (y_front, y_back) of all four pieces
    corner_volume: float  # m^3, all four
    block_min: tuple[float, float, float]  # (x, y, z)
    block_max: tuple[float, float, float]
    block_volume: float  # m^3

    @property
    def total_volume(self) -> float:
        return self.corner_volume + self.block_volume


def _mm(x: float) -> str:
    return f"{x * 1000:.1f} mm"


def polygon_area(pts: Sequence[tuple[float, float]]) -> float:
    n = len(pts)
    s = sum(pts[i][0] * pts[(i + 1) % n][1] - pts[(i + 1) % n][0] * pts[i][1] for i in range(n))
    return abs(s) / 2


def corner_profile(cx: float, cz: float, x_wall: float, z_wall: float,
                   r: float, trim: float) -> tuple[tuple[float, float], ...]:
    """One corner piece's cross-section: the rectangle between the trim lines and
    the two walls, less the disc of radius r about the coax centre."""
    sx, sz = math.copysign(1.0, x_wall - cx), math.copysign(1.0, z_wall - cz)
    X, Z = abs(x_wall - cx), abs(z_wall - cz)
    a = math.sqrt(r**2 - trim**2)
    if a > min(X, Z) - CORNER_MIN_T:
        raise ValueError(
            f"corner filler too thin: the coax sweep (radius {_mm(r)}) reaches {_mm(a)} from the "
            f"coax centre along the trim line, but the chamber walls are {_mm(X)} (side) and "
            f"{_mm(Z)} (top/bottom) away; at least {_mm(CORNER_MIN_T)} of wall is needed"
        )
    local = [(a, trim), (X, trim), (X, Z), (trim, Z), (trim, a)]
    t0, t1 = math.atan2(a, trim), math.atan2(trim, a)  # arc from (trim, a) to (a, trim)
    for k in range(1, ARC_SEGMENTS):
        t = t0 + (t1 - t0) * k / ARC_SEGMENTS
        local.append((r * math.cos(t), r * math.sin(t)))
    return tuple((cx + sx * u, cz + sz * v) for u, v in local)


def compute_fillers(p: Params, coax: Driver, *, coax_centre: tuple[float, float],
                    chamber_z: tuple[float, float], volume_to_fill: float) -> Fillers:
    c, v = p.cabinet, p.vents
    cx, cz = coax_centre
    z_lo, z_hi = chamber_z
    x_lo, x_hi = c.panel_t, c.width - c.panel_t
    r = coax.outer_d / 2 + CORNER_CLEARANCE

    y_front = c.baffle_t
    vent_y_start = v.centre_from_front - v.slot_w / 2
    y_back = vent_y_start - VENT_GAP
    if y_back <= y_front:
        raise ValueError(
            f"vent slots start at y={_mm(vent_y_start)}, so corner fillers would end at "
            f"y={_mm(y_back)} ({_mm(VENT_GAP)} before the vents), not behind the baffle rear "
            f"face at y={_mm(y_front)}; move the vents back"
        )
    walls = (("top_left", x_lo, z_hi), ("top_right", x_hi, z_hi),
             ("bottom_left", x_lo, z_lo), ("bottom_right", x_hi, z_lo))
    corners = tuple(
        (f"filler_corner_{n}", corner_profile(cx, cz, xw, zw, r, CORNER_TRIM)) for n, xw, zw in walls
    )
    corner_volume = sum(polygon_area(pts) for _, pts in corners) * (y_back - y_front)

    remaining = volume_to_fill - corner_volume
    if remaining < 0:
        raise ValueError(
            f"corner fillers alone remove {corner_volume * 1000:.3f} L, more than the "
            f"{volume_to_fill * 1000:.3f} L to fill (by {-remaining * 1000:.3f} L); "
            f"lower the target net volume"
        )
    bx0, bx1 = x_lo + SIDE_CHANNEL, x_hi - SIDE_CHANNEL
    coax_rear = c.baffle_t + coax.depth - coax.flange_t
    by0, by1 = max(y_back, coax_rear + MAGNET_GAP), c.depth - c.panel_t
    if bx1 <= bx0 or by1 <= by0:
        raise ValueError(
            f"rear block has no room: x {_mm(bx0)} to {_mm(bx1)}, y {_mm(by0)} to {_mm(by1)}"
        )
    h = remaining / ((bx1 - bx0) * (by1 - by0))
    if h > z_hi - z_lo:
        raise ValueError(
            f"rear block would need to be {_mm(h)} tall to fill {remaining * 1000:.2f} L, but the "
            f"coax chamber is {_mm(z_hi - z_lo)} tall; raise the target net volume"
        )
    z_mid = (z_lo + z_hi) / 2
    return Fillers(
        corners=corners,
        corner_y=(y_front, y_back),
        corner_volume=corner_volume,
        block_min=(bx0, by0, z_mid - h / 2),
        block_max=(bx1, by1, z_mid + h / 2),
        block_volume=remaining,
    )
