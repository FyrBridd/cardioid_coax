"""3D-printed coax chamber fillers as solids, from the polygons in ``speaker.fillers``.

Solids are in the global frame in millimetres, already placed (no Location).
"""
from __future__ import annotations

import cadquery as cq

from ..layout import Layout


def _mm(x: float) -> float:
    return x * 1000.0


def filler_solids(l: Layout) -> dict[str, cq.Workplane]:
    f = l.fillers
    y0, y1 = (_mm(y) for y in f.corner_y)
    # Plane with normal +y and xDir +x: its local y axis is global -z.
    plane = cq.Plane(origin=(0, y0, 0), xDir=(1, 0, 0), normal=(0, 1, 0))
    out = {
        name: cq.Workplane(plane).polyline([(_mm(x), -_mm(z)) for x, z in pts]).close().extrude(y1 - y0)
        for name, pts in f.corners
    }
    lo, hi = cq.Vector(*(_mm(v) for v in f.block_min)), cq.Vector(*(_mm(v) for v in f.block_max))
    out["filler_rear_block"] = cq.Workplane("XY").add(
        cq.Solid.makeBox(*(hi - lo).toTuple(), pnt=lo)
    )
    return out
