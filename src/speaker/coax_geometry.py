"""Coax flare and cavity geometry shared by the CAD and the BEM mesher. SI units (m).

Layout frame: y runs backwards from the front face at y=0. Profiles are (r, y)
polylines about the coax axis.
"""
from __future__ import annotations

import math

from speaker.drivers import Driver
from speaker.params import Params

FLARE_FIT_CLEARANCE = 0.2e-3  # m on diameter, push fit
FLARE_LIP = 3.0e-3  # m radial wall at the insert's front face


def flare_radii(p: Params, coax: Driver) -> tuple[float, float, float]:
    """(r_throat, r_mouth, r_out) of the flare insert."""
    r_out = (p.cabinet.flare_mouth_d - FLARE_FIT_CLEARANCE) / 2
    return coax.cutout_d / 2, r_out - FLARE_LIP, r_out


def cavity_profile(
    p: Params, coax: Driver, n_arc: int = 16
) -> tuple[tuple[str, tuple[tuple[float, float], ...]], ...]:
    """Labelled (r, y) polylines from the flare mouth inwards to the axis."""
    missing = [
        n for n in ("diaphragm_d", "cone_depth", "cone_inner_d", "hf_exit_d")
        if getattr(coax, n) is None
    ]
    if missing:
        raise ValueError(f"coax {coax.name!r} lacks cone geometry fields: {', '.join(missing)}")
    r_throat, r_mouth, _ = flare_radii(p, coax)
    if r_mouth <= r_throat:
        raise ValueError("flare mouth is not larger than the throat")
    t, bt = p.cabinet.panel_t, p.cabinet.baffle_t
    a = r_mouth - r_throat
    # Quarter-ellipse centred (r_mouth, t): tangent to the front face at the mouth, to the axis at the throat.
    flare = tuple(
        (r_mouth - a * math.sin(th), t - t * math.cos(th))
        for th in (math.pi / 2 * i / n_arc for i in range(n_arc + 1))
    )
    r_dia, r_ci, r_hf = coax.diaphragm_d / 2, coax.cone_inner_d / 2, coax.hf_exit_d / 2
    y_cone = bt + coax.cone_depth
    return (
        ("flare", flare),
        ("bore", ((r_throat, t), (r_throat, bt))),
        ("surround", ((r_throat, bt), (r_dia, bt))),
        ("cone", ((r_dia, bt), (r_ci, y_cone))),
        ("pole", ((r_ci, y_cone), (r_hf, y_cone))),
        ("hf_exit", ((r_hf, y_cone), (0.0, y_cone))),
    )
