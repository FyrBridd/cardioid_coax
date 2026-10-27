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


Profile = tuple[tuple[str, tuple[tuple[float, float], ...]], ...]


def cavity_profile(p: Params, coax: Driver, n_arc: int = 16) -> Profile:
    """Labelled (r, y) polylines from the flare mouth inwards to the axis.

    The flare runs through the outer baffle layer (y from 0 to outer_baffle_t), the bore at the
    throat radius through the inner layer (to baffle_t). flare_profile "ellipse" is a quarter
    ellipse of n_arc segments, "cone" one straight segment from the mouth to the throat.
    """
    _require_cone(coax)
    r_throat, r_mouth, _ = flare_radii(p, coax)
    if r_mouth <= r_throat:
        raise ValueError("flare mouth is not larger than the throat")
    t, bt = p.cabinet.outer_baffle_t, p.cabinet.baffle_t
    a = r_mouth - r_throat
    if p.cabinet.flare_profile == "ellipse":
        # Quarter-ellipse centred (r_mouth, t): tangent to the front face at the mouth, to the axis at the throat.
        flare = tuple(
            (r_mouth - a * math.sin(th), t - t * math.cos(th))
            for th in (math.pi / 2 * i / n_arc for i in range(n_arc + 1))
        )
    elif p.cabinet.flare_profile == "cone":
        flare = ((r_mouth, 0.0), (r_throat, t))
    else:
        raise ValueError(f"flare_profile {p.cabinet.flare_profile!r} is not 'ellipse' or 'cone'")
    return (
        ("flare", flare),
        ("bore", ((r_throat, t), (r_throat, bt))),
        *_cone_segments(coax, r_throat, bt),
    )


def recess_profile(coax: Driver) -> Profile:
    """The coax front-mounted on a flat baffle: its frame opening, then the cone in its own recess.

    The opening is cutout_d across and flange_t deep (the surround sits flange_t behind the
    mounting face at y=0); the surround, cone, pole and HF exit are cavity_profile's.
    """
    _require_cone(coax)
    r = coax.cutout_d / 2
    return (("bore", ((r, 0.0), (r, coax.flange_t))), *_cone_segments(coax, r, coax.flange_t))


def _require_cone(coax: Driver) -> None:
    missing = [
        n for n in ("diaphragm_d", "cone_depth", "cone_inner_d", "hf_exit_d")
        if getattr(coax, n) is None
    ]
    if missing:
        raise ValueError(f"coax {coax.name!r} lacks cone geometry fields: {', '.join(missing)}")


def _cone_segments(coax: Driver, r_open: float, y_surround: float) -> Profile:
    """surround (from the opening radius r_open in to the diaphragm), cone, pole and hf_exit, with
    the surround at depth y_surround."""
    r_dia, r_ci, r_hf = coax.diaphragm_d / 2, coax.cone_inner_d / 2, coax.hf_exit_d / 2
    y_cone = y_surround + coax.cone_depth
    return (
        ("surround", ((r_open, y_surround), (r_dia, y_surround))),
        ("cone", ((r_dia, y_surround), (r_ci, y_cone))),
        ("pole", ((r_ci, y_cone), (r_hf, y_cone))),
        ("hf_exit", ((r_hf, y_cone), (0.0, y_cone))),
    )
