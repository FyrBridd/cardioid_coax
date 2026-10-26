"""Cabinet layout, internal volumes and fit validation.

Frame: origin at bottom-left-front corner; x across width, y backwards, z up.
Front face at y=0. Interior spans y in [baffle_t, depth - panel_t]; floor
occupies z in [0, panel_t], top panel z in [height - panel_t, height].
All SI.
"""
from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass

from .cardioid import CardioidParams
from .drivers import Driver
from .fillers import Fillers, compute_fillers
from .params import Params

DIVIDER_CLEARANCE = 0.010  # m between coax frame bottom and divider top
BASS_NET_TARGET = 0.050  # m^3, soft target (warning)
BASS_NET_MIN = 0.030  # m^3, hard minimum (ValueError)
EAR_HEIGHT_TOLERANCE = 0.020  # m, coax centre vs seated ear height (warning)


@dataclass(frozen=True)
class Layout:
    height: float
    interior_w: float
    interior_d: float
    coax_centre: tuple[float, float]  # (x, z)
    bass_centres: tuple[tuple[float, float], tuple[float, float]]
    divider_z: float  # bottom face of the divider
    coax_chamber_gross: float
    coax_chamber_net: float
    filler_volume: float  # printed fillers (in the CAD), all five parts
    coax_chamber_empty: float  # gross less coax displacement, i.e. without fillers
    coax_chamber_as_built: float  # empty less fillers; fed to the cardioid model
    bass_gross: float
    bass_net: float
    vent_slots: tuple[tuple[float, float, float, float], ...]  # (y, z, w, h)
    warnings: tuple[str, ...] = ()
    fillers: Fillers | None = None  # None only inside compute_layout, before validation


def _displacement(d: Driver) -> float:
    return math.pi / 4 * d.basket_d**2 * (d.depth - d.flange_t)


def compute_layout(p: Params, coax: Driver, bass: Driver) -> Layout:
    c, v = p.cabinet, p.vents
    height = c.coax_centre_z + coax.outer_d / 2 + c.top_margin + c.panel_t
    interior_w = c.width - 2 * c.panel_t
    interior_d = c.depth - c.baffle_t - c.panel_t
    divider_top = c.coax_centre_z - coax.outer_d / 2 - DIVIDER_CLEARANCE
    divider_z = divider_top - c.divider_t
    top_underside = height - c.panel_t

    coax_h = top_underside - divider_top
    coax_gross = interior_w * interior_d * coax_h
    bass_gross = interior_w * interior_d * (divider_z - c.panel_t)
    bass_net = bass_gross - 2 * _displacement(bass)
    coax_net = c.coax_chamber_net_volume
    coax_empty = coax_gross - _displacement(coax)

    # Slots: one side's stack, centred vertically on the coax chamber.
    z_mid = (divider_top + top_underside) / 2
    n = v.slots_per_side
    vent_slots = tuple(
        (v.centre_from_front, z_mid + (i - (n - 1) / 2) * v.slot_pitch, v.slot_w, v.slot_h)
        for i in range(n)
    )

    layout = Layout(
        height=height,
        interior_w=interior_w,
        interior_d=interior_d,
        coax_centre=(c.width / 2, c.coax_centre_z),
        bass_centres=((c.width / 2, c.bass_centres_z[0]), (c.width / 2, c.bass_centres_z[1])),
        divider_z=divider_z,
        coax_chamber_gross=coax_gross,
        coax_chamber_net=coax_net,
        filler_volume=0.0,
        coax_chamber_empty=coax_empty,
        coax_chamber_as_built=coax_empty,
        bass_gross=bass_gross,
        bass_net=bass_net,
        vent_slots=vent_slots,
    )
    validate(layout, p, coax, bass)
    fillers = compute_fillers(
        p, coax, coax_centre=layout.coax_centre, chamber_z=(divider_top, top_underside),
        volume_to_fill=coax_empty - coax_net,
    )
    warns: list[str] = []
    if bass_net < BASS_NET_TARGET:
        warns.append(
            f"bass net volume {bass_net * 1000:.1f} L is below the "
            f"{BASS_NET_TARGET * 1000:.0f} L target (25 L per driver)"
        )
    ear = p.room.seat[2]
    if abs(c.coax_centre_z - ear) > EAR_HEIGHT_TOLERANCE:
        warns.append(
            f"coax centre {_mm(c.coax_centre_z)} differs from the seat ear height {_mm(ear)} "
            f"by more than {_mm(EAR_HEIGHT_TOLERANCE)}"
        )
    return dataclasses.replace(
        layout,
        fillers=fillers,
        filler_volume=fillers.total_volume,
        coax_chamber_as_built=coax_empty - fillers.total_volume,
        warnings=tuple(warns),
    )


def cardioid_params(p: Params, l: Layout, coax: Driver) -> CardioidParams:
    c, v = p.cabinet, p.vents
    zs = [s[1] for s in l.vent_slots]
    z_rel = (min(zs) + max(zs)) / 2 - c.coax_centre_z
    x = -(v.centre_from_front - c.baffle_t)
    half = c.width / 2
    return CardioidParams(
        chamber_volume=l.coax_chamber_as_built,
        vent_area=2 * v.slots_per_side * v.slot_w * v.slot_h,
        vent_length=c.panel_t,
        vent_resistance=v.specific_resistance,
        vent_positions=((x, half, z_rel), (x, -half, z_rel)),
        n_slots=2 * v.slots_per_side,
    )


def _mm(x: float) -> str:
    return f"{x * 1000:.1f} mm"


def validate(l: Layout, p: Params, coax: Driver, bass: Driver) -> None:
    c = p.cabinet
    if abs(c.baffle_t - 2 * c.panel_t) > 1e-9:
        raise ValueError(
            f"baffle_t {_mm(c.baffle_t)} must equal two panel layers (2 x panel_t {_mm(c.panel_t)} = "
            f"{_mm(2 * c.panel_t)}): the baffle is built as an outer and an inner panel_t layer"
        )
    if c.top_margin < 0:
        raise ValueError(
            f"top_margin {_mm(c.top_margin)} is negative: coax frame top at "
            f"{_mm(c.coax_centre_z + coax.outer_d / 2)} would cut into the top panel "
            f"(underside at {_mm(l.height - c.panel_t)})"
        )
    if p.vents.slots_per_side < 1:
        raise ValueError(f"slots_per_side must be at least 1, got {p.vents.slots_per_side}")
    if coax.outer_d > l.interior_w:
        raise ValueError(
            f"coax driver outer diameter {_mm(coax.outer_d)} exceeds interior width {_mm(l.interior_w)}"
        )
    if bass.outer_d > l.interior_w:
        raise ValueError(
            f"bass driver outer diameter {_mm(bass.outer_d)} exceeds interior width {_mm(l.interior_w)}"
        )
    if coax.depth - coax.flange_t > l.interior_d:
        raise ValueError(
            f"coax driver depth behind flange {_mm(coax.depth - coax.flange_t)} exceeds "
            f"interior depth {_mm(l.interior_d)}"
        )
    if bass.depth - bass.flange_t > l.interior_d:
        raise ValueError(
            f"bass driver depth behind flange {_mm(bass.depth - bass.flange_t)} exceeds "
            f"interior depth {_mm(l.interior_d)}"
        )
    z1, z2 = sorted(z for _, z in l.bass_centres)
    if z2 - z1 < bass.outer_d:
        raise ValueError(
            f"bass frames overlap each other: centres {_mm(z1)} and {_mm(z2)} are "
            f"{_mm(z2 - z1)} apart, frame diameter is {_mm(bass.outer_d)}"
        )
    if z1 - bass.outer_d / 2 < c.panel_t:
        raise ValueError(
            f"lower bass frame overlaps the floor panel: frame bottom at "
            f"{_mm(z1 - bass.outer_d / 2)}, floor top at {_mm(c.panel_t)}"
        )
    if z2 + bass.outer_d / 2 > l.divider_z:
        raise ValueError(
            f"upper bass frame overlaps the divider: frame top at "
            f"{_mm(z2 + bass.outer_d / 2)}, divider bottom face at {_mm(l.divider_z)}"
        )
    if l.coax_chamber_empty < l.coax_chamber_net:
        raise ValueError(
            f"coax chamber too small: gross {l.coax_chamber_gross * 1000:.2f} L less coax "
            f"displacement is {l.coax_chamber_empty * 1000:.2f} L "
            f"net-available, below target {l.coax_chamber_net * 1000:.2f} L"
        )
    back = c.depth - c.panel_t
    chamber_lo = l.divider_z + c.divider_t
    chamber_hi = l.height - c.panel_t
    for i, (y, z, w, h) in enumerate(l.vent_slots, start=1):
        if y + w / 2 > back + 1e-9:
            raise ValueError(
                f"vent slot overlap with back panel: slot extends to y={_mm(y + w / 2)}, "
                f"back panel front face at {_mm(back)}"
            )
        if y - w / 2 < c.baffle_t - 1e-9:
            raise ValueError(
                f"vent slot overlap with baffle: slot starts at y={_mm(y - w / 2)}, "
                f"baffle rear face at {_mm(c.baffle_t)}"
            )
        if z - h / 2 < chamber_lo - 1e-9 or z + h / 2 > chamber_hi + 1e-9:
            raise ValueError(
                f"vent slot {i} of {len(l.vent_slots)} spans z={_mm(z - h / 2)} to {_mm(z + h / 2)}, "
                f"outside the coax chamber (divider top face at {_mm(chamber_lo)}, "
                f"top panel underside at {_mm(chamber_hi)}); reduce slots_per_side or slot_pitch"
            )
    if l.bass_net < BASS_NET_MIN:
        raise ValueError(
            f"bass chamber net volume {l.bass_net * 1000:.1f} L is below the hard minimum "
            f"of {BASS_NET_MIN * 1000:.0f} L"
        )
