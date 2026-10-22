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
from .params import Params

DIVIDER_CLEARANCE = 0.010  # m between coax frame bottom and divider top
BASS_NET_TARGET = 0.050  # m^3, soft target (warning)
BASS_NET_MIN = 0.030  # m^3, hard minimum (ValueError)


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
    filler_volume: float
    bass_gross: float
    bass_net: float
    vent_slots: tuple[tuple[float, float, float, float], ...]  # (y, z, w, h)
    warnings: tuple[str, ...] = ()


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
    filler = coax_gross - _displacement(coax) - coax_net

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
        filler_volume=filler,
        bass_gross=bass_gross,
        bass_net=bass_net,
        vent_slots=vent_slots,
    )
    validate(layout, p, coax, bass)
    warns: list[str] = []
    if bass_net < BASS_NET_TARGET:
        warns.append(
            f"bass net volume {bass_net * 1000:.1f} L is below the "
            f"{BASS_NET_TARGET * 1000:.0f} L target (25 L per driver)"
        )
    return dataclasses.replace(layout, warnings=tuple(warns))


def cardioid_params(p: Params, l: Layout, coax: Driver) -> CardioidParams:
    c, v = p.cabinet, p.vents
    zs = [s[1] for s in l.vent_slots]
    z_rel = (min(zs) + max(zs)) / 2 - c.coax_centre_z
    x = -(v.centre_from_front - c.baffle_t)
    half = c.width / 2
    return CardioidParams(
        chamber_volume=l.coax_chamber_net,
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
    if l.filler_volume < 0:
        raise ValueError(
            f"coax chamber too small: gross {l.coax_chamber_gross * 1000:.2f} L less coax "
            f"displacement is {(l.coax_chamber_gross - l.coax_chamber_net - l.filler_volume) * 1000:.2f} L "
            f"net-available, below target {l.coax_chamber_net * 1000:.2f} L"
        )
    back = c.depth - c.panel_t
    for y, z, w, h in l.vent_slots:
        if y + w / 2 > back:
            raise ValueError(
                f"vent slot overlap with back panel: slot extends to y={_mm(y + w / 2)}, "
                f"back panel front face at {_mm(back)}"
            )
    if l.bass_net < BASS_NET_MIN:
        raise ValueError(
            f"bass chamber net volume {l.bass_net * 1000:.1f} L is below the hard minimum "
            f"of {BASS_NET_MIN * 1000:.0f} L"
        )
