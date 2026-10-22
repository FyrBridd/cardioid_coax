"""Parameter loading. The only place millimetres/litres are converted to SI."""
from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path

MAX_DEPTH = 0.350  # external depth limit, m

DEFAULT_PATH = Path(__file__).resolve().parents[2] / "params" / "default.toml"


@dataclass(frozen=True)
class Cabinet:
    width: float
    depth: float
    panel_t: float
    baffle_t: float
    coax_centre_z: float
    top_margin: float
    divider_t: float
    side_edge_radius: float
    flare_mouth_d: float
    bass_centres_z: tuple[float, float]
    coax_chamber_net_volume: float


@dataclass(frozen=True)
class Vents:
    slots_per_side: int
    slot_w: float
    slot_h: float
    slot_pitch: float
    centre_from_front: float
    specific_resistance: float  # Pa*s/m


@dataclass(frozen=True)
class Room:
    length_x: float
    width_y: float
    height_z: float
    back_gap: float
    speaker_y: tuple[float, float]
    seat: tuple[float, float, float]


@dataclass(frozen=True)
class Drive:
    amp_v_rms: float
    inductor_dcr: float
    active_v_rms: float


@dataclass(frozen=True)
class Params:
    cabinet: Cabinet
    vents: Vents
    room: Room
    drive: Drive
    drivers_measured: bool
    coax: str
    bass: str


def _mm(v):
    return v / 1000.0


def _mm_seq(vs):
    return tuple(_mm(v) for v in vs)


def load_params(path: Path | None = None) -> Params:
    path = DEFAULT_PATH if path is None else Path(path)
    with open(path, "rb") as fh:
        raw = tomllib.load(fh)
    c, v, r, d = raw["cabinet"], raw["vents"], raw["room"], raw["drive"]

    cabinet = Cabinet(
        width=_mm(c["width_mm"]),
        depth=_mm(c["depth_mm"]),
        panel_t=_mm(c["panel_t_mm"]),
        baffle_t=_mm(c["baffle_t_mm"]),
        coax_centre_z=_mm(c["coax_centre_z_mm"]),
        top_margin=_mm(c["top_margin_mm"]),
        divider_t=_mm(c["divider_t_mm"]),
        side_edge_radius=_mm(c["side_edge_radius_mm"]),
        flare_mouth_d=_mm(c["flare_mouth_d_mm"]),
        bass_centres_z=_mm_seq(c["bass_centres_z_mm"]),
        coax_chamber_net_volume=c["coax_chamber_net_volume_l"] / 1000.0,
    )
    if cabinet.depth > MAX_DEPTH + 1e-9:
        raise ValueError(
            f"cabinet depth {cabinet.depth * 1000:.0f} mm exceeds the "
            f"{MAX_DEPTH * 1000:.0f} mm limit"
        )
    if cabinet.side_edge_radius != 0:
        raise ValueError(
            "side_edge_radius is not implemented yet: the CAD builds square side edges, "
            f"set side_edge_radius_mm = 0 (got {cabinet.side_edge_radius * 1000:g} mm)"
        )
    vents = Vents(
        slots_per_side=int(v["slots_per_side"]),
        slot_w=_mm(v["slot_w_mm"]),
        slot_h=_mm(v["slot_h_mm"]),
        slot_pitch=_mm(v["slot_pitch_mm"]),
        centre_from_front=_mm(v["centre_from_front_mm"]),
        specific_resistance=float(v["specific_resistance_rayl"]),
    )
    room = Room(
        length_x=_mm(r["length_x_mm"]),
        width_y=_mm(r["width_y_mm"]),
        height_z=_mm(r["height_z_mm"]),
        back_gap=_mm(r["back_gap_mm"]),
        speaker_y=_mm_seq(r["speaker_y_mm"]),
        seat=_mm_seq(r["seat_mm"]),
    )
    drive = Drive(
        amp_v_rms=float(d["amp_v_rms"]),
        inductor_dcr=float(d["inductor_dcr_ohm"]),
        active_v_rms=float(d["active_v_rms"]),
    )
    return Params(
        cabinet=cabinet,
        vents=vents,
        room=room,
        drive=drive,
        drivers_measured=bool(raw["drivers_measured"]),
        coax=raw["coax"],
        bass=raw["bass"],
    )
