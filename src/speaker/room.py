"""Room model: first-order image sources and axial modes."""
from __future__ import annotations

from typing import Iterable, Literal

import numpy as np

from speaker.constants import C
from speaker.params import Room

Surface = Literal["front", "back", "left", "right", "floor", "ceiling"]
ALL: tuple[Surface, ...] = ("front", "back", "left", "right", "floor", "ceiling")


def _image(src, surface: str, room: Room) -> tuple[float, float, float]:
    x, y, z = src
    if surface == "front":
        return (-x, y, z)
    if surface == "back":
        return (2 * room.length_x - x, y, z)
    if surface == "left":
        return (x, -y, z)
    if surface == "right":
        return (x, 2 * room.width_y - y, z)
    if surface == "floor":
        return (x, y, -z)
    if surface == "ceiling":
        return (x, y, 2 * room.height_z - z)
    raise ValueError(f"unknown surface: {surface!r}")


def boundary_response(
    src: tuple[float, float, float],
    rcv: tuple[float, float, float],
    f: np.ndarray,
    room: Room,
    surfaces: Iterable[Surface] = ALL,
    reflection: float = 1.0,
) -> np.ndarray:
    """Complex pressure transfer: direct + first-order images, each e^{-jkr}/r."""
    f = np.atleast_1d(np.asarray(f, dtype=float))
    if np.any(f <= 0):
        raise ValueError("frequencies must be > 0")
    k = 2 * np.pi * f / C
    rcv_a = np.asarray(rcv, dtype=float)

    def term(pos) -> np.ndarray:
        r = float(np.linalg.norm(np.asarray(pos, dtype=float) - rcv_a))
        return np.exp(-1j * k * r) / r

    h = term(src)
    for s in surfaces:
        h = h + reflection * term(_image(src, s, room))
    return h


def sbir_null(d: float) -> float:
    """First speaker-boundary interference null for distance d to the wall."""
    return C / (4 * d)


def axial_modes(room: Room, n_max: int = 3) -> dict[str, list[float]]:
    dims = {"x": room.length_x, "y": room.width_y, "z": room.height_z}
    return {k: [n * C / (2 * L) for n in range(1, n_max + 1)] for k, L in dims.items()}
