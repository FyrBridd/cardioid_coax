"""CTA-2034 style spinorama from a total complex response, with a listening window
re-centred on the listener (+23° horizontal).

Inputs are complex PEAK pressures (Pa); levels are dB SPL re 20 µPa RMS. Power averages
are means of |p|² taken before the log.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from speaker.bem.run import Observation

P_REF_PEAK = 20e-6 * np.sqrt(2)  # Pa peak at 0 dB SPL

LW_H_ANGLES = (0, 10, 20, 30)  # horizontal, positive side only; plus the two az 23 el +-10 directions
FLOOR_V = (-20, -30, -40)
CEILING_V = (40, 50, 60)
FRONT_H = (0, 10, 20, 30)  # ±, 0 once
SIDE_H = (40, 50, 60, 70, 80)  # ±
REAR_H = (90, 180)  # ±90 and the 180 entry


@dataclass
class Spin:
    f: np.ndarray | None
    on_axis: np.ndarray
    at23: np.ndarray
    lw: np.ndarray
    er: np.ndarray
    sp: np.ndarray
    erdi: np.ndarray
    spdi: np.ndarray
    er_parts: dict[str, np.ndarray]


def _p2(total: np.ndarray, idx) -> np.ndarray:
    """Mean |p|² over the given direction indices, (nf,)."""
    return np.mean(np.abs(total[:, np.asarray(idx, dtype=int)]) ** 2, axis=1)


def _db(p2: np.ndarray) -> np.ndarray:
    return 10 * np.log10(p2 / P_REF_PEAK**2)


def _ring(angles: np.ndarray, idx: np.ndarray, wanted) -> list[int]:
    out = []
    for a in wanted:
        a = ((a + 180) % 360) - 180  # the ring covers -180..175, so 180 is -180
        hit = np.flatnonzero(np.isclose(angles, a))
        if not len(hit):
            raise ValueError(f"angle {a} is not on the ring")
        out.append(int(idx[hit[0]]))
    return out


def spinorama(total: np.ndarray, obs: Observation, *, f: np.ndarray | None = None) -> Spin:
    """Spinorama curves from total (nf, nd) complex peak Pa over obs.dirs."""
    h = lambda angs: _ring(obs.h_angles, obs.h_idx, angs)  # noqa: E731
    v = lambda angs: _ring(obs.v_angles, obs.v_idx, angs)  # noqa: E731
    pm = lambda angs: [*angs, *(-a for a in angs if a)]  # noqa: E731

    lw_idx = h(LW_H_ANGLES) + list(obs.lw_extra_idx)
    parts_idx = {
        "floor": v(FLOOR_V),
        "ceiling": v(CEILING_V),
        "front": h(pm(FRONT_H)),
        "side": h(pm(SIDE_H)),
        "rear": h((90, -90, 180)),
    }
    er_parts = {k: _db(_p2(total, ix)) for k, ix in parts_idx.items()}
    lw = _db(_p2(total, lw_idx))
    er = _db(np.mean([_p2(total, ix) for ix in parts_idx.values()], axis=0))
    sp = _db(np.sum(np.abs(total[:, obs.sphere_idx]) ** 2 * obs.sphere_w, axis=1))
    return Spin(
        f=f,
        on_axis=_db(_p2(total, h([0]))),
        at23=_db(_p2(total, [obs.at23_idx])),
        lw=lw, er=er, sp=sp, erdi=lw - er, spdi=lw - sp, er_parts=er_parts,
    )
