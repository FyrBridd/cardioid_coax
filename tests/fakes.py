"""Shared test doubles and design helpers."""
from __future__ import annotations

from collections.abc import Callable
from typing import NamedTuple

import numpy as np

from speaker.bem.mesh import CABINET_TAGS, PUCK_TAGS
from speaker.bem.run import AXIS_R, BAND_EDGES, FAR_R, OBS_R, SourceField, observation_set
from speaker.constants import C, RHO
from speaker.drivers import Driver, load_driver
from speaker.layout import Layout, compute_layout
from speaker.params import Params, load_params
from speaker.system import Ways


class DesignTuple(NamedTuple):
    p: Params
    coax: Driver
    bass: Driver

    @property
    def layout(self) -> Layout:
        return compute_layout(self.p, self.coax, self.bass)


def default_design() -> DesignTuple:
    """params/default.toml with the drivers it names."""
    p = load_params()
    return DesignTuple(p, load_driver(p.coax), load_driver(p.bass))


def _monopoles(f: np.ndarray, sources: np.ndarray, points: np.ndarray) -> np.ndarray:
    """(nf, nsrc, np) free-field pressure per unit volume velocity, jωρ e^{-jkr} / (4πr), e^{+jωt}."""
    r = np.linalg.norm(points[None, :, :] - sources[:, None, :], axis=2)  # (nsrc, np)
    w = 2 * np.pi * f[:, None, None]
    return 1j * w * RHO * np.exp(-1j * w / C * r[None]) / (4 * np.pi * r[None])


def _band_meta(f: np.ndarray) -> dict:
    edge = BAND_EDGES[min(int(np.searchsorted(BAND_EDGES, f[-1])), len(BAND_EDGES) - 1)]
    return {"bands": [{"edge": edge, "h": 0.05, "freqs": len(f), "nodes": 1, "triangles": 2,
                       "mesh_s": 0.0, "solve_s": 0.0}]}


def monopole_field(p: Params, coax: Driver, bass: Driver, f: np.ndarray) -> SourceField:
    """cabinet_field's stand-in: a free-field point source at each source's layout position.

    Sources are named as CABINET_TAGS: both coax tags at the coax centre on the front face, each vent
    at its side's slot-stack centre, the bass drivers at their centres (lower first). Observed on
    observation_set(OBS_R) about the coax centre, with meta["far_h0_h180"] at FAR_R on H 0° and 180°.
    """
    f = np.asarray(f, dtype=float)
    lay = compute_layout(p, coax, bass)
    cx, cz = lay.coax_centre
    zs = [z for _, z, _, _ in lay.vent_slots]
    vy, vz = p.vents.centre_from_front, (min(zs) + max(zs)) / 2
    lower, upper = sorted(lay.bass_centres, key=lambda c: c[1])
    sources = np.array([[cx, 0.0, cz], [cx, 0.0, cz], [0.0, vy, vz], [p.cabinet.width, vy, vz],
                        [lower[0], 0.0, lower[1]], [upper[0], 0.0, upper[1]]])
    obs = observation_set(OBS_R)
    far = FAR_R * np.array([[0.0, -1.0, 0.0], [0.0, 1.0, 0.0]])
    points = np.array([cx, 0.0, cz]) + np.concatenate([obs.r * obs.dirs, far])
    pr = _monopoles(f, sources, points)
    return SourceField(f, CABINET_TAGS, obs, pr[..., :-2], {**_band_meta(f), "far_h0_h180": pr[..., -2:]})


def monopole_puck(
    p: Params, coax: Driver, f: np.ndarray, *, hf_axis_db: Callable[[np.ndarray], np.ndarray] | None = None
) -> SourceField:
    """puck_field's stand-in: coax_lf and coax_hf as point sources at the puck's coax centre.

    meta["axis_1m"] is the (nf, 2) pressure AXIS_R on axis. With hf_axis_db, the coax_hf column of it has
    the magnitude 10^(hf_axis_db(f)/20), keeping the monopole's phase; the observation field is unchanged.
    """
    f = np.asarray(f, dtype=float)
    centre = np.array([p.cabinet.width / 2, 0.0, p.cabinet.coax_centre_z])
    obs = observation_set(OBS_R)
    points = centre + np.concatenate([obs.r * obs.dirs, [[0.0, -AXIS_R, 0.0]]])
    pr = _monopoles(f, np.array([centre, centre]), points)
    axis = pr[..., -1].copy()
    if hf_axis_db is not None:
        col = PUCK_TAGS.index("coax_hf")
        axis[:, col] = 10 ** (np.asarray(hf_axis_db(f), dtype=float) / 20) * np.exp(1j * np.angle(axis[:, col]))
    return SourceField(f, PUCK_TAGS, obs, pr[..., :-1], {**_band_meta(f), "axis_1m": axis})


# Hand-written values for every name `speaker.passive.value_names` can return (all options on), inside BOUNDS.
SAMPLE = {"L_b1": 6.5e-3, "C_b1": 100e-6, "L_b2": 2.0e-3, "R_bz": 6.8, "C_bz": 47e-6,
          "C_h1": 100e-6, "L_h1": 8.0e-3, "R_ls": 1.5, "R_lp": 15.0, "L_l1": 1.0e-3, "C_l1": 8.2e-6,
          "L_l2": 0.5e-3, "L_ln": 0.3e-3, "C_ln": 10e-6, "R_ln": 4.7, "R_lc": 10.0, "L_lc": 3.3e-3,
          "C_t1": 10e-6, "L_t1": 0.6e-3,
          "C_t2": 15e-6, "R_ts": 4.7, "R_tp": 3.3, "R_tr": 6.8, "L_tr": 0.4e-3, "C_tr": 4.7e-6}


def omni_ways(f: np.ndarray, level_db: float = 0.0) -> Ways:
    """Every way flat at level_db (re 1 Pa peak) in every direction of observation_set(2.0), no phase."""
    f = np.asarray(f, dtype=float)
    o = observation_set(2.0)
    a = 10 ** (level_db / 20)
    return Ways(f, o, {k: a * np.ones((len(f), len(o.dirs)), complex) for k in ("bass", "coax_lf", "coax_hf")}, {})
