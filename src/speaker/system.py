"""Source strengths and way responses: lumped models and datasheets driving the BEM fields.

Every complex quantity here is e^{+jωt}. Ways hold the complex PEAK pressure (Pa) at
obs.r for 2.83 V RMS into that way, unfiltered: the datasheet de-embedding gives peak
volume velocities directly, and the lumped RMS volume velocities are scaled by √2.

Coax LF volume velocity is the lumped cardioid model below LF_JOIN, the datasheet
(de-embedded with a baffled piston, minimum phase) above it, blended across it. Coax HF
is the datasheet de-embedded with the smoothed on-axis puck response. Fields are the
cabinet BEM up to trusted_limit/2, blended into the puck BEM up to trusted_limit, then
the puck alone; above PUCK_HOLD_F each direction keeps its 12 kHz level relative to the
datasheet on-axis level.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from speaker.bass import cone_velocity
from speaker.bem.analytic import baffled_piston_on_axis
from speaker.bem.run import Observation, SourceField
from speaker.cardioid import simulate
from speaker.constants import C
from speaker.curves import interp_db, load_curve, minimum_phase, smooth
from speaker.drivers import Driver
from speaker.layout import Layout, cardioid_params
from speaker.params import Params

LF_JOIN = (300.0, 500.0)  # Hz, lumped -> datasheet blend for the coax LF volume velocity
PUCK_HOLD_F = 12000.0  # Hz, top of the puck field; directivity is held above it
REF_SMOOTH = 1 / 3  # octaves, smoothing of the puck's on-axis reference
V_DRIVE = 2.83  # V RMS
P_PEAK_REF = 20e-6 * np.sqrt(2)  # Pa peak at 0 dB SPL
COAX_CURVES = "bms_10c262"  # data/curves stem of the coax datasheet traces
_AXIS_COLS = ("coax_lf", "coax_hf")  # column order of a puck's meta["axis_1m"]
_PPO = 48  # points per octave of the internal datasheet grid
_DS_RANGE = (20.0, 20000.0)  # Hz, minimum span of the internal datasheet grid


@dataclass
class Ways:
    f: np.ndarray
    obs: Observation
    ways: dict[str, np.ndarray]  # "bass", "coax_lf", "coax_hf": (nf, nd) complex peak Pa
    meta: dict


def coax_section(ways: Ways) -> np.ndarray:
    """The coax LF way (cone plus both vents), the section the cardioid gate scores."""
    return ways.ways["coax_lf"]


def _lerp_log(fs: np.ndarray, y: np.ndarray, ft: np.ndarray) -> np.ndarray:
    """y(fs) at ft, linear in log f along axis 0, held at the ends; exact where ft hits fs."""
    lfs, lft = np.log(fs), np.log(ft)
    i1 = np.clip(np.searchsorted(lfs, lft), 1, len(fs) - 1)
    i0 = i1 - 1
    w = np.clip((lft - lfs[i0]) / (lfs[i1] - lfs[i0]), 0.0, 1.0).reshape(-1, *[1] * (y.ndim - 1))
    out = y[i0] * (1 - w) + y[i1] * w
    hit = np.isin(ft, fs)
    out[hit] = y[np.searchsorted(fs, ft[hit])]
    return out


def _field(fld: SourceField, name: str, ft: np.ndarray) -> np.ndarray:
    """(len(ft), nd) field of one source on ft: dB and unwrapped phase, bulk e^{-jkr} removed."""
    p = fld.p[:, fld.names.index(name), :]
    if np.array_equal(fld.f, ft):
        return p.copy()
    r = fld.obs.r
    q = p * np.exp(1j * 2 * np.pi * fld.f / C * r)[:, None]
    db = 20 * np.log10(np.maximum(np.abs(q), 1e-300))
    ph = np.unwrap(np.angle(q), axis=0)
    out = 10 ** (_lerp_log(fld.f, db, ft) / 20) * np.exp(1j * _lerp_log(fld.f, ph, ft))
    out *= np.exp(-1j * 2 * np.pi * ft / C * r)[:, None]
    hit = np.isin(ft, fld.f)
    out[hit] = p[np.searchsorted(fld.f, ft[hit])]
    return out


def _blend(a: np.ndarray, b: np.ndarray, t: np.ndarray) -> np.ndarray:
    """a -> b: dB and phase linear in t, the phase step taken at its principal value."""
    t = t.reshape(-1, *[1] * (a.ndim - 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        mid = a * (b / a) ** t
    return np.where(t <= 0, a, np.where(t >= 1, b, mid))


def _log_ramp(f: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.clip(np.log(f / lo) / np.log(hi / lo), 0.0, 1.0)


def _piston_db(f: np.ndarray, d: float) -> np.ndarray:
    return 20 * np.log10(np.abs([baffled_piston_on_axis(x, d / 2, 1.0) for x in f]))


def _covers(fld: SourceField, lo: float, hi: float) -> bool:
    return fld.f[0] <= lo * (1 + 1e-9) and fld.f[-1] >= hi * (1 - 1e-9)


def way_responses(
    p: Params,
    coax: Driver,
    bass: Driver,
    layout: Layout,
    cab: SourceField,
    puck: SourceField | None,
    f: np.ndarray,
    *,
    trusted_limit: float,
) -> Ways:
    """Complex peak pressure of each way at cab.obs for 2.83 V into it, on the grid f."""
    f = np.asarray(f, dtype=float)
    tl = float(trusted_limit)
    if np.any(np.diff(f) <= 0):
        raise ValueError("f must be strictly increasing")
    if not _covers(cab, f[0], min(f[-1], tl)):
        raise ValueError(f"cabinet field covers {cab.f[0]:g}-{cab.f[-1]:g} Hz, "
                         f"need {f[0]:g}-{min(f[-1], tl):g} Hz")
    if puck is None:
        if f[-1] > tl:
            raise ValueError(f"without a puck field f must stop at trusted_limit, {tl:g} Hz")
    else:
        need = (tl / 2, min(f[-1], PUCK_HOLD_F))
        if need[0] <= need[1] and not _covers(puck, *need):
            raise ValueError(f"puck field covers {puck.f[0]:g}-{puck.f[-1]:g} Hz, "
                             f"need {need[0]:g}-{need[1]:g} Hz")
        if tl > PUCK_HOLD_F:
            raise ValueError(f"trusted_limit {tl:g} Hz is above the puck hold, {PUCK_HOLD_F:g} Hz")
        if puck.obs.r != cab.obs.r or not np.allclose(puck.obs.dirs, cab.obs.dirs):
            raise ValueError("puck and cabinet fields use different observation sets")

    # Internal grid: f up to the hold frequency, ending exactly on it when f goes past it.
    hold = puck is not None and f[-1] > PUCK_HOLD_F
    n_lo = int(np.searchsorted(f, PUCK_HOLD_F, side="right")) if hold else len(f)
    fi = f if not hold or f[n_lo - 1] == PUCK_HOLD_F else np.append(f[:n_lo], PUCK_HOLD_F)

    # Datasheet volume velocities on a dense grid, minimum phase, then onto fi.
    g_lo, g_hi = min(f[0], _DS_RANGE[0]), max(f[-1], _DS_RANGE[1])
    g = np.geomspace(g_lo, g_hi, int(np.ceil(_PPO * np.log2(g_hi / g_lo))) + 1)
    lf_spl, hf_spl = (interp_db(load_curve(f"{COAX_CURVES}_{w}_spl"), g) for w in ("lf", "hf"))
    hf_ref = _piston_db(g, coax.hf_exit_d)
    if puck is not None:
        col = puck.meta["axis_1m"][:, _AXIS_COLS.index("coax_hf")]
        sm = smooth(puck.f, 20 * np.log10(np.abs(col)), REF_SMOOTH)
        inside = g >= puck.f[0]
        hf_ref[inside] = np.interp(np.log(g[inside]), np.log(puck.f), sm)
    u_db = {"lf": lf_spl + 20 * np.log10(P_PEAK_REF) - _piston_db(g, coax.diaphragm_d),
            "hf": hf_spl + 20 * np.log10(P_PEAK_REF) - hf_ref}
    u_ph = {k: minimum_phase(g, v) for k, v in u_db.items()}

    def datasheet_u(key: str, x: np.ndarray) -> np.ndarray:
        return 10 ** (_lerp_log(g, u_db[key], x) / 20) * np.exp(1j * _lerp_log(g, u_ph[key], x))

    # Coax LF: lumped cone (RMS -> peak) below the join, datasheet above, blended across it.
    cp = cardioid_params(p, layout, coax)
    lump = simulate(coax, cp, fi, V_DRIVE)
    u_lump = np.sqrt(2) * lump.u_cone
    u_lf = _blend(u_lump, datasheet_u("lf", fi), _log_ramp(fi, *LF_JOIN))
    u_vent = lump.u_vent / lump.u_cone * u_lf / 2  # each of the two vents
    u_hf = datasheet_u("hf", fi)
    at = np.array([400.0])
    lump400 = np.sqrt(2) * simulate(coax, cp, at, V_DRIVE).u_cone
    mismatch = float(20 * np.log10(np.abs(lump400[0]) / np.abs(datasheet_u("lf", at)[0])))

    # Fields on fi: cabinet below tl, puck above tl/2 for the coax sources.
    def cab_on(name: str, x: np.ndarray) -> np.ndarray:
        out = np.zeros((len(x), len(cab.obs.dirs)), dtype=complex)
        m = x <= tl
        out[m] = _field(cab, name, x[m])
        return out

    def coax_on(name: str) -> np.ndarray:
        if puck is None:
            return cab_on(name, fi)
        out = cab_on(name, fi)
        m = fi > tl / 2
        out[m] = _blend(out[m], _field(puck, name, fi[m]), _log_ramp(fi[m], tl / 2, tl))
        return out

    vents = cab_on("vent_left", fi) + cab_on("vent_right", fi)
    lf_i = u_lf[:, None] * coax_on("coax_lf") + u_vent[:, None] * vents
    hf_i = u_hf[:, None] * coax_on("coax_hf")

    def extend(way: np.ndarray, key: str) -> np.ndarray:
        """way on f: fi below the hold, then the 12 kHz pattern on the datasheet on-axis curve."""
        if not hold:
            return way
        x = f[n_lo:]
        spl = lf_spl if key == "lf" else hf_spl
        ph = minimum_phase(g, spl)
        d_db = _lerp_log(g, spl, x) - _lerp_log(g, spl, np.array([PUCK_HOLD_F]))
        d_ph = _lerp_log(g, ph, x) - _lerp_log(g, ph, np.array([PUCK_HOLD_F]))
        d_k = 2 * np.pi * (x - PUCK_HOLD_F) / C * cab.obs.r
        ratio = 10 ** (d_db / 20) * np.exp(1j * (d_ph - d_k))
        return np.concatenate([way[:n_lo], way[-1][None, :] * ratio[:, None]])

    # Bass: two drivers in parallel at 2.83 V, each Sd·u (RMS -> peak), zero above tl.
    _, u = cone_velocity(bass, layout.bass_net / 2, f, V_DRIVE)
    u_bass = np.sqrt(2) * bass.sd * u
    bass_way = u_bass[:, None] * (cab_on("bass_lower", f) + cab_on("bass_upper", f))

    return Ways(
        f=f, obs=cab.obs,
        ways={"bass": bass_way, "coax_lf": extend(lf_i, "lf"), "coax_hf": extend(hf_i, "hf")},
        meta={"lf_join_mismatch_db": mismatch, "trusted_limit": tl,
              "hf_reference": "puck" if puck is not None else "piston"},
    )
