"""Source strengths and way responses: lumped models and datasheets driving the BEM fields.

Every complex quantity here is e^{+jωt}. Ways hold the complex PEAK pressure (Pa) at
obs.r for 2.83 V RMS into that way, unfiltered: the datasheet de-embedding gives peak
volume velocities directly, and the lumped RMS volume velocities are scaled by √2.

Coax LF volume velocity is the lumped cardioid model below LF_JOIN, the datasheet
(de-embedded with a baffled piston, minimum phase) above it, blended across it. Given the
recess ratio R (speaker.bem.run.recess_ratio), the LF reference is the piston times R: the
datasheet was measured with the cone in its own frame recess, which the 3D fields add again.
Coax HF is the datasheet de-embedded with the on-axis puck response, smoothed as an excess
over the HF exit's baffled piston. Above the top of R or of the puck, the reference continues
as the piston curve, shifted by the last smoothed excess. Fields are the
cabinet BEM alone up to trusted_limit, then blended into the puck BEM up to an octave
above it (or the top of the cabinet field, at least BLEND_MIN_OCT), then the puck alone;
the vents fade out across the blend, since the puck has none. Above PUCK_HOLD_F each
direction keeps its 12 kHz level relative to the datasheet on-axis level. When the cabinet
field carries far-field points (meta["far_h0_h180"]), the coax section is also composed
there, up to the trusted limit.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from speaker.bass import cone_velocity
from speaker.bem.analytic import baffled_piston_on_axis
from speaker.bem.run import FAR_R, Observation, SourceField
from speaker.cardioid import CardioidParams, simulate
from speaker.constants import C
from speaker.curves import interp_db, load_curve, minimum_phase, smooth
from speaker.drivers import Driver
from speaker.layout import Layout, cardioid_params
from speaker.params import Params

LF_JOIN = (300.0, 500.0)  # Hz, lumped -> datasheet blend for the coax LF volume velocity
PUCK_HOLD_F = 12000.0  # Hz, top of the puck field; directivity is held above it
BLEND_MIN_OCT = 1 / 3  # octaves, narrowest cabinet -> puck blend above the trusted limit
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
    # lf_join_mismatch_db, trusted_limit, hf_reference ("puck" | "piston"), lf_reference ("recess" | "piston");
    # u_lf: coax LF peak volume velocity up to the hold;
    # puck_blend: (lo, hi) Hz of the cabinet -> puck crossfade (lo = trusted limit), None without a puck
    # or when f stops at the trusted limit;
    # coax_far_h0_h180: (nf, 2) coax section at FAR_R on H 0° and 180°, NaN above the trusted
    # limit, None when the cabinet field has no far-field points
    meta: dict


def coax_section(ways: Ways) -> np.ndarray:
    """The coax LF way (cone plus both vents), the section the cardioid gate scores."""
    return ways.ways["coax_lf"]


def vent_scaled_section(
    ways: Ways, cab: SourceField, coax: Driver, cp: CardioidParams, scale: float
) -> np.ndarray:
    """The coax section with the vents re-strengthened for vent_resistance x scale, no new BEM run.

    Each vent's strength is the lumped u_vent/u_cone ratio times the coax LF strength u_lf
    (ways.meta["u_lf"]) at every frequency; u_lf stays as it is and only the ratio is taken at the
    scaled resistance, on the same cabinet vent fields, with the same fade across the
    cabinet-to-puck blend. Above the cabinet field (the trusted limit, or the blend's top) the
    section is unchanged.
    """
    f = ways.f
    fade = _vent_fade(f, ways.meta["trusted_limit"], ways.meta["puck_blend"])
    m = fade > 0
    out = coax_section(ways).copy()
    if scale == 1.0 or not m.any():
        return out
    nom = simulate(coax, cp, f[m], V_DRIVE)
    alt = simulate(coax, replace(cp, vent_resistance=cp.vent_resistance * scale), f[m], V_DRIVE)
    u_lf = ways.meta["u_lf"][: int(m.sum())]  # on way_responses' internal grid, equal to f up to the hold
    d_vent = (alt.u_vent / alt.u_cone - nom.u_vent / nom.u_cone) * u_lf / 2
    vents = _field(cab, "vent_left", f[m]) + _field(cab, "vent_right", f[m])
    out[m] += (fade[m] * d_vent)[:, None] * vents
    return out


def _vent_fade(f: np.ndarray, tl: float, blend: tuple[float, float] | None) -> np.ndarray:
    """Weight of the cabinet vent fields: 1 up to tl, fading to 0 across the puck blend."""
    if blend is None:
        return (f <= tl).astype(float)
    return np.where(f < blend[1], 1.0 - _log_ramp(f, *blend), 0.0)


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
    return _interp(fld.f, fld.p[:, fld.names.index(name), :], fld.obs.r, ft)


def _interp(fs: np.ndarray, p: np.ndarray, r: float, ft: np.ndarray) -> np.ndarray:
    """p (len(fs), n) at distance r, onto ft: dB and unwrapped phase, bulk e^{-jkr} removed."""
    if np.array_equal(fs, ft):
        return p.copy()
    q = p * np.exp(1j * 2 * np.pi * fs / C * r)[:, None]
    db = 20 * np.log10(np.maximum(np.abs(q), 1e-300))
    ph = np.unwrap(np.angle(q), axis=0)
    out = 10 ** (_lerp_log(fs, db, ft) / 20) * np.exp(1j * _lerp_log(fs, ph, ft))
    out *= np.exp(-1j * 2 * np.pi * ft / C * r)[:, None]
    hit = np.isin(ft, fs)
    out[hit] = p[np.searchsorted(fs, ft[hit])]
    return out


def _blend(a: np.ndarray, b: np.ndarray, t: np.ndarray) -> np.ndarray:
    """a -> b along axis 0 (frequency): dB and phase linear in t.

    The phase difference angle(b/a) is unwrapped along frequency across the band
    0 < t < 1, anchored to its principal value at the band's low end.
    """
    out = np.where((t <= 0).reshape(-1, *[1] * (a.ndim - 1)), a, b)
    band = (t > 0) & (t < 1)
    if band.any():
        r = b[band] / a[band]
        tb = t[band].reshape(-1, *[1] * (a.ndim - 1))
        out[band] = a[band] * np.abs(r) ** tb * np.exp(1j * tb * np.unwrap(np.angle(r), axis=0))
    return out


def _log_ramp(f: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.clip(np.log(f / lo) / np.log(hi / lo), 0.0, 1.0)


def _piston_db(f: np.ndarray, d: float) -> np.ndarray:
    return 20 * np.log10(np.abs([baffled_piston_on_axis(x, d / 2, 1.0) for x in f]))


def recess_db(lf_recess: tuple[np.ndarray, np.ndarray], x: np.ndarray) -> np.ndarray:
    """The LF reference's recess term on the increasing grid x (Hz): 20·log10|R|, 1/3-octave smoothed.

    20·log10|R| is interpolated log-linearly onto the points of x inside R's range (and R's end
    frequencies) and smoothed there, so the window is even whatever R's own spacing. Outside R's
    range it holds its end values: the reference continues as the piston curve, shifted.
    """
    fr, r = (np.asarray(a) for a in lf_recess)
    fr, x = fr.astype(float), np.asarray(x, dtype=float)
    db = 20 * np.log10(np.abs(r))
    xs = np.union1d(x[(x >= fr[0]) & (x <= fr[-1])], fr[[0, -1]])
    d = np.interp(np.log(xs), np.log(fr), db)
    sm = smooth(xs, d, REF_SMOOTH) if len(xs) > 1 else d
    return np.interp(np.log(x), np.log(xs), sm)


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
    lf_recess: tuple[np.ndarray, np.ndarray] | None = None,
) -> Ways:
    """Complex peak pressure of each way at cab.obs for 2.83 V into it, on the grid f.

    lf_recess: (f, R) from recess_ratio, the coax cone's own recess relative to a flat piston. When
    given, the coax LF datasheet is de-embedded with the baffled piston times R (see recess_db).
    """
    f = np.asarray(f, dtype=float)
    tl = float(trusted_limit)
    if np.any(np.diff(f) <= 0):
        raise ValueError("f must be strictly increasing")
    if not _covers(cab, f[0], min(f[-1], tl)):
        raise ValueError(f"cabinet field covers {cab.f[0]:g}-{cab.f[-1]:g} Hz, "
                         f"need {f[0]:g}-{min(f[-1], tl):g} Hz")
    blend = None
    if puck is None:
        if f[-1] > tl:
            raise ValueError(f"without a puck field f must stop at trusted_limit, {tl:g} Hz")
    else:
        if f[-1] > tl:
            blend = (tl, float(min(2 * tl, cab.f[-1], PUCK_HOLD_F)))
        if blend is not None and np.log2(blend[1] / tl) < BLEND_MIN_OCT - 1e-9:
            raise ValueError(f"cabinet field stops at {cab.f[-1]:g} Hz: the cabinet-to-puck blend above the "
                             f"trusted limit {tl:g} Hz needs it to reach {tl * 2 ** BLEND_MIN_OCT:g} Hz")
        need = (tl, min(f[-1], PUCK_HOLD_F))
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
    if not hold:
        fi = f
    elif f[n_lo - 1] == PUCK_HOLD_F:
        fi = f[:n_lo]
    else:
        fi = np.append(f[:n_lo], PUCK_HOLD_F)

    # Datasheet volume velocities on a dense grid, minimum phase, then onto fi.
    g_lo, g_hi = min(f[0], _DS_RANGE[0]), max(f[-1], _DS_RANGE[1])
    g = np.geomspace(g_lo, g_hi, int(np.ceil(_PPO * np.log2(g_hi / g_lo))) + 1)
    lf_spl, hf_spl = (interp_db(load_curve(f"{COAX_CURVES}_{w}_spl"), g) for w in ("lf", "hf"))
    hf_ref = _piston_db(g, coax.hf_exit_d)
    if puck is not None:
        # The puck's on-axis level is smoothed as an excess over the HF exit's baffled piston (puck dB
        # minus piston dB), and the piston added back. Above the puck's top the excess holds its last
        # smoothed value, so the reference continues as the piston curve, shifted to meet the puck; below
        # its bottom the reference is the piston alone. Smoothing the excess, not the absolute level,
        # avoids the bias of a window cut short at the puck's ends on a rising curve. Against round 0b,
        # which smoothed the absolute level, the reference moves by about 0.05 dB mid-band and up to
        # about 0.5 dB at the puck's 500 Hz and 12 kHz ends.
        col = puck.meta["axis_1m"][:, _AXIS_COLS.index("coax_hf")]
        excess = smooth(puck.f, 20 * np.log10(np.abs(col)) - _piston_db(puck.f, coax.hf_exit_d), REF_SMOOTH)
        inside = g >= puck.f[0]
        hf_ref[inside] += np.interp(np.log(g[inside]), np.log(puck.f), excess)
    lf_ref = _piston_db(g, coax.diaphragm_d)
    if lf_recess is not None:
        lf_ref += recess_db(lf_recess, g)
    u_db = {"lf": lf_spl + 20 * np.log10(P_PEAK_REF) - lf_ref,
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

    # Fields on fi: cabinet alone up to tl; for the coax sources, cabinet -> puck across the blend.
    def cab_on(name: str, x: np.ndarray, top: float = tl) -> np.ndarray:
        out = np.zeros((len(x), len(cab.obs.dirs)), dtype=complex)
        m = x <= top
        out[m] = _field(cab, name, x[m])
        return out

    def coax_on(name: str) -> np.ndarray:
        if blend is None:
            return cab_on(name, fi)
        out = cab_on(name, fi, blend[1])
        m = fi > tl
        out[m] = _blend(out[m], _field(puck, name, fi[m]), _log_ramp(fi[m], *blend))
        return out

    fade = _vent_fade(fi, tl, blend)
    top = tl if blend is None else blend[1]
    vents = fade[:, None] * (cab_on("vent_left", fi, top) + cab_on("vent_right", fi, top))
    lf_i = u_lf[:, None] * coax_on("coax_lf") + u_vent[:, None] * vents

    # Coax section at the far-field points, from the cabinet field alone (the puck has none).
    far = None
    if "far_h0_h180" in cab.meta:
        far = np.full((len(f), 2), np.nan, dtype=complex)
        n = int(np.count_nonzero(fi <= tl))  # fi equals f below the hold, and tl <= the hold
        fp = cab.meta["far_h0_h180"]
        src = {k: _interp(cab.f, fp[:, cab.names.index(k), :], FAR_R, fi[:n])
               for k in ("coax_lf", "vent_left", "vent_right")}
        far[:n] = u_lf[:n, None] * src["coax_lf"] + u_vent[:n, None] * (src["vent_left"] + src["vent_right"])
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
        meta={"lf_join_mismatch_db": mismatch, "trusted_limit": tl, "u_lf": u_lf,
              "hf_reference": "puck" if puck is not None else "piston",
              "lf_reference": "recess" if lf_recess is not None else "piston",
              "puck_blend": blend, "coax_far_h0_h180": far},
    )
