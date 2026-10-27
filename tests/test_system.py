import numpy as np, pytest
from speaker.bem import analytic
from speaker.bem.run import SourceField, observation_set
from speaker.cardioid import simulate
from speaker.drivers import load_driver
from speaker.layout import cardioid_params, compute_layout
from speaker.params import load_params
from speaker.system import way_responses, coax_section, LF_JOIN
from speaker.curves import load_curve, interp_db

P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
L = compute_layout(P, CX, B)
F = np.geomspace(100, 16000, 120)
NAMES = ("coax_lf", "coax_hf", "vent_left", "vent_right", "bass_lower", "bass_upper")

def _piston_field(f, r=2.0, scale=1.0):
    o = observation_set(r)
    pp = np.array([analytic.baffled_piston_on_axis(x, CX.diaphragm_d / 2, r) for x in f]) * scale
    p = np.broadcast_to(pp[:, None, None], (len(f), len(NAMES), len(o.dirs))).copy()
    p[:, 2:4, :] = 0.0   # vents silent, so coax_lf is the cone alone
    return SourceField(f, NAMES, o, p, {"axis_1m": np.stack([np.array([analytic.baffled_piston_on_axis(x, CX.hf_exit_d / 2, 1.0) for x in f])] * 2, 1)})

def test_deembedding_reproduces_datasheet_at_2m():
    fld = _piston_field(F)
    w = way_responses(P, CX, B, L, fld, fld, F, trusted_limit=3000.0)
    m = (F > 600) & (F < 2500)
    got = 20 * np.log10(np.abs(w.ways["coax_lf"][m, w.obs.h_idx[w.obs.h_angles == 0][0]]) / 20e-6 / np.sqrt(2))
    # Vents are silent in this field: the cone must match the datasheet minus 6 dB (1 m -> 2 m) within 1 dB.
    want = interp_db(load_curve("bms_10c262_lf_spl"), F[m]) - 20 * np.log10(2)
    assert np.abs(got - want).max() < 1.0

def test_blend_moves_monotonically_from_cabinet_to_puck():
    cab, puck = _piston_field(F), _piston_field(F, scale=2.0)   # puck 6 dB hotter everywhere
    w0 = way_responses(P, CX, B, L, cab, cab, F, trusted_limit=3000.0)
    w1 = way_responses(P, CX, B, L, cab, puck, F, trusted_limit=3000.0)
    i0 = w0.obs.h_idx[w0.obs.h_angles == 0][0]
    d = 20 * np.log10(np.abs(w1.ways["coax_hf"][:, i0]) / np.abs(w0.ways["coax_hf"][:, i0]))
    # The trusted range is the cabinet alone; the blend runs an octave above the limit.
    assert np.all(d[F <= 3000] == 0) and np.allclose(d[(F > 6000) & (F < 12000)], 6.02, atol=0.05)
    m = (F >= 3000) & (F <= 6000)
    assert np.all(np.diff(d[m]) >= -1e-9) and 0.5 < d[m][len(d[m]) // 2] < 5.5
    assert np.all(w1.ways["bass"][F > 3000] == 0)

def test_lf_join_mismatch_reported():
    fld = _piston_field(F)
    w = way_responses(P, CX, B, L, fld, fld, F, trusted_limit=3000.0)
    assert abs(w.meta["lf_join_mismatch_db"]) < 6.0

def test_hold_anchors_on_12k_when_the_grid_contains_it():
    def hf_axis(f):
        fld = _piston_field(f)
        w = way_responses(P, CX, B, L, fld, fld, f, trusted_limit=3000.0)
        return 20 * np.log10(np.abs(w.ways["coax_hf"][-2:, w.obs.h_idx[w.obs.h_angles == 0][0]]))
    exact = hf_axis(np.array([1000.0, 5000.0, 11000.0, 12000.0, 14000.0, 16000.0]))
    near = hf_axis(np.array([1000.0, 5000.0, 11000.0, 11999.999, 14000.0, 16000.0]))
    assert np.allclose(exact, near, atol=0.01)

def test_cab_to_puck_blend_phase_is_continuous_through_pi():
    f = np.geomspace(100, 16000, 400)   # ~54 points per octave: the smooth blend moves < 4° per step
    cab = _piston_field(f)
    # The puck leads the cabinet by a phase that crosses pi inside the 3-6 kHz blend band.
    ramp = np.clip(np.log(f / 3000.0) / np.log(2.0), 0.0, 1.0)
    puck = _piston_field(f)
    puck.p[:] *= np.exp(1j * np.pi * (0.8 + 0.4 * ramp))[:, None, None]
    w0 = way_responses(P, CX, B, L, cab, cab, f, trusted_limit=3000.0)
    w1 = way_responses(P, CX, B, L, cab, puck, f, trusted_limit=3000.0)
    i0 = w0.obs.h_idx[w0.obs.h_angles == 0][0]
    r = w1.ways["coax_hf"][:, i0] / w0.ways["coax_hf"][:, i0]
    assert np.degrees(np.abs(np.angle(r[1:] / r[:-1]))).max() < 10.0

# Measured at R=190 Pa·s/m: BEM coax-section front/back 5.9/5.7 dB against lumped
# front_back 22.3/24.0 dB at 150/200 Hz. The cabinet shadows the cone from the rear and the
# vents from the front, which the free-field lumped model cannot see. The 3D model is the
# authority on the cardioid (spec §4/§8), so this test pins the vent sign and the shadowing
# instead of agreement with the lumped model.
@pytest.mark.bem
def test_bem_vent_sign_and_cabinet_shadowing(tmp_path):
    from speaker.bem.run import cabinet_field
    from speaker.params import DEFAULT_PATH
    f = np.array([150.0, 200.0])
    cab = cabinet_field(P, CX, B, f, h_max=0.06, cache_dir=tmp_path)
    w = way_responses(P, CX, B, L, cab, None, f, trusted_limit=3000.0)
    o = w.obs
    front, back = o.h_idx[o.h_angles == 0][0], o.h_idx[o.h_angles == -180][0]
    fb = lambda s: 20 * np.log10(np.abs(s[:, front]) / np.abs(s[:, back]))
    # Below LF_JOIN the coax LF strengths are the lumped ones, RMS -> peak.
    res = simulate(CX, cardioid_params(P, L, CX), f, 2.83)
    cone, vl, vr = (cab.p[:, cab.names.index(n), :] for n in ("coax_lf", "vent_left", "vent_right"))
    vents = np.sqrt(2) * res.u_vent[:, None] / 2 * (vl + vr)
    s = coax_section(w)
    assert np.allclose(s, np.sqrt(2) * res.u_cone[:, None] * cone + vents)
    flipped = np.sqrt(2) * res.u_cone[:, None] * cone - vents
    assert np.all(fb(s) >= fb(flipped) + 1.0)
    # Regression guard for the hard gate: the measured 5.9 / 5.7 dB at 150 / 200 Hz.
    assert fb(s) == pytest.approx([5.9, 5.7], abs=1.5)
    # Shadowing at 150 Hz, per unit volume velocity.
    assert abs(cone[0, front]) > abs(cone[0, back])
    for v in (vl, vr):
        assert abs(v[0, back]) > abs(v[0, front])


def test_vent_scaled_section_keeps_the_cone_and_rescales_the_vents():
    import dataclasses
    from speaker.system import vent_scaled_section
    fld = _piston_field(F)
    o = fld.obs
    for k, n in ((2, "vent_left"), (3, "vent_right")):  # direction-dependent vent fields
        fld.p[:, k, :] = fld.p[:, 0, :] * (0.5 + o.dirs[:, 1] * (0.3 if k == 2 else -0.2))[None, :]
    cp = cardioid_params(P, L, CX)
    tl = 3000.0
    w = way_responses(P, CX, B, L, fld, fld, F, trusted_limit=tl)
    assert np.array_equal(vent_scaled_section(w, fld, CX, cp, 1.0), coax_section(w))
    s = 1.3
    got = vent_scaled_section(w, fld, CX, cp, s)
    # Oracle: a full way_responses at R x s, with its coax LF strength swapped back to the nominal one.
    # The vents follow the lumped ratio wherever the cabinet field is used: up to the top of the
    # cabinet-to-puck blend, fading out across it as the puck (which has no vents) takes over.
    p_alt = dataclasses.replace(P, vents=dataclasses.replace(P.vents, specific_resistance=P.vents.specific_resistance * s))
    w_alt = way_responses(p_alt, CX, B, L, fld, fld, F, trusted_limit=tl)
    lo, hi = w.meta["puck_blend"]
    m = F < hi
    n = int(m.sum())
    u_nom, u_alt = w.meta["u_lf"][:n], w_alt.meta["u_lf"][:n]
    assert np.allclose(got[m], coax_section(w_alt)[m] * (u_nom / u_alt)[:, None], rtol=1e-9)
    j = (F > 500) & m   # datasheet strength above the join: the oracle needs no swap there
    assert np.allclose(u_nom[j[:n]], u_alt[j[:n]]) and not np.allclose(got[j], coax_section(w)[j])
    assert not np.allclose(got[F <= LF_JOIN[0]], coax_section(w)[F <= LF_JOIN[0]])
    assert not np.allclose(got[(F > tl) & m], coax_section(w)[(F > tl) & m])   # vents still there above tl
    assert np.array_equal(got[~m], coax_section(w)[~m])


def test_far_pair_is_composed_like_the_coax_section():
    fld = _piston_field(F)
    o = fld.obs
    for k, s in ((2, 0.3), (3, -0.2)):   # direction-dependent vents, louder to the rear
        fld.p[:, k, :] = -fld.p[:, 0, :] * (0.5 + s * o.dirs[:, 1])[None, :]
    i0, i180 = o.h_idx[o.h_angles == 0][0], o.h_idx[o.h_angles == -180][0]
    fld.meta["far_h0_h180"] = fld.p[:, :, [i0, i180]]   # stand-in: the 2 m H 0/180 values
    tl = 3000.0
    w = way_responses(P, CX, B, L, fld, fld, F, trusted_limit=tl)
    far = w.meta["coax_far_h0_h180"]
    assert far.shape == (len(F), 2)
    m = F <= tl
    assert np.allclose(far[m], coax_section(w)[m][:, [i0, i180]], rtol=1e-12)
    assert np.all(np.isnan(far[~m]))   # the cabinet field stops at the trusted limit
    # Without the far points (synthetic fields, older caches) there is no pair.
    del fld.meta["far_h0_h180"]
    assert way_responses(P, CX, B, L, fld, fld, F, trusted_limit=tl).meta["coax_far_h0_h180"] is None


def test_puck_blend_range_in_meta():
    fld = _piston_field(F)
    # The blend starts at the trusted limit and runs up to an octave above it ...
    assert way_responses(P, CX, B, L, fld, fld, F, trusted_limit=3000.0).meta["puck_blend"] == (3000.0, 6000.0)
    # ... or to the top of the cabinet field when that comes first.
    f4 = F[F <= 4000.0]
    cab4 = SourceField(f4, fld.names, fld.obs, fld.p[: len(f4)], {"axis_1m": fld.meta["axis_1m"][: len(f4)]})
    assert way_responses(P, CX, B, L, cab4, fld, F, trusted_limit=3000.0).meta["puck_blend"] == (3000.0, f4[-1])
    f = F[F <= 3000.0]
    assert way_responses(P, CX, B, L, fld, None, f, trusted_limit=3000.0).meta["puck_blend"] is None


def test_vents_fade_out_across_the_puck_blend():
    from speaker.system import _log_ramp
    silent = _piston_field(F)
    loud = _piston_field(F)
    o = loud.obs
    for k, sg in ((2, 0.3), (3, -0.2)):   # direction-dependent vents
        loud.p[:, k, :] = loud.p[:, 0, :] * (0.5 + sg * o.dirs[:, 1])[None, :]
    tl = 3000.0
    w_s = way_responses(P, CX, B, L, silent, silent, F, trusted_limit=tl)
    w_l = way_responses(P, CX, B, L, loud, silent, F, trusted_limit=tl)   # same puck, so only the vents differ
    lo, hi = w_l.meta["puck_blend"]
    assert (lo, hi) == (tl, 2 * tl)
    n = len(w_l.meta["u_lf"])   # grid points up to the 12 kHz hold
    f, vent = F[:n], (coax_section(w_l) - coax_section(w_s))[:n]
    # The full vent term, from the cabinet field at full strength: u_vent per vent times both vent fields.
    lump = simulate(CX, cardioid_params(P, L, CX), f, 2.83)
    u_vent = lump.u_vent / lump.u_cone * w_l.meta["u_lf"] / 2
    full = u_vent[:, None] * (loud.p[:n, 2, :] + loud.p[:n, 3, :])
    m = f <= tl
    assert np.allclose(vent[m], full[m], rtol=1e-9, atol=0)
    j = (f > lo) & (f < hi)
    assert j.sum() > 10
    assert np.allclose(vent[j], (1 - _log_ramp(f[j], lo, hi))[:, None] * full[j], rtol=1e-9, atol=0)
    assert np.all(vent[f >= hi] == 0)


def test_no_blend_when_the_grid_stops_at_the_trusted_limit():
    fld = _piston_field(F)
    f = F[F <= 3000.0]
    w = way_responses(P, CX, B, L, fld, fld, f, trusted_limit=3000.0)
    assert w.meta["puck_blend"] is None
    assert w.meta["hf_reference"] == "puck"


def test_puck_blend_stops_at_the_hold():
    fld = _piston_field(F)   # cabinet field to 16 kHz
    w = way_responses(P, CX, B, L, fld, fld, F, trusted_limit=8000.0)
    assert w.meta["puck_blend"] == (8000.0, 12000.0)


def test_puck_blend_narrower_than_a_third_octave_is_refused():
    fld = _piston_field(F)
    f32 = F[F <= 3200.0]
    cab = SourceField(f32, fld.names, fld.obs, fld.p[: len(f32)], {"axis_1m": fld.meta["axis_1m"][: len(f32)]})
    with pytest.raises(ValueError, match="blend"):
        way_responses(P, CX, B, L, cab, fld, F, trusted_limit=3000.0)
