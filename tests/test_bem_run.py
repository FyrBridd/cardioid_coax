import sys
from dataclasses import replace

import numpy as np, pytest
from speaker.bem import analytic
from speaker.bem.mesh import box_piston_mesh
from fakes import default_design
from speaker.bem.run import freq_grid, observation_set, cabinet_field, code_version
from speaker.bem.solver import solve_band
from speaker.constants import C, RHO
from speaker.drivers import load_driver
from speaker.layout import compute_layout
from speaker.params import load_params, DEFAULT_PATH
import speaker.bem.run as run

def test_observation_set():
    o = observation_set(2.0)
    assert len(o.h_idx) == 72 and len(o.v_idx) == 72
    assert np.allclose(o.dirs[o.h_idx[o.h_angles == 0][0]], [0, -1, 0])
    assert np.allclose(o.dirs[o.h_idx[o.h_angles == 90][0]], [1, 0, 0])
    assert o.sphere_w.sum() == pytest.approx(1.0)
    # Sphere weights integrate cos^2 (dipole power) to 1/3.
    d = o.dirs[o.sphere_idx]
    assert (o.sphere_w * d[:, 1] ** 2).sum() == pytest.approx(1 / 3, abs=0.01)
    a, e = np.radians(23), np.radians(10)
    assert np.allclose(o.dirs[o.at23_idx], [np.sin(a), -np.cos(a), 0])
    assert o.at23_idx not in o.sphere_idx and o.at23_idx not in o.h_idx
    assert np.allclose(o.dirs[list(o.lw_extra_idx)],
                       [[np.sin(a) * np.cos(e), -np.cos(a) * np.cos(e), np.sin(e)],
                        [np.sin(a) * np.cos(e), -np.cos(a) * np.cos(e), -np.sin(e)]])

def test_observation_set_conventions():
    o = observation_set(2.0)
    assert np.allclose(np.linalg.norm(o.dirs, axis=1), 1)
    assert np.allclose(o.dirs[o.v_idx[o.v_angles == 90][0]], [0, 0, 1])   # positive φ is up
    assert len(o.sphere_idx) == 36 * 17 + 2 and len(set(o.sphere_idx)) == len(o.sphere_idx)

def test_freq_grid():
    f = freq_grid(50, 20000, 12)
    assert f[0] == 50 and f[-1] <= 20000 and np.allclose(np.diff(np.log2(f)), 1 / 12)

def _fake_solver(calls, value=1j):
    def fake(mesh, freqs, points, **kw):
        calls.append((len(mesh.triangles), tuple(freqs)))
        return np.full((len(freqs), len(mesh.tag_names), points.shape[1]), value)
    return fake

def test_banding_and_stored_phase_convention(monkeypatch):
    P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
    calls = []
    monkeypatch.setattr(run, "solve_band", _fake_solver(calls, 1j))   # +90° in solve_band's e^{-iωt}
    f = np.array([100.0, 250.0, 260.0, 1500.0, 12000.0])
    fld = cabinet_field(P, CX, B, f, per_wavelength=1.0, h_max=0.2, cache_dir=None)
    assert [c[1] for c in calls] == [(100.0, 250.0), (260.0,), (1500.0,), (12000.0,)]
    bands = fld.meta["bands"]
    assert [b["edge"] for b in bands] == [250.0, 500.0, 2000.0, 12000.0]
    assert [b["h"] for b in bands] == pytest.approx([0.2, 0.2, C / 2000, C / 12000])
    assert [b["triangles"] for b in bands] == [c[0] for c in calls]
    assert all(b["mesh_s"] >= 0 and b["solve_s"] >= 0 for b in bands)
    assert fld.names == ("coax_lf", "coax_hf", "vent_left", "vent_right", "bass_lower", "bass_upper")
    assert fld.p.shape == (5, 6, len(fld.obs.dirs))
    assert np.all(fld.p == -1j)   # stored e^{+jωt}: the conjugate of solve_band
    with pytest.raises(ValueError, match="12000"):
        cabinet_field(P, CX, B, np.array([12001.0]), cache_dir=None)

def test_cache_key_covers_inputs_and_code_version(tmp_path, monkeypatch):
    import dataclasses
    P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
    calls = []
    monkeypatch.setattr(run, "solve_band", _fake_solver(calls))
    kw = dict(h_max=0.1, cache_dir=tmp_path)
    cabinet_field(P, CX, B, np.array([60.0]), **kw)
    cabinet_field(P, CX, B, np.array([60.0]), **kw)
    assert len(calls) == 1
    cabinet_field(P, CX, B, np.array([70.0]), **kw)
    cabinet_field(P, CX, B, np.array([60.0]), **{**kw, "h_max": 0.09})
    cabinet_field(P, CX, B, np.array([60.0]), per_wavelength=5.0, **kw)
    cabinet_field(P, CX, load_driver("bms_10c262"), np.array([60.0]), **kw)   # bass driver changed
    cabinet_field(P, CX, B, np.array([60.0]), h_coax=0.008, **kw)
    cabinet_field(P, CX, B, np.array([60.0]), h_vent=0.015, **kw)
    # Geometry enters the key through the meshes.
    P2 = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, coax_centre_z=P.cabinet.coax_centre_z + 0.01))
    cabinet_field(P2, CX, B, np.array([60.0]), **kw)
    monkeypatch.setattr(run, "code_version", lambda: "other")
    cabinet_field(P, CX, B, np.array([60.0]), **kw)
    assert len(calls) == 9


def test_code_version_includes_the_bempp_version(monkeypatch):
    import importlib.metadata
    base = run.code_version()
    real = importlib.metadata.version
    monkeypatch.setattr(importlib.metadata, "version", lambda d: "0.0" if d == "bempp-cl" else real(d))
    assert run.code_version() != base
    monkeypatch.setattr(importlib.metadata, "version", real)
    assert run.code_version() == base


@pytest.fixture
def counting_solver(monkeypatch):
    class Counter:
        calls = 0

    def fake(mesh, f, points, **kw):
        Counter.calls += 1
        return np.zeros((len(f), len(mesh.tag_names), points.shape[1]), complex)
    monkeypatch.setattr(run, "solve_band", fake)
    return Counter


F = np.array([100.0])


def test_cache_ignores_non_geometry_params(counting_solver, tmp_path):
    p, coax, bass = default_design()
    cabinet_field(p, coax, bass, F, h_max=0.2, cache_dir=tmp_path)
    q = replace(p, vents=replace(p.vents, specific_resistance=400.0),
                drive=replace(p.drive, amp_v_rms=10.0),
                cabinet=replace(p.cabinet, coax_chamber_net_volume=0.012))
    cabinet_field(q, coax, bass, F, h_max=0.2, cache_dir=tmp_path)
    assert counting_solver.calls == 1


def test_cache_resolves_when_geometry_moves(counting_solver, tmp_path):
    p, coax, bass = default_design()
    cabinet_field(p, coax, bass, F, h_max=0.2, cache_dir=tmp_path)
    q = replace(p, vents=replace(p.vents, centre_from_front=0.240))
    cabinet_field(q, coax, bass, F, h_max=0.2, cache_dir=tmp_path)
    assert counting_solver.calls == 2


def test_code_version_covers_only_the_solver(monkeypatch, tmp_path):
    before = code_version()
    fake = tmp_path / "layout.py"; fake.write_text("# changed\n")
    monkeypatch.setattr(sys.modules["speaker.layout"], "__file__", str(fake))
    assert code_version() == before
    fake2 = tmp_path / "solver.py"; fake2.write_text("# changed\n")
    monkeypatch.setattr(sys.modules["speaker.bem.solver"], "__file__", str(fake2))
    assert code_version() != before


def test_cabinet_field_far_points_on_h0_and_h180(monkeypatch):
    P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
    seen = []

    def fake(mesh, freqs, points, **kw):   # e^{-iωt} value i * (y of the point)
        seen.append(points)
        return np.broadcast_to(1j * points[1], (len(freqs), len(mesh.tag_names), points.shape[1])).copy()
    monkeypatch.setattr(run, "solve_band", fake)
    f = np.array([150.0, 300.0])
    fld = cabinet_field(P, CX, B, f, h_max=0.2, cache_dir=None)
    o, (cx, cz) = fld.obs, compute_layout(P, CX, B).coax_centre
    assert fld.p.shape == (2, 6, len(o.dirs))
    far = fld.meta["far_h0_h180"]
    assert far.shape == (2, 6, 2)
    pts = seen[0][:, -2:].T
    assert np.allclose(pts, [[cx, -run.FAR_R, cz], [cx, run.FAR_R, cz]])
    assert np.allclose(far[..., 0], 1j * run.FAR_R) and np.allclose(far[..., 1], -1j * run.FAR_R)  # conjugated
    assert np.allclose(fld.p, np.conj(1j * seen[0][1, :-2]))


def test_cabinet_field_mesh_sizes_reach_the_mesher(monkeypatch):
    P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
    calls = []
    monkeypatch.setattr(run, "solve_band", _fake_solver(calls))
    for kw in ({}, {"h_coax": 0.006}, {"h_vent": 0.01}):
        cabinet_field(P, CX, B, np.array([100.0]), h_max=0.05, cache_dir=None, **kw)
    assert calls[1][0] > calls[0][0] and calls[2][0] > calls[0][0]   # finer regions, more triangles


@pytest.mark.bem
def test_cabinet_radiates_as_monopole_at_low_frequency(tmp_path):
    P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
    fld = cabinet_field(P, CX, B, np.array([20.0]),
                        h_max=0.08, cache_dir=tmp_path)   # ka ~ 0.2 for the whole cabinet
    # Each source against a monopole at its own position. The bass drivers sit up to 0.75 m from the
    # reference point, so at 2 m the 1/r spread alone is +4/-3 dB; compare per source position instead.
    L = compute_layout(P, CX, B)
    (cx, cz), w = L.coax_centre, 2 * np.pi * 20
    zv = np.mean([z for _, z, _, _ in L.vent_slots])
    zb = sorted(z for _, z in L.bass_centres)
    src = [(cx, 0, cz), (cx, 0, cz), (0, P.vents.centre_from_front, zv),
           (P.cabinet.width, P.vents.centre_from_front, zv), (cx, 0, zb[0]), (cx, 0, zb[1])]
    o = fld.obs
    pts = np.array([cx, 0, cz]) + o.r * o.dirs
    for s, x in enumerate(src):
        r = np.linalg.norm(pts - np.array(x), axis=1)
        mono = 1j * w * RHO * np.exp(-1j * w / C * r) / (4 * np.pi * r)   # e^{+jωt}, unit volume velocity
        q = fld.p[0, s] / mono
        # Power over the sphere pins the unit-volume-velocity normalisation (Review Focus 2).
        sw, si = o.sphere_w, o.sphere_idx
        power = (sw * np.abs(fld.p[0, s, si]) ** 2).sum() / (sw * np.abs(mono[si]) ** 2).sum()
        assert 10 * np.log10(power) == pytest.approx(0, abs=0.1), fld.names[s]
        # Per direction: at kr = 0.73 the flow round the body adds a near-field dipole of about
        # ±0.9 dB at 2 m (±0.08 dB at 20 m), so the bound is looser than the power check.
        assert np.all(np.abs(20 * np.log10(np.abs(q))) < 1.5), fld.names[s]
        # Phase: stored e^{+jωt}; the unconjugated e^{-iωt} field would sit ~2(90° - kr) = 96° away.
        assert np.all(np.abs(np.degrees(np.angle(q))) < 10), fld.names[s]

@pytest.mark.bem
def test_small_box_baffle_step():
    m = box_piston_mesh((0.2, 0.15, 0.3), 0.015, h=0.012)
    f = np.array([100.0, 2000.0, 3000.0, 4000.0])
    p = np.abs(solve_band(m, f, np.array([[0.0], [-2.0], [0.0]]) + np.array([[0.1], [0], [0.15]]))[:, 0, 0])
    step = 20 * np.log10(p[1:].mean() / p[0]) - 20 * np.log10(f[1:].mean() / f[0])  # remove the ω rise
    assert step == pytest.approx(6.0, abs=1.5)

@pytest.mark.bem
def test_cache_hit_skips_solver(tmp_path, monkeypatch):
    P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
    a = cabinet_field(P, CX, B, np.array([60.0]), h_max=0.1, cache_dir=tmp_path)
    import speaker.bem.run as run
    monkeypatch.setattr(run, "solve_band", lambda *a, **k: (_ for _ in ()).throw(AssertionError("solver called")))
    b = cabinet_field(P, CX, B, np.array([60.0]), h_max=0.1, cache_dir=tmp_path)
    assert np.array_equal(a.p, b.p)
    import dataclasses
    P2 = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, coax_centre_z=P.cabinet.coax_centre_z + 0.01))
    with pytest.raises(AssertionError, match="solver called"):
        cabinet_field(P2, CX, B, np.array([60.0]), h_max=0.1, cache_dir=tmp_path)
