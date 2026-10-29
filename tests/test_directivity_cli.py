import json
from pathlib import Path
import pytest
from speaker.cli import main

@pytest.mark.bem
def test_directivity_fast_end_to_end(tmp_path):
    assert main(["directivity", "--fast", "--out", str(tmp_path)]) == 0
    html = (tmp_path / "report.html").read_text()
    for s in ("Round 0", "unverified", "cardioid_hard", "Listening window", "Confidence"):
        assert s in html
    keys = {d["key"] for d in json.loads((tmp_path / "scores.json").read_text())}
    assert {"cardioid_hard", "cardioid_aim", "window_spread", "rear_hemisphere", "whole_rear"} <= keys
    for png in ("polar_h", "polar_v", "spinorama", "window_family", "rear_rejection"):
        assert (tmp_path / "plots" / f"{png}.png").stat().st_size > 0

def test_directivity_rejects_bad_params(tmp_path, capsys):
    text = Path("params/default.toml").read_text().replace("depth_mm = 350", "depth_mm = 351")
    pf = tmp_path / "p.toml"; pf.write_text(text)
    assert main(["directivity", "--params", str(pf), "--out", str(tmp_path / "o")]) == 2
    assert capsys.readouterr().err.startswith("error: ")


def _az(d):
    import numpy as np
    return np.degrees(np.arctan2(d[:, 0], -d[:, 1]))


def _fake_cabinet_field(calls, hot=(1100.0, 1900.0), at=None, gain_db=3.0, null=None):
    """Analytic stand-in for cabinet_field: the 8-per-wavelength mesh is gain_db hot over `hot` (default
    1.1-1.9 kHz, a range outside every crossover level-matching band; None: converged everywhere), in
    every direction or only at H angle `at`. null = (at, coarse_db, fine_db) also cuts H angle `at` by
    coarse_db / fine_db on the 6 / 8 per-wavelength mesh, over `hot`."""
    import numpy as np
    from speaker.bem import analytic
    from speaker.bem.run import BAND_EDGES, CAB_H_COAX, CAB_H_VENT, SourceField, observation_set
    from speaker.bem.mesh import CABINET_TAGS

    def fake(p, coax, bass, f, *, per_wavelength, h_max, cache_dir,
             h_coax=CAB_H_COAX, h_vent=CAB_H_VENT):
        calls.append((per_wavelength, float(f[0]), float(f[-1]), h_max, h_coax, h_vent, tuple(f)))
        o = observation_set(2.0)
        pp = np.array([analytic.baffled_piston_on_axis(x, 0.1, 2.0) for x in f])
        pr = np.broadcast_to(pp[:, None, None], (len(f), len(CABINET_TAGS), len(o.dirs))).copy()
        if per_wavelength == 8.0 and hot is not None:
            fm = (f >= hot[0]) & (f <= hot[1])
            dm = np.ones(len(o.dirs), bool) if at is None else np.abs(_az(o.dirs) - at) < 0.5
            pr[np.ix_(fm, np.arange(len(CABINET_TAGS)), dm)] *= 10 ** (gain_db / 20)
        if null is not None:
            fm, dm = (f >= hot[0]) & (f <= hot[1]), np.abs(_az(o.dirs) - null[0]) < 0.5
            pr[np.ix_(fm, np.arange(len(CABINET_TAGS)), dm)] *= 10 ** (-null[1 if per_wavelength == 6.0 else 2] / 20)
        edge = BAND_EDGES[int(np.searchsorted(BAND_EDGES, f[0]))]
        band = {"edge": edge, "h": 0.05, "freqs": len(f), "nodes": 100, "triangles": 196, "mesh_s": 0.1, "solve_s": 1.0}
        return SourceField(np.asarray(f, float), CABINET_TAGS, o, pr, {"bands": [band]})
    return fake


def _small_dofs(p, coax, bass, h, h_coax, h_vent):
    return 1000


def test_mesh_check_limit_is_the_last_band_within_half_a_db(tmp_path, monkeypatch):
    import speaker.directivity_report as dr
    calls = []
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field(calls))
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)
    out = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    assert out["per_wavelength"] == [6.0, 8.0] and out["capped_by_memory"] is False
    assert out["trusted_limit_hz"] == pytest.approx(1000 * 2 ** (1 / 6))
    dev = out["band_max_dev_db"]
    assert dev["1000"] < 0.01 and dev["1260"] > 0.5
    # Every band is solved at both densities, one call per band and density.
    assert set(out["elements"]) == {"250", "500", "1000", "2000", "3000", "4000"}
    assert out["elements"]["2000"]["fine"]["nodes"] == 100 and "solve_s" in out["elements"]["2000"]["coarse"]
    assert sorted(c[0] for c in calls) == [6.0] * 6 + [8.0] * 6


def test_mesh_check_refines_every_mesh_size(tmp_path, monkeypatch):
    """The finer run scales h_max, h_coax and h_vent by 6/8, so even the bands h_max caps are refined."""
    import speaker.directivity_report as dr
    from speaker.bem.run import CAB_H_COAX, CAB_H_VENT
    calls = []
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field(calls, hot=None))
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)
    out = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    coarse = {c[3:6] for c in calls if c[0] == 6.0}
    fine = {c[3:6] for c in calls if c[0] == 8.0}
    assert coarse == {(dr.CAB_H_MAX, CAB_H_COAX, CAB_H_VENT)}   # the full run's meshes, so its cache is shared
    assert len(fine) == 1 and list(fine)[0] == pytest.approx((dr.CAB_H_MAX * 0.75, CAB_H_COAX * 0.75, CAB_H_VENT * 0.75))
    assert out["mesh_sizes_m"]["fine"] == pytest.approx({"h_max": 0.0375, "h_coax": 0.0075, "h_vent": 0.015})
    assert "rear rejection" in out["compared"] and "40-90" in out["compared"]


def test_fine_meshes_have_more_nodes_in_every_band():
    """gmsh only (no BEM): the finer mesh of every mesh-check band up to 4 kHz has more nodes."""
    import speaker.directivity_report as dr
    d = dr.load_design(None)
    for edge in (250.0, 500.0, 1000.0, 2000.0, 3000.0, 4000.0):
        n = [dr._fine_dofs(d.p, d.coax, d.bass, *dr._mesh_sizes(edge, pw)) for pw in (dr.COARSE, dr.FINE)]
        assert n[1] > 1.3 * n[0], (edge, n)


@pytest.mark.parametrize("at, hot, limit", [
    (125.0, (200.0, 300.0), 1000 * 2 ** (-7.5 / 3)),  # coax rear rejection: H 125 in the 150-500 Hz window
    (125.0, (1100.0, 1900.0), None),                  # ... but not outside it
    (-70.0, (1100.0, 1900.0), 1000 * 2 ** (1 / 6)),   # H -40...-90 side
])
def test_mesh_check_compares_cardioid_and_off_axis_curves(tmp_path, monkeypatch, at, hot, limit):
    import speaker.directivity_report as dr
    from speaker.bem.run import freq_grid
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field([], hot=hot, at=at, gain_db=6.0))
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)
    out = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    assert out["trusted_limit_hz"] == pytest.approx(limit or float(freq_grid(50, 4000, 12)[-1]))


@pytest.mark.parametrize("coarse_db, fine_db, limit", [
    (58.0, 55.7, None),                         # a rear null deeper than the cap reads 30 dB on both meshes
    (20.0, 23.0, 1000 * 2 ** (-7.5 / 3)),       # below the cap the rejection is compared as is
])
def test_mesh_check_caps_the_rear_rejection(tmp_path, monkeypatch, coarse_db, fine_db, limit):
    import speaker.directivity_report as dr
    from speaker.bem.run import freq_grid
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field([], hot=(200.0, 300.0), gain_db=0.0,
                                                                  null=(165.0, coarse_db, fine_db)))
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)
    out = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    assert out["rejection_cap_db"] == dr.REJ_CAP_DB == 30.0
    assert "capped at 30 dB" in out["compared"]
    assert out["trusted_limit_hz"] == pytest.approx(limit or float(freq_grid(50, 4000, 12)[-1]))


def test_mesh_check_stops_before_a_fine_mesh_over_the_memory_budget(tmp_path, monkeypatch):
    import speaker.directivity_report as dr
    calls, asked = [], []
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field(calls))

    def dofs(p, coax, bass, h, h_coax, h_vent):
        asked.append((h, h_coax, h_vent))
        return 10**6 if h < 0.03 else 1000
    monkeypatch.setattr(dr, "_fine_dofs", dofs)
    out = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    assert out["capped_by_memory"] is True
    assert out["trusted_limit_hz"] == pytest.approx(50 * 2 ** (51 / 12))  # capped at the top solved frequency
    assert set(out["elements"]) == {"250", "500", "1000"}
    assert max(c[2] for c in calls) <= 1000.0
    assert all(a[1:] == pytest.approx((0.0075, 0.015)) for a in asked)   # the guard sizes the real fine mesh


def test_memory_cap_at_the_4k_fine_band_keeps_the_3k_band(tmp_path, monkeypatch):
    """The real case: only the <= 4 kHz fine mesh is over budget. The check verifies through the <= 3 kHz
    band, its limit is the top solved grid point there, and a full run accepts that limit."""
    import numpy as np
    import speaker.directivity_report as dr
    calls = []
    _stub_fields(monkeypatch, [])
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field(calls, hot=None))
    monkeypatch.setattr(dr, "_fine_dofs", lambda p, coax, bass, h, h_coax, h_vent: 19370 if h < 0.012 else 7000)
    mc = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    g = dr.grid(*dr.FULL_GRID)
    top3k = float(g[g <= 3000.0][-1])   # 2850.9 Hz
    assert mc["capped_by_memory"] is True and "3000" in mc["elements"] and "4000" not in mc["elements"]
    assert max(c[2] for c in calls) == pytest.approx(top3k)
    assert mc["trusted_limit_hz"] == pytest.approx(top3k)
    _stub_fields(monkeypatch, [])   # the full run's fields (the mesh check's fake is replaced)
    page = dr.run_directivity(None, tmp_path).read_text()
    assert f"{top3k:.0f} Hz" in page and "unverified" not in page
    cab_top = float(g[g <= 4000.0][-1])
    assert f"{top3k:.0f} Hz to {cab_top:.0f} Hz" in page   # the puck blend sits above the verified limit
    assert "above the trusted limit" in page
    assert f"The cabinet fields in the blend above {top3k:.0f} Hz are not mesh-checked" in page
    gains = json.loads((tmp_path / "gains.json").read_text())
    assert set(gains) == {"bass", "coax_lf", "coax_hf"} and all(isinstance(v, float) for v in gains.values())
    trusted = [s for s in json.loads((tmp_path / "scores.json").read_text()) if s["confidence"] == "trusted"]
    assert all(s["band"][1] <= top3k + 1e-6 for s in trusted)


def test_output_grids_hold_the_score_band_edges():
    import numpy as np
    import speaker.directivity_report as dr
    from speaker.bem.run import freq_grid
    full, fast = dr.grid(*dr.FULL_GRID), dr.grid(*dr.FAST_GRID)
    for g, spec in ((full, dr.FULL_GRID), (fast, dr.FAST_GRID)):
        assert np.all(np.diff(g) > 0) and set(freq_grid(*spec)) <= set(g)
        assert g[0] == spec[0] and g[-1] <= spec[1]
    assert {150.0, 300.0, 500.0, 2000.0, 5000.0} <= set(full)
    assert {150.0, 300.0, 500.0} <= set(fast) and fast[-1] <= 1000.0


def test_mesh_check_grid_is_the_full_run_cabinet_grid(tmp_path, monkeypatch):
    import numpy as np
    import speaker.directivity_report as dr
    calls = []
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field(calls, hot=None))
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)
    dr.mesh_check(None, tmp_path, f_max=4000.0)
    g = dr.grid(*dr.FULL_GRID)
    want = [tuple(g[(g <= 4000.0) & (np.searchsorted(dr.BAND_EDGES, g) == b)]) for b in range(6)]
    assert [c[6] for c in calls if c[0] == 6.0] == want


def test_converged_mesh_check_feeds_a_usable_limit_to_the_run(tmp_path, monkeypatch):
    import speaker.directivity_report as dr
    from speaker.bem.run import freq_grid
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field([], hot=None))
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)
    mc = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    top = float(freq_grid(50, 4000, 12)[-1])  # 3805.5 Hz: the last band's upper edge is past the grid
    assert mc["trusted_limit_hz"] == pytest.approx(top)
    # The full run's cabinet grid is the same, so the limit fits inside it.
    f = dr.grid(*dr.FULL_GRID)
    cab_top = float(f[f <= 4000][-1])
    tl, text = dr._trusted_limit(mc, cab_top=cab_top, f_top=cab_top, puck_lo=None)
    assert tl == pytest.approx(top) and "unverified" not in text and "capped" not in text
    # With the puck, the cabinet field must reach 1/3 octave past the limit for the blend above it:
    # a limit at the cabinet top is capped there, and says so.
    tl, text = dr._trusted_limit(mc, cab_top=cab_top, f_top=float(f[-1]), puck_lo=500.0)
    assert tl == pytest.approx(cab_top / 2 ** (1 / 3)) and "capped" in text and "blend" in text
    # A limit with room for the blend is kept as it is.
    tl, text = dr._trusted_limit({**mc, "trusted_limit_hz": 2850.9}, cab_top=cab_top, f_top=20000.0, puck_lo=500.0)
    assert tl == 2850.9 and "capped" not in text
    # An older check that reached past the cabinet field is capped too.
    tl, text = dr._trusted_limit({**mc, "trusted_limit_hz": 4489.8}, cab_top=top, f_top=20000.0, puck_lo=500.0)
    assert tl == pytest.approx(top / 2 ** (1 / 3)) and "capped" in text
    # A limit below the puck's lowest frequency (where the blend starts) is an error, not a crash.
    with pytest.raises(ValueError, match="trusted 3D limit"):
        dr._trusted_limit({**mc, "trusted_limit_hz": 480.0}, cab_top=top, f_top=20000.0, puck_lo=503.97)
    assert dr._trusted_limit(None, cab_top=898.0, f_top=898.0, puck_lo=None) == (3000.0, "3000 Hz (unverified)")


def test_mesh_check_json_is_used_only_by_a_matching_run(tmp_path, monkeypatch):
    import hashlib
    import speaker.directivity_report as dr
    import speaker.bem.run as bem_run
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field([], hot=None))
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)
    mc = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    d = dr.load_design(None)
    assert mc["params_sha256"] == hashlib.sha256(d.path.read_bytes()).hexdigest()
    assert mc["code_version"] == bem_run.code_version()
    assert mc["per_wavelength"] == [6.0, 8.0] and mc["h_max_m"] == dr.CAB_H_MAX
    full = dict(h_max=dr.CAB_H_MAX, per_wavelength=6.0)
    assert dr._usable_mesh_check(mc, d, **full) == (mc, "")
    assert dr._usable_mesh_check(None, d, **full)[0] is None
    for bad, kw, word in (({"params_sha256": "0" * 64}, full, "params"),
                          ({}, dict(h_max=dr.FAST_H_MAX, per_wavelength=6.0), "h_max"),
                          ({}, dict(h_max=dr.CAB_H_MAX, per_wavelength=5.0), "per_wavelength"),
                          ({"per_wavelength": [6.0, 10.0]}, full, "per_wavelength"),
                          ({"rejection_cap_db": None}, full, "rear rejection cap None vs 30.0")):
        got, why = dr._usable_mesh_check({**mc, **bad}, d, **kw)
        assert got is None and word in why
    monkeypatch.setattr(bem_run, "code_version", lambda: "other")
    got, why = dr._usable_mesh_check(mc, d, **full)
    assert got is None and "code" in why


def test_mesh_check_with_no_converged_band_is_an_error(tmp_path, monkeypatch, capsys):
    import speaker.directivity_report as dr
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field([], hot=(0.0, 1e9)))
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)
    assert main(["directivity", "--mesh-check", "--out", str(tmp_path)]) == 2
    err = capsys.readouterr().err
    assert "error: " in err and "no 1/3-octave band" in err


def _stub_fields(monkeypatch, puck_calls):
    """Analytic stand-ins for cabinet_field and puck_field on whatever grid the run asks for (no BEM).

    A forward-leaning coax, rear-leaning vents of opposite sign, far-field points and a puck on-axis curve.
    """
    import numpy as np
    import speaker.directivity_report as dr
    from speaker.bem import analytic
    from speaker.bem.mesh import CABINET_TAGS, PUCK_TAGS
    from speaker.bem.run import BAND_EDGES, CAB_H_COAX, CAB_H_VENT, SourceField, freq_grid, observation_set
    o = observation_set(2.0)
    i0, i180 = o.h_idx[o.h_angles == 0][0], o.h_idx[o.h_angles == -180][0]

    def make(f, names):
        f = np.asarray(f, float)
        pp = np.array([analytic.baffled_piston_on_axis(x, 0.1, 2.0) for x in f]).conj()
        p = (pp[:, None, None] * (0.6 - 0.4 * o.dirs[:, 1])[None, None, :] * np.ones((1, len(names), 1))).astype(complex)
        if len(names) > 2:
            p[:, 2:4, :] *= -0.3 * (1 + 0.5 * o.dirs[None, :, 1])
        edge = BAND_EDGES[int(np.searchsorted(BAND_EDGES, f[0]))]
        band = {"edge": edge, "h": 0.05, "freqs": len(f), "nodes": 1, "triangles": 2, "mesh_s": 0.0, "solve_s": 0.0}
        return f, p, {"bands": [band]}

    def cab(p, coax, bass, f, *, per_wavelength, h_max, cache_dir, h_coax=CAB_H_COAX, h_vent=CAB_H_VENT):
        f, pr, meta = make(f, CABINET_TAGS)
        return SourceField(f, CABINET_TAGS, o, pr, {**meta, "far_h0_h180": pr[:, :, [i0, i180]] * 2 / 50})

    def puck(p, coax, f, *, per_wavelength, cache_dir):
        puck_calls.append(np.asarray(f, float))
        f, pr, meta = make(f, PUCK_TAGS)
        ax = np.array([analytic.baffled_piston_on_axis(x, coax.hf_exit_d / 2, 1.0) for x in f]).conj()
        return SourceField(f, PUCK_TAGS, o, pr, {**meta, "axis_1m": np.stack([ax, ax], 1)})

    def recess(coax, f=None, *, per_wavelength=6.0, cache_dir=None):   # a recess about 1 dB hot near 1 kHz
        f = freq_grid(300.0, 4000.0, 6) if f is None else np.asarray(f, float)
        return f, 10 ** (np.exp(-np.log2(f / 1000.0) ** 2) / 20) + 0j
    monkeypatch.setattr(dr, "cabinet_field", cab)
    monkeypatch.setattr(dr, "puck_field", puck)
    monkeypatch.setattr(dr, "recess_ratio", recess)
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)


def test_full_run_end_to_end_on_stub_fields(tmp_path, monkeypatch):
    """The full (not fast) run on the real full grid, BEM fields replaced by analytic stand-ins."""
    import numpy as np
    import speaker.directivity_report as dr
    puck_calls = []
    _stub_fields(monkeypatch, puck_calls)
    mc = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    tl = mc["trusted_limit_hz"]
    path = dr.run_directivity(None, tmp_path)
    page = path.read_text()
    fp = np.concatenate(puck_calls)
    assert fp[0] == 500.0 and fp[-1] == 12000.0   # the puck reaches the hold frequency
    assert f"{tl:.0f} Hz" in page and "unverified" not in page
    scores = {(s["key"], s["confidence"]): s for s in json.loads((tmp_path / "scores.json").read_text())}
    assert {("presence_lift", "trusted"), ("presence_lift", "low"), ("cardioid_hard", "trusted"),
            ("cardioid_hard_far", "trusted"), ("window_spread", "low")} <= set(scores)
    assert scores[("cardioid_hard_far", "trusted")]["status"] == "info"
    assert 150.0 <= scores[("cardioid_hard", "trusted")]["detail"]["worst_hz"] <= 300.0
    # The trusted range is the cabinet alone; the blend sits in the low-confidence range above it.
    ws, ws_low = scores[("window_spread", "trusted")], scores[("window_spread", "low")]
    assert ws["detail"]["puck_fraction"] == 0.0 and ws["detail"]["puck_mean_weight"] == 0.0
    assert ws_low["detail"]["puck_fraction"] == 1.0 and 0 < ws_low["detail"]["puck_mean_weight"] < 1
    assert scores[("cardioid_hard", "trusted")]["detail"]["puck_fraction"] == 0.0
    # Confidence notes: the blend (I3), the far field (I5), de-embedding and estimated volumes (I6).
    run_tl = ws["band"][1]
    assert f"{run_tl:.0f} Hz to " in page and "above the trusted limit" in page and "round cylinder" in page
    assert "The gate is scored at 2 m (spec §4); far field is shown for reference." in page
    assert "infinite-baffle piston" in page and "double-count" in page and "puck's on-axis response" in page
    assert "filler volumes" in page and "§7.1.1" in page
    assert "rear rejection" in page   # the mesh check note says what was compared
    assert "times the recess ratio R" in page and "front-mounted on a large baffle" in page
    assert "R runs 300-3810 Hz" in page and "plots/lf_reference.png" in page
    for png in ("polar_h", "polar_v", "spinorama", "window_family", "rear_rejection", "whole_rear", "lf_reference",
                "crossover_23"):
        assert (tmp_path / "plots" / f"{png}.png").stat().st_size > 0
    cx = json.loads((tmp_path / "crossover.json").read_text())
    assert f"{cx['bass_coax']:g} Hz" in page and f"{cx['lf_hf']:g} Hz" in page and "acoustic LR4 target" in page
    assert "No EQ, no delays" not in page and "plots/crossover_23.png" in page
    assert "crossover.txt" in page and "vituixcad/README.txt" in page and "Passive-network warnings" in page
    vc = tmp_path / "vituixcad"
    assert (vc / "README.txt").exists()
    assert len(list(vc.glob("*.frd"))) == 5 * 143 and len(list(vc.glob("*.zma"))) == 4
    # A fast run in the same directory does not trust the full-run mesh check, and says why.
    page = dr.run_directivity(None, tmp_path, fast=True).read_text()
    assert "3000 Hz (unverified)" in page and "h_max" in page
    assert "No EQ, no delays" in page and "crossover_23" not in page
    assert "piston of the HF exit" in page
    assert "flat piston (fast run)" in page and "recess ratio" not in page


def test_trusted_limit_override_in_header(tmp_path, monkeypatch):
    _stub_fields(monkeypatch, [])
    assert main(["directivity", "--fast", "--trusted-limit", "900", "--out", str(tmp_path)]) == 0
    page = (tmp_path / "report.html").read_text()
    assert "900 Hz (given on the command line)" in page and "unverified" not in page
    trusted = [s for s in json.loads((tmp_path / "scores.json").read_text()) if s["confidence"] == "trusted"]
    assert trusted and all(s["band"][1] <= 900.0 for s in trusted)


def test_trusted_limit_override_is_still_capped_for_the_blend():
    import speaker.directivity_report as dr
    tl, text = dr._trusted_limit(None, cab_top=3805.5, f_top=20000.0, puck_lo=500.0, given=(4000.0, "from round 1"))
    assert tl == pytest.approx(3805.5 / 2 ** (1 / 3)) and text.startswith(f"{tl:.0f} Hz (from round 1)")
    assert "capped" in text and "unverified" not in text
    assert dr._trusted_limit({"trusted_limit_hz": 1000.0}, cab_top=3805.5, f_top=20000.0, puck_lo=500.0,
                             given=(2000.0, "given")) == (2000.0, "2000 Hz (given)")


def test_frequency_axis_shades_the_puck_blend():
    import matplotlib.pyplot as plt
    import numpy as np
    from speaker.directivity_plots import _freq_axis
    fig, ax = plt.subplots()
    _freq_axis(ax, np.geomspace(50, 20000, 50), 3805.0, "3805 Hz", blend=(1902.7, 3805.0))
    spans = [p for p in ax.patches if "puck" in p.get_label()]
    assert len(spans) == 1 and "1903-3805 Hz" in spans[0].get_label()
    x0, x1 = spans[0].get_x(), spans[0].get_x() + spans[0].get_width()
    assert (x0, x1) == pytest.approx((1902.7, 3805.0))
    plt.close(fig)


def test_full_run_writes_the_crossover_files(tmp_path, monkeypatch):
    import numpy as np
    from fakes import monopole_field, monopole_puck
    import speaker.directivity_report as dr
    from speaker.bem.run import freq_grid
    from speaker.params import DEFAULT_PATH
    from speaker.xo_fit import BASS_COAX_GRID, LF_HF_GRID

    def cab(p, coax, bass, f, *, per_wavelength, h_max, cache_dir, **kw):
        return monopole_field(p, coax, bass, f)

    def puck(p, coax, f, *, per_wavelength, cache_dir):
        return monopole_puck(p, coax, f)

    def recess(coax, f=None, *, per_wavelength=6.0, cache_dir=None):
        f = freq_grid(300.0, 4000.0, 6) if f is None else np.asarray(f, float)
        return f, np.ones(len(f), complex)

    monkeypatch.setattr(dr, "cabinet_field", cab)
    monkeypatch.setattr(dr, "puck_field", puck)
    monkeypatch.setattr(dr, "recess_ratio", recess)
    dr.run_directivity(DEFAULT_PATH, tmp_path, trusted_limit=(2851.0, "test"))
    xo = json.loads((tmp_path / "crossover.json").read_text())
    assert xo["bass_coax"] in BASS_COAX_GRID and xo["lf_hf"] in LF_HF_GRID
    txt = (tmp_path / "crossover.txt").read_text()
    assert "# fs = 48000" in txt and "# fs = 96000" in txt
    assert sum(l.split("\t")[0] in ("bass", "coax_lf", "coax_hf") and "delay_ms" in l
               for l in txt.splitlines()) == 3
    assert (tmp_path / "gains.json").exists()


def test_crossover_plot_tolerates_a_way_without_data(tmp_path):
    """The bass is zero above the trusted limit: its level is NaN or -inf there and must not reach log10."""
    import warnings
    import numpy as np
    from speaker import directivity_plots as plots
    f = np.geomspace(50, 4000, 120)
    bass = np.where(f < 1000, 80.0, -np.inf)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        plots.crossover_23(tmp_path / "x.png", f, {"bass": bass, "coax_lf": np.full(len(f), 80.0),
                                                   "sum": np.full(len(f), 82.0)}, 2851.0, "2851 Hz")
    assert (tmp_path / "x.png").stat().st_size > 0
