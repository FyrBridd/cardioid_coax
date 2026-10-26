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


def _fake_cabinet_field(calls, hot=(1100.0, 1900.0), at=None, gain_db=3.0):
    """Analytic stand-in for cabinet_field: the 8-per-wavelength mesh is gain_db hot over `hot` (default
    1.1-1.9 kHz, a range outside every crossover level-matching band; None: converged everywhere), in
    every direction or only at H angle `at`."""
    import numpy as np
    from speaker.bem import analytic
    from speaker.bem.run import BAND_EDGES, CAB_H_COAX, CAB_H_VENT, SourceField, observation_set
    from speaker.bem.mesh import CABINET_TAGS

    def fake(p, coax, bass, f, *, params_text, per_wavelength, h_max, cache_dir,
             h_coax=CAB_H_COAX, h_vent=CAB_H_VENT):
        calls.append((per_wavelength, float(f[0]), float(f[-1]), h_max, h_coax, h_vent, tuple(f)))
        o = observation_set(2.0)
        pp = np.array([analytic.baffled_piston_on_axis(x, 0.1, 2.0) for x in f])
        pr = np.broadcast_to(pp[:, None, None], (len(f), len(CABINET_TAGS), len(o.dirs))).copy()
        if per_wavelength == 8.0 and hot is not None:
            fm = (f >= hot[0]) & (f <= hot[1])
            dm = np.ones(len(o.dirs), bool) if at is None else np.abs(_az(o.dirs) - at) < 0.5
            pr[np.ix_(fm, np.arange(len(CABINET_TAGS)), dm)] *= 10 ** (gain_db / 20)
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
    assert set(out["elements"]) == {"250", "500", "1000", "2000", "4000"}
    assert out["elements"]["2000"]["fine"]["nodes"] == 100 and "solve_s" in out["elements"]["2000"]["coarse"]
    assert sorted(c[0] for c in calls) == [6.0] * 5 + [8.0] * 5


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
    d = dr._design(None)
    for edge in (250.0, 500.0, 1000.0, 2000.0, 4000.0):
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


def test_output_grids_hold_the_score_band_edges():
    import numpy as np
    import speaker.directivity_report as dr
    from speaker.bem.run import freq_grid
    full, fast = dr._grid(dr.FULL_GRID), dr._grid(dr.FAST_GRID)
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
    g = dr._grid(dr.FULL_GRID)
    want = [tuple(g[(g <= 4000.0) & (np.searchsorted(dr.BAND_EDGES, g) == b)]) for b in range(5)]
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
    f = dr._grid(dr.FULL_GRID)
    tl, text = dr._trusted_limit(mc, cab_top=float(f[f <= 4000][-1]), f_top=float(f[-1]), puck_lo=500.0)
    assert tl == pytest.approx(top) and "unverified" not in text
    # An older check that reached past the cabinet field is capped at its top, and says so.
    tl, text = dr._trusted_limit({**mc, "trusted_limit_hz": 4489.8}, cab_top=top, f_top=20000.0, puck_lo=500.0)
    assert tl == top and "capped" in text
    # A limit too low for the puck blend (tl/2 >= the puck's lowest frequency) is an error, not a crash.
    with pytest.raises(ValueError, match="trusted 3D limit"):
        dr._trusted_limit({**mc, "trusted_limit_hz": 1000.0}, cab_top=top, f_top=20000.0, puck_lo=503.97)
    assert dr._trusted_limit(None, cab_top=898.0, f_top=898.0, puck_lo=None) == (3000.0, "3000 Hz (unverified)")


def test_mesh_check_json_is_used_only_by_a_matching_run(tmp_path, monkeypatch):
    import hashlib
    import speaker.directivity_report as dr
    import speaker.bem.run as bem_run
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field([], hot=None))
    monkeypatch.setattr(dr, "_fine_dofs", _small_dofs)
    mc = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    d = dr._design(None)
    assert mc["params_sha256"] == hashlib.sha256(d.path.read_bytes()).hexdigest()
    assert mc["code_version"] == bem_run.code_version()
    assert mc["per_wavelength"] == [6.0, 8.0] and mc["h_max_m"] == dr.CAB_H_MAX
    full = dict(h_max=dr.CAB_H_MAX, per_wavelength=6.0)
    assert dr._usable_mesh_check(mc, d, **full) == (mc, "")
    assert dr._usable_mesh_check(None, d, **full)[0] is None
    for bad, kw, word in (({"params_sha256": "0" * 64}, full, "params"),
                          ({}, dict(h_max=dr.FAST_H_MAX, per_wavelength=6.0), "h_max"),
                          ({}, dict(h_max=dr.CAB_H_MAX, per_wavelength=5.0), "per_wavelength"),
                          ({"per_wavelength": [6.0, 10.0]}, full, "per_wavelength")):
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
    from speaker.bem.run import BAND_EDGES, CAB_H_COAX, CAB_H_VENT, SourceField, observation_set
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

    def cab(p, coax, bass, f, *, params_text, per_wavelength, h_max, cache_dir, h_coax=CAB_H_COAX, h_vent=CAB_H_VENT):
        f, pr, meta = make(f, CABINET_TAGS)
        return SourceField(f, CABINET_TAGS, o, pr, {**meta, "far_h0_h180": pr[:, :, [i0, i180]] * 2 / 50})

    def puck(p, coax, f, *, params_text, per_wavelength, cache_dir):
        puck_calls.append(np.asarray(f, float))
        f, pr, meta = make(f, PUCK_TAGS)
        ax = np.array([analytic.baffled_piston_on_axis(x, coax.hf_exit_d / 2, 1.0) for x in f]).conj()
        return SourceField(f, PUCK_TAGS, o, pr, {**meta, "axis_1m": np.stack([ax, ax], 1)})
    monkeypatch.setattr(dr, "cabinet_field", cab)
    monkeypatch.setattr(dr, "puck_field", puck)
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
    ws = scores[("window_spread", "trusted")]["detail"]
    assert 0 < ws["puck_fraction"] < 1 and 0 < ws["puck_mean_weight"] < ws["puck_fraction"]
    assert scores[("cardioid_hard", "trusted")]["detail"]["puck_fraction"] == 0.0
    # Confidence notes: the blend (I3), the far field (I5), de-embedding and estimated volumes (I6).
    assert f"{tl / 2:.0f} Hz to {tl:.0f} Hz" in page and "round cylinder" in page
    assert "The gate is scored at 2 m (spec §4); far field is shown for reference." in page
    assert "infinite-baffle piston" in page and "double-count" in page and "puck's on-axis response" in page
    assert "filler volumes" in page and "§7.1.1" in page
    assert "rear rejection" in page   # the mesh check note says what was compared
    for png in ("polar_h", "polar_v", "spinorama", "window_family", "rear_rejection", "whole_rear"):
        assert (tmp_path / "plots" / f"{png}.png").stat().st_size > 0
    # A fast run in the same directory does not trust the full-run mesh check, and says why.
    page = dr.run_directivity(None, tmp_path, fast=True).read_text()
    assert "3000 Hz (unverified)" in page and "h_max" in page
    assert "piston of the HF exit" in page


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
