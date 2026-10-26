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


def _fake_cabinet_field(calls, hot=(1100.0, 1900.0)):
    """Analytic stand-in for cabinet_field: the 8-per-wavelength mesh is 3 dB hot over `hot` (default
    1.1-1.9 kHz, a range outside every crossover level-matching band; None: converged everywhere)."""
    import numpy as np
    from speaker.bem import analytic
    from speaker.bem.run import BAND_EDGES, SourceField, observation_set
    from speaker.bem.mesh import CABINET_TAGS

    def fake(p, coax, bass, f, *, params_text, per_wavelength, h_max, cache_dir):
        calls.append((per_wavelength, float(f[0]), float(f[-1])))
        o = observation_set(2.0)
        pp = np.array([analytic.baffled_piston_on_axis(x, 0.1, 2.0) for x in f])
        if per_wavelength == 8.0 and hot is not None:
            pp = pp * np.where((f >= hot[0]) & (f <= hot[1]), 10 ** (3 / 20), 1.0)
        pr = np.broadcast_to(pp[:, None, None], (len(f), len(CABINET_TAGS), len(o.dirs))).copy()
        edge = BAND_EDGES[int(np.searchsorted(BAND_EDGES, f[0]))]
        band = {"edge": edge, "h": 0.05, "freqs": len(f), "nodes": 100, "triangles": 196, "mesh_s": 0.1, "solve_s": 1.0}
        return SourceField(np.asarray(f, float), CABINET_TAGS, o, pr, {"bands": [band]})
    return fake


def test_mesh_check_limit_is_the_last_band_within_half_a_db(tmp_path, monkeypatch):
    import speaker.directivity_report as dr
    calls = []
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field(calls))
    monkeypatch.setattr(dr, "_fine_dofs", lambda p, coax, bass, h: 1000)
    out = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    assert out["per_wavelength"] == [6.0, 8.0] and out["capped_by_memory"] is False
    assert out["trusted_limit_hz"] == pytest.approx(1000 * 2 ** (1 / 6))
    dev = out["band_max_dev_db"]
    assert dev["1000"] < 0.01 and dev["1260"] > 0.5
    # Bands where h_max sets both meshes are solved once; the others at both densities, one call per band.
    assert set(out["elements"]) == {"250", "500", "1000", "2000", "4000"}
    assert out["elements"]["250"]["same_mesh"] is True and out["elements"]["2000"]["same_mesh"] is False
    assert out["elements"]["2000"]["fine"]["nodes"] == 100 and "solve_s" in out["elements"]["2000"]["coarse"]
    assert sorted(c[0] for c in calls) == [6.0] * 5 + [8.0] * 3


def test_mesh_check_stops_before_a_fine_mesh_over_the_memory_budget(tmp_path, monkeypatch):
    import speaker.directivity_report as dr
    calls = []
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field(calls))
    monkeypatch.setattr(dr, "_fine_dofs", lambda p, coax, bass, h: 10**6 if h < 0.03 else 1000)
    out = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    assert out["capped_by_memory"] is True
    assert out["trusted_limit_hz"] == pytest.approx(50 * 2 ** (51 / 12))  # capped at the top solved frequency
    assert set(out["elements"]) == {"250", "500", "1000"}
    assert max(c[2] for c in calls) < 1000.0


def test_converged_mesh_check_feeds_a_usable_limit_to_the_run(tmp_path, monkeypatch):
    import speaker.directivity_report as dr
    from speaker.bem.run import freq_grid
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field([], hot=None))
    monkeypatch.setattr(dr, "_fine_dofs", lambda p, coax, bass, h: 1000)
    mc = json.loads(dr.mesh_check(None, tmp_path, f_max=4000.0).read_text())
    top = float(freq_grid(50, 4000, 12)[-1])  # 3805.5 Hz: the last band's upper edge is past the grid
    assert mc["trusted_limit_hz"] == pytest.approx(top)
    # The full run's cabinet grid is the same, so the limit fits inside it.
    f = freq_grid(50, 20000, 12)
    tl, text = dr._trusted_limit(mc, cab_top=float(f[f <= 4000][-1]), f_top=float(f[-1]), puck=True)
    assert tl == pytest.approx(top) and "unverified" not in text
    # An older check that reached past the cabinet field is capped at its top, and says so.
    tl, text = dr._trusted_limit({**mc, "trusted_limit_hz": 4489.8}, cab_top=top, f_top=20000.0, puck=True)
    assert tl == top and "capped" in text
    # A limit too low for the puck blend (tl/2 >= 500 Hz) is an error, not a crash.
    with pytest.raises(ValueError, match="trusted 3D limit"):
        dr._trusted_limit({**mc, "trusted_limit_hz": 800.0}, cab_top=top, f_top=20000.0, puck=True)
    assert dr._trusted_limit(None, cab_top=898.0, f_top=898.0, puck=False) == (3000.0, "3000 Hz (unverified)")


def test_mesh_check_with_no_converged_band_is_an_error(tmp_path, monkeypatch, capsys):
    import speaker.directivity_report as dr
    monkeypatch.setattr(dr, "cabinet_field", _fake_cabinet_field([], hot=(0.0, 1e9)))
    monkeypatch.setattr(dr, "_fine_dofs", lambda p, coax, bass, h: 1000)
    assert main(["directivity", "--mesh-check", "--out", str(tmp_path)]) == 2
    err = capsys.readouterr().err
    assert "error: " in err and "no 1/3-octave band" in err
