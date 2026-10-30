"""``speaker passive``: power ratings, the network sheet, parts list and VituixCAD netlist, and the report run."""
from __future__ import annotations

import csv
import io
import json
import warnings

import numpy as np
import pytest
from fakes import SAMPLE, omni_ways

from speaker.network import Part
from speaker.passive import PassiveFit, PassiveResult, Tolerance, Topology, build, value_names

TL = 2851.0


def _resistive(f, bass=8.0):
    f = np.asarray(f, dtype=float)
    return {"coax_lf": np.full(len(f), 8.0 + 0j), "coax_hf": np.full(len(f), 8.0 + 0j),
            "bass_pair": np.full(len(f), bass + 0j), "bass_driver": np.full(len(f), 2 * bass + 0j)}


def _result(topo=Topology(), z_ok=True, series=None, **values):
    v = {n: SAMPLE[n] for n in value_names(topo, "passive")} | values
    fit = PassiveFit(200.0, 1600.0, topo, v, dict(v), 3.5, 3.9 if z_ok else 3.2, z_ok, 0.0)
    tol = Tolerance(3.6, 1.20, 1.45, (5.0, 6.0), True)
    return PassiveResult(fit, tol, series, ["Chose 200 Hz bass to coax and 1600 Hz coax LF to HF."])


def test_power_rating_rounds_up_to_the_standard_series():
    from speaker.passive_report import power_ratings
    parts = [Part("R_ts", "R", "in", "o", 4.0), Part("coax_hf", "Z", "o", "0", 0.0)]
    f = np.geomspace(20, 20000, 50)
    r = power_ratings(parts, {"coax_hf": np.full(50, 4 + 0j)}, f, 25.0)["R_ts"]
    # 25/3 V across 8 ohm -> 1.04 A -> 4.34 W in R_ts; full voltage 39.1 W
    assert r[0] == pytest.approx(4.34, abs=0.05) and r[1] == pytest.approx(39.06, abs=0.1) and r[2] == 5
    r = power_ratings(parts, {"coax_hf": np.full(50, 4 + 0j)}, f, 60.0)["R_ts"]   # 25 W at a third
    assert r[2] == 30
    assert set(power_ratings(parts, {"coax_hf": np.full(50, 4 + 0j)}, f, 25.0)) == {"R_ts"}


def test_network_text_shows_floor_flatness_links_and_both_networks(monkeypatch):
    import speaker.passive_report as pr
    from speaker.evaluate import load_design
    monkeypatch.setattr(pr, "impedances", lambda d, f: _resistive(f))
    d = load_design(None)
    series = _result(Topology(bass_wiring="series", bass_inverted=True))
    res = _result(z_ok=False, series=series)
    text = pr.network_text(res, d, TL)
    assert "3.5 Ω" in text and "flat23 1.20 dB" in text
    assert "parallel" in text and "series" in text
    for k in ("K1", "K2", "K3", "K4", "K5"):
        assert k in text
    assert "Passive" in text and "Semi-active" in text and "Fully active" in text
    assert "cored inductor" in text and "polypropylene film capacitor" in text and "5 %" in text
    assert "Chose 200 Hz" in text and "tolerance" in text.lower()
    for p in build(res.fit.values, res.fit.topo, "passive"):
        assert p.name in text
    # semi-active: one netlist, the parallel network's, fed at cx; the series block only points to it
    assert text.count("-- Semi-active mode (") == 1 and "cx (COAX)-l1" in text
    series_block = text.split("== Alternative: the woofers in series")[1]
    pointer = "Semi-active mode uses the parallel network's coax boards; semi_active.txt is fitted to those."
    assert pointer in series_block
    assert "cx (COAX)" not in series_block


def test_network_text_states_the_sensitivity(monkeypatch):
    import speaker.passive_report as pr
    from speaker.evaluate import load_design
    monkeypatch.setattr(pr, "impedances", lambda d, f: _resistive(f))
    w = omni_ways(np.geomspace(50.0, 16000.0, 150))
    res = _result(z_ok=False, series=_result(Topology(bass_wiring="series")))
    text = pr.network_text(res, load_design(None), TL, ways=w, active_total=w.ways["coax_lf"])
    lines = [x for x in text.splitlines() if "Sensitivity" in x]
    assert len(lines) == 3   # active, parallel, series
    assert all("dB at 2 m for 2.83 V" in x and "what the woofers alone give" in x for x in lines)
    # the active total here is one way at the woofers' level: 0.0 dB below them
    assert "0.0 dB below what the woofers alone give (woofers in parallel)" in lines[0]
    assert "woofers in series" in lines[2]


def test_network_text_flags_a_value_at_a_search_limit(monkeypatch):
    import speaker.passive_report as pr
    from speaker.evaluate import load_design
    monkeypatch.setattr(pr, "impedances", lambda d, f: _resistive(f))
    text = pr.network_text(_result(L_b1=20e-3), load_design(None), TL)
    assert "L_b1 is at the 20 mH limit of the search" in text and "bass loss" in text
    assert "limit of the search" not in pr.network_text(_result(), load_design(None), TL)
    text = pr.network_text(_result(R_tp=100.0), load_design(None), TL)
    assert "R_tp is at the 100 Ω limit of the search" in text and "barely acts" in text


def test_parts_csv_lists_a_pair_with_types_and_ratings(monkeypatch):
    import speaker.passive_report as pr
    from speaker.evaluate import load_design
    monkeypatch.setattr(pr, "impedances", lambda d, f: _resistive(f))
    series = _result(Topology(bass_wiring="series"))
    res = _result(Topology(zobel=True), z_ok=False, series=series)
    rows = list(csv.DictReader(io.StringIO(pr.parts_csv(res, load_design(None)))))
    assert list(rows[0]) == ["part", "value", "unit", "type", "tolerance", "dcr_ohm", "peak_a", "rating_w", "qty_pair",
                             "network"]
    par = {r["part"]: r for r in rows if r["network"] == "parallel"}
    ser = {r["part"]: r for r in rows if r["network"] == "series"}
    assert set(par) == set(value_names(Topology(zobel=True), "passive"))
    assert set(ser) == set(value_names(Topology(), "passive"))
    assert all(r["qty_pair"] == "2" and r["tolerance"] == "5 %" for r in rows)
    assert par["L_b1"]["type"] == "cored inductor" and par["L_h1"]["type"] == "air-core inductor"
    assert float(par["L_b1"]["dcr_ohm"]) == pytest.approx(0.14 * np.sqrt(6.5), abs=0.005)
    assert par["C_b1"]["type"] == "bipolar electrolytic capacitor"   # 100 µF in the bass network
    assert par["C_bz"]["type"] == "bipolar electrolytic capacitor"   # 47 µF Zobel
    assert par["C_h1"]["type"] == "polypropylene film capacitor"     # 100 µF, but in the coax network
    assert par["C_l1"]["type"] == "polypropylene film capacitor"
    assert par["R_ts"]["type"].startswith("resistor (") and par["R_ts"]["type"].endswith(" W)")
    assert int(par["R_ts"]["rating_w"]) in (5, 10, 15, 20, 30, 50) and par["R_ts"]["dcr_ohm"] == ""
    assert par["L_b1"]["value"] == "6.5" and par["L_b1"]["unit"] == "mH"
    assert par["C_b1"]["value"] == "100" and par["C_b1"]["unit"] == "uF"
    assert par["R_ls"]["value"] == "1.5" and par["R_ls"]["unit"] == "ohm"
    # the series alternative is rated in passive mode only (semi-active uses the parallel boards)
    d = load_design(None)
    passive_only = pr.power_ratings(build(series.fit.values, series.fit.topo, "passive"),
                                    pr.driver_loads(_resistive(pr.Z_GRID), series.fit.topo), pr.Z_GRID,
                                    d.p.drive.amp_v_rms)
    assert all(int(ser[n]["rating_w"]) == passive_only[n][2] for n in passive_only)


def test_vituixcad_netlist_is_in_build_order_with_inverted_drivers_swapped():
    from speaker.passive_report import vituixcad_netlist
    res = _result(Topology(bass_inverted=True))
    text = vituixcad_netlist(res)
    lines = [x for x in text.splitlines() if x and not x.startswith("#")]
    parts = build(res.fit.values, res.fit.topo, "passive")
    assert [x.split()[0] for x in lines] == [p.name for p in parts]
    by = {x.split()[0]: x.split() for x in lines}
    assert by["L_b1"][1:5] == ["L", "6.5mH", "in", "b1"]
    assert by["bass"][3:5] == ["0", "b1"]           # inverted: + on the network's ground side
    assert by["coax_lf"][3:5] == ["lp", "0"]
    assert "cx" in text and "# " in text and "in " in text


def test_tolerance_keeps_each_draw():
    from speaker.passive import fit_passive_at, tolerance
    w = omni_ways(np.geomspace(100, 16000, 160))
    z = _resistive
    fit = fit_passive_at(w, z, TL, 200.0, 1600.0, Topology())
    tol = tolerance(w, z, TL, fit, n=8)
    assert len(tol.draws_z_min) == 8 and len(tol.draws_flat23) == 8
    assert min(tol.draws_z_min) == pytest.approx(tol.z_min)


def _synthetic(monkeypatch):
    """design_ways on omni ways (the bass exactly zero above the limit, a little directivity on the coax so the
    polar maps have contours), 8 Ω resistive drivers, and one-point grids."""
    import speaker.passive_report as pr
    from speaker.directivity_report import DesignWays
    from speaker.evaluate import load_design
    f = np.geomspace(50.0, 16000.0, 150)
    w = omni_ways(f)
    tilt = 1 + 0.4 * w.obs.dirs[:, 1] * (f[:, None] / 16000.0)   # -y is the front
    for k in ("coax_lf", "coax_hf"):
        w.ways[k] = w.ways[k] * (1 - 0.4 * w.obs.dirs[None, :, 1]) * tilt
    w.ways["bass"][f > TL] = 0
    w.meta.update(puck_blend=None, coax_far_h0_h180=None)
    d = load_design(None)

    def fake(params_path, *, trusted_limit, **kw):
        return DesignWays(d, f, w, TL, f"{TL:.0f} Hz (test)", None, None, None)
    monkeypatch.setattr(pr, "design_ways", fake)
    monkeypatch.setattr(pr, "impedances", lambda d, f: _resistive(f))
    monkeypatch.setattr("speaker.passive.BASS_COAX_GRID", (200.0,))
    monkeypatch.setattr("speaker.passive.LF_HF_GRID", (1600.0,))
    monkeypatch.setattr("speaker.passive.OPTIONS", ())
    monkeypatch.setattr("speaker.passive.TOL_DRAWS", 10)
    monkeypatch.setattr("speaker.semi_active.BASS_COAX_GRID", (200.0,))
    monkeypatch.setattr("speaker.xo_fit.BASS_COAX_GRID", (200.0,))
    monkeypatch.setattr("speaker.xo_fit.LF_HF_GRID", (1600.0,))


def test_run_passive_writes_every_output(tmp_path, monkeypatch):
    from speaker.passive_report import run_passive
    _synthetic(monkeypatch)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        out = run_passive(None, tmp_path, trusted_limit=(TL, "test"))
    assert out == tmp_path / "report.html"
    for name in ("network.txt", "parts.csv", "semi_active.txt", "vituixcad_netlist.txt", "scores.json",
                 "schematic_bass.svg", "schematic_coax.svg", "links.svg", "report.html"):
        assert (tmp_path / name).exists()
    for png in ("crossover_23_passive", "crossover_23_semi_active", "crossover_23_active", "impedance",
                "polar_h_passive", "polar_v_passive", "whole_rear_compare", "tolerance", "electrical"):
        assert (tmp_path / "plots" / f"{png}.png").stat().st_size > 0
    scores = json.loads((tmp_path / "scores.json").read_text())
    assert set(scores) == {"passive", "semi_active", "active"}
    assert all(any(s["key"] == "flat23" and s["confidence"] == "trusted" for s in v) for v in scores.values())
    assert "3.5" in (tmp_path / "network.txt").read_text()
    assert (tmp_path / "network.txt").read_text().count("what the woofers alone give") == 2   # active, parallel
    page = out.read_text()
    for s in ("Where passive falls short of active", "Tolerance", "bare woofer", "schematic_bass.svg",
              "links.svg", "plots/impedance.png", "plots/whole_rear_compare.png", "semi_active.txt",
              "Passive", "Semi-active", "Active", "8.00 Ω", "Passive behind active", "dB at 2 m for 2.83 V",
              "plots/electrical.png", "minimum is"):
        assert s in page, s
    assert "series" not in (tmp_path / "parts.csv").read_text().split("\n", 1)[1]


def test_run_passive_reports_the_series_fallback(tmp_path, monkeypatch):
    import speaker.passive_report as pr
    _synthetic(monkeypatch)
    monkeypatch.setattr(pr, "impedances", lambda d, f: _resistive(f, bass=1.0))   # no parallel network clears 3.5 Ω
    calls, rows = [], pr._part_rows

    def spy(fit, d, modes=pr.BOTH_MODES):
        calls.append((fit.topo.bass_wiring, tuple(modes)))
        return rows(fit, d, modes)
    monkeypatch.setattr(pr, "_part_rows", spy)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        page = pr.run_passive(None, tmp_path, trusted_limit=(TL, "test")).read_text()
    assert "Series fallback" in page and "schematic_bass_series.svg" in page
    assert (tmp_path / "schematic_bass_series.svg").exists() and (tmp_path / "schematic_coax_series.svg").exists()
    nets = {r["network"] for r in csv.DictReader(io.StringIO((tmp_path / "parts.csv").read_text()))}
    assert nets == {"parallel", "series"}
    text = (tmp_path / "network.txt").read_text()
    assert "Alternative: the woofers in series" in text and "does not clear 3.5 Ω" in text
    assert "for the series alternative passive mode only" in text
    # the series network is rated in passive mode only, everywhere (semi-active uses the parallel boards)
    assert ("series", ("passive",)) in calls and all(m == ("passive",) for w, m in calls if w == "series")
    assert "1.00 Ω" in page   # the bare woofer pair the semi-active bass channel drives


def test_link_states_are_the_single_source_of_every_mode_header(monkeypatch):
    import speaker.passive_report as pr
    from speaker.evaluate import load_design
    from speaker.schematic import link_states, links_text, terminal_text
    monkeypatch.setattr(pr, "impedances", lambda d, f: _resistive(f))
    res = _result(Topology(bass_inverted=True))
    text = pr.network_text(res, load_design(None), TL)
    m = link_states(res.fit.topo)
    semi = next(x for x in text.splitlines() if x.startswith("-- Semi-active mode ("))
    assert links_text(m["Semi-active"]) in semi and all(terminal_text(t) in semi for t in m["Semi-active"]["terminals"])
    full = next(x for x in text.splitlines() if x.startswith("-- Fully active mode"))
    assert links_text(m["Fully active"]) in full
    assert "two-pole" in text and "relative to each driver's own + terminal" in text
    from speaker.crossover import WayFilter
    from speaker.semi_active import SemiActive
    page = pr._semi_html(SemiActive(150.0, WayFilter(0.0, 0.0, -1, ()), WayFilter(0.0, 0.0, 1, ())), res.fit.topo)
    assert links_text(m["Semi-active"]) in page and "K1-K3" not in page


def test_inductor_peak_current_is_the_sine_peak_at_full_voltage():
    """1 mH (air core, DCR 0.5 ohm) in series with an 8 ohm load, 25 V rms: the most current flows at the lowest
    frequency, 25 / |8.5 + j 2 pi 20 Hz 1 mH| = 2.941 A rms, so 4.159 A peak."""
    from speaker.network import inductor_dcr
    from speaker.passive_report import peak_currents
    f = np.geomspace(20, 20000, 50)
    parts = [Part("L_l1", "L", "in", "o", 1e-3, dcr=inductor_dcr(1e-3, cored=False)),
             Part("coax_lf", "Z", "o", "0", 0.0)]
    peak = peak_currents(parts, {"coax_lf": np.full(50, 8 + 0j)}, f, 25.0)
    assert set(peak) == {"L_l1"}
    assert peak["L_l1"] == pytest.approx(np.sqrt(2) * 25.0 / abs(8.5 + 2j * np.pi * 20.0 * 1e-3), rel=1e-6)
    assert peak["L_l1"] == pytest.approx(4.159, abs=0.001)


def test_peak_currents_are_in_the_parts_list_and_the_network_sheet(monkeypatch):
    import speaker.passive_report as pr
    from speaker.evaluate import load_design
    monkeypatch.setattr(pr, "impedances", lambda d, f: _resistive(f))
    d = load_design(None)
    res = _result()
    rows = {r["part"]: r for r in csv.DictReader(io.StringIO(pr.parts_csv(res, d)))}
    v = d.p.drive.amp_v_rms
    want = {m: pr.peak_currents(build(res.fit.values, res.fit.topo, m), pr.driver_loads(_resistive(pr.Z_GRID),
                                                                                         res.fit.topo), pr.Z_GRID, v)
            for m in ("passive", "semi_active")}
    for n in ("L_b1", "L_h1", "L_l1", "L_t1"):
        assert float(rows[n]["peak_a"]) == pytest.approx(max(w[n] for w in want.values() if n in w), abs=0.005)
    assert rows["C_b1"]["peak_a"] == "" and rows["R_ls"]["peak_a"] == ""
    text = pr.network_text(res, d, TL)
    line = next(x for x in text.splitlines() if x.startswith("L_b1 "))
    assert f"{max(w['L_b1'] for w in want.values() if 'L_b1' in w):.2f} A" in line
    assert "cored: saturation current at least this" in line
    assert "cored: saturation" not in next(x for x in text.splitlines() if x.startswith("L_l1 "))


def test_network_sheet_gives_the_bass_band_minimum_impedance_and_the_stopbands(monkeypatch):
    import speaker.passive_report as pr
    from speaker.evaluate import load_design
    from speaker.passive import STOP_DB, stopbands
    monkeypatch.setattr(pr, "impedances", lambda d, f: _resistive(f))
    text = pr.network_text(_result(), load_design(None), TL)
    assert "minimum impedance in the bass band (20-500 Hz):" in text
    lines = [x for x in text.splitlines() if x.startswith("electrical stopband")]
    assert len(lines) == len(stopbands(200.0, 1600.0)) and all(f"limit {STOP_DB:g} dB" in x for x in lines)
