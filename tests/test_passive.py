"""The passive topology (value names, netlists in both modes, driver loads, polarity) and the value fit."""
from __future__ import annotations

import warnings
from dataclasses import replace

import numpy as np
import pytest
from fakes import SAMPLE, omni_ways

from speaker.crossover import Peq, peq
from speaker.curves import smooth
from speaker.network import snap
from speaker.passive import (
    BOUNDS, SHUNT_Q_MAX, Z_GRID, PassiveFit, Topology, bass_scale, build, driver_loads, fit_passive_at, passive_total,
    start_values, transfers, value_names,
)

F = np.geomspace(20, 20000, 300)
R8 = {k: np.full(len(F), 8.0 + 0j) for k in ("bass", "coax_lf", "coax_hf")}
FULL = Topology(bass_order=3, zobel=True, lf_order=3, lf_notch=True, hf_shunt=True)


def _values(t, mode="passive"):
    return {n: SAMPLE[n] for n in value_names(t, mode)}


def test_second_order_bass_low_pass_falls_40_db_per_decade():
    t = Topology()
    h, _ = transfers(build(_values(t), t, "passive"), F, R8, t)
    db = 20 * np.log10(np.abs(h["bass"]))
    i1, i2 = np.searchsorted(F, [2000, 20000])
    assert db[i1] - db[i2] == pytest.approx(40, abs=4)


def test_inverting_a_branch_flips_only_its_sign():
    t, ti = Topology(), Topology(bass_inverted=True, hf_inverted=True)
    h, _ = transfers(build(_values(t), t, "passive"), F, R8, t)
    hi, _ = transfers(build(_values(ti), ti, "passive"), F, R8, ti)
    np.testing.assert_allclose(hi["coax_hf"], -h["coax_hf"])
    np.testing.assert_allclose(hi["bass"], -h["bass"])
    np.testing.assert_allclose(hi["coax_lf"], h["coax_lf"])


def test_semi_active_netlist_keeps_every_coax_value_and_drops_the_coax_high_pass():
    t = Topology(lf_notch=True, hf_shunt=True)
    v = _values(t)
    p, s = build(v, t, "passive"), build(v, t, "semi_active")
    names = {x.name for x in s}
    assert not names & {"C_h1", "L_h1", "L_b1", "C_b1", "bass"}
    for x in s:
        if x.kind != "Z":
            assert x.value == next(y.value for y in p if y.name == x.name)
    h, _ = transfers(s, F, R8, t)
    assert "bass" not in h


def test_series_woofers_quadruple_the_load_and_halve_the_drive():
    z = {"bass_pair": np.full(len(F), 3.0 + 0j), "coax_lf": np.full(len(F), 5.0 + 0j),
         "coax_hf": np.full(len(F), 6.0 + 0j)}
    par, ser = driver_loads(z, Topology()), driver_loads(z, Topology(bass_wiring="series"))
    np.testing.assert_allclose(par["bass"], z["bass_pair"])
    np.testing.assert_allclose(ser["bass"], 4 * z["bass_pair"])
    assert set(par) == {"bass", "coax_lf", "coax_hf"}
    assert bass_scale(Topology()) == 1.0 and bass_scale(Topology(bass_wiring="series")) == 0.5


def test_value_names_order_and_options():
    assert value_names(Topology(), "passive") == [
        "L_b1", "C_b1", "C_h1", "L_h1", "R_ls", "R_lp", "L_l1", "C_l1",
        "C_t1", "L_t1", "C_t2", "R_ts", "R_tp"]
    names = value_names(FULL, "passive")
    assert set(names) == set(SAMPLE) and len(names) == len(SAMPLE)
    assert names.index("L_b2") == 2 and names.index("R_bz") == 3
    semi = value_names(FULL, "semi_active")
    assert not set(semi) & {"L_b1", "C_b1", "L_b2", "R_bz", "C_bz", "C_h1", "L_h1"}


def test_inductor_dcr_cored_for_bass_only():
    parts = build(_values(FULL), FULL, "passive")
    ls = {p.name: p for p in parts if p.kind == "L"}
    assert ls["L_b1"].dcr == pytest.approx(0.14 * np.sqrt(6.5))
    assert ls["L_l1"].dcr == pytest.approx(0.5 * np.sqrt(1.0))


def test_every_name_has_bounds_and_sample_is_inside():
    for n, v in SAMPLE.items():
        lo, hi = BOUNDS[n]
        assert lo <= v <= hi
    assert BOUNDS["R_lp"] == (0.5, 100.0) and BOUNDS["L_b1"] == (0.05e-3, 20e-3)


def test_full_topology_transfers_are_finite():
    h, z = transfers(build(_values(FULL), FULL, "passive"), F, R8, FULL)
    assert set(h) == {"bass", "coax_lf", "coax_hf"}
    assert all(np.all(np.isfinite(x)) for x in h.values()) and np.all(np.isfinite(z))


def _resistive(bass=8.0, lf=8.0, hf=8.0):
    return lambda f: {"coax_lf": np.full(len(f), lf + 0j), "coax_hf": np.full(len(f), hf + 0j),
                      "bass_pair": np.full(len(f), bass + 0j)}


FIT_F = np.geomspace(100, 16000, 160)


def test_start_values_are_the_textbook_butterworth_values():
    t = Topology(bass_order=3)
    v = start_values(t, 200.0, 1600.0, driver_loads(_resistive(bass=4.0)(Z_GRID), t))
    wb, wh = 2 * np.pi * 200.0, 2 * np.pi * 1600.0
    assert v["L_b1"] == pytest.approx(1.5 * 4.0 / wb) and v["C_b1"] == pytest.approx(4 / 3 / (wb * 4.0))
    assert v["L_b2"] == pytest.approx(0.5 * 4.0 / wb)
    assert v["C_h1"] == pytest.approx(1 / (np.sqrt(2) * wb * 8.0))
    assert v["L_h1"] == pytest.approx(np.sqrt(2) * 8.0 / wb)
    assert v["L_l1"] == pytest.approx(np.sqrt(2) * 8.0 / wh)
    assert v["C_l1"] == pytest.approx(1 / (np.sqrt(2) * wh * 8.0))
    assert v["C_t1"] == pytest.approx(2 / 3 / (wh * 8.0)) and v["L_t1"] == pytest.approx(0.75 * 8.0 / wh)
    assert v["C_t2"] == pytest.approx(2 / (wh * 8.0))
    assert set(v) == set(value_names(t, "passive"))


def test_start_values_pad_down_by_the_level_difference():
    t = Topology()
    v = start_values(t, 200.0, 1600.0, driver_loads(_resistive()(Z_GRID), t),
                     pad_db={"coax_lf": 6.0, "coax_hf": 10.0})
    for s, p, db in (("R_ls", "R_lp", 6.0), ("R_ts", "R_tp", 10.0)):
        rs, rp = v[s], v[p]
        assert rp * 8 / (rp + 8) / (rs + rp * 8 / (rp + 8)) == pytest.approx(10 ** (-db / 20))  # attenuation
        assert rs + rp * 8 / (rp + 8) == pytest.approx(8.0)  # the filter still sees 8 ohm


def test_fit_on_flat_resistive_ways_is_flat_and_safe():
    w = omni_ways(FIT_F)
    z = _resistive()
    fit = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, Topology())
    total = passive_total(w, fit, z)
    db = 20 * np.log10(np.abs(total[:, w.obs.at23_idx]))
    m = (w.f >= 150) & (w.f <= 2851)
    assert np.ptp(smooth(w.f, db, 1 / 3)[m]) < 3.0
    assert fit.z_ok and fit.z_min >= 3.5
    assert all(v == pytest.approx(snap(n[0], v)) for n, v in fit.values.items())
    assert list(fit.values) == value_names(Topology(), "passive") == list(fit.continuous)
    assert (fit.bass_coax, fit.lf_hf, fit.topo, fit.z_floor) == (200.0, 1600.0, Topology(), 3.5)
    assert np.isfinite(fit.cost) and fit.cost >= 0


def test_floor_is_respected_or_flagged():
    w = omni_ways(FIT_F)
    z = _resistive(bass=3.0)
    t = Topology()
    fit = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, t)
    _, z_in = transfers(build(fit.values, t, "passive"), Z_GRID, driver_loads(z(Z_GRID), t), t)
    assert fit.z_min == pytest.approx(np.min(np.abs(z_in)))
    assert fit.z_ok == (fit.z_min >= 3.5)


def test_fit_ignores_exact_zero_bass_above_the_limit():
    w = omni_ways(FIT_F)
    w.ways["bass"][w.f > 2851] = 0
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        fit = fit_passive_at(w, _resistive(), 2851.0, 200.0, 1600.0, Topology())
        total = passive_total(w, fit, _resistive())
    assert np.isfinite(fit.cost) and np.all(np.isfinite(total))


def test_semi_active_total_has_no_bass():
    w = omni_ways(FIT_F)
    w.ways["bass"] *= 1e6  # would dominate if it were summed
    t = Topology()
    fit = PassiveFit(200.0, 1600.0, t, {n: SAMPLE[n] for n in value_names(t, "passive")}, {}, 3.5, 4.0, True, 0.0)
    h, _ = transfers(build(fit.values, t, "semi_active"), w.f, driver_loads(_resistive()(w.f), t), t)
    np.testing.assert_allclose(passive_total(w, fit, _resistive(), mode="semi_active"),
                               h["coax_lf"][:, None] * w.ways["coax_lf"] + h["coax_hf"][:, None] * w.ways["coax_hf"])


def test_fit_keeps_the_crossovers_at_the_points():
    # on flat ways a flat sum alone does not pin the points (the sections could stagger anywhere)
    w, z, t = omni_ways(FIT_F), _resistive(), Topology()
    fit = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, t)
    h, _ = transfers(build(fit.values, t, "passive"), FIT_F, driver_loads(z(FIT_F), t), t)
    a = {k: np.abs(v) for k, v in h.items()}
    lo, hi = FIT_F < 800, FIT_F > 800
    xb = FIT_F[lo][np.argmin(np.abs(np.log(a["bass"] / a["coax_lf"]))[lo])]
    xh = FIT_F[hi][np.argmin(np.abs(np.log(a["coax_hf"] / a["coax_lf"]))[hi])]
    assert xb == pytest.approx(200.0, rel=0.1) and xh == pytest.approx(1600.0, rel=0.1)


def test_shunt_q_is_held_at_the_limit():
    w = omni_ways(FIT_F)
    w.ways["coax_hf"] = w.ways["coax_hf"] * peq(FIT_F, Peq(2600.0, 8.0, 6.0))[:, None]  # a narrow hump
    t = Topology(hf_shunt=True, bass_inverted=True, hf_inverted=True)
    fit = fit_passive_at(w, _resistive(), 16000.0, 200.0, 1600.0, t)
    for v in (fit.continuous, fit.values):
        assert np.sqrt(v["L_tr"] / v["C_tr"]) / v["R_tr"] <= SHUNT_Q_MAX * 1.1


def test_fit_smoothing_is_curves_smooth():
    from speaker.passive import FIT_SMOOTH, _smoother
    db = 10 * np.sin(np.linspace(0, 20, len(FIT_F)))
    np.testing.assert_allclose(10 * np.log10(_smoother(FIT_F, FIT_SMOOTH) @ 10 ** (db / 10)),
                               smooth(FIT_F, db, FIT_SMOOTH))


def _small_search(monkeypatch):
    # two points each; the 1000 Hz point sits below lf_hf_min=1200 and must be floored out
    monkeypatch.setattr("speaker.passive.BASS_COAX_GRID", (170.0, 210.0))
    monkeypatch.setattr("speaker.passive.LF_HF_GRID", (1000.0, 1600.0))
    monkeypatch.setattr("speaker.passive.TOL_DRAWS", 20)


def test_search_returns_a_grid_point_at_or_above_the_floor(monkeypatch):
    from speaker.passive import OPTIONS, PassiveResult, Tolerance, search_passive
    _small_search(monkeypatch)
    res = search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=1200.0)
    assert isinstance(res, PassiveResult) and isinstance(res.tolerance, Tolerance)
    assert res.fit.bass_coax in (170.0, 210.0) and res.fit.lf_hf == 1600.0
    assert res.fit.z_ok and res.series is None
    assert res.fit.z_floor >= 3.5
    assert len(res.notes) >= 1 + len(OPTIONS) + 1  # grid choice, each option, the tolerance check
    assert all(isinstance(n, str) and n.endswith(".") for n in res.notes)


def test_search_with_no_lf_hf_point_above_the_driver_minimum_raises(monkeypatch):
    from speaker.passive import search_passive
    _small_search(monkeypatch)
    with pytest.raises(ValueError, match="2500 Hz"):
        search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=2500.0)


def test_search_falls_back_to_series_woofers_when_no_parallel_network_clears_the_floor(monkeypatch):
    # 1 ohm, 2nd-order bass: even the largest bass inductor (20 mH, with its DCR) leaves |Z_in| at 20 Hz under
    # 3.5 ohm. No options: a 3rd-order bass adds a second series inductor and can clear it (and would be kept).
    from speaker.passive import search_passive
    _small_search(monkeypatch)
    monkeypatch.setattr("speaker.passive.OPTIONS", ())
    res = search_passive(omni_ways(FIT_F), _resistive(bass=1.0), 2851.0, lf_hf_min=1200.0)
    assert res.fit.z_ok is False
    assert res.series is not None and res.series.fit.topo.bass_wiring == "series"
    assert res.series.series is None
    assert any("series" in n for n in res.notes)
    # the parallel network misses its floor already, so the floor is not raised for it
    assert res.fit.z_floor == 3.5 and not any(n.startswith("Raised") for n in res.notes)


def test_tolerance_is_deterministic():
    from speaker.passive import passive_scores, tolerance
    w, z = omni_ways(FIT_F), _resistive()
    fit = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, Topology())
    a, b = tolerance(w, z, 2851.0, fit, n=20), tolerance(w, z, 2851.0, fit, n=20)
    assert a == b
    assert a.z_min < fit.z_min * 1.2 and a.whole_rear[0] <= a.whole_rear[1]
    flat23 = [s.value for s in passive_scores(w, fit, z, 2851.0) if s.key == "flat23" and s.confidence == "trusted"]
    assert a.flat23_nominal == pytest.approx(flat23[0])


def _near_floor_fit(z_target):
    """A 2nd-order network whose minimum |Z_in| sits at 20 Hz on the bass inductor (20 mH, so its reactance
    is large), with the bass load chosen so that minimum is z_target ohm at nominal values."""
    from scipy.optimize import brentq
    t = Topology()
    v = {n: SAMPLE[n] for n in value_names(t, "passive")} | {"L_b1": 20e-3, "C_b1": 10e-6}

    def z_min(rb):
        _, z_in = transfers(build(v, t, "passive"), Z_GRID, driver_loads(_resistive(bass=rb)(Z_GRID), t), t)
        return float(np.min(np.abs(z_in)))
    rb = brentq(lambda r: z_min(r) - z_target, 0.5, 8.0)
    return PassiveFit(200.0, 1600.0, t, v, v, 3.3, z_min(rb), True, 0.0), _resistive(bass=rb)


def test_tolerance_fails_below_3_3_ohm():
    from speaker.passive import tolerance
    fit, z = _near_floor_fit(3.35)
    assert fit.z_min == pytest.approx(3.35, abs=1e-6)
    tol = tolerance(omni_ways(FIT_F), z, 2851.0, fit, n=20)
    assert tol.z_min < 3.3 and tol.passed is False


def test_tolerance_passes_well_clear_of_the_floor():
    from speaker.passive import tolerance
    fit, z = _near_floor_fit(4.0)
    tol = tolerance(omni_ways(FIT_F), z, 2851.0, fit, n=20)
    assert 3.3 < tol.z_min < 4.0
    assert tol.passed == (tol.flat23_p95 <= tol.flat23_nominal + 0.5)


def test_an_option_is_kept_only_for_a_real_flat23_gain_at_little_cost_elsewhere():
    from speaker.passive import _option_verdict
    old = dict(window_spread=1.0, presence_lift=1.0, di_smooth_bass=0.5, di_smooth_lfhf=0.5, flat23=1.0,
               whole_rear=10.0)
    assert _option_verdict(old, old | {"flat23": 0.9}) == ""
    assert "0.1 dB" in _option_verdict(old, old | {"flat23": 0.95})
    assert "no better" in _option_verdict(old, old | {"flat23": 1.2})
    assert "window_spread" in _option_verdict(old, old | {"flat23": 0.5, "window_spread": 1.3})
    assert _option_verdict(old, old | {"flat23": 0.5, "window_spread": 1.2}) == ""
    assert "whole_rear" in _option_verdict(old, old | {"flat23": 0.5, "whole_rear": 9.7})  # higher is better
    assert _option_verdict(old, old | {"flat23": 0.5, "whole_rear": 12.0}) == ""


def test_passive_scores_are_score_all_on_the_passive_total():
    from speaker.passive import passive_scores
    from speaker.scores import score_all
    from speaker.spinorama import spinorama
    w, z, t = omni_ways(FIT_F), _resistive(), Topology()
    fit = PassiveFit(200.0, 1600.0, t, {n: SAMPLE[n] for n in value_names(t, "passive")}, {}, 3.5, 4.0, True, 0.0)
    total = passive_total(w, fit, z)
    want = score_all(w.f, total, w.ways["coax_lf"], w.obs, spinorama(total, w.obs, f=w.f), trusted_limit=2851.0)
    got = passive_scores(w, fit, z, 2851.0)
    assert [(s.key, s.confidence, s.value) for s in got] == [(s.key, s.confidence, s.value) for s in want]


def _one_point_search(monkeypatch, tolerances):
    """A one-pair search on flat 8 ohm ways, the tolerance check replaced by `tolerances` in turn (the last
    repeats); returns the floors fit_passive_at was called with."""
    import speaker.passive as P
    monkeypatch.setattr(P, "BASS_COAX_GRID", (170.0,))
    monkeypatch.setattr(P, "LF_HF_GRID", (1600.0,))
    queue = list(tolerances)

    def fake_tolerance(*a, **k):
        return queue.pop(0) if len(queue) > 1 else queue[0]
    floors = []
    real = P.fit_passive_at

    def spy(*a, **k):
        floors.append(a[6] if len(a) > 6 else k.get("z_floor", 3.5))
        return real(*a, **k)
    monkeypatch.setattr(P, "tolerance", fake_tolerance)
    monkeypatch.setattr(P, "fit_passive_at", spy)
    return floors


def test_a_raised_floor_the_refit_cannot_meet_keeps_the_last_network_that_met_its_floor(monkeypatch):
    import speaker.passive as P
    failed = P.Tolerance(3.2, 0.5, 0.7, (0.0, 0.0), False)
    floors = _one_point_search(monkeypatch, [failed])
    real = P.fit_passive_at

    def missing_raised_floors(*a, **k):
        fit = real(*a, **k)
        return replace(fit, z_ok=False, z_min=fit.z_floor - 0.05) if fit.z_floor > 3.5 else fit
    monkeypatch.setattr(P, "fit_passive_at", missing_raised_floors)
    res = P.search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=None)
    assert res.fit.z_floor == 3.5 and res.fit.z_ok
    assert res.tolerance is failed
    assert res.series is None  # it clears 3.5 ohm: no fallback for a missed raised floor
    assert any("still fails" in n and "3.20 Ω" in n and "+0.20 dB" in n for n in res.notes)
    assert not any(n.startswith("Raised") for n in res.notes)
    assert max(floors) == pytest.approx(3.6)  # one refit at 3.6 ohm missed it; 3.7 was never tried


def test_the_floor_is_raised_while_the_tolerance_check_fails_and_the_refit_meets_it(monkeypatch):
    import speaker.passive as P
    failed = P.Tolerance(3.2, 0.5, 0.7, (0.0, 0.0), False)
    passed = P.Tolerance(3.4, 0.5, 0.6, (0.0, 0.0), True)
    floors = _one_point_search(monkeypatch, [failed, passed])
    res = P.search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=None)
    assert res.fit.z_floor == pytest.approx(3.6) and res.fit.z_ok and res.tolerance is passed
    assert floors[-1] == pytest.approx(3.6) and max(floors) == pytest.approx(3.6)
    assert res.series is None
    assert any(n.startswith("Raised the impedance floor to 3.6 Ω") for n in res.notes)


_SCORES = dict(window_spread=1.0, presence_lift=1.0, di_smooth_bass=0.5, di_smooth_lfhf=0.5, flat23=1.0,
               whole_rear=10.0)


def _fake_search(monkeypatch, fit_of, scores_of=lambda fit: {}, grid=((170.0,), (1600.0,))):
    """search_passive with fit_passive_at, passive_scores and tolerance faked: fit_of(fb, topo) gives
    (z_min, cost), scores_of(fit) overrides _SCORES, and the tolerance check always passes. Returns the
    topologies fit_passive_at was asked for."""
    import speaker.passive as P
    from speaker.scores import Score
    monkeypatch.setattr(P, "BASS_COAX_GRID", grid[0])
    monkeypatch.setattr(P, "LF_HF_GRID", grid[1])
    calls = []

    def fit(ways, z, tl, fb, fh, topo, z_floor=3.5):
        calls.append(topo)
        z_min, cost = fit_of(fb, topo)
        return PassiveFit(fb, fh, topo, {n: SAMPLE[n] for n in value_names(topo, "passive")}, {}, z_floor, z_min,
                          z_min >= z_floor, cost)

    def scores(ways, fit, z, tl, mode="passive"):
        return [Score(k, k, v, "dB", "info", (100.0, tl), "trusted", "") for k, v in (_SCORES | scores_of(fit)).items()]
    monkeypatch.setattr(P, "fit_passive_at", fit)
    monkeypatch.setattr(P, "passive_scores", scores)
    monkeypatch.setattr(P, "tolerance", lambda *a, **k: P.Tolerance(5.0, 1.0, 1.0, (10.0, 10.0), True))
    return calls


def test_polarity_choice_puts_the_floor_before_the_cost(monkeypatch):
    from speaker.passive import search_passive

    def fit_of(fb, t):  # the cheapest polarity pair misses the floor, a dearer one clears it
        return {(False, False): (3.0, 1.0), (True, False): (4.0, 2.0)}.get((t.bass_inverted, t.hf_inverted),
                                                                            (3.0, 5.0))
    _fake_search(monkeypatch, fit_of)
    res = search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=None)
    assert (res.fit.topo.bass_inverted, res.fit.topo.hf_inverted) == (True, False)
    assert res.fit.z_ok and res.series is None


def test_the_pick_prefers_a_grid_point_that_clears_the_floor_and_then_needs_no_fallback(monkeypatch):
    from speaker.passive import search_passive
    calls = _fake_search(monkeypatch, lambda fb, t: (3.0, 1.0) if fb == 170.0 else (4.0, 2.0),
                         lambda fit: {"whole_rear": 20.0} if fit.bass_coax == 170.0 else {},
                         grid=((170.0, 210.0), (1600.0,)))
    res = search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=None)
    assert res.fit.bass_coax == 210.0 and res.fit.z_ok  # 170 Hz scores far better but misses 3.5 ohm
    assert res.series is None and all(t.bass_wiring == "parallel" for t in calls)


def test_the_fallback_note_says_no_network_in_the_search_clears_the_floor(monkeypatch):
    from speaker.passive import search_passive
    _fake_search(monkeypatch, lambda fb, t: (3.0, 1.0) if t.bass_wiring == "parallel" else (4.0, 1.0))
    res = search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=None)
    assert res.series is not None and res.series.fit.z_ok
    assert any(n.startswith("No parallel-woofer network in the search clears 3.5 Ω") for n in res.notes)


def test_an_option_that_lifts_the_network_over_the_floor_is_kept_whatever_flat23_does(monkeypatch):
    from speaker.passive import search_passive
    _fake_search(monkeypatch, lambda fb, t: (4.0, 1.0) if t.zobel else (3.0, 1.0),
                 lambda fit: {"flat23": 1.5} if fit.topo.zobel else {})
    res = search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=None)
    assert res.fit.topo.zobel and res.fit.z_ok and res.series is None
    assert any(n.startswith("Kept the Zobel") and "3.5 Ω" in n for n in res.notes)


def _fit_level(w, fit, z, tl=2851.0):
    """The level fit_passive_at anchors: the weighted mean of the smoothed 23° dB over [FIT_LO, tl] (and
    ABOVE_WEIGHT above it, to FIT_TOP)."""
    from speaker.passive import ABOVE_WEIGHT, FIT_LO, FIT_SMOOTH, FIT_TOP, _smoother
    t = passive_total(w, fit, z)[:, w.obs.at23_idx]
    db = 10 * np.log10(_smoother(w.f, FIT_SMOOTH) @ np.abs(t) ** 2)
    wt = np.where((w.f >= FIT_LO) & (w.f <= tl), 1.0, np.where((w.f > tl) & (w.f <= FIT_TOP), ABOVE_WEIGHT, 0.0))
    return float(np.sum(wt * db) / np.sum(wt))


def test_fit_holds_the_level_near_the_woofers_own(monkeypatch):
    """On equal omni ways the coax L-pads can throw level away for free (the level is otherwise not fitted); the
    level term keeps the sum within LEVEL_SLACK_DB (and a small margin) of the woofers' own level."""
    import speaker.passive as P
    w, z = omni_ways(FIT_F), _resistive()
    target = P.woofer_level_db(w, Topology())
    assert target == pytest.approx(0.0, abs=0.01)
    monkeypatch.setattr(P, "LEVEL_WEIGHT", 0.0)
    free = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, Topology())
    assert _fit_level(w, free, z) < target - P.LEVEL_SLACK_DB - 2.0   # without the term, level is thrown away
    monkeypatch.setattr(P, "LEVEL_WEIGHT", 3.0)
    held = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, Topology())
    assert _fit_level(w, held, z) >= target - P.LEVEL_SLACK_DB - 0.5


def test_no_level_term_without_woofer_data_in_its_band():
    from speaker.passive import woofer_level_db
    w = omni_ways(np.geomspace(200, 16000, 120))   # no point in 80-150 Hz
    assert woofer_level_db(w, Topology()) is None
    fit = fit_passive_at(w, _resistive(), 2851.0, 200.0, 1600.0, Topology())
    assert np.isfinite(fit.cost)
