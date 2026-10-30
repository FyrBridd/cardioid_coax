"""The passive topology (value names, netlists in both modes, driver loads, polarity) and the value fit."""
from __future__ import annotations

import warnings
from dataclasses import replace

import numpy as np
import pytest
from fakes import SAMPLE, omni_ways

from speaker.crossover import Peq, peq
from speaker.curves import smooth
from speaker.network import inductor_dcr, snap
from speaker.passive import (
    BOUNDS, SHUNT_Q_MAX, Z_GRID, PassiveFit, Topology, bass_scale, build, driver_loads, fit_passive_at, passive_total,
    start_values, transfers, value_names,
)

F = np.geomspace(20, 20000, 300)
R8 = {k: np.full(len(F), 8.0 + 0j) for k in ("bass", "coax_lf", "coax_hf")}
FULL = Topology(bass_order=3, zobel=True, bass_comp=True, lf_order=3, lf_notch=True, lf_shelf=True, hf_shunt=True)


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
    monkeypatch.setattr(P, "LEAK_WEIGHT", 0.0)   # the level term alone (the bass-under term also pulls on level)
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


def test_the_contour_is_a_resistor_with_an_inductor_across_it_before_the_l_pad():
    t = Topology(lf_shelf=True, lf_notch=True)
    parts = {p.name: p for p in build(_values(t), t, "passive")}
    rc, lc, rs = parts["R_lc"], parts["L_lc"], parts["R_ls"]
    assert (rc.a, rc.b) == (lc.a, lc.b) == ("ln", "lc") and rs.a == "lc"
    assert lc.dcr == pytest.approx(inductor_dcr(lc.value, cored=False))
    assert value_names(t, "passive").index("R_lc") == value_names(t, "passive").index("R_ln") + 1


def test_the_contour_shelves_from_no_loss_to_the_resistive_loss():
    """A pure R parallel L in series with an 8 ohm load: about 0 dB well below its corner R / (2 pi L), and the
    loss 8 / (R + 8) well above it."""
    from speaker.network import Part
    r, ind = 10.0, 3.3e-3
    parts = [Part("R_lc", "R", "in", "lc", r), Part("L_lc", "L", "in", "lc", ind), Part("coax_lf", "Z", "lc", "0", 1.0)]
    f = np.array([5.0, 50000.0])
    h, _ = transfers(parts, f, {"coax_lf": np.full(2, 8.0 + 0j)}, Topology())
    db = 20 * np.log10(np.abs(h["coax_lf"]))
    assert db[0] == pytest.approx(0.0, abs=0.05)
    assert db[1] == pytest.approx(-20 * np.log10((r + 8) / 8), abs=0.05)


def test_the_bass_sits_under_the_coax_lf_above_twice_its_point(monkeypatch):
    """On equal omni ways flatness alone lets the bass run close to the coax LF well above fb; the bass-under term
    pulls it down towards BASS_UNDER_DB under the coax LF over [2 fb, 600 Hz]."""
    import speaker.passive as P
    from speaker.passive import FIT_SMOOTH, _smoother
    w, z, t = omni_ways(FIT_F), _resistive(), Topology()
    sm, i = _smoother(FIT_F, FIT_SMOOTH), w.obs.at23_idx
    m = (FIT_F >= 400.0) & (FIT_F <= 600.0)

    def worst(fit):
        h, _ = transfers(build(fit.values, t, "passive"), FIT_F, driver_loads(z(FIT_F), t), t)
        b, c = (sm @ np.abs(h[k] * w.ways[k][:, i]) ** 2 for k in ("bass", "coax_lf"))
        return float((10 * np.log10(b / c))[m].max())
    monkeypatch.setattr(P, "LEAK_WEIGHT", 0.0)
    free = worst(fit_passive_at(w, z, 2851.0, 200.0, 1600.0, t))
    monkeypatch.undo()   # the module's LEAK_WEIGHT
    held = worst(fit_passive_at(w, z, 2851.0, 200.0, 1600.0, t))
    assert free > -P.BASS_UNDER_DB + 3.0 and held < free - 2.0 and held <= -P.BASS_UNDER_DB + 1.0


def test_the_bass_under_term_skips_points_where_a_way_is_exactly_zero(monkeypatch):
    """With the bass exactly zero from 2 fb up the term has nothing to count: the fit is the same with it off and
    with it made severe."""
    import speaker.passive as P
    w, z = omni_ways(FIT_F), _resistive()
    w.ways["bass"][FIT_F >= 400.0] = 0
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        monkeypatch.setattr(P, "LEAK_WEIGHT", 0.0)
        off = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, Topology())
        monkeypatch.setattr(P, "LEAK_WEIGHT", 10.0)
        monkeypatch.setattr(P, "BASS_UNDER_DB", 100.0)
        severe = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, Topology())
    assert severe.values == off.values and severe.cost == pytest.approx(off.cost)


# ---- the electrical stopband term (C1)

def _inductive_bass(f, r=3.3, le=3e-3):
    """The coax ways 8 ohm; the woofer pair r ohm plus a voice-coil inductance le (H)."""
    f = np.asarray(f, dtype=float)
    return {"coax_lf": np.full(len(f), 8.0 + 0j), "coax_hf": np.full(len(f), 8.0 + 0j),
            "bass_pair": r + 2j * np.pi * f * le}


def _stop_h(v, t, z, fb, fh, mode="passive"):
    from speaker.passive import stop_grid
    g = stop_grid(fb, fh)
    h, _ = transfers(build(v, t, mode), g, driver_loads(z(g), t), t)
    return h


def test_a_bass_low_pass_resonating_with_the_voice_coil_is_penalised():
    """8 mH with 3.3 µF on a 3.3 ohm + 3 mH load resonates near 1.5 kHz: the bass rises over its level at fb there."""
    from speaker.passive import STOP_DB, stop_excess_db, stopbands
    t, fb, fh = Topology(), 150.0, 1200.0
    v = _values(t) | {"L_b1": 8e-3, "C_b1": 3.3e-6}
    h = _stop_h(v, t, _inductive_bass, fb, fh)
    bands = stopbands(fb, fh)
    ex = stop_excess_db(h, fb, fh)
    assert len(ex) == sum(int(m.sum()) for _, _, m in bands)
    j = [k for k, _, _ in bands].index("bass")
    start = sum(int(m.sum()) for _, _, m in bands[:j])
    bass = ex[start:start + int(bands[j][2].sum())]
    assert bass.max() > 20.0   # the resonance sits more than 20 dB over the STOP_DB line
    # the same values on a resistive load are a plain low-pass: a much smaller excess
    h8 = _stop_h(v, t, lambda f: _inductive_bass(f) | {"bass_pair": np.full(len(f), 3.3 + 0j)}, fb, fh)
    assert stop_excess_db(h8, fb, fh)[start:start + len(bass)].max() < bass.max() - 10.0
    assert STOP_DB == -12.0


def test_the_fit_keeps_the_bass_out_of_its_far_stopband_on_an_inductive_load(monkeypatch):
    """A 2.5 ohm + 1 mH woofer pair under the 3.5 ohm floor, its way exactly zero above the trusted limit (as the BEM
    bass field is): the floor pushes the bass low-pass towards a small C_b1 that resonates with the voice coil where
    the acoustic data cannot see it. The stopband term holds it STOP_DB under its level at fb from 4 fb up."""
    import speaker.passive as P
    w, t, fb, fh = omni_ways(FIT_F), Topology(), 150.0, 1200.0
    w.ways["bass"][FIT_F > 2851.0] = 0
    m = Z_GRID >= 4 * fb

    def z(f):
        return _inductive_bass(f, 2.5, 1e-3)

    def worst(fit):
        h = _stop_h(fit.values, t, z, fb, fh)
        a = 20 * np.log10(np.abs(h["bass"]))
        return float(a[:len(Z_GRID)][m].max() - a[len(Z_GRID)])
    monkeypatch.setattr(P, "STOP_WEIGHT", 0.0)
    free = worst(fit_passive_at(w, z, 2851.0, fb, fh, t))
    monkeypatch.undo()
    held = worst(fit_passive_at(w, z, 2851.0, fb, fh, t))
    assert free > P.STOP_DB + 3.0   # without the term the fit lets the bass leak (this load invites it)
    assert held <= P.STOP_DB + 0.5


def test_the_stopband_term_covers_the_coax_hf_below_a_quarter_of_its_point():
    from speaker.passive import STOP_DB, stop_excess_db, stopbands
    t, fb, fh = Topology(), 150.0, 1200.0
    leaky = _values(t) | {"C_t1": 200e-6, "C_t2": 200e-6, "L_t1": 20e-3}   # the HF high-pass barely acts
    tight = _values(t) | {"C_t1": 6.8e-6, "C_t2": 10e-6, "L_t1": 0.6e-3}
    for mode in ("passive", "semi_active"):
        bands = stopbands(fb, fh, mode)
        names = [k for k, _, _ in bands]
        assert names.count("coax_hf") == 1
        j = names.index("coax_hf")
        _, fc, m = bands[j]
        assert fc == fh and Z_GRID[m].max() <= fh / 4 and Z_GRID[m].min() == Z_GRID[0]
        start = sum(int(x.sum()) for _, _, x in bands[:j])
        sl = slice(start, start + int(m.sum()))
        assert stop_excess_db(_stop_h(leaky, t, _resistive(), fb, fh, mode), fb, fh, mode)[sl].max() > 5.0
        h = _stop_h(tight, t, _resistive(), fb, fh, mode)
        assert stop_excess_db(h, fb, fh, mode)[sl].max() == 0.0
        a = 20 * np.log10(np.abs(h["coax_hf"]))
        assert a[:len(Z_GRID)][m].max() - a[len(Z_GRID) + 1] <= STOP_DB
    # the coax LF's lower stopband exists in passive mode only (in semi-active the DSP owns that edge)
    assert [k for k, _, _ in stopbands(fb, fh, "passive")].count("coax_lf") == 2
    assert [k for k, _, _ in stopbands(fb, fh, "semi_active")].count("coax_lf") == 1
    assert "bass" not in [k for k, _, _ in stopbands(fb, fh, "semi_active")]


def test_an_option_that_lowers_the_impedance_further_under_the_floor_is_left_out(monkeypatch):
    """Neither network clears the floor: an option may not buy flat23 with an even lower impedance."""
    from speaker.passive import search_passive
    _fake_search(monkeypatch, lambda fb, t: (2.9 if t.zobel else 3.0, 1.0),
                 lambda fit: {"flat23": 0.5} if fit.topo.zobel else {})
    res = search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=None)
    assert not res.fit.topo.zobel
    assert any(n.startswith("Left out the Zobel") and "3.00 Ω to 2.90 Ω" in n for n in res.notes)


def test_a_failing_series_search_keeps_the_parallel_result_with_a_note(monkeypatch):
    import speaker.passive as P
    _fake_search(monkeypatch, lambda fb, t: (3.0, 1.0))
    real = P._search

    def boom(ways, z, tl, lf_hf, base, z_floor):
        if base.bass_wiring == "series":
            raise ValueError("no passive crossover candidate has all 6 trusted search scores")
        return real(ways, z, tl, lf_hf, base, z_floor)
    monkeypatch.setattr(P, "_search", boom)
    res = P.search_passive(omni_ways(FIT_F), _resistive(), 2851.0, lf_hf_min=None)
    assert res.series is None and res.fit.topo.bass_wiring == "parallel"
    assert any("woofers-in-series" in n and "failed" in n and "no passive crossover candidate" in n
               for n in res.notes)


# ---- no bass peak below the fitted band (follow-up: the +4.7 dB resonance at 54 Hz)

LOW_F = np.geomspace(40.0, 16000.0, 180)   # reaches under FIT_LO, as the BEM grid does


def _sealed_bass(f, re=3.3, res=15.0, fc=35.0, qm=2.0, le=0.5e-3):
    """The coax ways 8 ohm; the woofer pair a sealed-box driver: re + j w le plus its motional impedance (a parallel
    RLC peaking at res ohm at fc). Above fc the pair looks capacitive, and the bass low-pass's series inductor
    resonates with it under FIT_LO, as on the real woofers."""
    w, w0 = 2 * np.pi * np.asarray(f, dtype=float), 2 * np.pi * fc
    zm = 1 / (1 / res + 1 / (1j * w * res / (qm * w0)) + 1j * w * qm / (res * w0))
    return {"coax_lf": np.full(len(w), 8.0 + 0j), "coax_hf": np.full(len(w), 8.0 + 0j),
            "bass_pair": re + 1j * w * le + zm}


def _max_gain_db(fit, z, branch="bass"):
    t = fit.topo
    h, _ = transfers(build(fit.values, t, "passive"), Z_GRID, driver_loads(z(Z_GRID), t), t)
    return float(20 * np.log10(np.abs(h[branch])).max())


def _low_rise_db(w, fit, z, tl=2851.0):
    """The largest rise of the smoothed 23° sum above the fitted level below FIT_LO (independent of passive.py's)."""
    from speaker.passive import FIT_LO, FIT_SMOOTH, _smoother
    t = passive_total(w, fit, z)[:, w.obs.at23_idx]
    db = 10 * np.log10(_smoother(w.f, FIT_SMOOTH) @ np.abs(t) ** 2)
    return float(db[w.f < FIT_LO].max() - _fit_level(w, fit, z, tl))


def test_electrical_gain_cap_penalises_a_resonant_branch(monkeypatch):
    """1 mH with 390 µF into 3.3 ohm is a low-pass of Q ~2: about +5 dB near 240 Hz, which the cap penalises; the
    Butterworth start values on the same load do not exceed 0 dB and cost nothing. A fit on that load stays under
    the cap, and on a sealed woofer's impedance (which invites the resonance under FIT_LO, where the residual cannot
    see it) the cap is what holds it there."""
    import speaker.passive as P
    from speaker.passive import ELEC_GAIN_MAX_DB, gain_excess_db
    assert ELEC_GAIN_MAX_DB == 1.0
    t, z = Topology(), _resistive(bass=3.3)
    loads = driver_loads(z(Z_GRID), t)
    v = _values(t) | {"L_b1": 1e-3, "C_b1": 390e-6}
    h, _ = transfers(build(v, t, "passive"), Z_GRID, loads, t)
    assert 20 * np.log10(np.abs(h["bass"])).max() > 3.0
    ex = gain_excess_db(h)
    assert len(ex) == 3 * len(Z_GRID) and ex.max() > 2.0 and np.all(ex >= 0)
    calm = _values(t) | start_values(t, 170.0, 1200.0, loads)
    hc, _ = transfers(build(calm, t, "passive"), Z_GRID, loads, t)
    assert gain_excess_db({"bass": hc["bass"]}).max() == 0.0
    w = omni_ways(LOW_F)
    assert _max_gain_db(fit_passive_at(w, z, 2851.0, 170.0, 1200.0, t), z) <= ELEC_GAIN_MAX_DB + 0.2
    weight = P.PEAK_WEIGHT
    monkeypatch.setattr(P, "PEAK_TOL_DB", 100.0)   # the electrical cap alone
    monkeypatch.setattr(P, "PEAK_WEIGHT", 0.0)
    free = _max_gain_db(fit_passive_at(w, _sealed_bass, 2851.0, 170.0, 1200.0, t), _sealed_bass)
    monkeypatch.setattr(P, "PEAK_WEIGHT", weight)
    held = _max_gain_db(fit_passive_at(w, _sealed_bass, 2851.0, 170.0, 1200.0, t), _sealed_bass)
    assert free > ELEC_GAIN_MAX_DB + 2.0 and held <= ELEC_GAIN_MAX_DB + 0.2


def test_no_acoustic_peak_below_the_fitted_band(monkeypatch):
    """On omni ways from 40 Hz with a sealed woofer's impedance, the fit without the guard lets the sum rise well over
    its level under FIT_LO; with it the rise stays within PEAK_TOL_DB. A sum that rolls off under FIT_LO costs
    nothing: the fit is the same with the guard as without it."""
    import speaker.passive as P
    from speaker.passive import PEAK_TOL_DB
    assert PEAK_TOL_DB == 1.0
    t, w = Topology(), omni_ways(LOW_F)
    monkeypatch.setattr(P, "ELEC_GAIN_MAX_DB", 100.0)   # the acoustic guard alone
    monkeypatch.setattr(P, "PEAK_TOL_DB", 100.0)
    free = _low_rise_db(w, fit_passive_at(w, _sealed_bass, 2851.0, 170.0, 1200.0, t), _sealed_bass)
    monkeypatch.setattr(P, "PEAK_TOL_DB", 1.0)
    held = _low_rise_db(w, fit_passive_at(w, _sealed_bass, 2851.0, 170.0, 1200.0, t), _sealed_bass)
    assert free > PEAK_TOL_DB + 3.0 and held <= PEAK_TOL_DB + 0.2
    # the woofers roll off steeply under 100 Hz (4th order at 100 Hz) on a resistive load: nothing to penalise
    s = 1j * LOW_F / 100.0
    rolled = omni_ways(LOW_F)
    rolled.ways["bass"] *= ((s ** 2 / (s ** 2 + s / 0.707 + 1)) ** 2)[:, None]
    z = _resistive(bass=4.0)
    held = fit_passive_at(rolled, z, 2851.0, 170.0, 1200.0, t)
    assert _low_rise_db(rolled, held, z) < 0.0
    monkeypatch.setattr(P, "PEAK_TOL_DB", 100.0)
    free = fit_passive_at(rolled, z, 2851.0, 170.0, 1200.0, t)
    assert held.values == free.values and held.cost == pytest.approx(free.cost)


# ---- the impedance-peak compensator across the woofers (follow-up: the owner chose it, 2026-10-30)

def _round3_pair(f):
    """The real woofer pair (two RSS265HO-44, parallel) in round 3's sealed box: the lumped impedance, and the
    driver and the volume each woofer sees."""
    from pathlib import Path

    from speaker.bass import sealed_impedance
    from speaker.drivers import load_driver
    from speaker.layout import compute_layout
    from speaker.params import load_params
    p = load_params(Path(__file__).resolve().parents[1] / "params" / "round-3.toml")
    d = load_driver(p.bass)
    vb = compute_layout(p, load_driver(p.coax), d).bass_net / 2
    return sealed_impedance(d, vb, f) / 2, d, vb


def _textbook_comp(re, fc, qec, qmc):
    """The series R-L-C that cancels a sealed woofer's box resonance: (R, L, C)."""
    return re * (1 + qec / qmc), re * qec / (2 * np.pi * fc), 1 / (2 * np.pi * fc * re * qec)


def _rlc(f, r, l, c):
    w = 2 * np.pi * np.asarray(f, dtype=float)
    return r + 1j * w * l + 1 / (1j * w * c)


def test_bass_comp_flattens_the_woofer_impedance_peak():
    """Across the real pair in round 3's box, the textbook compensator (from the sealed alignment) takes the 20-200 Hz
    peak-to-minimum |Z| ratio from ~9 to under 1.5, and start_values derives the same parts from the impedance."""
    from speaker.bass import sealed_alignment
    f = np.geomspace(20.0, 200.0, 400)
    z, d, vb = _round3_pair(f)
    a = sealed_alignment(d, vb)
    r, l, c = _textbook_comp(d.re / 2, a.fc, d.qes * a.fc / d.fs, d.qms * a.fc / d.fs)
    zc = _rlc(f, r, l, c)
    both = z * zc / (z + zc)
    assert np.abs(z).max() / np.abs(z).min() > 8.0
    assert np.abs(both).max() / np.abs(both).min() < 1.5
    t = Topology(bass_comp=True)
    v = start_values(t, 170.0, 1600.0, driver_loads({"bass_pair": _round3_pair(Z_GRID)[0],
                                                     "coax_lf": np.full(len(Z_GRID), 8.0 + 0j),
                                                     "coax_hf": np.full(len(Z_GRID), 8.0 + 0j)}, t))
    assert v["R_bc"] == pytest.approx(r, rel=1e-3)
    assert v["L_bc"] == pytest.approx(l, rel=1e-3)
    assert v["C_bc"] == pytest.approx(c, rel=1e-3)


def test_bass_comp_values_and_netlist():
    from speaker.passive import OPTION_TEXT, OPTIONS
    t = Topology(zobel=True, bass_comp=True)
    names = value_names(t, "passive")
    assert names[:7] == ["L_b1", "C_b1", "R_bz", "C_bz", "R_bc", "L_bc", "C_bc"]
    assert value_names(Topology(bass_comp=True), "passive")[:5] == ["L_b1", "C_b1", "R_bc", "L_bc", "C_bc"]
    assert not {"R_bc", "L_bc", "C_bc"} & set(value_names(t, "semi_active"))
    assert "R_bc" not in value_names(Topology(), "passive")
    assert BOUNDS["C_bc"] == (100e-6, 3000e-6) and BOUNDS["L_bc"] == (1e-3, 50e-3) and BOUNDS["R_bc"] == (1.0, 20.0)
    assert OPTIONS[0] == "bass_comp" and OPTIONS.index("lf_shelf") == 1
    assert OPTION_TEXT["bass_comp"] == "an impedance-peak compensator across the woofers"
    # a series R-L-C from the woofers' node to ground, on the network side (the woofers hang on the same node)
    for topo, end in ((t, "b1"), (Topology(bass_order=3, bass_comp=True), "b2")):
        parts = {p.name: p for p in build(_values(topo), topo, "passive")}
        assert (parts["R_bc"].a, parts["R_bc"].b) == (end, "bc1")
        assert (parts["L_bc"].a, parts["L_bc"].b) == ("bc1", "bc2")
        assert (parts["C_bc"].a, parts["C_bc"].b) == ("bc2", "0")
        assert (parts["bass"].a, parts["bass"].b) == (end, "0")
        assert parts["L_bc"].dcr == pytest.approx(inductor_dcr(SAMPLE["L_bc"], cored=True))
    assert not {"R_bc", "L_bc", "C_bc"} & {p.name for p in build(_values(t), t, "semi_active")}
    # start values: the textbook compensator of a hand-made woofer, Re 3.2 ohm with a box resonance at 40 Hz
    re, fc, qec, qmc, le = 3.2, 40.0, 0.75, 6.6, 1.5e-3
    w, w0 = 2 * np.pi * Z_GRID, 2 * np.pi * fc
    res = re * qmc / qec
    zm = 1 / (1 / res + 1 / (1j * w * res / (qmc * w0)) + 1j * w * qmc / (res * w0))
    z = {"bass_pair": re + 1j * w * le + zm, "coax_lf": np.full(len(w), 8.0 + 0j),
         "coax_hf": np.full(len(w), 8.0 + 0j)}
    tc = Topology(bass_comp=True)
    v = start_values(tc, 170.0, 1600.0, driver_loads(z, tc))
    r, l, c = _textbook_comp(re, fc, qec, qmc)
    assert v["R_bc"] == pytest.approx(r, rel=1e-3) and v["L_bc"] == pytest.approx(l, rel=1e-3)
    assert v["C_bc"] == pytest.approx(c, rel=1e-3)
    # the same woofers in series: four times the impedance, so R and L four times, C a quarter
    ts = Topology(bass_comp=True, bass_wiring="series")
    vs = start_values(ts, 170.0, 1600.0, driver_loads(z, ts))
    assert vs["R_bc"] == pytest.approx(4 * r, rel=1e-3) and vs["L_bc"] == pytest.approx(4 * l, rel=1e-3)
    assert vs["C_bc"] == pytest.approx(c / 4, rel=1e-3)
    # a load with no resonance (resistive) still gives values inside the bounds
    vr = start_values(tc, 170.0, 1600.0, driver_loads(_resistive(bass=4.0)(Z_GRID), tc))
    for n in ("R_bc", "L_bc", "C_bc"):
        assert BOUNDS[n][0] < vr[n] < BOUNDS[n][1]
