"""The passive topology (value names, netlists in both modes, driver loads, polarity) and the value fit."""
from __future__ import annotations

import warnings

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
