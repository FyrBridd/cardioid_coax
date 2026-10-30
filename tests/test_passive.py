"""The passive topology: value names, netlists in both modes, driver loads and polarity."""
from __future__ import annotations

import numpy as np
import pytest
from fakes import SAMPLE

from speaker.passive import BOUNDS, Topology, bass_scale, build, driver_loads, transfers, value_names

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
