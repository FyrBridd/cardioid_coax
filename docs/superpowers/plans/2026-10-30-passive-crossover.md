# Passive Crossover Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Design, score and report the passive and semi-active crossovers of the current design, with
component values, parts list, schematics and plots, from a new `speaker passive` command.

**Architecture:**
- A small nodal circuit solver (`network.py`) turns a netlist of R, L (with DCR), C and driver
  impedances into each driver's terminal voltage and the speaker's input impedance.
- `passive.py` builds the fixed topology from named values and fits them to the existing simulated
  way responses (`speaker.system.Ways`, pressure for 2.83 V at the terminals). It searches the same
  crossover grid as `xo_fit` and picks with `xo_fit.pick`.
- `semi_active.py` fits a two-channel DSP stage in front of the passive coax network.
  `passive_report.py` writes the outputs, and `schematic.py` draws them with `schemdraw`.

**Tech Stack:** Python, numpy, scipy (`least_squares`), matplotlib, schemdraw (new), pytest, uv.

**Spec:** `docs/superpowers/specs/2026-10-30-passive-crossover-design.md`. The main spec is
`docs/superpowers/specs/2026-10-22-cardioid-coax-floorstander-design.md` §7.2.

## Global Constraints

- Run everything through `uv run` (there is no `python` on PATH). The full suite is `uv run pytest -q`
  and takes about 2.5 min. Baseline: 372 passed, 13 deselected.
- Do not edit `src/speaker/bem/run.py`, `src/speaker/bem/solver.py` or `src/speaker/constants.py`.
  They are the BEM cache key.
- Do not change `speaker.evaluate.SCORING_VERSION` (3). Passive work does not change round scoring.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. Never push or merge.
- Complex convention is e^{+jωt} throughout (`s = j f / fc`), as in `speaker.crossover`.
- Minimum speaker impedance floor: **3.5 Ω** on the modelled |Z_in| at nominal values, 20 Hz–20 kHz.
- Tolerance: **200** draws, **±5 %** uniform on every L and C, fixed seed **0**. The network fails if
  any draw's min |Z_in| < **3.3 Ω**, or if the 95th-percentile flat23 is more than **0.5 dB** above the
  nominal flat23. On failure, refit with the floor raised by **0.1 Ω**, up to **3.8 Ω**.
- Standard values: capacitors E12, resistors E24, inductors in **0.1 mH** steps below 1 mH and
  **0.5 mH** steps from 1 mH. No part snaps to zero.
- Inductor DCR model *(spec §5.5)*: air core **0.5 Ω·√(L / 1 mH)**, cored **0.14 Ω·√(L / 1 mH)**.
  Bass low-pass inductors are cored; all others are air core.
- Capacitor types: polypropylene film in the coax network; bipolar electrolytic allowed above **30 µF**.
- Resistor rating: the largest sine dissipation over 20 Hz–20 kHz at **amp_v_rms / 3** (params
  `[drive] amp_v_rms`, 25 V), rounded up to the next of 5, 10, 15, 20, 30, 50 W. The full-voltage figure
  is listed beside it.
- Crossover grid and pick rule are reused unchanged: `xo_fit.BASS_COAX_GRID`, `xo_fit.LF_HF_GRID`
  floored at `Driver.hf_xo_min`, and `xo_fit.pick`.
- Polarity: bass and coax HF may each be inverted. Coax LF is the reference and is never inverted.

## Review Focus

1. **A bass pair below the floor on its own** (3.26 Ω at 147 Hz, under 3.5 Ω). The search must still
   return a result: either a parallel network that clears the floor, or the woofers-in-series
   alternative with a plain flag. It must not raise. Pinned by
   `test_search_falls_back_to_series_woofers_when_no_parallel_network_clears_the_floor` (Task 4).
2. **Bass data that is exactly zero above the trusted limit** (`Ways` stops the bass field at tl).
   Transfers, sums, dB plots and the fit residual must not produce inf or NaN warnings that leak into
   scores. Pinned by `test_fit_ignores_exact_zero_bass_above_the_limit` (Task 3).
3. **Semi-active must reuse the passive coax values exactly** (one set of boards). Pinned by
   `test_semi_active_netlist_keeps_every_coax_value_and_drops_the_coax_high_pass` (Task 2).
4. **Snapping that removes a part or leaves a non-standard value.** Pinned by
   `test_snap_never_returns_zero_and_lands_on_the_series` (Task 1).
5. **The tolerance run must be repeatable,** so reports don't change between runs. Pinned by
   `test_tolerance_is_deterministic` (Task 4).

---

### Task 1: Circuit solver, standard values and inductor losses

**Files:**
- Create: `src/speaker/network.py`
- Test: `tests/test_network.py`

**Interfaces:**
- Produces:
  - `Part` (frozen dataclass): `name: str`, `kind: Literal["R", "L", "C", "Z"]`, `a: str`, `b: str`,
    `value: float` (ohm, H or F; ignored for "Z"), `dcr: float = 0.0` (ohm, series loss of an L).
    Kind "Z" is a driver load whose impedance comes from the `loads` argument, keyed by `name`.
  - `solve(parts: list[Part], f: np.ndarray, loads: dict[str, np.ndarray]) -> Solution`. A 1 V source
    drives node `"in"` against ground `"0"`.
  - `Solution` (dataclass): `v: dict[str, np.ndarray]` (complex node voltages, "0" included) and
    `z_in: np.ndarray` (complex input impedance), with `across(p: Part) -> np.ndarray` returning
    `v[p.a] - v[p.b]`.
  - `snap(kind: Literal["R", "L", "C"], value: float) -> float`, `neighbours(kind, value) -> tuple[float, float]`
    (the standard values just below and just above a standard value),
    `inductor_dcr(value_h: float, cored: bool) -> float`.
  - Constants: `E12`, `E24` (mantissas), `AIR_DCR = 0.5`, `CORED_DCR = 0.14` (ohm at 1 mH).

- [ ] **Step 1: Write the failing tests**

```python
def test_rc_low_pass_matches_the_analytic_transfer():
    f = np.geomspace(20, 20000, 200)
    parts = [Part("R1", "R", "in", "o", 8.0), Part("C1", "C", "o", "0", 10e-6)]
    s = solve(parts, f, {})
    expect = 1 / (1 + 2j * np.pi * f * 8.0 * 10e-6)
    assert np.allclose(s.v["o"], expect, rtol=1e-9)
    assert np.allclose(s.z_in, 8.0 + 1 / (2j * np.pi * f * 10e-6), rtol=1e-9)

def test_inductor_dcr_is_in_series_with_the_inductance():
    f = np.array([1.0, 1000.0])
    s = solve([Part("L1", "L", "in", "o", 1e-3, dcr=0.5), Part("D", "Z", "o", "0", 0.0)], f,
              {"D": np.full(2, 8.0 + 0j)})
    assert np.allclose(s.z_in, 8.5 + 2j * np.pi * f * 1e-3)
    assert np.allclose(s.across(Part("D", "Z", "o", "0", 0.0)), 8.0 / s.z_in)

def test_snap_never_returns_zero_and_lands_on_the_series():
    assert snap("C", 9.4e-6) == pytest.approx(10e-6)
    assert snap("C", 1e-9) > 0
    assert snap("R", 3.05) == pytest.approx(3.0)
    assert snap("L", 0.43e-3) == pytest.approx(0.4e-3)
    assert snap("L", 2.3e-3) == pytest.approx(2.5e-3)
    assert snap("L", 1e-7) == pytest.approx(0.1e-3)
    lo, hi = neighbours("C", 10e-6)
    assert lo == pytest.approx(8.2e-6) and hi == pytest.approx(12e-6)

def test_inductor_dcr_model():
    assert inductor_dcr(1e-3, cored=False) == pytest.approx(0.5)
    assert inductor_dcr(4e-3, cored=True) == pytest.approx(0.28)
```

- [ ] **Step 2: Run them to see them fail**

Run: `uv run pytest tests/test_network.py -q`
Expected: FAIL (ModuleNotFoundError: speaker.network).

- [ ] **Step 3: Implement `network.py`**

Use modified nodal analysis, vectorised over frequency. Build an (nf, n, n) admittance matrix. Each
L with dcr > 0 is the impedance `dcr + jωL`, not an internal node. Add one extra row for the source.
Solve with `np.linalg.solve`. An unknown load name raises `ValueError` naming it.

- [ ] **Step 4: Run the tests and see them pass**

Run: `uv run pytest tests/test_network.py -q`
Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/network.py tests/test_network.py
git commit -m "feat: nodal circuit solver, standard part values and inductor losses"
```

---

### Task 2: The passive topology and its transfers

**Files:**
- Create: `src/speaker/passive.py`
- Test: `tests/test_passive.py`

**Interfaces:**
- Consumes: `network.Part`, `network.solve`, `network.inductor_dcr`.
- Produces:
  - `Topology` (frozen dataclass), with these fields and defaults: `bass_order: int = 2` (2 or 3),
    `zobel: bool = False`, `lf_order: int = 2` (2 or 3), `lf_notch: bool = False`,
    `hf_shunt: bool = False`, `bass_inverted: bool = False`, `hf_inverted: bool = False`,
    `bass_wiring: Literal["parallel", "series"] = "parallel"`.
  - `value_names(topo: Topology, mode: Literal["passive", "semi_active"]) -> list[str]`. The names, in
    this order where present:
    - bass: `L_b1, C_b1, L_b2` (3rd order), `R_bz, C_bz` (Zobel);
    - coax high-pass: `C_h1, L_h1`;
    - coax LF: `R_ls, R_lp` (L-pad), `L_l1, C_l1, L_l2` (3rd order), `L_ln, C_ln, R_ln` (notch: a
      parallel LCR in series with the LF);
    - coax HF: `C_t1, L_t1, C_t2`, `R_ts, R_tp` (L-pad), `R_tr, L_tr, C_tr` (series RLC across the HF).
    
    semi_active has no bass and no `C_h1, L_h1`.
  - `build(values: dict[str, float], topo: Topology, mode) -> list[Part]`. The node layout:
    - `in` → bass network → load `"bass"`;
    - `in` → (passive only: series `C_h1`, shunt `L_h1`) → node `cx`;
    - `cx` → LF L-pad and low-pass → load `"coax_lf"`;
    - `cx` → `C_t1`, shunt `L_t1`, `C_t2` → HF L-pad → load `"coax_hf"`.
    
    Inductors named `L_b*` get cored DCR and all others air-core DCR (`inductor_dcr`). Loads are Parts
    of kind "Z" named `bass`, `coax_lf`, `coax_hf`.
  - `driver_loads(z: dict[str, np.ndarray], topo: Topology) -> dict[str, np.ndarray]`. It maps
    `speaker.vituixcad.impedances` keys to load names: bass = `bass_pair` (parallel) or
    `4 * bass_pair` (series).
  - `transfers(parts, f, loads, topo) -> tuple[dict[str, np.ndarray], np.ndarray]`. The first item is
    each way's terminal voltage per input volt, polarity applied. In semi_active there is no "bass" key.
    The second item is z_in.
  - `bass_scale(topo) -> float`: 1.0 for parallel, 0.5 for series. It scales `Ways.ways["bass"]`,
    because in series each woofer sees half the voltage.
  - `BOUNDS: dict[str, tuple[float, float]]` per value name prefix: L 0.05–20 mH, C 0.5–400 µF,
    R 0.1–50 Ω; `R_lp` and `R_tp` 0.5–100 Ω.

- [ ] **Step 1: Write the failing tests**

```python
F = np.geomspace(20, 20000, 300)
R8 = {k: np.full(len(F), 8.0 + 0j) for k in ("bass", "coax_lf", "coax_hf")}

# Hand-written values for every name value_names can return (all options on), inside BOUNDS.
# Task 6's tests reuse it: put it in tests/test_passive.py at module level.
SAMPLE = {"L_b1": 6.5e-3, "C_b1": 100e-6, "L_b2": 2.0e-3, "R_bz": 6.8, "C_bz": 47e-6,
          "C_h1": 100e-6, "L_h1": 8.0e-3, "R_ls": 1.5, "R_lp": 15.0, "L_l1": 1.0e-3, "C_l1": 8.2e-6,
          "L_l2": 0.5e-3, "L_ln": 0.3e-3, "C_ln": 10e-6, "R_ln": 4.7, "C_t1": 10e-6, "L_t1": 0.6e-3,
          "C_t2": 15e-6, "R_ts": 4.7, "R_tp": 3.3, "R_tr": 6.8, "L_tr": 0.4e-3, "C_tr": 4.7e-6}

def test_second_order_bass_low_pass_falls_40_db_per_decade():
    t = Topology()
    v = {n: SAMPLE[n] for n in value_names(t, "passive")}
    h, _ = transfers(build(v, t, "passive"), F, R8, t)
    db = 20 * np.log10(np.abs(h["bass"]))
    i1, i2 = np.searchsorted(F, [2000, 20000])
    assert db[i1] - db[i2] == pytest.approx(40, abs=4)

def test_inverting_a_branch_flips_only_its_sign(): ...  # h_inv["coax_hf"] == -h["coax_hf"], others equal

def test_semi_active_netlist_keeps_every_coax_value_and_drops_the_coax_high_pass():
    t = Topology(lf_notch=True, hf_shunt=True)
    v = {n: SAMPLE[n] for n in value_names(t, "passive")}
    p, s = build(v, t, "passive"), build(v, t, "semi_active")
    names = {x.name for x in s}
    assert not names & {"C_h1", "L_h1", "L_b1", "C_b1", "bass"}
    for x in s:
        if x.kind != "Z":
            assert x.value == next(y.value for y in p if y.name == x.name)

def test_series_woofers_quadruple_the_load_and_halve_the_drive(): ...
    # driver_loads(... series)["bass"] == 4 * z["bass_pair"]; bass_scale == 0.5
```

Fill in the two `...` bodies with the asserts in their comments.

- [ ] **Step 2: Run them to see them fail**

Run: `uv run pytest tests/test_passive.py -q`
Expected: FAIL (ModuleNotFoundError: speaker.passive).

- [ ] **Step 3: Implement the topology, `build`, `driver_loads`, `transfers`, `bass_scale`, `BOUNDS`**

- [ ] **Step 4: Run the tests and see them pass**

Run: `uv run pytest tests/test_passive.py tests/test_network.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/passive.py tests/test_passive.py
git commit -m "feat: passive crossover topology, netlists for passive and semi-active modes"
```

---

### Task 3: Fitting the values at one pair of crossover points

**Files:**
- Modify: `src/speaker/passive.py`
- Test: `tests/test_passive.py`

**Interfaces:**
- Consumes: Task 2's interfaces; `speaker.curves.smooth`; `speaker.system.Ways`.
- Produces:
  - `Z_GRID = np.geomspace(20.0, 20000.0, 300)`: where the impedance floor is checked.
  - `FIT_LO = 100.0`, `FIT_TOP = 10000.0`, `ABOVE_WEIGHT = 0.25`, `FIT_SMOOTH = 1/3` octave,
    `SHUNT_Q_MAX = 2.0`.
  - `PassiveFit` (frozen dataclass): `bass_coax: float`, `lf_hf: float`, `topo: Topology`,
    `values: dict[str, float]` (snapped), `continuous: dict[str, float]`, `z_floor: float`,
    `z_min: float` (min |z_in| on Z_GRID at the snapped values), `z_ok: bool` (z_min ≥ z_floor),
    `cost: float`.
  - `start_values(topo, fb, fh, z_at: dict[str, np.ndarray]) -> dict[str, float]`. The textbook
    Butterworth values at fb and fh on each load's |Z| at that frequency, the L-pads from the level
    difference, and the notch at 2.7 kHz.
  - `passive_total(ways: Ways, fit: PassiveFit, z: Callable[[np.ndarray], dict[str, np.ndarray]], mode="passive")
    -> np.ndarray`: the (nf, nd) complex sum, using `bass_scale`. In semi_active there is no bass.
  - `fit_passive_at(ways, z, tl, fb, fh, topo, z_floor=3.5) -> PassiveFit`. Here `z` is a callable
    `f -> {"coax_lf", "coax_hf", "bass_pair", ...}` complex impedances (in production
    `lambda f: impedances(d, f)`).

The objective, which the signature does not fix:
- Residual: the `FIT_SMOOTH`-smoothed dB of the summed 23° response (`ways.obs.at23_idx`) minus its
  weighted mean. The level is free, so only the shape counts.
- Weights: 1 on [FIT_LO, tl] and ABOVE_WEIGHT on (tl, FIT_TOP]. Points where the sum is exactly 0 get
  weight 0.
- Penalty: `10 * max(0, z_floor - |z_in|)` on Z_GRID, appended to the residual.
- Variables are the logs of the values, inside `BOUNDS`, fitted with `scipy.optimize.least_squares`
  from `start_values`.
- The RLC shunt's Q, `sqrt(L/C)/R`, is held at SHUNT_Q_MAX or less by a penalty on its excess.
- Snap: snap every value, then make coordinate-descent passes. Each pass tries each value's two
  `neighbours` and keeps a move that lowers the cost without breaking the floor. Stop when a pass
  changes nothing.

- [ ] **Step 1: Write the failing tests**

```python
def omni(f, level_db=None):  # every way flat 0 dB at 2 m in every direction, like test_xo_fit.omni_ways
    ...

def test_fit_on_flat_resistive_ways_is_flat_and_safe():
    w = omni(np.geomspace(100, 16000, 160))
    z = lambda f: {"coax_lf": np.full(len(f), 8 + 0j), "coax_hf": np.full(len(f), 8 + 0j),
                   "bass_pair": np.full(len(f), 8 + 0j)}
    fit = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, Topology())
    total = passive_total(w, fit, z)
    db = 20 * np.log10(np.abs(total[:, w.obs.at23_idx]))
    m = (w.f >= 150) & (w.f <= 2851)
    assert np.ptp(smooth(w.f, db, 1/3)[m]) < 3.0
    assert fit.z_ok and fit.z_min >= 3.5
    assert all(v == pytest.approx(snap(n[0], v)) for n, v in fit.values.items())

def test_floor_is_respected_or_flagged():
    # bass_pair 3.0 ohm resistive: the fit either clears 3.5 (z_ok) or reports z_ok False with z_min < 3.5
    ...

def test_fit_ignores_exact_zero_bass_above_the_limit():
    w = omni(...); w.ways["bass"][w.f > 2851] = 0
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        fit = fit_passive_at(w, z, 2851.0, 200.0, 1600.0, Topology())
    assert np.isfinite(fit.cost)
```

The `omni` helper builds `Ways` the way `tests/test_xo_fit.py::omni_ways` does, at `observation_set`'s
2 m. `snap(n[0], v)` relies on the first letter of each value name being its kind.

- [ ] **Step 2: Run them to see them fail**

Run: `uv run pytest tests/test_passive.py -q -k "fit or floor"`
Expected: FAIL (ImportError: fit_passive_at).

- [ ] **Step 3: Implement `start_values`, `passive_total`, `fit_passive_at` and the constants above**

- [ ] **Step 4: Run the tests and see them pass**

Run: `uv run pytest tests/test_passive.py -q`
Expected: all pass, in under 60 s.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/passive.py tests/test_passive.py
git commit -m "feat: fit passive component values at a crossover pair, snapped to standard values"
```

---

### Task 4: The search, the tolerance check and the series fallback

**Files:**
- Modify: `src/speaker/passive.py`
- Test: `tests/test_passive.py`

**Interfaces:**
- Consumes: Task 3; `xo_fit.BASS_COAX_GRID`, `xo_fit.LF_HF_GRID`, `xo_fit.SEARCH_KEYS`, `xo_fit.pick`;
  `scores.score_all`, `spinorama.spinorama`, `system.coax_section`.
- Produces:
  - `passive_scores(ways, fit, z, tl, mode="passive") -> list[Score]`. These are `score_all` on the
    passive total, with the same arguments `xo_fit.fit_active` passes.
  - `OPTIONS = ("zobel", "bass_order", "lf_order", "lf_notch", "hf_shunt")`. For the two order fields,
    "on" means 3.
  - `Tolerance` (frozen dataclass): `z_min: float` (worst draw), `flat23_nominal: float`,
    `flat23_p95: float`, `whole_rear: tuple[float, float]` (5th, 95th percentile), `passed: bool`.
  - `tolerance(ways, z, tl, fit, n=200, spread=0.05, seed=0) -> Tolerance`.
  - `PassiveResult` (frozen dataclass): `fit: PassiveFit`, `tolerance: Tolerance`,
    `series: PassiveResult | None` (the series-woofer alternative, set only when the parallel fit is
    not z_ok), `notes: list[str]` (one plain sentence per decision: options kept or dropped, floor
    raised, fallback).
  - `search_passive(ways, z, tl, lf_hf_min: float | None, z_floor=3.5) -> PassiveResult`.

The search order, which the signature does not fix:
1. Grid stage, base Topology. For every grid pair, try every polarity pair (bass ±, HF ±) and keep the
   one with the lowest cost. Score each kept fit; candidates without all SEARCH_KEYS trusted scores
   are dropped. Choose the points with `xo_fit.pick`.
2. Options, at the chosen points, in `OPTIONS` order. Refit with the option on. Keep it if trusted
   flat23 improves by at least 0.1 dB and no other SEARCH_KEYS score gets worse by more than 0.25 dB
   (higher whole_rear is better; every other key, lower is better).
3. Tolerance. While it fails and the floor is below 3.8 Ω, raise the floor by 0.1 and refit at the
   same points and topology.
4. Fallback. If the final fit is not z_ok, run steps 1–3 again with `bass_wiring="series"` and attach
   the result as `series`.

Each step adds a sentence to `notes`.

- [ ] **Step 1: Write the failing tests**

```python
def test_search_returns_a_grid_point_at_or_above_the_floor(): ...
    # omni ways, 8 ohm loads, lf_hf_min=1200: fit.bass_coax in BASS_COAX_GRID, fit.lf_hf in LF_HF_GRID and >= 1200

def test_search_falls_back_to_series_woofers_when_no_parallel_network_clears_the_floor():
    # bass_pair 2.0 ohm resistive: result.fit.z_ok is False, result.series is not None,
    # result.series.fit.topo.bass_wiring == "series", and a note mentions "series"
    ...

def test_tolerance_is_deterministic():
    a, b = tolerance(w, z, 2851.0, fit), tolerance(w, z, 2851.0, fit)
    assert a == b

def test_tolerance_fails_below_3_3_ohm(): ...
    # a fit whose nominal z_min is 3.35 ohm with a 5 % spread: passed is False
```

Speed: these tests may monkeypatch `xo_fit.BASS_COAX_GRID` and `LF_HF_GRID` down to two points each,
and use `n=20` draws, to keep each test under 60 s. Leave the production constants alone.

- [ ] **Step 2: Run them to see them fail**

Run: `uv run pytest tests/test_passive.py -q -k "search or tolerance"`
Expected: FAIL (ImportError: search_passive).

- [ ] **Step 3: Implement `passive_scores`, `OPTIONS`, `Tolerance`, `tolerance`, `PassiveResult`, `search_passive`**

- [ ] **Step 4: Run the tests and see them pass**

Run: `uv run pytest tests/test_passive.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/passive.py tests/test_passive.py
git commit -m "feat: passive crossover search, tolerance check and series-woofer fallback"
```

---

### Task 5: Semi-active mode: the DSP stage in front of the passive coax boards

**Files:**
- Create: `src/speaker/semi_active.py`
- Test: `tests/test_semi_active.py`

**Interfaces:**
- Consumes:
  - from Tasks 2–4: `passive.build`, `passive.transfers`, `passive.driver_loads`, `PassiveFit`;
  - from `crossover`: `lr4`, `peq`, `Peq`, `WayFilter`, `biquad_response`, and the private
    `_rbj_butterworth` and `_rbj_peaking` (make these two public as `rbj_butterworth` and
    `rbj_peaking`, and update their callers in `crossover.py`);
  - from `xo_fit`: `fit_peq`, `_align` (make it public as `align` and update its caller),
    `BASS_COAX_GRID`, `BASS_DELAYS`, `N_PEQ["bass"]`, `eq_band`, `pick`, `SEARCH_KEYS`, `EXPORT_FS`.
- Produces:
  - `SemiActive` (frozen dataclass): `bass_coax: float`, `bass: WayFilter` (gain, delay, polarity, up to
    2 PEQs), `coax: WayFilter` (gain only; delay 0, polarity +1, no PEQ).
  - `semi_ways(ways, z, fit: PassiveFit) -> Ways`. This is a two-way `Ways` with keys `"bass"`
    (unfiltered, with `bass_scale` applied) and `"coax"`, the LF and HF through the semi_active
    network, summed.
  - `semi_total(sw: Ways, sa: SemiActive) -> np.ndarray`.
  - `fit_semi_active(ways, z, tl, fit: PassiveFit) -> SemiActive`. It covers `BASS_COAX_GRID`. At each
    point:
    1. Set the coax gain to 0 dB and the bass gain to match the coax.
    2. Fit the bass PEQs with `fit_peq` over `eq_band("bass", fb, ·)` towards the coax level over
       300–1000 Hz.
    3. Align the bass delay and polarity with `align`.
    4. Score.
    
    Choose with `pick` on SEARCH_KEYS.
  - `semi_active_text(sa: SemiActive, tl, *, design: str, sha: str) -> str`. A DSP sheet in the style
    of `xo_fit.crossover_text`: '#' header lines, rows of gain, delay and polarity per channel, the
    bass PEQs, then the biquads for each of `EXPORT_FS`. Its header says the coax channel feeds the
    passive coax boards with the coax high-pass links fitted.

- [ ] **Step 1: Write the failing tests**

```python
def test_semi_active_recovers_a_known_bass_delay():
    # omni ways with the bass way delayed by 0.8 ms; the fit's bass delay is the compensating shift
    # (coax delayed or bass advanced) within 0.05 ms, relative delays as in xo_fit
    ...

def test_semi_active_text_has_two_channels_and_biquads_at_both_rates():
    txt = semi_active_text(sa, 2851.0, design="round-3.toml", sha="abcd1234")
    assert "bass\tdelay_ms" in txt and "coax\tgain" in txt and "# fs = 48000" in txt and "# fs = 96000" in txt
    assert "high-pass links" in txt
```

- [ ] **Step 2: Run them to see them fail**

Run: `uv run pytest tests/test_semi_active.py -q`
Expected: FAIL (ModuleNotFoundError).

- [ ] **Step 3: Rename the three private helpers, then implement `semi_active.py`**

- [ ] **Step 4: Run the semi-active, crossover and xo_fit tests and see them pass**

Run: `uv run pytest tests/test_semi_active.py tests/test_crossover.py tests/test_xo_fit.py -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/semi_active.py src/speaker/crossover.py src/speaker/xo_fit.py tests/test_semi_active.py
git commit -m "feat: semi-active mode: DSP bass-to-coax stage ahead of the passive coax boards"
```

---

### Task 6: Schematics with schemdraw

**Files:**
- Modify: `pyproject.toml` (add `"schemdraw>=0.19"` to `dependencies`), `uv.lock` (`uv sync`)
- Create: `src/speaker/schematic.py`
- Test: `tests/test_schematic.py`

**Interfaces:**
- Consumes: `network.Part`, `passive.build`, `passive.Topology`.
- Produces:
  - `draw_network(parts: list[Part], title: str, path: Path) -> Path`. An SVG of one sub-network.
    Every part is labelled with its name and value in engineering units (µF, mH, Ω) and, for an
    inductor, its DCR. Driver loads are drawn as speakers.
  - `draw_links(topo: Topology, path: Path) -> Path`. An SVG table-diagram of the three modes:
    passive, 1 terminal pair; semi-active, 2 (coax, bass), with the coax high-pass links fitted;
    fully active, 3 (HF, LF, bass), with every network bypassed.
  - `fmt_value(kind: str, value: float) -> str`, e.g. `fmt_value("C", 6.8e-6) == "6.8 µF"`,
    `fmt_value("L", 2.5e-3) == "2.5 mH"`, `fmt_value("R", 3.3) == "3.3 Ω"`.

- [ ] **Step 1: Write the failing tests**

```python
def test_network_svg_names_every_part(tmp_path):
    t = Topology(hf_shunt=True)
    parts = build({n: v for n, v in SAMPLE.items()}, t, "passive")  # SAMPLE: a hand-written value dict
    svg = draw_network([p for p in parts if p.name.startswith(("C_t", "L_t", "R_t")) or p.name == "coax_hf"],
                       "Coax HF", tmp_path / "hf.svg").read_text()
    for p in parts:
        if p.name.startswith(("C_t", "L_t", "R_t")):
            assert p.name.replace("_", "") in svg.replace("_", "")

def test_links_svg_lists_three_modes(tmp_path):
    svg = draw_links(Topology(), tmp_path / "links.svg").read_text()
    assert "Passive" in svg and "Semi-active" in svg and "Fully active" in svg

def test_fmt_value(): ...  # the three examples above
```

schemdraw writes text as SVG `<text>` when `schemdraw.config(svg2=True)` (its default backend for
SVG). If labels come out as paths, set `schemdraw.use("svg")` so the test can find the text.

- [ ] **Step 2: Run them to see them fail**

Run: `uv sync && uv run pytest tests/test_schematic.py -q`
Expected: FAIL (ModuleNotFoundError: speaker.schematic).

- [ ] **Step 3: Implement `schematic.py`**

The layout is a ladder: series parts along the top, shunt parts down to a ground rail, and the load
at the right.

- [ ] **Step 4: Run them and see them pass**

Run: `uv run pytest tests/test_schematic.py -q`
Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml uv.lock src/speaker/schematic.py tests/test_schematic.py
git commit -m "feat: passive crossover schematics and the three-mode link diagram"
```

---

### Task 7: The `speaker passive` command and its report

**Files:**
- Modify: `src/speaker/directivity_report.py`. Extract the field and ways preparation at the top of
  `run_directivity` (from `load_design` to `rec`) into `design_ways(params_path, *, trusted_limit) -> DesignWays`.
  `run_directivity` then calls it, with no change in behaviour.
- Modify: `src/speaker/cli.py`. Add the `passive` sub-command with `--params`, `--out` (default
  `out/passive/`) and `--trusted-limit HZ`, the same as `directivity`.
- Create: `src/speaker/passive_report.py`
- Test: `tests/test_passive_report.py`, `tests/test_cli.py`

**Interfaces:**
- Consumes: Tasks 1–6; `xo_fit.fit_active`; `vituixcad.impedances`; `directivity_plots` (polar_map,
  crossover_23, whole_rear); `evaluate.ways_scores`.
- Produces:
  - `DesignWays` (dataclass): `d: Design`, `f`, `ways: Ways`, `tl: float`, `tl_text: str`, `cab`,
    `puck`, `rec`.
  - `power_ratings(parts, z, f, v_rms) -> dict[str, tuple[float, float, int]]`: per resistor, the
    dissipation at v_rms/3, the dissipation at v_rms, and the chosen standard rating (W).
  - `network_text(result: PassiveResult, d: Design, tl) -> str`, `parts_csv(result, d) -> str`
    (quantities for a pair) and `vituixcad_netlist(result) -> str`.
  - `run_passive(params_path, out_dir, *, trusted_limit: tuple[float, str] | None) -> Path`. It writes
    the files below and returns the report path.

Files written to `out_dir`:
- `network.txt`:
  - per mode, the netlist: part, value, type, tolerance, DCR, and power rating at 1/3 and at full
    voltage;
  - which links to fit for each mode;
  - the notes, the tolerance result and the floor used;
  - with a series fallback, both networks.
- `parts.csv`: columns `part,value,unit,type,tolerance,dcr_ohm,rating_w,qty_pair`.
- `semi_active.txt`: from `semi_active_text`.
- `vituixcad_netlist.txt`: one line per part, `name kind value node_a node_b`, in build order, with a
  header that explains the node names.
- `scores.json`: `{"passive": [...], "semi_active": [...], "active": [...]}` (Score dicts).
- `schematic_bass.svg`, `schematic_coax.svg`, `links.svg`.
- `plots/`:
  - `crossover_23_<mode>.png` for each mode;
  - `impedance.png`: |Z_in| for passive with a 3.5 Ω line, plus the semi-active coax channel;
  - `polar_h_passive.png` and `polar_v_passive.png`;
  - `whole_rear_compare.png`: the three modes on one axis;
  - `tolerance.png`: a flat23 histogram and the min |Z| per draw.
- `report.html`:
  - the three modes' trusted scores side by side;
  - the plots and the schematics, inlined as `<img>`;
  - "Where passive falls short of active", plain sentences computed from the score differences of
    0.25 dB or more, and from `xo_fit.passive_warnings` of the active fit;
  - the tolerance result;
  - the series fallback, when present;
  - a warning that in semi-active mode the bass channel drives the bare woofer pair (its minimum |Z| from
    `driver_loads`, 3.26 Ω for the current design), so that amp must handle it. This is a warning, not a
    pass or fail.

  It uses the existing report CSS.

- [ ] **Step 1: Write the failing tests**

```python
def test_power_rating_rounds_up_to_the_standard_series():
    parts = [Part("R_ts", "R", "in", "o", 4.0), Part("coax_hf", "Z", "o", "0", 0.0)]
    f = np.geomspace(20, 20000, 50)
    r = power_ratings(parts, {"coax_hf": np.full(50, 4 + 0j)}, f, 25.0)["R_ts"]
    # 25/3 V across 8 ohm -> 1.04 A -> 4.34 W in R_ts; full voltage 39.1 W
    assert r[0] == pytest.approx(4.34, abs=0.05) and r[1] == pytest.approx(39.06, abs=0.1) and r[2] == 5

def test_run_passive_writes_every_output(tmp_path, monkeypatch):
    # monkeypatch passive_report.design_ways to return synthetic omni DesignWays (as in tests/test_passive.py)
    # and impedances to 8 ohm resistive; shrink the grids as in Task 4's tests
    out = run_passive(None, tmp_path, trusted_limit=(2851.0, "test"))
    for name in ("network.txt", "parts.csv", "semi_active.txt", "vituixcad_netlist.txt", "scores.json",
                 "schematic_bass.svg", "schematic_coax.svg", "links.svg", "report.html"):
        assert (tmp_path / name).exists()
    assert set(json.loads((tmp_path / "scores.json").read_text())) == {"passive", "semi_active", "active"}
    assert "3.5" in (tmp_path / "network.txt").read_text()

def test_cli_passive_subcommand_parses(monkeypatch): ...  # main(["passive", "--trusted-limit", "2851"]) calls run_passive with (2851.0, "given on the command line")
```

- [ ] **Step 2: Run them to see them fail**

Run: `uv run pytest tests/test_passive_report.py tests/test_cli.py -q`
Expected: FAIL (ModuleNotFoundError: speaker.passive_report).

- [ ] **Step 3: Extract `design_ways`, then implement `passive_report.py` and the CLI sub-command**

- [ ] **Step 4: Run the full suite**

Run: `uv run pytest -q`
Expected: every test passes (the 372 baseline plus the new ones), 13 deselected.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/directivity_report.py src/speaker/cli.py src/speaker/passive_report.py tests/test_passive_report.py tests/test_cli.py
git commit -m "feat: speaker passive: passive and semi-active crossovers, parts list, schematics and report"
```

---

### Task 8: Run it on the design and write round 3d

**Files:**
- Create: `docs/rounds/round-3d.md`

- [ ] **Step 1: Run the command on the current design**

Run: `uv run speaker passive --params params/round-3.toml --out out/passive/round-3-design --trusted-limit 2851 2>&1 | tee out/passive/run.log`
Expected: `wrote out/passive/round-3-design/report.html`, with the BEM fields read from the cache
(no new solves; minutes, not hours).

- [ ] **Step 2: Check the outputs by eye**

Open the files and confirm that:
- every SVG renders, with readable labels;
- `impedance.png` shows the 3.5 Ω line;
- the three modes in `scores.json` have every trusted SEARCH_KEYS score;
- `network.txt` values are standard values.

Write anything wrong into the report under "Problems found", and fix it before continuing if it's a
code defect.

- [ ] **Step 3: Write `docs/rounds/round-3d.md`**

Follow the style of `docs/rounds/round-3c.md`:
- what was run (commit, command, time);
- the passive network (points, topology, every value, notes);
- the semi-active DSP settings;
- a trusted-score table of passive, semi-active and active;
- the minimum impedance and the tolerance result;
- the series fallback, if any;
- where passive falls short and why;
- what stays open.

Outputs are not committed.

The owner's report page (the claude.ai artifact) gets its "Passive" option from the controller after
this task. It is not part of this task.

- [ ] **Step 4: Commit**

```bash
git add docs/rounds/round-3d.md
git commit -m "docs: round 3d: the passive and semi-active crossovers of bass260-580"
```
