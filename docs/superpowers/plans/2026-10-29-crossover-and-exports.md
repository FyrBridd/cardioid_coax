# Crossover and Exports Implementation Plan (plan 3 of 3)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the fixed 170 Hz / 1.4 kHz crossover with an active crossover fitted to each variant (points, LR4 slopes, delays, polarity, parametric EQ), export it as DSP data, and write VituixCAD FRD/ZMA sets for any full run.

**Architecture:** `crossover.py` gains the filter blocks and an `ActiveCrossover` that `apply_active` sums; a new `xo_fit.py` fits one per variant (EQ per way, then delay and polarity per crossover, then a grid search of the two points ranked by the existing directivity key). `evaluate` (screen stages) and `run_directivity` (finals, single runs) use it; the vents stage keeps the fixed crossover. A new `vituixcad.py` writes FRD sets from `Ways` and ZMA files from the lumped models and the datasheet traces; `run_directivity` calls it. Last, round 3 is re-scored and the chosen design (`params/round-3.toml`) is exported.

**Tech Stack:** Python 3 via `uv run`, numpy, scipy (already a dependency). No new dependencies.

**Spec:** `docs/superpowers/specs/2026-10-26-directivity-optimisation-design.md` (§6.3, §6.4, §7 item 5, §9 item 3), with the main spec `docs/superpowers/specs/2026-10-22-cardioid-coax-floorstander-design.md` §7.2 (active-mode deliverables, passive-mode requirements).

## Global Constraints

- Complex quantities are e^{+jωt}; way pressures are complex PEAK Pa at 2 m for 2.83 V RMS into the way (`speaker.system` docstring).
- Crossover point ranges: bass-to-coax 150–250 Hz, coax LF-to-HF 1.0–2.0 kHz (spec §5).
- Acoustic target slopes LR4, plus a delay per way, plus parametric EQ fitted to make the 23° response flat (spec §6.3).
- Ranking and every choice the fit makes use the trusted band only (≤ the trusted limit); scores above it stay low confidence (spec §4).
- Passive-mode flags: EQ needing more than about 6 dB of boost, or what a passive network cannot copy (spec §6.3).
- FRD sets per way, horizontal and vertical, 5° steps; the coax LF set includes the vents; bass each driver and the pair; ZMA per way including the chamber's effect on the coax LF; a header marking each file simulated with the confidence band; FRDs at 2 m, free space, with a note on VituixCAD settings (spec §6.4).
- Do not edit `src/speaker/bem/run.py`, `src/speaker/bem/solver.py` or `src/speaker/constants.py`: they make up the BEM cache key (`SOLVER_MODULES`), and every cached field must stay valid.
- `params/default.toml` stays the round 0 design (Ruling below).
- Tests run with `uv run pytest -q`; BEM tests stay deselected by the `bem` marker.

**Rulings made while planning (the owner may overrule at review):**

- Ruling: `params/default.toml` is not overwritten with the chosen design, against spec §2's wording — 24 test files pin round 0 values from it, and `speaker build --params params/round-3.toml` already builds the CAD, DXF and STEP from the chosen design — cost if wrong: one task to swap the files and move the tests to `params/round-0.toml`.
- Ruling: filters are evaluated as analog prototypes in the simulation, and exported as digital biquads (bilinear, pre-warped) at 48 and 96 kHz — the scores then do not depend on a DSP's sample rate — cost if wrong: small (< 0.1 dB below 5 kHz, tested).
- Ruling: no EQ boost below a way's own high-pass point (bass: below 80 Hz) — boosting a rolled-off driver costs excursion that the model does not check — cost if wrong: slightly less flat acoustic slopes, which the delay search partly absorbs.
- Ruling: the VituixCAD FRDs cover full circles (−175…180°), a superset of the spec's 0–180°, since the data exist and the cabinet is not symmetric vertically.
- Ruling: the FRD grid is the run's own grid (12 points per octave, 50 Hz–20 kHz); the bass sets stop at the trusted limit, where the bass field stops.

## Review Focus

1. Re-running an existing round directory after this change must not reuse its fixed-crossover scores: a stored result without the new scoring version is re-evaluated (Task 4, `test_old_scoring_results_are_not_reused`).
2. A way whose EQ band holds less than 1/3 octave of grid points (the screen grid stops at the trusted limit, so a 2 kHz LF-HF point leaves the HF band short) gets no EQ rather than an error (Task 2, `test_short_band_gets_no_eq`).
3. The bass field is zero above the trusted limit: no FRD may contain `inf` or `nan`, and the bass sets stop at the limit (Task 7, `test_no_non_finite_values_and_bass_stops_at_the_limit`).
4. Round directories from rounds 1–3 have no `crossover.json`; `--report-only` must still render their crossover section from `gains.json` (Task 5, `test_round_report_without_crossover_json_falls_back_to_gains`).
5. Angle naming: −180° and 180° are the same direction and must give one file, `hor 180`; `ver 0` is not written (VituixCAD ignores it) (Task 7, `test_frd_names_are_unique_and_parseable`).

---

### Task 1: Filter blocks and the active crossover

**Files:**
- Modify: `src/speaker/crossover.py`
- Test: `tests/test_crossover.py`

**Interfaces:**
- Consumes: `Ways` (`speaker.system`), existing `lr4`.
- Produces:
  - `Peq(f0: float, gain_db: float, q: float)` (frozen dataclass).
  - `peq(f: np.ndarray, p: Peq) -> np.ndarray` — analog RBJ peaking: `H(s) = (s² + s·A/Q + 1) / (s² + s/(A·Q) + 1)`, `s = j f/f0`, `A = 10^(gain_db/40)`.
  - `WayFilter(gain_db: float, delay_s: float, polarity: int, peq: tuple[Peq, ...])` (frozen).
  - `ActiveCrossover(bass_coax: float, lf_hf: float, ways: dict[str, WayFilter], fit_top: float)` (frozen; `fit_top` is the top of the data the EQ was fitted on). Property `gains_db -> dict[str, float]`.
  - `way_transfer(f, xo: ActiveCrossover, way: str) -> np.ndarray` — gain × polarity × LR4 sections × Π peq × `exp(-j2πf·delay_s)`. Topology: `bass` low-pass at `bass_coax`; `coax_lf` high-pass at `bass_coax` and low-pass at `lf_hf`; `coax_hf` high-pass at `lf_hf`.
  - `apply_active(w: Ways, xo: ActiveCrossover) -> np.ndarray` — (nf, nd) total.
  - `biquads(xo: ActiveCrossover, fs: float) -> dict[str, list[tuple[str, tuple[float, float, float, float, float]]]]` — per way, `(label, (b0, b1, b2, a1, a2))` with a0 = 1 and `H(z) = (b0 + b1 z⁻¹ + b2 z⁻²) / (1 + a1 z⁻¹ + a2 z⁻²)`. Each LR4 edge is two identical RBJ 2nd-order Butterworth sections (Q = 1/√2), labels `"lr4_hp <fc> Hz (1/2)"` etc.; each PEQ is an RBJ peakingEQ, label `"peq <f0> Hz <gain> dB Q<q>"`. Gain, delay and polarity are not biquads.
  - `biquad_response(c: tuple[float, ...], f: np.ndarray, fs: float) -> np.ndarray`.

- [ ] **Step 1: Write the failing tests**

```python
def test_peq_has_its_gain_at_f0_and_none_far_away():
    f = np.array([20.0, 2000.0, 20000.0])
    h = peq(f, Peq(2000.0, -6.0, 2.0))
    assert 20 * np.log10(abs(h[1])) == pytest.approx(-6.0, abs=1e-9)
    assert abs(20 * np.log10(abs(h[0]))) < 0.05 and abs(20 * np.log10(abs(h[2]))) < 0.05


def _flat_xo(**kw):
    w = {k: WayFilter(0.0, 0.0, 1, ()) for k in ("bass", "coax_lf", "coax_hf")}
    w.update(kw)
    return ActiveCrossover(170.0, 1400.0, w, 20000.0)


def test_active_without_eq_or_delay_matches_the_fixed_crossover():
    o = observation_set(2.0)
    f = np.geomspace(20, 20000, 200)
    w = Ways(f, o, {k: np.ones((len(f), len(o.dirs)), complex) for k in ("bass", "coax_lf", "coax_hf")}, {})
    fixed, _ = apply_crossover(w, FixedCrossover())
    assert np.allclose(apply_active(w, _flat_xo()), fixed, atol=1e-12)


def test_delay_and_polarity():
    f = np.array([1000.0])
    xo = _flat_xo(coax_hf=WayFilter(-3.0, 1e-4, -1, ()))
    h = way_transfer(f, xo, "coax_hf") / lr4(f, 1400.0, "high")
    assert abs(h[0]) == pytest.approx(10 ** (-3 / 20))
    assert np.angle(h[0]) == pytest.approx(np.angle(-np.exp(-2j * np.pi * 1000.0 * 1e-4)))


@pytest.mark.parametrize("fs", [48000.0, 96000.0])
def test_biquads_match_the_analog_filters_below_5k(fs):
    xo = _flat_xo(coax_hf=WayFilter(0.0, 0.0, 1, (Peq(2000.0, 6.0, 2.0),)))
    f = np.geomspace(100, 5000, 50)
    for way in ("bass", "coax_lf", "coax_hf"):
        digital = np.prod([biquad_response(c, f, fs) for _, c in biquads(xo, fs)[way]], axis=0)
        analog = way_transfer(f, xo, way)
        assert np.max(np.abs(20 * np.log10(np.abs(digital / analog)))) < 0.1
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_crossover.py -q`
Expected: FAIL (ImportError: `Peq`, `ActiveCrossover`, … not defined).

- [ ] **Step 3: Implement the Interfaces above in `src/speaker/crossover.py`**

Keep `FixedCrossover` and `apply_crossover` unchanged (the vents stage and the mesh check still use them). RBJ cookbook formulas for the digital sections, with `w0 = 2π f/fs`.

- [ ] **Step 4: Run them to verify they pass**

Run: `uv run pytest tests/test_crossover.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/crossover.py tests/test_crossover.py
git commit -m "feat: active crossover blocks: PEQ, delay, polarity, biquad export"
```

### Task 2: Per-way EQ fit

**Files:**
- Create: `src/speaker/xo_fit.py`
- Test: `tests/test_xo_fit.py`

**Interfaces:**
- Consumes: `Peq`, `peq` (Task 1); `smooth` (`speaker.curves`).
- Produces:
  - Constants: `N_PEQ = {"bass": 2, "coax_lf": 4, "coax_hf": 5}`, `GAIN_RANGE = (-15.0, 12.0)` dB, `Q_RANGE = (0.5, 6.0)`, `EQ_SMOOTH = 1 / 3` octave, `MIN_GAIN = 0.5` dB, `MIN_BAND_OCT = 1 / 3`, `EQ_TOP = 16000.0` Hz, `BASS_EQ_LO = 80.0` Hz.
  - `eq_band(way: str, fb: float, fh: float) -> tuple[float, float]` — `bass`: (80, 2·fb); `coax_lf`: (fb, 2·fh); `coax_hf`: (fh, 16000). No band reaches below the way's own high-pass point (Ruling).
  - `fit_peq(f: np.ndarray, db: np.ndarray, ref: float, band: tuple[float, float], n: int) -> tuple[Peq, ...]` — `db` is the way's 23° level after its gain; the fit flattens `smooth(f, db, EQ_SMOOTH)` towards `ref` over `band` clipped to `f`. Returns `()` when the clipped band spans less than `MIN_BAND_OCT`.

Algorithm (the signature does not determine it): greedy start, then a joint refine.
1. Residual `r = smoothed − ref` on the band points. For each of `n` filters: put `f0` at the band point of largest `|r|`, `gain = clip(−r, GAIN_RANGE)`, `Q = 2`; subtract that filter's dB from `r`.
2. `scipy.optimize.least_squares` over all filters jointly, in variables `(log f0, gain, log Q)` with bounds (`f0` inside the band, `GAIN_RANGE`, `Q_RANGE`), residual `smoothed + Σ peq dB − ref` on the band points.
3. Drop filters with `|gain| < MIN_GAIN`; return the rest sorted by `f0`.

- [ ] **Step 1: Write the failing tests**

```python
def test_fit_removes_a_single_bump():
    f = np.geomspace(100, 10000, 160)
    bump = Peq(1500.0, 5.0, 2.0)
    db = 80 + 20 * np.log10(np.abs(peq(f, bump)))
    eq = fit_peq(f, db, 80.0, (500.0, 5000.0), 3)
    after = db + sum(20 * np.log10(np.abs(peq(f, p))) for p in eq)
    m = (f >= 500) & (f <= 5000)
    assert np.max(np.abs(after[m] - 80)) < 0.5
    assert any(abs(p.f0 - 1500) < 150 and p.gain_db == pytest.approx(-5.0, abs=0.7) for p in eq)


def test_boost_and_cut_stay_in_range():
    f = np.geomspace(100, 10000, 160)
    db = np.where((f > 1000) & (f < 2000), 50.0, 80.0)  # a 30 dB hole
    for p in fit_peq(f, db, 80.0, (500.0, 5000.0), 4):
        assert GAIN_RANGE[0] <= p.gain_db <= GAIN_RANGE[1] and Q_RANGE[0] <= p.q <= Q_RANGE[1]


def test_short_band_gets_no_eq():
    f = np.geomspace(100, 2851, 60)
    assert fit_peq(f, np.full(len(f), 70.0), 80.0, (2600.0, 16000.0), 5) == ()


def test_eq_bands_never_reach_below_a_high_pass():
    assert eq_band("coax_lf", 170, 1400) == (170, 2800)
    assert eq_band("coax_hf", 170, 1400) == (1400, 16000.0)
    assert eq_band("bass", 170, 1400) == (80.0, 340)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_xo_fit.py -q`
Expected: FAIL (ModuleNotFoundError: `speaker.xo_fit`).

- [ ] **Step 3: Implement `eq_band` and `fit_peq` in `src/speaker/xo_fit.py`**

- [ ] **Step 4: Run them to verify they pass**

Run: `uv run pytest tests/test_xo_fit.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/xo_fit.py tests/test_xo_fit.py
git commit -m "feat: per-way parametric EQ fit towards a flat 23-degree response"
```

### Task 3: Delays, polarity and the crossover-point search

**Files:**
- Modify: `src/speaker/xo_fit.py`, `src/speaker/ranking.py` (rename `_directivity_key` to public `directivity_key`; update its one caller)
- Test: `tests/test_xo_fit.py`

**Interfaces:**
- Consumes: Task 1 (`ActiveCrossover`, `WayFilter`, `way_transfer`, `apply_active`), Task 2 (`fit_peq`, `eq_band`, `N_PEQ`); `MATCH_BANDS`, `REF_WAY` (`speaker.crossover`); `score_all` (`speaker.scores`); `spinorama`; `coax_section` (`speaker.system`); `directivity_key` (`speaker.ranking`).
- Produces:
  - Constants: `BASS_COAX_GRID = (150.0, 170.0, 190.0, 210.0, 230.0, 250.0)`, `LF_HF_GRID = (1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0)`, `BASS_DELAYS = np.arange(-3.0, 3.0001, 0.05) * 1e-3` s, `HF_DELAYS = np.arange(-1.0, 1.0001, 0.01) * 1e-3` s, `ALIGN_SMOOTH = 1 / 6` octave.
  - `fit_at(w: Ways, tl: float, fb: float, fh: float, fallback_db: dict[str, float] | None = None) -> ActiveCrossover` — for fixed points: gains as `apply_crossover` computes them (same `MATCH_BANDS`, same fallback); EQ per way with `fit_peq` on the 23° level (`w.obs.at23_idx`) after the gain, `ref` = the coax LF band level, band from `eq_band`, data clipped to `f ≤ min(f[-1], tl)` for `bass` and `f ≤ f[-1]` for the coax ways; then the bass (against coax LF) and the coax HF (against coax LF) each get the delay from its grid and polarity ±1 that minimise the RMS of `smooth(23° dB of the two filtered ways' sum, ALIGN_SMOOTH) − ref` over `[fc/2, 2fc] ∩ [f[0], tl]`; finally all three delays are shifted so the smallest is 0. `fit_top = f[-1]`.
  - `fit_active(w: Ways, tl: float, fallback_db: dict[str, float] | None = None) -> ActiveCrossover` — `fit_at` on every (fb, fh) of the two grids; each candidate's total → `spinorama` → `score_all(..., trusted_limit=tl, xo=candidate, puck_blend=w.meta.get("puck_blend"), coax_far=w.meta.get("coax_far_h0_h180"))`; the candidate with the smallest `directivity_key` over its trusted `window_spread`, `presence_lift`, `di_smooth_bass`, `di_smooth_lfhf`, `flat23` wins; ties keep grid order (fb outer, fh inner). A candidate missing any of the five trusted scores is skipped; if all are skipped, `ValueError("no crossover candidate has all five trusted directivity scores")`.

- [ ] **Step 1: Write the failing tests**

```python
def omni_ways():
    o = observation_set(2.0)
    f = np.geomspace(50, 16000, 200)
    return Ways(f, o, {k: np.ones((len(f), len(o.dirs)), complex) for k in ("bass", "coax_lf", "coax_hf")}, {})


def test_alignment_recovers_a_known_delay_and_inverted_polarity():
    w = omni_ways()
    w.ways["coax_hf"] = -w.ways["coax_hf"] * np.exp(-2j * np.pi * w.f * 0.3e-3)[:, None]
    xo = fit_at(w, 16000.0, 170.0, 1400.0)
    lf, hf = xo.ways["coax_lf"], xo.ways["coax_hf"]
    assert hf.polarity == -1
    assert (lf.delay_s - hf.delay_s) == pytest.approx(0.3e-3, abs=0.02e-3)
    assert min(x.delay_s for x in xo.ways.values()) == 0.0


def test_flat_omni_ways_need_no_eq():
    xo = fit_at(omni_ways(), 16000.0, 170.0, 1400.0)
    assert all(x.peq == () for x in xo.ways.values())


def test_search_stays_in_the_spec_ranges_and_returns_a_grid_point():
    xo = fit_active(omni_ways(), 2851.0)
    assert xo.bass_coax in BASS_COAX_GRID and xo.lf_hf in LF_HF_GRID


def test_directivity_key_is_public():
    from speaker.ranking import directivity_key
    v = dict(window_spread=1.0, presence_lift=1.0, di_smooth_bass=0.5, di_smooth_lfhf=0.5, flat23=2.0)
    assert directivity_key(v) == (8, 1.0, 2.0)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_xo_fit.py -q`
Expected: FAIL (ImportError: `fit_at`, `fit_active`, `directivity_key`).

- [ ] **Step 3: Implement `fit_at` and `fit_active`; rename the ranking key**

- [ ] **Step 4: Run the tests and the ranking tests**

Run: `uv run pytest tests/test_xo_fit.py tests/test_ranking.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/xo_fit.py src/speaker/ranking.py tests/test_xo_fit.py
git commit -m "feat: fit the active crossover: delays, polarity and the point search"
```

### Task 4: The fitted crossover inside the loop

**Files:**
- Modify: `src/speaker/evaluate.py`, `src/speaker/scores.py` (type hint only), `src/speaker/ranking.py` (`VariantResult`), `src/speaker/round.py` (reuse rule), `src/speaker/directivity_report.py` (`run_directivity`, `_check_curves` call)
- Test: `tests/test_evaluate.py`, `tests/test_round.py`, `tests/test_directivity_cli.py`

**Interfaces:**
- Consumes: `fit_active` (Task 3), `ActiveCrossover`, `apply_active`, `biquads` (Task 1).
- Produces:
  - `EvalSettings.fit_crossover: bool = False` (new last field); `SCREEN` sets it `True`, `VENTS` keeps `False`.
  - `responses(d, cab, puck, f, tl, xo: FixedCrossover | None, lf_recess=None) -> tuple[Ways, np.ndarray, dict[str, float], Spin, FixedCrossover | ActiveCrossover]` — `xo=None` fits (`fit_active`, with `_datasheet_fallback()`); the fifth element is the crossover used. The gains dict is `xo.gains_db` when fitted.
  - `Evaluation.xo` (new field) holds the crossover used.
  - `SCORING_VERSION = 2` in `speaker.evaluate`; `VariantResult.scoring: int = 1` (new last field, so stored results without it load as 1). `_run_variant` reuses a stored result only if `prev.scoring == SCORING_VERSION`, and stores `scoring=SCORING_VERSION`.
  - `run_directivity` (not `--fast`) fits the crossover and writes `crossover.json` (`dataclasses.asdict(xo)` plus `"trusted_limit"`) and `crossover.txt` (below), and still writes `gains.json`. `--fast` keeps `FixedCrossover()`.
  - `crossover_text(xo: ActiveCrossover, tl: float) -> str` in `speaker.xo_fit`: `#` comment lines (points, LR4 acoustic target, "delays are relative, the smallest is 0", "PEQ above <tl> Hz rests on low-confidence data"), then per way the rows `gain`, `delay_ms`, `polarity`, `peq` (f0, gain, Q, and `low` when f0 > tl), then one block per sample rate (48000, 96000) of `biquads` rows: `way  label  b0  b1  b2  a1  a2`, tab-separated, 10 significant digits.
  - `score_all`'s `xo` parameter is typed `FixedCrossover | ActiveCrossover` (it reads only `bass_coax` and `lf_hf`).

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_evaluate.py (its autouse fake_fields fixture supplies the fields)
def test_screen_fits_the_crossover_and_vents_does_not(tmp_path):
    assert SCREEN.fit_crossover and not VENTS.fit_crossover
    ev = evaluate(DEFAULT_PATH, SCREEN, trusted_limit=2851.0, tune=False, cache_dir=tmp_path)
    assert isinstance(ev.xo, ActiveCrossover)
    assert ev.xo.bass_coax in BASS_COAX_GRID and ev.xo.lf_hf in LF_HF_GRID
    ev = evaluate(DEFAULT_PATH, VENTS, trusted_limit=2851.0, tune=False, cache_dir=tmp_path)
    assert isinstance(ev.xo, FixedCrossover)


# tests/test_round.py, with the fake_eval fixture, after test_round_reruns_on_a_new_code_version's pattern
def test_old_scoring_results_are_not_reused(tmp_path, fake_eval):
    # run the round once; strip "scoring" from one variant's result.json; run again:
    # that variant's log line has no "(reused)", and its new result.json has scoring == SCORING_VERSION;
    # an untouched variant is still reused.


# tests/test_directivity_cli.py
def test_full_run_writes_the_crossover_files(tmp_path, monkeypatch):
    # monkeypatch speaker.directivity_report's cabinet_field, puck_field and recess_ratio with
    # fakes.monopole_field / fakes.monopole_puck / a unit recess ratio, as tests/test_evaluate.py's
    # fake_fields does; then run_directivity(DEFAULT_PATH, tmp_path, trusted_limit=(2851.0, "test")):
    xo = json.loads((tmp_path / "crossover.json").read_text())
    assert xo["bass_coax"] in BASS_COAX_GRID and xo["lf_hf"] in LF_HF_GRID
    txt = (tmp_path / "crossover.txt").read_text()
    assert "# fs = 48000" in txt and "# fs = 96000" in txt
    assert sum(l.split("\t")[0] in ("bass", "coax_lf", "coax_hf") and "delay_ms" in l
               for l in txt.splitlines()) == 3
    assert (tmp_path / "gains.json").exists()
```

`test_presets` changes to the new `EvalSettings` field.

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_evaluate.py tests/test_round.py tests/test_directivity_cli.py -q`
Expected: the three new tests FAIL.

- [ ] **Step 3: Implement the Interfaces above**

`_check_curves` (mesh check) keeps `FixedCrossover()`: it compares two meshes, and a crossover fitted to each would move with the mesh.

- [ ] **Step 4: Run the whole suite**

Run: `uv run pytest -q`
Expected: PASS (existing tests that unpack `responses` or construct `EvalSettings` are updated to the new arity).

- [ ] **Step 5: Commit**

```bash
git add src/speaker tests
git commit -m "feat: fit the active crossover for screens and full runs; re-score old round results"
```

### Task 5: Crossover in the reports, with passive-mode warnings

**Files:**
- Modify: `src/speaker/xo_fit.py`, `src/speaker/directivity_report.py`, `src/speaker/directivity_plots.py`, `src/speaker/round_report.py`
- Test: `tests/test_xo_fit.py`, `tests/test_directivity_cli.py`, `tests/test_round_report.py`

**Interfaces:**
- Consumes: `ActiveCrossover`, `way_transfer`, `peq` (Task 1); `crossover.json` (Task 4).
- Produces:
  - `passive_warnings(xo: ActiveCrossover, tl: float) -> list[str]` in `speaker.xo_fit`, one sentence per finding:
    - a way whose summed EQ (Σ peq dB) exceeds `PASSIVE_BOOST_DB = 6.0` anywhere in its `eq_band` (names the way, the boost and the frequency);
    - `lf_hf < HF_RESONANCE_OCT_HZ = 1500.0` (within an octave of the 750 Hz HF resonance, main spec §6.1);
    - `|delay(coax_hf) − delay(coax_lf)| > PASSIVE_DELAY_S = 0.1e-3` (a passive coax network cannot add that delay);
    - any PEQ with `f0 > tl` (rests on low-confidence data).
  - `plots.crossover_23(path, f, curves: dict[str, np.ndarray], tl, tl_text, blend=None)` — each way's 23° level after its filter, and the sum, 1/6-octave smoothed, the trusted limit marked as the other directivity plots do.
  - `run_directivity`'s report: the Crossover section shows the points, "acoustic LR4 target", and per way gain, delay (ms), polarity and its PEQs (f0, gain, Q, "low confidence" above the limit), the passive warnings as a list (or "none"), the `crossover_23.png` plot, and links to `crossover.txt` and `vituixcad/README.txt`. `--fast` keeps today's fixed-crossover text.
  - `round_report._crossover`: per final, read `crossover.json`; one row per final with the points, the delays (ms), the PEQ count per way, the largest boost and the warning count. A final without `crossover.json` falls back to today's `gains.json` table and text.

- [ ] **Step 1: Write the failing tests**

```python
def test_passive_warnings():
    eq = (Peq(3000.0, 8.0, 1.0),)
    xo = ActiveCrossover(170.0, 1200.0, {"bass": WayFilter(0, 0, 1, ()), "coax_lf": WayFilter(0, 0.4e-3, 1, ()),
                                         "coax_hf": WayFilter(-5, 0.0, 1, eq)}, 20000.0)
    w = " ".join(passive_warnings(xo, 2851.0))
    assert "coax_hf" in w and "8.0 dB" in w      # boost over 6 dB
    assert "750 Hz" in w                          # 1.2 kHz is within an octave of the HF resonance
    assert "0.40 ms" in w                         # LF-HF delay a passive network can't add
    assert "3000 Hz" in w                         # PEQ above the trusted limit
    flat = {k: WayFilter(0.0, 0.0, 1, ()) for k in ("bass", "coax_lf", "coax_hf")}
    assert passive_warnings(ActiveCrossover(170.0, 1600.0, flat, 20000.0), 2851.0) == []


# tests/test_round_report.py, with its round_dir fixture (finals carry gains.json only)
def test_round_report_without_crossover_json_falls_back_to_gains(round_dir):
    html = html_of(round_dir)
    # the Crossover section still has the bass / coax_lf / coax_hf gains table, and no traceback or
    # "crossover.json" error


def test_round_report_shows_the_fitted_crossover(round_dir):
    # put(<first final dir>/"crossover.json", dataclasses.asdict(xo) | {"trusted_limit": 2851.0}) for an
    # ActiveCrossover at 190 Hz / 1600 Hz with a 0.25 ms coax LF delay; re-render with html_of:
    # the Crossover section shows "190 Hz", "1600 Hz" and "0.25 ms"
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_xo_fit.py tests/test_round_report.py tests/test_directivity_cli.py -q`
Expected: the new tests FAIL.

- [ ] **Step 3: Implement the Interfaces above**

- [ ] **Step 4: Run the whole suite**

Run: `uv run pytest -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/speaker tests
git commit -m "feat: report the fitted crossover and where a passive network would struggle"
```

### Task 6: Electrical impedances

**Files:**
- Modify: `src/speaker/cardioid.py`, `src/speaker/bass.py`
- Create: `src/speaker/vituixcad.py`
- Test: `tests/test_cardioid.py`, `tests/test_bass.py`, `tests/test_vituixcad.py`

**Interfaces:**
- Consumes: `CardioidParams`, `cardioid_params(p, layout, coax)`, `load_curve`, `interp_db`, `minimum_phase`, `LF_JOIN` and `_log_ramp` (`speaker.system`), `Design` (`speaker.evaluate`).
- Produces:
  - `cardioid.electrical_impedance(d: Driver, p: CardioidParams, f) -> np.ndarray` — `Re + jωLe + Bl² / zm`, with `zm = jωMms + Rms + 1/(jωCms) + Sd²·zbox` and `zbox` as in `simulate`.
  - `bass.sealed_impedance(d: Driver, vb: float, f) -> np.ndarray` — the same with a sealed box `vb`.
  - `vituixcad.ZMA_GRID = (10.0, 20000.0, 48)` (f_lo, f_hi, points per octave) and `vituixcad.impedances(d: Design, f: np.ndarray) -> dict[str, np.ndarray]` with keys `coax_lf`, `coax_hf`, `bass_driver`, `bass_pair`:
    - `coax_lf`: |Z| from `electrical_impedance` (with the design's `cardioid_params`) below `LF_JOIN[0]`, from the datasheet trace `bms_10c262_lf_imp` above `LF_JOIN[1]`, blended in dB across it with `_log_ramp`; phase `minimum_phase(f, joined dB)`.
    - `coax_hf`: the trace `bms_10c262_hf_imp`, phase `minimum_phase`.
    - `bass_driver`: `sealed_impedance(bass, layout.bass_net / 2, f)` (the lumped complex value); `bass_pair` = `bass_driver / 2` (two drivers in parallel).

- [ ] **Step 1: Write the failing tests**

```python
def test_sealed_impedance_peaks_at_fc():
    d = load_driver("dayton_rss265ho_44")
    f = np.geomspace(10, 1000, 2000)
    z = sealed_impedance(d, 0.028, f)
    assert f[np.argmax(np.abs(z))] == pytest.approx(sealed_alignment(d, 0.028).fc, rel=0.02)
    assert d.re < abs(z[0]) < 1.5 * d.re  # well below fc the motional part is small


def test_minimum_phase_reproduces_the_lumped_phase():
    # validates the method used on the datasheet traces: impedance is minimum phase
    d = load_driver("dayton_rss265ho_44")
    f = np.geomspace(10, 20000, 48 * 11)
    z = sealed_impedance(d, 0.028, f)
    ph = minimum_phase(f, 20 * np.log10(np.abs(z)))
    m = (f >= 30) & (f <= 300)
    assert np.max(np.abs(np.degrees(ph[m] - np.angle(z[m])))) < 5.0


def test_closed_chamber_matches_a_sealed_box():
    d = load_driver("bms_10c262")
    f = np.geomspace(20, 2000, 100)
    p = CardioidParams(0.016, 0.0168, 0.018, math.inf, ((0.0, 0.1, 0.0),), 12)
    assert np.allclose(electrical_impedance(d, p, f), sealed_impedance(d, 0.016, f), rtol=1e-9)


def test_impedances_join_the_model_and_the_datasheet():
    d = load_design(DEFAULT_PATH)  # no BEM: load_design builds only the layout
    f = np.geomspace(*ZMA_GRID[:2], 400)
    z = impedances(d, f)
    assert set(z) == {"coax_lf", "coax_hf", "bass_driver", "bass_pair"}
    # coax_lf |Z| equals the lumped model below 300 Hz and the datasheet trace above 500 Hz (within 0.01 dB)
    assert np.allclose(z["bass_pair"], z["bass_driver"] / 2)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_cardioid.py tests/test_bass.py tests/test_vituixcad.py -q`
Expected: the new tests FAIL.

- [ ] **Step 3: Implement the Interfaces above**

- [ ] **Step 4: Run them to verify they pass**

Run: `uv run pytest tests/test_cardioid.py tests/test_bass.py tests/test_vituixcad.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/speaker tests
git commit -m "feat: electrical impedances of every way, with the coax chamber's loading"
```

### Task 7: VituixCAD FRD and ZMA export

The format follows the import test set sent to the owner on 2026-10-29 (`out/vituixcad_test.zip`). **If the owner's import check reported a problem, change `FRD_COMMENT`, the name pattern or the header to what they found before Step 1.**

**Files:**
- Modify: `src/speaker/system.py` (`way_responses` meta only), `src/speaker/vituixcad.py`, `src/speaker/directivity_report.py`
- Test: `tests/test_vituixcad.py`, `tests/test_system.py`

**Interfaces:**
- Consumes: `Ways`, `impedances` (Task 6), `P_REF_PEAK`, `C`.
- Produces:
  - `Ways.meta["bass_each"] = {"bass_lower": (nf, nd), "bass_upper": (nf, nd)}` — each driver at 2.83 V, the terms `way_responses` already sums into `bass`.
  - `vituixcad.FRD_SETS = ("coax_lf", "coax_hf", "bass_pair", "bass_lower", "bass_upper")`, `FRD_COMMENT = "*"`.
  - `vituixcad.export_vituixcad(out_dir: Path, ways: Ways, d: Design, tl: float, tl_text: str) -> Path` — writes `out_dir/vituixcad/` and returns it:
    - per set and plane, one file per ring angle named `f"{set} {plane} {angle}.frd"` with `plane` in `hor`, `ver`; angles from the rings (−175…175 and −180 written as `180`), `ver 0` skipped: 72 `hor` and 71 `ver` files per set;
    - rows `f"{f:.4f}\t{dB:.3f}\t{deg:.2f}"`: SPL `20·log10(|p| / P_REF_PEAK)` and phase in degrees of `p·exp(+j2πf·r/C)` (the common 2 m flight time removed); the bass sets only up to `tl`;
    - five header lines, each starting with `FRD_COMMENT`: what the file is and that it is SIMULATED; the design file and its sha256 (first 8); "SPL dB re 20 uPa at 2 m for 2.83 V RMS into the way, unfiltered, free space"; "phase deg, common 2 m flight time removed; rotation centre: the coax centre on the front face"; "trusted below <tl> Hz (<tl_text>); low confidence above";
    - `coax_lf.zma`, `coax_hf.zma`, `bass_driver.zma`, `bass_pair.zma` on `ZMA_GRID`, rows `f, ohm, deg`, with a header saying how each was made (Task 6);
    - `README.txt`: how to import (each set as one driver, "hor"/"ver" keywords, angle multiplier 1); set every driver's X, Y and Z to 0, since the positions are already in the phase; the data are at 2 m in free space, so do not add VituixCAD's own baffle or diffraction simulation; vertical angles are positive above the axis; the bass sets stop at the trusted limit; the coax LF set includes both vents; the bass pair is the two drivers in parallel at 2.83 V; the confidence band.
  - `run_directivity` (not `--fast`) calls `export_vituixcad(out_dir, ways, d, tl, tl_text)`.

- [ ] **Step 1: Write the failing tests**

Fixture for this file: `d = load_design(DEFAULT_PATH)`; `ways` on `observation_set(2.0)` and `f = np.geomspace(50, 20000, 120)`, every way (and both `bass_each` entries, the bass ones zero above 2851 Hz as in `way_responses`) equal to `P_REF_PEAK * 10**(80/20) * exp(-j 2π f r / C)` in every direction.

```python
def test_frd_names_are_unique_and_parseable(tmp_path, ways, d):
    out = export_vituixcad(tmp_path, ways, d, 2851.0, "test")
    names = sorted(p.name for p in out.glob("coax_lf *.frd"))
    assert len(names) == 72 + 71 == len(set(names))
    assert "coax_lf hor 180.frd" in names and "coax_lf hor -175.frd" in names
    assert "coax_lf hor -180.frd" not in names and "coax_lf ver 0.frd" not in names
    assert all(re.fullmatch(r"coax_lf (hor|ver) -?\d+\.frd", n) for n in names)


def test_no_non_finite_values_and_bass_stops_at_the_limit(tmp_path, ways, d):
    out = export_vituixcad(tmp_path, ways, d, 2851.0, "test")
    for p in out.glob("*.frd"):
        rows = [l.split("\t") for l in p.read_text().splitlines() if not l.startswith("*")]
        vals = np.array(rows, dtype=float)
        assert np.isfinite(vals).all()
        if p.name.startswith("bass"):
            assert vals[-1, 0] <= 2851.0


def test_level_and_phase(tmp_path, ways, d):
    out = export_vituixcad(tmp_path, ways, d, 2851.0, "test")
    rows = [l.split("\t") for l in (out / "coax_hf hor 30.frd").read_text().splitlines() if not l.startswith("*")]
    vals = np.array(rows, dtype=float)
    assert np.allclose(vals[:, 1], 80.0, atol=1e-3) and np.allclose(vals[:, 2], 0.0, atol=0.01)


def test_every_file_is_marked_simulated(tmp_path, ways, d):
    # every .frd and .zma file's first line starts with "*" and contains "SIMULATED";
    # the four .zma files and README.txt exist


# tests/test_system.py
def test_bass_each_sums_to_the_bass_way():  # build w as the existing way_responses tests in this file do
    each = w.meta["bass_each"]
    assert np.allclose(each["bass_lower"] + each["bass_upper"], w.ways["bass"])
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_vituixcad.py tests/test_system.py -q`
Expected: the new tests FAIL.

- [ ] **Step 3: Implement the Interfaces above**

- [ ] **Step 4: Run the whole suite**

Run: `uv run pytest -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/speaker tests
git commit -m "feat: VituixCAD FRD and ZMA export of every full run"
```

### Task 8: Re-score round 3 and export the chosen design

No new code. Run, check, write up. Every BEM field is cached, so these runs should not start a new solve; a log line `cabinet: band <= … Hz, … s` taking more than a few seconds means a cache miss — stop and find out why before going on.

**Files:**
- Create: `docs/rounds/round-3b.md`

- [ ] **Step 1: Full run of the chosen design**

Run: `uv run speaker directivity --params params/round-3.toml --out out/plan3/round-3-design --trusted-limit 2851`
Expected: `wrote out/plan3/round-3-design/report.html`; the directory holds `crossover.json`, `crossover.txt` and `vituixcad/` (5 × 143 FRDs, 4 ZMAs, README.txt). The points lie in the spec ranges. Trusted window spread, presence lift and the cardioid rows are close to the fixed-crossover final of `out/rounds/3/finals/bass260-580/scores.json`; flat23 should fall, since the EQ now flattens 23°.

- [ ] **Step 2: Re-score round 3 into a new directory**

Run: `uv run speaker round rounds/round-3.toml --out out/rounds/3b`
Expected: exit 0; `out/rounds/3b/report.html`. A new directory, so nothing from `out/rounds/3` is reused.

- [ ] **Step 3: Build the CAD of the chosen design**

Run: `uv run speaker build --params params/round-3.toml --out out/plan3/round-3-build`
Expected: exit 0 (the first CAD build of the round 2/3 geometry: 70 mm slots at 88 mm, 320 mm flare). If it fails, stop and report: that is a finding for the owner, not part of this task.

- [ ] **Step 4: Write `docs/rounds/round-3b.md`**

In the style of `docs/rounds/round-1.md`: what was run (commit, commands, times); the fitted crossover of the chosen design (points, delays, polarity, EQ, passive warnings); its scores against the fixed-crossover final; whether the round 3 ranking changed with the fitted crossover, and if so what that means for the owner's choice of bass260-580; the CAD build result; where the VituixCAD set is.

- [ ] **Step 5: Commit**

```bash
git add docs/rounds/round-3b.md
git commit -m "docs: round 3b: round 3 re-scored with the fitted crossover"
```
