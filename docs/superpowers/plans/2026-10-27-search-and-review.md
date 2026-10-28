# Search and Review Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the round-0 simulation into a repeatable round: generate design variants from a round file, tune each variant's cardioid vents, screen them cheaply, run the best at full resolution, rank them, and write the spec §7 review report. Then run round 1 and recommend a design.

**Architecture:**
- **The BEM cache is keyed on the mesh itself**, not on the params file. Changing the vent resistance, the chamber volume or a comment then costs no new solve. A geometry change produces a new mesh and so a new solve.
- **New shape parameters** (spec §5): front-edge family (square, roundover, chamfer, angled upper facets), coax horizontal offset, a thicker outer baffle layer, and the flare profile. They live in `params/*.toml`, with buildability rules in one module (`front_shape.py`), and the BEM mesher builds them. The CAD for them is deferred (see Global Constraints).
- **The coax LF de-embedding gets a recess correction.** Round 0's open question 1 asks whether the +13.5 dB coax LF peak at ~570 Hz is real or counted twice. Today the datasheet is de-embedded with a flat piston, while the 3D model adds the cone's own recess again. A BEM ratio fixes this: the cone in its own frame recess on a disc baffle, divided by a flat piston on the same disc.
- **Each variant is tuned for the cardioid before it is scored.** The vent resistance and the chamber net volume change no BEM field, so a grid search over both runs in seconds per variant.
- **A round runner** (`speaker round rounds/round-N.toml`) evaluates its stages in order. Each stage's best variant is the base of the next. The top few of the last stage are run at full resolution, and an HTML report is written. Everything resumes after an interruption.

**Tech Stack:** Python 3.12 (uv), numpy, scipy, matplotlib, bempp-cl 0.4.x (OpenCL via pocl), gmsh ≥ 4.13 (PyPI, 4.15.2 installed), and the new `tomli-w` for writing variant TOML files.

**Spec:** `docs/superpowers/specs/2026-10-26-directivity-optimisation-design.md` (this plan is its §9 item 2). Read §4 (scores and ranking), §5 (constraints and parameters), §6.2 (search) and §7 (report). Round 0's results and open questions are in `docs/rounds/round-0.md` and `docs/rounds/round-0b.md`.

## Global Constraints

- **Units and frames as Plan 1.**
  - SI inside the code. mm, litres and degrees appear only in TOML, converted by `params.py`.
  - Layout frame: x across the width, y backwards (front face at y=0), z up.
  - Directions: H angle θ is `(sin θ, −cos θ, 0)`, and the listener is at **+23°**, towards +x.
  - Engineering e^{+jωt} convention.
- **Hard limits** (spec §5): external width ≤ 330 mm, depth ≤ 350 mm. The coax centre stays at ear height (`coax_centre_z_mm`). Bass net volume ≥ 25 L per driver, so ≥ 50 L for the pair: a variant below it is rejected (Task 3 rule (k)). Top and bottom edges stay square. No coax section narrower than the bass section.
- **Cardioid gate** (spec §4): coax section alone, H 0° minus H 180°.
  - **Hard:** ≥ 10 dB over 150–300 Hz.
  - **Aim:** ≥ 6 dB over 300–500 Hz.
  - Variants that fail the hard gate are out of the ranking, unless every variant fails (then see Task 8).
- **Ranking** (spec §4, my formula): trusted scores only.
  - **Primary:** `window_spread + presence_lift`, bucketed to 0.25 dB.
  - **Then:** `di_smooth_bass + di_smooth_lfhf`.
  - **Then:** `flat23`.
- **Trusted limit:** 2851 Hz, inherited from round 0's mesh check. Variants are not mesh-checked; the round report says so. Only the recommended design gets its own mesh check (Task 10).
- **Crossover:** fixed ideal LR4 at 170 Hz and 1.4 kHz, as in Plan 1. Plan 3 replaces it (spec §9).
- **Round budget:** one unattended round ≤ 12 h on this machine (spec §3). The BEM memory budget stays 22 GB.
- **CAD is deferred (ruling, to confirm with the owner).** `speaker build` (CAD, DXF, STEP) refuses the new shape parameters at non-default values, with a clear message.
  - The CAD for the shape the owner accepts is drawn after the acceptance, in a short follow-on plan.
  - **Why:** the families' build details (corner strips, facet pieces) are only worth drawing once one is chosen.
  - **Cost if wrong:** the follow-on plan is needed before CNC release anyway.
- **Do not edit** `speaker.bem.solver`, `speaker.bem.run` or `speaker.constants` while a BEM job runs. After Task 1 these are the only source files in the cache key.
- **Tests:** BEM-running tests carry `@pytest.mark.bem`. The default `uv run pytest` skips them.

## Review Focus

1. **A round killed mid-variant and re-run.** Expected: it resumes without redoing finished variants, and a half-written result file is never trusted. Test: `test_round_resumes_and_ignores_partial_results` (Task 8).
2. **An owner suggestion with a misspelt key** (`cabinet.widht_mm`). Expected: `speaker round` stops before any BEM, naming the key and the variant. Test: `test_unknown_override_key_names_it` (Task 7).
3. **Every variant fails the hard cardioid gate.** Expected: the round still ranks (by cardioid_hard first), and the report says the gate is unmet. It must not produce an empty ranking. Test: `test_rank_when_every_variant_fails_the_gate` (Task 8).
4. **Two variants that differ only in non-geometry values** (resistance, chamber volume, comments). Expected: the second reuses the cached BEM fields. Test: `test_cache_ignores_non_geometry_params` (Task 1).
5. **A front treatment that cuts into the chamber, or squeezes the flare mouth or a bass cutout.** Expected: the variant is rejected with the limit it broke, and the round continues. Tests: Task 3's validation tests, plus `test_round_records_rejected_variants` (Task 8).

---

## File Structure

| File | Responsibility |
|---|---|
| `src/speaker/bem/run.py` | (modify) cache key from the meshes; `code_version()` over the solver modules only; `recess_ratio` |
| `src/speaker/bem/mesh.py` | (modify) front-edge families, facets, coax offset in `cabinet_mesh`; new `disc_baffle_mesh` |
| `src/speaker/coax_geometry.py` | (modify) outer baffle thickness, flare profile, `recess_profile` |
| `src/speaker/params.py` | (modify) new `[cabinet]` keys |
| `src/speaker/front_shape.py` | (create) front cross-section, flat widths, buildability rules |
| `src/speaker/layout.py` | (modify) coax offset, calls `validate_front` |
| `src/speaker/system.py` | (modify) LF recess correction, HF reference continued above the puck |
| `src/speaker/vent_tune.py` | (create) resistance × chamber-volume grid search per variant |
| `src/speaker/evaluate.py` | (create) `Design`, `load_design`, `banded`, `grid`, `EvalSettings`, `evaluate` |
| `src/speaker/directivity_report.py` | (modify) uses `evaluate.py` helpers; trusted-limit override; LF-reference note and plot |
| `src/speaker/variants.py` | (create) round files, overrides, variant TOML writing |
| `src/speaker/ranking.py` | (create) stage ranking rules |
| `src/speaker/round.py` | (create) the round runner, persistence, resume, finals |
| `src/speaker/round_report.py` | (create) the spec §7 HTML report |
| `src/speaker/shape_drawing.py` | (create) front and plan drawings of designs, for the report |
| `src/speaker/report.py` | (modify) `speaker build` refuses BEM-only geometry |
| `src/speaker/cli.py` | (modify) `speaker round`; `speaker directivity --trusted-limit` |
| `params/default.toml` | (modify) new keys at today's values |
| `rounds/round-1.toml` | (create, Task 10) round 1 |
| `tests/fakes.py` | (create, Task 5) free-field monopole `SourceField` for fast tests |

---

### Task 1: BEM cache keyed on the mesh

**Files:**
- Modify: `src/speaker/bem/run.py` (`cabinet_field`, `puck_field`, `_banded_field`, `_cache_key`, `code_version`; drop `_imports`)
- Modify: `src/speaker/directivity_report.py` and every test that passes `params_text=`
- Test: `tests/test_bem_run.py`

**Interfaces:**
- Produces:
  - `cabinet_field(p, coax, bass, f, *, per_wavelength=6.0, h_max=0.05, h_coax=CAB_H_COAX, h_vent=CAB_H_VENT, cache_dir=...) -> SourceField`. The `params_text` argument is removed.
  - `puck_field(p, coax, f, *, per_wavelength=6.0, cache_dir=...) -> SourceField`, likewise.
  - `code_version() -> str`: sha256 over the source of `SOLVER_MODULES = ("speaker.bem.run", "speaker.bem.solver", "speaker.constants")` and the installed bempp-cl version.

`_banded_field` builds every band's mesh first. Its key is the sha256 of:
- the kind, `per_wavelength`, `h_max` and `code_version()`;
- `f` and `points` as float64 bytes;
- for each band in order: the edge, `h`, `vertices` (float64), `triangles` (int64), `tags` (int64) and `tag_names`.

On a miss, it solves those same meshes. It does not mesh twice. Geometry enters the key only through the meshes and the observation points. gmsh's version no longer needs to be in the key: the mesh bytes cover it.

- [ ] **Step 1: Write the failing tests.** In `tests/test_bem_run.py`, add a fixture that monkeypatches `speaker.bem.run.solve_band` with a counter returning `np.zeros((len(f), len(mesh.tag_names), points.shape[1]), complex)`. Call `cabinet_field(..., f=np.array([100.0]), h_max=0.2, cache_dir=tmp_path)` with:

```python
def test_cache_ignores_non_geometry_params(counting_solver, tmp_path):
    p, coax, bass = default_design()
    cabinet_field(p, coax, bass, F, h_max=0.2, cache_dir=tmp_path)
    q = replace(p, vents=replace(p.vents, specific_resistance=400.0),
                drive=replace(p.drive, amp_v_rms=10.0),
                cabinet=replace(p.cabinet, coax_chamber_net_volume=0.012))
    cabinet_field(q, coax, bass, F, h_max=0.2, cache_dir=tmp_path)
    assert counting_solver.calls == 1

def test_cache_resolves_when_geometry_moves(counting_solver, tmp_path):
    p, coax, bass = default_design()
    cabinet_field(p, coax, bass, F, h_max=0.2, cache_dir=tmp_path)
    q = replace(p, vents=replace(p.vents, centre_from_front=0.240))
    cabinet_field(q, coax, bass, F, h_max=0.2, cache_dir=tmp_path)
    assert counting_solver.calls == 2

def test_code_version_covers_only_the_solver(monkeypatch, tmp_path):
    before = code_version()
    fake = tmp_path / "layout.py"; fake.write_text("# changed\n")
    monkeypatch.setattr(sys.modules["speaker.layout"], "__file__", str(fake))
    assert code_version() == before
    fake2 = tmp_path / "solver.py"; fake2.write_text("# changed\n")
    monkeypatch.setattr(sys.modules["speaker.bem.solver"], "__file__", str(fake2))
    assert code_version() != before
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_bem_run.py -k "cache_ignores or cache_resolves or covers_only" -v`
Expected: FAIL. `cabinet_field` still requires `params_text`, and `code_version` still follows the import closure.

- [ ] **Step 3: Implement** the key, the mesh-first banding and `code_version` as described under Interfaces. Remove `params_text` from both field functions and from every caller (`directivity_report.py`'s `kw`, `mesh_check`, tests). Update the module docstring's cache paragraph.

- [ ] **Step 4: Run the default suite and the BEM suite**

Run: `uv run pytest -q` then `uv run pytest -m bem -q`
Expected: all pass. The BEM suite takes about 8 min. The cache directory of a real run is not touched by tests.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/bem/run.py src/speaker/directivity_report.py tests
git commit -m "feat: key the BEM cache on the meshes, so only geometry changes cost a solve"
```

---

### Task 2: Coax LF de-embedding with the recess ratio

Round 0's open question 1 (`docs/rounds/round-0.md`) leads here. The datasheet LF curve was measured with the driver front-mounted, so it already contains the cone's own recess. The 3D model adds that recess again. The new reference removes it once:

`u_lf = datasheet ÷ (baffled piston at 1 m × R(f))`

- R(f) is the BEM on-axis pressure at 1 m of the coax cone in its own frame recess, front-mounted on a disc baffle, divided by a flat piston of the same diameter on the same disc.
- Both are per unit volume velocity. The ratio cancels most of the disc's own edge diffraction.
- **Assumption (state it in the report):** the datasheet was measured front-mounted on a large baffle. The frame opening is `cutout_d`, and the surround sits `flange_t` behind the mounting face.

This task also fixes the HF reference above the top of a puck field. Above that top it continues as the baffled piston curve, shifted to meet the puck's last smoothed value. Today it holds the puck's absolute end value, which would be wrong for a puck that stops at 4 kHz (Task 6). The LF ratio follows the same rule: above its top it holds its last value, so the reference continues as the piston curve shifted.

**Files:**
- Modify: `src/speaker/coax_geometry.py` (shared cone segments; `recess_profile`)
- Modify: `src/speaker/bem/mesh.py` (`disc_baffle_mesh`; `_cavity` takes a profile)
- Modify: `src/speaker/bem/run.py` (`recess_ratio`)
- Modify: `src/speaker/system.py` (`way_responses(..., lf_recess=None)`; HF reference above the puck)
- Modify: `src/speaker/directivity_report.py` (passes `lf_recess`; report note; `plots/lf_reference.png`)
- Create: `tests/fakes.py` (`monopole_field`, see Interfaces; Tasks 5, 6 and 8 reuse it)
- Test: `tests/test_coax_geometry.py`, `tests/test_bem_mesh.py`, `tests/test_system.py`, `tests/test_bem_run.py`

**Interfaces:**
- Produces:
  - `recess_profile(coax: Driver) -> tuple[tuple[str, tuple[tuple[float, float], ...]], ...]`: segments `bore` (r = cutout_d/2, y from 0 to flange_t), `surround`, `cone`, `pole` and `hf_exit`, with the cone starting at y = flange_t.
  - `disc_baffle_mesh(coax, *, h, h_coax, recess: bool) -> SurfaceMesh`:
    - A cylinder of radius `DISC_R = 0.5` m along +y, from y=0 to `flange_t + cone_depth + PUCK_EXTRA_DEPTH`, centred on the y axis.
    - `recess=True`: the revolved recess cut in, with the cone tagged `coax_lf` and the HF exit rigid.
    - `recess=False`: a flat disc of `diaphragm_d` imprinted on the front face, tagged `coax_lf`.
  - `recess_ratio(coax, f=None, *, per_wavelength=6.0, cache_dir=...) -> tuple[np.ndarray, np.ndarray]`:
    - Returns `(f, R)`, with R complex, on axis at `AXIS_R` = 1 m.
    - The default f is `freq_grid(300.0, 4000.0, 6)`.
    - It uses `_banded_field` with kinds `"recess_disc"` and `"flat_disc"`, so it is cached.
  - `way_responses(..., lf_recess: tuple[np.ndarray, np.ndarray] | None = None)`:
    - When given, `u_db["lf"] -= ` the 1/3-octave-smoothed `20·log10|R|`, interpolated on the internal grid and held at its end values outside its range.
    - `meta["lf_reference"]` is `"recess"` or `"piston"`.
  - `tests/fakes.py`: `monopole_field(p, coax, bass, f) -> SourceField`.
    - Free-field point sources at each source's layout position (coax centre, the two vent-stack centres, the bass centres), named as `CABINET_TAGS`, on `observation_set(OBS_R)` about the coax centre.
    - Pressure per unit volume velocity is `jωρ e^{−jkr} / (4πr)`.
    - It carries `meta["bands"]` (one band entry) and `meta["far_h0_h180"]` like `cabinet_field`.
    - `monopole_puck(p, coax, f, *, hf_axis_db=None)` does the same for `coax_lf` and `coax_hf`, with `meta["axis_1m"]`. When `hf_axis_db(f)` is given, the coax_hf axis column has that magnitude.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_system.py
def test_lf_recess_lowers_the_datasheet_strength(design):
    f = np.geomspace(120.0, 2000.0, 40)
    cab = monopole_field(design.p, design.coax, design.bass, f)
    flat = way_responses(design.p, design.coax, design.bass, design.layout, cab, None, f,
                         trusted_limit=2000.0, lf_recess=None)
    plus6 = way_responses(design.p, design.coax, design.bass, design.layout, cab, None, f, trusted_limit=2000.0,
                          lf_recess=(np.array([300.0, 4000.0]), np.full(2, 10 ** (6 / 20) + 0j)))
    i1k, i150 = np.argmin(abs(f - 1000)), np.argmin(abs(f - 150))
    assert abs(plus6.meta["u_lf"][i1k]) == pytest.approx(abs(flat.meta["u_lf"][i1k]) / 10 ** (6 / 20), rel=1e-3)
    assert plus6.meta["u_lf"][i150] == pytest.approx(flat.meta["u_lf"][i150])  # lumped below the join
    assert plus6.meta["lf_reference"] == "recess"

def test_hf_reference_continues_as_piston_above_the_puck(design):
    # monopole_puck with hf_axis_db = the baffled-piston curve of hf_exit_d + 3 dB, from 500 Hz to 4 kHz.
    # Run way_responses on f up to 12 kHz with that puck and a cabinet field to its blend top.
    # At 8 kHz, |u_hf| must equal the piston-only de-embedding (hf_reference "piston") minus 3 dB,
    # within 0.1 dB, not the level implied by holding the puck's 4 kHz value.
```

```python
# tests/test_bem_mesh.py
@pytest.mark.parametrize("recess", [True, False])
def test_disc_baffle_mesh_closed_and_tagged(recess):
    m = disc_baffle_mesh(load_driver("bms_10c262"), h=0.08, h_coax=0.02, recess=recess)
    assert m.tag_names == ("coax_lf",)
    assert np.count_nonzero(m.tags == 1) > 0  # edge_check already ran inside the mesher
```

```python
# tests/test_bem_run.py
@pytest.mark.bem
def test_recess_ratio_near_unity_at_300_hz(tmp_path):
    f, r = recess_ratio(load_driver("bms_10c262"), np.array([300.0]), cache_dir=tmp_path)
    assert abs(20 * np.log10(abs(r[0]))) < 1.0  # a ~75 mm recess is small against λ = 1.14 m
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_system.py tests/test_bem_mesh.py -k "recess or hf_reference or disc_baffle" -v`
Expected: FAIL (`lf_recess` and `disc_baffle_mesh` don't exist).

- [ ] **Step 3: Implement.**
  - Refactor `cavity_profile` so it and `recess_profile` share the surround, cone, pole and hf_exit segments, built from a surround depth.
  - `_cavity` takes the profile as an argument.
  - In `run_directivity`, unless `fast`: call `recess_ratio(d.coax, cache_dir=CACHE_DIR)` and pass it to `way_responses`.
  - Add a report note stating the method and the assumption above.
  - Add `plots/lf_reference.png`: 20·log10|R| from 300 to 4000 Hz, with the smoothed curve that is applied.
  - In fast mode the note says "flat piston (fast run)".

- [ ] **Step 4: Run the tests**

Run: `uv run pytest -q` then `uv run pytest -m bem -q -k recess`
Expected: all pass.

- [ ] **Step 5: Measure R once, and write down what it shows.** Run `uv run python -c "from speaker.bem.run import recess_ratio; from speaker.drivers import load_driver; import numpy as np; f,r=recess_ratio(load_driver('bms_10c262')); [print(f'{a:7.0f} {20*np.log10(abs(b)):6.2f}') for a,b in zip(f,r)]"` (about 40 min, cached afterwards).
  - Record in the ledger the dB of R at 300, 570, 1000, 1850 and 4000 Hz.
  - Record whether R accounts for the 570 Hz excess of round 0's open question 1 (+13.5 dB over a point source against a baffle's ~+6 dB). Task 10's round note uses this.

- [ ] **Step 6: Commit**

```bash
git add src/speaker tests
git commit -m "feat: de-embed the coax LF datasheet with the cone's own recess, from a disc-baffle BEM ratio"
```

---

### Task 3: Front-shape parameters and buildability rules

**Files:**
- Modify: `params/default.toml`, `src/speaker/params.py`, `src/speaker/layout.py`, `src/speaker/report.py`
- Create: `src/speaker/front_shape.py`
- Test: `tests/test_front_shape.py` (create), `tests/test_params.py`, `tests/test_layout.py`

The new `[cabinet]` keys in `params/default.toml` replace `baffle_t_mm` and `side_edge_radius_mm`. Their values keep today's design:

```toml
outer_baffle_t_mm = 18     # the baffle is this outer layer plus one panel_t inner layer
front_edge = "square"      # square | roundover | chamfer | facets
edge_size_mm = 0           # roundover radius, or the 45° chamfer leg; 0 for square and facets
facet_angle_deg = 0        # facets: angle back from the front plane
facet_width_mm = 0         # facets: width of each facet, from the side edge inwards
facet_base_z_mm = 0        # facets: height where the facet starts (the bass section below stays flat)
facet_ramp_mm = 0          # facets: 0 = a horizontal step at the base; > 0 = a planar ramp this tall
coax_offset_x_mm = 0       # + towards +x, the listener's side; mirrored on the other cabinet
flare_profile = "ellipse"  # ellipse | cone
```

**Interfaces:**
- Produces:
  - `Cabinet` fields: `outer_baffle_t`, `front_edge: str`, `edge_size`, `facet_angle` (radians), `facet_width`, `facet_base_z`, `facet_ramp`, `coax_offset_x`, `flare_profile: str`. `baffle_t` stays a field, set to `outer_baffle_t + panel_t`. `side_edge_radius` is removed.
  - `Layout.coax_centre == (width / 2 + coax_offset_x, coax_centre_z)`.
  - `front_shape.FRONT_LAND = 0.005` m, `CHAMFER_MAX = 0.045` m, `FAMILIES = ("square", "roundover", "chamfer", "facets")`.
  - `front_shape.cut_depth(c: Cabinet) -> float`: how far back the treatment reaches on the side face. This is `edge_size` for roundover and chamfer, and `facet_width · tan(facet_angle)` for facets.
  - `front_shape.flat_width(c: Cabinet, z: float) -> float`: the flat front width at height z. Facets apply in full at z ≥ base + ramp, linearly across the ramp, and not below the base.
  - `front_shape.section(c: Cabinet, z: float, n_arc: int = 8) -> list[tuple[float, float]]`: the outer horizontal cross-section at z, as a counter-clockwise (x, y) polygon.
  - `front_shape.validate_front(p: Params, coax: Driver, bass: Driver, l: Layout) -> None`: raises ValueError, called from `layout.validate`.

Rules in `validate_front`. Each rule is one ValueError whose message names the values and the limit:
- **(a)** `front_edge` must be in FAMILIES.
- **(b)** For square and facets, `edge_size` must be 0. For roundover and chamfer it must be > 0. A chamfer's `edge_size` must be ≤ CHAMFER_MAX.
- **(c)** For facets: 0 < angle ≤ 45°, width > 0 and ramp ≥ 0. For other families the four facet values must all be 0.
- **(d)** The treatment must not reach the interior. The interior front corners (panel_t, baffle_t) and (width − panel_t, baffle_t) must lie inside `section` at every treated height.
- **(e)** At the coax centre height, `flat_width ≥ flare_mouth_d + 2·|coax_offset_x| + 2·FRONT_LAND`.
- **(f)** At each bass centre height, `flat_width ≥ bass.cutout_d + 2·FRONT_LAND`.
- **(g)** `|coax_offset_x| + coax.outer_d / 2 ≤ interior_w / 2`, so the coax frame clears the side panels.
- **(h)** For facets, two height limits:
  - `facet_base_z ≥` the top of the upper bass cutout + FRONT_LAND;
  - `facet_base_z + facet_ramp ≤ coax_centre_z − flare_mouth_d / 2`.
- **(i)** For facets, the vent slots must start behind the facet cut: `centre_from_front − slot_w / 2 ≥ cut_depth + FRONT_LAND`.
  The facets may start below the divider, but only above the upper bass driver: the part of the front that carries the bass drivers stays flat (spec §5).
- **(j)** `flare_profile` must be in ("ellipse", "cone"). `panel_t ≤ outer_baffle_t ≤ 3·panel_t`, a laminate of standard sheets.
- **(k)** `bass_net ≥ BASS_NET_TARGET` (50 L: spec §5's 25 L per driver). This replaces `compute_layout`'s warning for it; `BASS_NET_MIN` stays as `validate`'s floor for direct callers.

`cardioid_params` places the vents laterally at `width − cx` and `−cx` from the coax centre (today `±width/2`).

`report.build` (`speaker build`) raises ValueError if any of these differ from the default: front_edge != "square", coax_offset_x != 0, outer_baffle_t != panel_t, flare_profile != "ellipse". The message names them and says: "the CAD does not model these yet; `speaker directivity` and `speaker round` accept them".

- [ ] **Step 1: Write the failing tests.**
  - `tests/test_front_shape.py`: one test per rule (a)–(k), each loading `params/default.toml` text with the keys replaced and asserting `pytest.raises(ValueError, match=...)` with a phrase from that rule's message.
  - Plus section-area tests:

```python
def test_section_areas():
    sq = cab(front_edge="square")
    assert polygon_area(section(sq, 0.5)) == pytest.approx(W * D)
    ch = cab(front_edge="chamfer", edge_size=0.03)
    assert polygon_area(section(ch, 0.5)) == pytest.approx(W * D - 0.03**2)          # two s²/2 corners
    ro = cab(front_edge="roundover", edge_size=0.03)
    assert polygon_area(section(ro, 0.5), ) == pytest.approx(W * D - 2 * (1 - math.pi / 4) * 0.03**2, rel=1e-3)
    fa = cab(front_edge="facets", facet_angle=math.radians(20), facet_width=0.04, facet_base_z=0.7, facet_ramp=0.1)
    d = 0.04 * math.tan(math.radians(20))
    assert polygon_area(section(fa, 0.69)) == pytest.approx(W * D)
    assert polygon_area(section(fa, 0.9)) == pytest.approx(W * D - 0.04 * d)        # two fw·d/2 triangles
```

  - `tests/test_layout.py`:
    - `test_default_layout_unchanged`: compare against values captured before the change: height 1.118 m, bass_net 54.09 L, coax_chamber_empty 18.83 L.
    - `test_coax_offset_moves_the_centre`.
  - `tests/test_params.py`: replace `test_side_edge_radius_not_implemented` with `test_facet_angle_converted_to_radians`.
  - `tests/test_report.py` (or `tests/test_cli.py`): `test_build_refuses_bem_only_geometry`.

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_front_shape.py tests/test_layout.py tests/test_params.py -v`
Expected: FAIL (module and keys missing).

- [ ] **Step 3: Implement** `front_shape.py`, the params keys, the layout offset and the call to `validate_front`, the `cardioid_params` lateral positions, and the `speaker build` guard.
  - The old rule "baffle_t must equal two panel layers" is replaced by rule (j).
  - The flare insert's CAD keeps using `panel_t` as its depth, which stays correct because `speaker build` refuses `outer_baffle_t != panel_t`.

- [ ] **Step 4: Run the suite**

Run: `uv run pytest -q`
Expected: all pass, with `speaker build` on the default params unchanged.

- [ ] **Step 5: Commit**

```bash
git add params/default.toml src/speaker tests
git commit -m "feat: front-edge families, coax offset, outer baffle thickness and flare profile as params, with buildability rules"
```

---

### Task 4: BEM meshes of the new shapes

**Files:**
- Modify: `src/speaker/coax_geometry.py` (flare depth = `outer_baffle_t`; `flare_profile`)
- Modify: `src/speaker/bem/mesh.py` (`cabinet_mesh` body)
- Test: `tests/test_bem_mesh.py`, `tests/test_coax_geometry.py`

**Interfaces:**
- Consumes: the Task 3 `Cabinet` fields, `Layout.coax_centre` and `front_shape.cut_depth`.
- Produces: `cabinet_mesh` and `puck_mesh` build every family. The signatures are unchanged.

Geometry:
- **Flare:** the flare runs over the outer layer, `t = outer_baffle_t`, and the bore over the inner layer, to `baffle_t`.
  - `"ellipse"` keeps today's quarter ellipse.
  - `"cone"` is the straight segment from (r_mouth, 0) to (r_throat, t).
- **Roundover:** `occ.fillet` on the two front vertical edges (x=0 and x=width at y=0) with radius `edge_size`.
- **Chamfer:** `occ.chamfer` on the same edges, with distance `edge_size` (45°).
- **Facets:** cut one wedge per side, then mirror it for the side at x=width. With fw = facet_width, d = fw·tan(angle), zb = facet_base_z and zr = zb + facet_ramp, the left side is:
  - **Prism:** the triangle (0,0), (fw,0), (0,d), from z = zr to above the top.
  - **Tetrahedron (only when ramp > 0):** A1=(0,0,zr), B1=(fw,0,zr), C1=(0,d,zr) and A0=(0,0,zb). Its face B1-C1-A0 is the planar ramp.
- The cavity is cut at `Layout.coax_centre`, which includes the offset. Vents and bass discs are imprinted as today.

- [ ] **Step 1: Write the failing tests** (`h=0.04`, `h_coax=0.02` for speed). `mesh_volume(m)` is the signed volume from the mesh module.

```python
CUTS = {  # analytic volume removed by both sides together
    "roundover": lambda c, H: 2 * (1 - math.pi / 4) * c.edge_size**2 * H,
    "chamfer":   lambda c, H: c.edge_size**2 * H,
    "facets_step": lambda c, H: c.facet_width * fd(c) * (H - c.facet_base_z),
    "facets_ramp": lambda c, H: c.facet_width * fd(c) * (H - c.facet_base_z - c.facet_ramp)
                                + c.facet_width * fd(c) * c.facet_ramp / 3,
}

@pytest.mark.parametrize("family", list(CUTS))
def test_front_edge_cut_volumes(family):
    p_x = shaped_params(family)                     # width 330, flare mouth 240, so every rule passes
    p_sq = square_like(p_x)                         # the same params with front_edge "square", no facets
    v_sq = mesh_volume(cabinet_mesh(p_sq, coax, bass, h=0.04, h_coax=0.02))
    v_x = mesh_volume(cabinet_mesh(p_x, coax, bass, h=0.04, h_coax=0.02))
    H = compute_layout(p_x, coax, bass).height
    assert v_sq - v_x == pytest.approx(CUTS[family](p_x.cabinet, H), rel=0.03)

def test_coax_offset_moves_the_cone():  # coax_lf triangle centroids' mean x == width/2 + offset ± 1 mm
def test_thicker_outer_baffle_deepens_the_cavity():  # coax_lf centroid y grows by the extra thickness ± 1 mm
def test_flare_profile_changes_cavity_volume_as_revolved():
    # mesh-volume(ellipse) - mesh-volume(cone) == revolved-volume difference of the two flare
    # polylines (Pappus over each closed (r, y) region), rel 3 %
def test_puck_mesh_builds_each_flare():  # puck_mesh closes for ellipse/cone and outer_baffle_t 18/30 mm

def test_default_cabinet_mesh_unchanged():
    # triangle count and volume of the default params' cabinet_mesh(h=0.04, h_coax=0.02),
    # captured before editing mesh.py and pasted here as constants
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_bem_mesh.py tests/test_coax_geometry.py -v`
Expected: FAIL (families not built; cone profile missing).

- [ ] **Step 3: Implement** as described under Geometry.

- [ ] **Step 4: Run the tests and the BEM suite**

Run: `uv run pytest -q` then `uv run pytest -m bem -q`
Expected: all pass, including `test_default_cabinet_mesh_unchanged`.

- [ ] **Step 5: Commit**

```bash
git add src/speaker tests
git commit -m "feat: mesh roundovers, chamfers, angled facets, coax offset, thicker outer baffle and cone flares"
```

---

### Task 5: Vent tuning per variant

Neither the vent resistance nor the chamber net volume changes a BEM field: the fillers are inside the box, and Task 1 keys the cache on the exterior mesh. So each variant can pick both from a grid. The tuning goal is robust: maximise the worst cardioid margin over the felt spread.

`margin(R, V) = min over s ∈ {0.7, 1.0, 1.3} of min(hard_s − 10, aim_s − 6)`

**Files:**
- Create: `src/speaker/vent_tune.py`
- Test: `tests/test_vent_tune.py` (uses Task 2's `tests/fakes.py`)

**Interfaces:**
- Consumes: `way_responses`, `vent_scaled_section`, `coax_section`, `compute_layout`, `cardioid_params`.
- Produces:
  - `R_GRID = np.geomspace(80.0, 1000.0, 15)` (Pa·s/m), `V_STEP = 0.5e-3` (m³), `V_MIN = 4e-3` (m³), `CARD_BAND = (150.0, 500.0)`.
  - `@dataclass VentTune`:
    - `resistance: float`, `net_volume: float`, `hard: float`, `aim: float`, `margin: float`;
    - `spread: dict[str, tuple[float, float]]` (keys "0.7" and "1.3", values (hard, aim));
    - `base: dict[str, float]` (hard, aim and margin at the params' own R and V).
  - `feasible_volumes(p, coax, bass) -> np.ndarray`: V_MIN up to the empty chamber in V_STEP steps, plus the params' own value, keeping each value where `compute_layout` succeeds.
  - `tune_vents(p, coax, bass, cab: SourceField, *, lf_recess=None) -> VentTune`:
    - Evaluates R_GRID ∪ {p's R} × `feasible_volumes`.
    - Uses `way_responses` on `cab.f` restricted to CARD_BAND (`trusted_limit` = the band top, no puck), and `vent_scaled_section` for 0.7 and 1.3.
    - hard and aim are the min over their bands of SPL(H 0°) − SPL(H 180°) of the section, the same formula as `scores.py`'s cardioid rows.
    - It picks the largest margin.
    - It raises ValueError if `cab.f` doesn't cover CARD_BAND.

- [ ] **Step 1: Write the failing tests**

```python
def test_tuned_margin_not_below_the_base(design):
    cab = monopole_field(*design, F_120_600)
    t = tune_vents(*design, cab)
    assert t.margin >= t.base["margin"] - 1e-9

def test_tuned_hard_matches_the_score(design):
    cab = monopole_field(*design, F_120_600)
    t = tune_vents(*design, cab)
    q = with_vents(design.p, t.resistance, t.net_volume)
    ways = way_responses(q, design.coax, design.bass, compute_layout(q, design.coax, design.bass), cab, None,
                         F_120_600, trusted_limit=600.0)
    sec = coax_section(ways)  # also passed as "total": only the cardioid rows are compared
    sc = score_all(F_120_600, sec, sec, cab.obs, spinorama(sec, cab.obs, f=F_120_600),
                   trusted_limit=600.0, xo=FixedCrossover())
    assert next(s.value for s in sc if s.key == "cardioid_hard") == pytest.approx(t.hard, abs=1e-6)

def test_feasible_volumes_all_lay_out(design):
    for v in feasible_volumes(*design):
        compute_layout(with_vents(design.p, design.p.vents.specific_resistance, v), design.coax, design.bass)

def test_tune_needs_the_cardioid_band(design):
    with pytest.raises(ValueError, match="150"):
        tune_vents(*design, monopole_field(*design, np.geomspace(200, 600, 10)))
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_vent_tune.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement** `vent_tune.py`. `with_vents(p, R, V)` goes in `vent_tune.py` as a public helper using `dataclasses.replace`.

- [ ] **Step 4: Run the tests**

Run: `uv run pytest tests/test_vent_tune.py -q && uv run pytest -q`
Expected: all pass. Report how long one `tune_vents` takes on the fake field: it must stay under 60 s.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/vent_tune.py tests/test_vent_tune.py
git commit -m "feat: tune vent resistance and chamber volume per variant for the worst-case cardioid margin"
```

---

### Task 6: One evaluation API for screens and full runs

**Files:**
- Create: `src/speaker/evaluate.py`
- Modify: `src/speaker/directivity_report.py` (imports the moved helpers; `trusted_limit` override), `src/speaker/cli.py` (`--trusted-limit`)
- Test: `tests/test_evaluate.py`, `tests/test_directivity_cli.py`

**Interfaces:**
- Consumes: Task 1's field functions, Task 2's `recess_ratio`, Task 5's `tune_vents` and `with_vents`.
- Produces:
  - Moved from `directivity_report.py` and made public: `Design` (was `_Design`), `load_design(path) -> Design`, `banded(solve, f, label) -> SourceField`, `grid(f_lo, f_hi, ppo) -> np.ndarray` (with the SCORE_EDGES inside its range), `CACHE_DIR`. `directivity_report.py` imports them from here.
  - `@dataclass(frozen=True) EvalSettings`:
    - `name: str`, `f_lo: float`, `f_hi: float | None` (None means the trusted limit), `ppo: int`;
    - `per_wavelength: float`, `h_max: float`;
    - `puck_top: float | None` (Hz; a puck from 500 Hz to here, used only for the HF reference);
    - `recess: bool`.
  - Presets, as `SETTINGS = {"vents": ..., "screen": ...}`:
    - `VENTS = EvalSettings("vents", 120.0, 600.0, 12, 6.0, 0.05, None, True)`;
    - `SCREEN = EvalSettings("screen", 120.0, None, 12, 5.0, 0.05, 4000.0, True)`.
  - `@dataclass Evaluation`: `design: Design` (tuned when `tuned`), `settings`, `f`, `ways`, `total`, `spin`, `scores: list[Score]`, `gains: dict`, `tune: VentTune`, `tuned: bool`, `trusted_limit: float`, `bands: list[dict]`, `wall_s: float`.
  - `evaluate(path: Path, settings: EvalSettings, *, trusted_limit: float, tune: bool, cache_dir: Path = CACHE_DIR) -> Evaluation`. The steps:
    1. Fields: cabinet over the grid; puck if `puck_top`; recess ratio if `recess`.
    2. Always call `tune_vents` (it costs seconds). When `tune`, rebuild the design with the tuned R and V. When not, the result is kept for its `base` margin, which ranks an untuned variant on the same worst-case terms.
    3. `way_responses` with `lf_recess` and `trusted_limit`.
    4. `apply_crossover(..., FixedCrossover())`, `spinorama`, `score_all` (no blend).
  - `run_directivity(..., trusted_limit: tuple[float, str] | None = None)`. When given, it replaces the mesh check's limit and its header text, and is still capped by `_trusted_limit`'s blend rule. The CLI gets `speaker directivity --trusted-limit HZ`, with the note "given on the command line".

- [ ] **Step 1: Write the failing tests** (`tests/test_evaluate.py`, with `cabinet_field`, `puck_field` and `recess_ratio` monkeypatched to `tests/fakes.py` fields)

```python
def test_screen_scores_are_all_trusted(monkeypatch, tmp_path):
    ev = evaluate(DEFAULT_PATH, SCREEN, trusted_limit=2851.0, tune=False, cache_dir=tmp_path)
    assert 2851.0 / 2 ** (1 / 12) < ev.f[-1] <= 2851.0
    assert not ev.tuned and ev.design.p == load_params(DEFAULT_PATH)
    assert {s.confidence for s in ev.scores} == {"trusted"}
    assert any(s.key == "presence_lift" for s in ev.scores)

def test_tuned_design_carries_the_tuned_vents(monkeypatch, tmp_path):
    ev = evaluate(DEFAULT_PATH, VENTS, trusted_limit=2851.0, tune=True, cache_dir=tmp_path)
    assert ev.design.p.vents.specific_resistance == ev.tune.resistance
    assert ev.design.p.cabinet.coax_chamber_net_volume == ev.tune.net_volume

def test_vents_preset_has_no_puck_and_stops_at_600(monkeypatch, tmp_path):  # puck_field never called
```

  In `tests/test_directivity_cli.py`, add `test_trusted_limit_override_in_header`: a fast run with `--trusted-limit 900` shows "900 Hz (given on the command line)".

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_evaluate.py tests/test_directivity_cli.py -v`
Expected: FAIL.

- [ ] **Step 3: Implement**, moving the helpers without changing their behaviour. The existing CLI tests must pass unchanged, apart from import paths.

- [ ] **Step 4: Run the tests**

Run: `uv run pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/speaker tests
git commit -m "feat: evaluate() runs a design at screen or vent settings, tuned or as it stands"
```

---

### Task 7: Round files and variant params

**Files:**
- Modify: `pyproject.toml` (add `"tomli-w"` to dependencies; `uv lock`)
- Create: `src/speaker/variants.py`
- Test: `tests/test_variants.py`

A round file (the format; Task 10 writes `rounds/round-1.toml`):

```toml
number = 1
title = "Cardioid first, then the front shape"
base = "params/default.toml"
trusted_limit_hz = 2851
trusted_note = "inherited from round 0's mesh check (2026-10-27)"
finals = 2

[[compare]]
label = "Round 0 (round 0b run)"
scores = "out/round0b/scores.json"

[[stage]]
name = "vents"
settings = "vents"          # a key of evaluate.SETTINGS
rank_by = "cardioid"        # cardioid | directivity
include_round0 = true       # also score the base file as it stands, untuned, as "round0"
grid = { "vents.centre_from_front_mm" = [140, 190, 240, 280], "vents.slots_per_side" = [6, 9, 12] }

[[stage]]
name = "shape"
settings = "screen"
rank_by = "directivity"
include_round0 = true
[[stage.variant]]
name = "w330-ro30-m260"
set = { "cabinet.width_mm" = 330, "cabinet.front_edge" = "roundover", "cabinet.edge_size_mm" = 30, "cabinet.flare_mouth_d_mm" = 260 }
```

Every stage also scores its own base, tuned, as `stage-base`. Override keys are dotted paths into the base TOML.

**Interfaces:**
- Produces:
  - `@dataclass(frozen=True) VariantSpec(name: str, overrides: dict[str, Any])`.
  - `StageSpec(name, settings: str, rank_by: str, include_round0: bool, variants: tuple[VariantSpec, ...])`.
  - `Compare(label: str, scores: Path)`.
  - `RoundSpec(number: int, title: str, base: Path, trusted_limit: float, trusted_note: str, finals: int, compare: tuple[Compare, ...], stages: tuple[StageSpec, ...])`.
  - `load_round(path: Path) -> RoundSpec`:
    - Paths are relative to the repo root.
    - A grid expands as a cartesian product in sorted key order. Its variants are named by joining `"{last key part without _mm}={value}"` with `"-"`, for example `centre_from_front=140-slots_per_side=6`.
    - Explicit variants follow the grid's.
  - `read_raw(path) -> dict`, `apply_overrides(raw: dict, overrides: dict) -> dict` (deep copy), `write_params(raw: dict, path: Path) -> None`.
  - `with_tuned(raw: dict, resistance: float, net_volume: float) -> dict`: sets `vents.specific_resistance_rayl` and `cabinet.coax_chamber_net_volume_l` (litres).
  - `RESERVED = ("stage-base", "round0")`.

`load_round` rejects these with a ValueError naming the stage, the variant and the problem:
- an unknown settings name or `rank_by`;
- a duplicate or reserved variant name, or a name outside `[A-Za-z0-9._=-]+`;
- an override key absent from the base TOML;
- `finals < 0`;
- a stage with no variants.

- [ ] **Step 1: Write the failing tests**

```python
def test_unknown_override_key_names_it(tmp_path):
    spec = ROUND.replace('"cabinet.width_mm" = 330', '"cabinet.widht_mm" = 330')
    with pytest.raises(ValueError, match=r"cabinet\.widht_mm.*w330-ro30-m260"):
        load_round(write(tmp_path, spec))

def test_grid_expands_in_sorted_key_order(tmp_path):
    r = load_round(write(tmp_path, ROUND))
    names = [v.name for v in r.stages[0].variants]
    assert len(names) == 12 and names[0] == "centre_from_front=140-slots_per_side=6"

def test_written_variant_loads_with_the_override(tmp_path):
    raw = apply_overrides(read_raw(DEFAULT_PATH), {"cabinet.width_mm": 330, "cabinet.front_edge": "chamfer"})
    write_params(raw, tmp_path / "v.toml")
    p = load_params(tmp_path / "v.toml")
    assert p.cabinet.width == pytest.approx(0.330) and p.cabinet.front_edge == "chamfer"

def test_with_tuned_writes_rayl_and_litres(tmp_path):  # loads back to R and V in SI
def test_reserved_and_duplicate_names_rejected(tmp_path):
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_variants.py -v`
Expected: FAIL (module missing).

- [ ] **Step 3: Implement** `variants.py` and add the dependency.

- [ ] **Step 4: Run the tests**

Run: `uv run pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml uv.lock src/speaker/variants.py tests/test_variants.py
git commit -m "feat: round files: stages of variants as overrides of a base params file"
```

---

### Task 8: Ranking and the round runner

**Files:**
- Create: `src/speaker/ranking.py`, `src/speaker/round.py`
- Modify: `src/speaker/cli.py` (`speaker round ROUND_TOML [--out DIR]`)
- Test: `tests/test_ranking.py`, `tests/test_round.py`

**Interfaces:**
- Consumes: `evaluate`, `SETTINGS`, `run_directivity(..., trusted_limit=)`, and the `variants.py` functions.
- Produces:
  - `@dataclass VariantResult`:
    - `name`, `stage`, `status: Literal["ok", "rejected", "error"]`, `reason: str`, `overrides: dict`;
    - `params_sha256` (of the pre-tune `params.toml`), `settings: str`, `code_version: str`, `commit: str`;
    - `scores: list[dict]` (`dataclasses.asdict(Score)`), `tune: dict | None` (`asdict(VentTune)`), `tuned: bool`, `wall_s: float`.
    - It is stored as `result.json`.
  - `ranking.rank(results: list[VariantResult], rank_by: str) -> Ranking`, where `Ranking(order: list[str], out: dict[str, str], gate_unmet: bool)`:
    - **cardioid:** order by the worst-case margin, descending: `tune["margin"]` if `tuned`, else `tune["base"]["margin"]`. Any status that is not ok is out.
    - **directivity:** a status that is not ok is out. A trusted cardioid_hard with status "fail" is out, with the reason "fails the hard cardioid gate (x.x dB)". The rest are ordered ascending on the key `(floor((window_spread + presence_lift) / 0.25), di_smooth_bass + di_smooth_lfhf, flat23)`, using trusted values.
    - **If every ok variant fails the gate:** `gate_unmet = True`, and all ok variants are ranked by cardioid_hard (descending) and then by the key.
  - `round.run_round(spec_path: Path, out_dir: Path | None = None) -> Path`. Out defaults to `out/rounds/<number>`. It returns `out_dir / "results.json"`; Task 9 changes this to the report path.
    1. **For each stage in order.** The base raw params are the round's base for stage 1, else the previous stage's best `params_tuned.toml`. The best is the first in the order, excluding `round0`.
    2. **Variants to run:** `stage-base` (no overrides), then `round0` if `include_round0` (the round's base file, `tune=False`), then the stage's variants.
    3. **Each variant** gets its own directory, `out_dir/<stage>/<name>/`:
       - Write `params.toml`.
       - If `result.json` exists with the same `params_sha256`, `settings` and `code_version`, reuse it.
       - Otherwise call `evaluate(..., tune=(name != "round0"))` and write `params_tuned.toml` and `scores.json`.
       - Write `result.json` last, atomically (a temp file, then `os.replace`).
       - A ValueError means `rejected`, with `str(e)` as the reason. Any other exception means `error`, with the traceback's last 20 lines. The round continues either way.
       - Print one progress line to stderr per variant: `[stage i/n] name: status, wall s`.
    4. **After each stage**, write `out_dir/<stage>/ranking.json`.
    5. **Finals:** the top `finals` of the last stage's order (excluding `round0`). Each gets `run_directivity(params_tuned, out_dir/"finals"/name, trusted_limit=(spec.trusted_limit, spec.trusted_note))`, and `out_dir/"finals.json"` lists them in rank order.
    6. **Finally** write `out_dir/results.json`: every stage's results and ranking, plus the spec path and its sha256.

- [ ] **Step 1: Write the failing tests.** Use canned scores (`tests/test_ranking.py`), and monkeypatch `speaker.round.evaluate` and `speaker.round.run_directivity` with fakes (`tests/test_round.py`).

```python
def test_rank_excludes_gate_failures():
def test_rank_when_every_variant_fails_the_gate():
    r = rank([res("a", hard=6.0, ws=2, pl=2), res("b", hard=8.0, ws=3, pl=3)], "directivity")
    assert r.gate_unmet and r.order == ["b", "a"]
def test_rank_bucket_ties_go_to_di_smoothness():
    # window+presence 4.05 vs 4.20 (same 0.25 dB bucket); the lower DI sum wins
def test_round_carries_the_best_vents_into_the_next_stage(tmp_path, fake_eval):
    # stage 2's stage-base params.toml has stage 1's winner's tuned rayl and litres
def test_round_resumes_and_ignores_partial_results(tmp_path, fake_eval):
    # run once; delete one result.json and leave a stray result.json.<pid>.tmp; run again:
    # evaluate is called exactly once more, for that variant only
def test_round_records_rejected_variants(tmp_path, fake_eval):
    # an override breaking rule (e) gives status "rejected" with the flare-mouth message; the others are ok
def test_round_records_errors_and_continues(tmp_path, fake_eval):
def test_finals_are_the_top_of_the_last_stage(tmp_path, fake_eval):
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_ranking.py tests/test_round.py -v`
Expected: FAIL (modules missing).

- [ ] **Step 3: Implement** `ranking.py`, `round.py` and the CLI subcommand.

- [ ] **Step 4: Run the tests**

Run: `uv run pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/speaker tests
git commit -m "feat: speaker round runs a round's stages, tunes and ranks each variant, then the finals"
```

---

### Task 9: The round report

**Files:**
- Create: `src/speaker/round_report.py`, `src/speaker/shape_drawing.py`
- Modify: `src/speaker/round.py` (writes the report; `--report-only`), `src/speaker/cli.py`
- Test: `tests/test_round_report.py`, `tests/test_shape_drawing.py`

**Interfaces:**
- Consumes: the `results.json`, `ranking.json` and `finals.json` written by Task 8, each final's `report.html`, `scores.json` and `plots/`, and the `compare` scores files.
- Produces:
  - `shape_drawing.draw_designs(path: Path, designs: list[tuple[str, Params]]) -> None`. One column per design (the recommendation, then round 0), at the same scale, with two views:
    - **Front elevation:** the outline, with the edge treatment's inner lines, the flare mouth, the bass cutouts and the vent stack drawn dashed.
    - **Plan section** at the coax centre height: `front_shape.section`.
  - `round_report.write_report(out_dir: Path, spec: RoundSpec) -> Path`, which writes `out_dir/report.html`. Its sections are in spec §7's order, with these headings:
    1. **"Recommendation".** The text under `## Recommendation` in `out_dir/notes.md` (paragraphs and `- ` lists, HTML-escaped), or else "Draft: the recommendation is written in the next session." Then a key-dimensions table: rank 1 against round 0, covering width, front edge and size, flare mouth, profile and outer layer, coax offset, top margin, vent position, slots, R and V, and bass heights. Then `shape.png`.
    2. **"Curves".** Rank 1's `polar_h`, `polar_v`, `spinorama` (both DIs), `window_family` and `rear_rejection` plots, copied into `out_dir/plots/`.
    3. **"Scores".** Rows: the recommendation, the other finals, then each `compare` entry. Columns: every trusted score key, plus cardioid hard and aim, with pass, warn and fail styling. The low-confidence values go in a second, greyed table. If `gate_unmet` is set, a red line above the table says the hard cardioid gate is unmet by every variant.
    4. **"Confidence".**
       - The trusted limit and its note, and that variants are not mesh-checked.
       - The screen settings (5 elements per wavelength, to the trusted limit) against the finals (6, with the puck).
       - The coax LF de-embedding (Task 2).
       - The felt spread.
    5. **"Crossover".** The fixed LR4 points and each final's way gains. Plan 3 replaces it.
    6. **"What I'd try next".** From `## What I'd try next` in `notes.md`, or a placeholder.
    7. **Appendix "All variants".** One table per stage: rank, name, status, the reason if out, the key scores, and the tuned R (rayl) and V (L), and what tuning bought (margin − base margin). For a stage ranked by cardioid, the key scores are cardioid_hard, cardioid_aim and rear_hemisphere only: its grid stops at 600 Hz, so its directivity scores are partial and are not shown.
  - `speaker round ROUND_TOML --report-only` re-renders the report from the stored results without evaluating anything.

- [ ] **Step 1: Write the failing tests.** Build a round directory from canned JSON and two tiny fake plot PNGs.

```python
def test_report_sections_in_spec_order(round_dir):
    html = write_report(round_dir, SPEC).read_text()
    heads = ["Recommendation", "Curves", "Scores", "Confidence", "Crossover", "What I'd try next", "All variants"]
    assert [html.index(f">{h}<") for h in heads] == sorted(html.index(f">{h}<") for h in heads)

def test_notes_fill_the_recommendation(round_dir):  # notes.md text appears, escaped
def test_placeholder_without_notes(round_dir):      # "Draft: the recommendation is written in the next session."
def test_rejected_variants_listed_with_reason(round_dir):
def test_gate_unmet_banner(round_dir):
def test_report_only_does_not_evaluate(round_dir, monkeypatch):  # evaluate patched to raise; report still written
```

  In `tests/test_shape_drawing.py`, add `test_draws_each_family`: the PNG is written and non-empty for square, roundover, chamfer and facets.

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_round_report.py tests/test_shape_drawing.py -v`
Expected: FAIL.

- [ ] **Step 3: Implement** both modules, following `directivity_report.py`'s CSS and `_row` conventions. `run_round` now returns the report path.

- [ ] **Step 4: Run the tests**

Run: `uv run pytest -q`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/speaker tests
git commit -m "feat: round report in spec section 7 order, with design drawings against round 0"
```

---

### Task 10: Round 1: run it and recommend

This task is mostly running jobs and writing up. It has no TDD cycle. The final whole-branch review should come **before** Step 4, so that review fixes don't invalidate an 11-hour run.

**Files:**
- Create: `rounds/round-1.toml`, `docs/rounds/round-1.md`
- Not committed: `out/rounds/1/` (including `notes.md`)

- [ ] **Step 1: Write `rounds/round-1.toml`.**
  - **Header** as in Task 7. `compare` points at `out/round0b/scores.json`, already copied there from the round 0b run.
  - **Stage "vents"** (`settings = "vents"`, `rank_by = "cardioid"`, `include_round0 = true`) is a grid of `vents.centre_from_front_mm` [120, 160, 200, 240, 280] × `vents.slots_per_side` [6, 9, 12]. Infeasible combinations are rejected by the fillers' rules at no cost.
  - **Stage "shape"** (`settings = "screen"`, `rank_by = "directivity"`, `include_round0 = true`) has these variants, all on the vents winner (mm). W means `cabinet.width_mm`; m means `cabinet.flare_mouth_d_mm`.

| Name | Overrides |
|---|---|
| w330-sq-m290 | W 330 |
| w320-ro15-m280 | roundover 15, m 280 |
| w330-ro20-m280 | W 330, roundover 20, m 280 |
| w330-ro30-m260 | W 330, roundover 30, m 260 |
| w330-ro40-m240 | W 330, roundover 40, m 240 |
| w320-ro35-m240 | roundover 35, m 240 |
| w330-ch20-m280 | W 330, chamfer 20, m 280 |
| w330-ch30-m260 | W 330, chamfer 30, m 260 |
| w330-ch40-m240 | W 330, chamfer 40, m 240 |
| w320-ch35-m240 | chamfer 35, m 240 |
| w330-fa15-step-m260 | W 330, facets 15°, width 30, base 800, ramp 0, m 260 |
| w330-fa25-step-m260 | W 330, facets 25°, width 30, base 800, ramp 0, m 260 |
| w330-fa15-ramp-m240 | W 330, facets 15°, width 40, base 700, ramp 100, m 240 |
| w330-fa25-ramp-m240 | W 330, facets 25°, width 40, base 700, ramp 100, m 240 |
| w330-ro30-m260-ob30 | as w330-ro30-m260, `outer_baffle_t_mm` 30 |
| w330-ro30-m260-cone | as w330-ro30-m260, `flare_profile` "cone" |
| w330-ro20-m260-in10 | W 330, roundover 20, m 260, `coax_offset_x_mm` +10 |
| w330-ro20-m260-out10 | the same, offset −10 |
| w330-ro30-m260-top60 | as w330-ro30-m260, `top_margin_mm` 60 |
| w330-ch30-m260-bass | as w330-ch30-m260, `bass_centres_z_mm` [260, 580] |

  - Run `uv run python -c "from speaker.variants import load_round; load_round('rounds/round-1.toml')"`. Then apply each variant's overrides and call `compute_layout`, listing the ones that are rejected and why. A rejected variant is fine; record the list in the ledger.

- [ ] **Step 2: Measure before committing the night.**
  - Run `uv run pytest -m bem -q` (all pass).
  - Time one SCREEN evaluation of `w330-ro30-m260` and one VENTS evaluation of the default, cold cache, with `uv run python -c ...evaluate(...)`.
  - Estimate the total: vents × VENTS time + shape × SCREEN time + `finals` × (the full run, about 40 min of cabinet plus about 75 min of puck per distinct width and flare).
  - If it is over 12 h, drop shape variants from the bottom of the table until it fits, and record the ruling.
  - Commit the round file.

```bash
git add rounds/round-1.toml
git commit -m "chore: round 1: vents grid, then front-shape variants"
```

- [ ] **Step 3: Final whole-branch review** (per the executing skill), and its fixes, before the run.

- [ ] **Step 4: Run round 1 in the background**

Run: `mkdir -p out/rounds/1 && nohup uv run speaker round rounds/round-1.toml > out/rounds/1/run.log 2>&1 &`
Expected: it finishes within the budget, `out/rounds/1/report.html` exists, and the "All variants" tables list every variant with a status.

- [ ] **Step 5: Mesh-check the recommendation.** The recommendation is rank 1 of the finals, unless the curves argue otherwise; if so, say why in the notes.

Run: `uv run speaker directivity --mesh-check --params out/rounds/1/finals/<name>/params_tuned.toml --out out/rounds/1/finals/<name>`, then the full run again in the same directory (cache hits).
Expected: its own `trusted_limit_hz`, which is reported next to the inherited 2851 Hz.

- [ ] **Step 6: Write `out/rounds/1/notes.md`** (`## Recommendation`, `## What I'd try next`). Argue from the curves, with the dB each change bought against the stage base and round 0 (spec §2 step 5 and §7). Then re-render: `uv run speaker round rounds/round-1.toml --report-only`.

- [ ] **Step 7: Write `docs/rounds/round-1.md`**, in the style of `docs/rounds/round-0.md`. Include:
  - what ran, the commit, timings and the rejected variants;
  - the cardioid result after tuning, and whether the hard gate now passes;
  - the shape ranking and the recommendation;
  - the Task 2 recess finding on open question 1;
  - what remains open.

```bash
git add docs/rounds/round-1.md
git commit -m "docs: round 1 results and recommendation"
```

  Then present the report to the owner for review (spec §2 step 6).
