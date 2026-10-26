# Coax Chamber Fillers and Vents Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Model the coax chamber's printed fillers (four corner pieces plus a
rear block) as real geometry and CAD parts. Double the vent area. Make the
cardioid report use the as-built chamber.

**Architecture:**
- A new pure-geometry module, `speaker/fillers.py` (SI units, no CadQuery),
  computes:
  - the corner-piece profiles, as polygons;
  - the rear block box, solved to hit the target net volume;
  - the volumes of both.
- `compute_layout` calls it and stores the result on `Layout`.
- A new `speaker/cad/fillers.py` extrudes exactly those polygons and boxes,
  so the CAD volumes and the reported volumes come from one source.
- `report.py` tunes and reports the as-built chamber against a hard band and
  an aim band.

**Tech Stack:** Python 3.12, uv, numpy, scipy, CadQuery 2.8, ezdxf 1.4,
matplotlib, pytest.

**Spec:** `docs/superpowers/specs/2026-10-22-cardioid-coax-floorstander-design.md`.
This plan implements §7.1.1 and the revised cardioid target in §7.1.

## Global Constraints

- **Units:**
  - `layout.py` and the new `fillers.py` are SI.
  - `cad/` is millimetres, and converts with `_mm`.
  - Only `params.py` and `drivers.py` convert units from the TOML files.
- **Frame:** origin at the bottom-left-front corner; x across the width, y
  backwards (front face at y = 0), z up. "Left" means smaller x.
- **Cardioid hard requirement:** at least 10 dB minimum front/back over
  150–300 Hz, as built.
- **Cardioid aim:** at least 6 dB minimum over 300–500 Hz. A miss is a
  warning, never traded for 150–300 Hz margin.
- **Vents:** 12 slots per side, 100 mm long (in y) × 10 mm high, 22 mm pitch,
  centred 270 mm from the front.
- **Corner pieces** run from the baffle's rear face (y = `baffle_t`) to
  15 mm before the first vent slot edge.
- **Rear block:**
  - at least 10 mm behind the coax's rear;
  - at least 25 mm from each side panel;
  - net volume equals `coax_chamber_net_volume_l` within ±0.05 L.
- **Fillers are solid for volume purposes.** The parts are sealed, so the
  full outer envelope is subtracted.
- **Filler STEP/STL files** carry the `_PROVISIONAL` suffix while
  `drivers_measured = false`.
- **Rear block screw holes and the back-panel split/gasket are out of scope.**
  They belong to Plan 3. The back panel stays one piece here.
- Every error is a `ValueError` whose message names the part and the numbers
  involved (the existing project convention).

## Review Focus

1. A different coax with a larger frame leaves corner pieces too thin or
   impossible. Expect a clear `ValueError`, not a malformed solid.
   Test: `test_large_coax_frame_rejected_for_corner_fillers` (Task 1).
2. Vents moved forward leave no room for corner pieces. Expect a clear
   `ValueError`. Test: `test_vents_too_far_forward_rejected` (Task 1).
3. A target volume that can't be reached (block taller than the chamber, or
   the corners alone overshoot) must fail loudly, never quietly miss.
   Tests: `test_unreachable_small_target_rejected` and
   `test_corners_alone_overshoot_rejected` (Task 1).
4. Changing ear height moves the coax centre. The fillers must follow it and
   still hit the target.
   Test: `test_fillers_follow_ear_height` (Task 1).
5. Fewer or smaller vents make the aim fail. The summary must say so.
   Test: `test_cardioid_warnings` (Task 3).

---

### Task 1: Filler geometry, layout integration and new vent defaults

**Files:**
- Create: `src/speaker/fillers.py`
- Modify: `src/speaker/layout.py` (the `Layout` fields, `compute_layout`,
  `cardioid_params`, and the volume check in `validate`)
- Modify: `params/default.toml` (`[vents]`)
- Test: `tests/test_fillers.py` (new) and `tests/test_layout.py`

**Interfaces:**
- Consumes:
  - `Driver.outer_d`, `.depth`, `.flange_t` (SI);
  - `Params.cabinet`, `Params.vents`.
- Produces, in `speaker.fillers`:
  ```python
  CORNER_CLEARANCE = 0.003   # m, radial, beyond coax.outer_d / 2
  CORNER_TRIM = 0.040        # m, pieces keep |dx| >= and |dz| >= this from the coax centre lines
  CORNER_MIN_T = 0.005       # m, thinnest allowed wall between arc end and chamber wall
  VENT_GAP = 0.015           # m, corner-piece rear face to first vent slot edge
  MAGNET_GAP = 0.010         # m, coax rear to rear-block front face
  SIDE_CHANNEL = 0.025       # m, rear block to each side panel
  ARC_SEGMENTS = 64          # straight segments per corner arc
  VOLUME_TOL = 0.05e-3       # m^3

  @dataclass(frozen=True)
  class Fillers:
      corners: tuple[tuple[str, tuple[tuple[float, float], ...]], ...]
      # (name, closed polygon of (x, z) vertices in metres, global frame, no repeated first vertex).
      # Names: "filler_corner_top_left", "filler_corner_top_right",
      #        "filler_corner_bottom_left", "filler_corner_bottom_right".
      corner_y: tuple[float, float]               # (y_front, y_back) of all four pieces
      corner_volume: float                        # m^3, all four
      block_min: tuple[float, float, float]       # (x, y, z) m
      block_max: tuple[float, float, float]
      block_volume: float                         # m^3
      @property
      def total_volume(self) -> float: ...

  def polygon_area(pts: Sequence[tuple[float, float]]) -> float   # shoelace, always >= 0
  def corner_profile(cx: float, cz: float, x_wall: float, z_wall: float,
                     r: float, trim: float) -> tuple[tuple[float, float], ...]
  def compute_fillers(p: Params, coax: Driver, *, coax_centre: tuple[float, float],
                      chamber_z: tuple[float, float], volume_to_fill: float) -> Fillers
  ```
- Produces, in `speaker.layout.Layout`:
  - **Changed:**
    - `filler_volume` is now `fillers.total_volume`, and it is in the CAD.
    - `coax_chamber_as_built` is now gross less coax displacement less
      fillers.
  - **New:**
    - `coax_chamber_empty: float` is gross less coax displacement, i.e. the
      old `as_built`.
    - `fillers: Fillers | None = None`. It is `None` only inside
      `compute_layout` before validation.
  - `cardioid_params` now uses `chamber_volume=l.coax_chamber_as_built`.

**Geometry (decided, not in the spec):**
- **Corner profile.** Let `dx = x − cx` and `dz = z − cz`. The corner profile
  is the rectangle `|dx| ∈ [trim, |x_wall − cx|]`,
  `|dz| ∈ [trim, |z_wall − cz|]` minus the disc of radius `r`, where
  `r = coax.outer_d / 2 + CORNER_CLEARANCE`. The quadrant is given by the
  signs of `x_wall − cx` and `z_wall − cz`.
- Its boundary, in quadrant-local `(|dx|, |dz|)`, runs:
  - from `(a, trim)` to `(X, trim)`, then to `(X, Z)`, then to `(trim, Z)`,
    then to `(trim, a)`;
  - then along the arc of radius `r` back to `(a, trim)`, in
    `ARC_SEGMENTS` segments.

  Here `a = sqrt(r² − trim²)`, `X = |x_wall − cx|` and `Z = |z_wall − cz|`.
  Arc vertices lie exactly on the circle. The chords therefore cut slightly
  inside it, and the clearance absorbs this.
- `corner_profile` raises `ValueError` matching "corner filler" if
  `a > min(X, Z) − CORNER_MIN_T`.
- **Chamber walls.**
  - x runs from `panel_t` to `width − panel_t`.
  - z runs over `chamber_z = (divider top face, top-panel underside)`.
- **Corner pieces in y.**
  - `y_front = baffle_t`.
  - `y_back = vent_y_start − VENT_GAP`, where
    `vent_y_start = vents.centre_from_front − vents.slot_w / 2`.
  - If `y_back <= y_front`, raise `ValueError` matching "vent".
- **Rear block.**
  - x runs over `[panel_t + SIDE_CHANNEL, width − panel_t − SIDE_CHANNEL]`.
  - y runs from `max(y_back_corners, baffle_t + coax.depth − coax.flange_t + MAGNET_GAP)`
    to `depth − panel_t`.
  - z is centred on the chamber mid-height.
  - Height is `(volume_to_fill − corner_volume) / (w · d)`.
- **Rear block errors.**
  - If `volume_to_fill − corner_volume < 0`, raise `ValueError` matching
    "corner fillers alone".
  - If the height exceeds the chamber height, raise `ValueError` matching
    "rear block".

- [ ] **Step 1: Write the failing tests** in `tests/test_fillers.py`.
  Use the default params `P`, coax `CX` and bass `B`, as `test_layout.py`
  does. Build helpers like `_cab`/`_vents` in `test_layout.py`.

```python
def test_polygon_area_unit_square():
    assert polygon_area([(0, 0), (1, 0), (1, 1), (0, 1)]) == pytest.approx(1.0)
    assert polygon_area([(0, 1), (1, 1), (1, 0), (0, 0)]) == pytest.approx(1.0)

def test_corner_profile_stays_outside_coax_sweep_and_inside_walls():
    pts = corner_profile(0.160, 0.950, 0.302, 1.100, 0.133, 0.040)
    for x, z in pts:
        assert math.hypot(x - 0.160, z - 0.950) >= 0.133 - 1e-9
        assert 0.160 + 0.040 - 1e-9 <= x <= 0.302 + 1e-9
        assert 0.950 + 0.040 - 1e-9 <= z <= 1.100 + 1e-9

def test_default_fillers():
    l = compute_layout(P, CX, B)
    f = l.fillers
    assert [n for n, _ in f.corners] == ["filler_corner_top_left", "filler_corner_top_right",
                                         "filler_corner_bottom_left", "filler_corner_bottom_right"]
    assert f.corner_y == pytest.approx((0.036, 0.205))
    assert f.corner_volume == pytest.approx(3.85e-3, abs=0.05e-3)
    assert f.block_min[0] == pytest.approx(0.043) and f.block_max[0] == pytest.approx(0.277)
    assert f.block_min[1] == pytest.approx(0.205) and f.block_max[1] == pytest.approx(0.332)
    assert f.block_max[2] - f.block_min[2] == pytest.approx(0.167, abs=0.002)
    assert l.coax_chamber_empty == pytest.approx(0.01883, abs=2e-4)
    assert l.coax_chamber_as_built == pytest.approx(0.010, abs=VOLUME_TOL)
    assert l.filler_volume == pytest.approx(f.total_volume)

def test_large_coax_frame_rejected_for_corner_fillers():
    with pytest.raises(ValueError, match="corner filler"):
        compute_layout(P, dataclasses.replace(CX, outer_d=0.284), B)  # a = 139.4 > 137 mm

def test_vents_too_far_forward_rejected():
    with pytest.raises(ValueError, match="vent"):
        compute_layout(_vents(centre_from_front=0.100), CX, B)

def test_unreachable_small_target_rejected():
    with pytest.raises(ValueError, match="rear block"):
        compute_layout(_cab(coax_chamber_net_volume=0.004), CX, B)

def test_corners_alone_overshoot_rejected():
    with pytest.raises(ValueError, match="corner fillers alone"):
        compute_layout(_cab(coax_chamber_net_volume=0.016), CX, B)

@pytest.mark.parametrize("ear", [0.85, 1.05])
def test_fillers_follow_ear_height(ear):
    l = compute_layout(_cab(coax_centre_z=ear), CX, B)
    assert l.coax_chamber_as_built == pytest.approx(0.010, abs=VOLUME_TOL)
    for _, pts in l.fillers.corners:
        assert all(math.hypot(x - 0.160, z - ear) >= CX.outer_d / 2 for x, z in pts)
```

  In `tests/test_layout.py`:
  - `test_default_layout`: change `len(l.vent_slots) == 6` to `== 12`.
  - Replace `test_as_built_coax_chamber_excludes_filler` with
    `test_as_built_coax_chamber_includes_fillers`. It asserts
    `l.coax_chamber_as_built == approx(l.coax_chamber_empty - l.filler_volume)`
    and `l.coax_chamber_empty == approx(0.01883, abs=2e-4)`.
  - `test_ear_height_moves_everything_consistently` (0.85) may now hit the
    bass-volume warning only. Leave its assertions as they are.

- [ ] **Step 2: Run** `uv run pytest tests/test_fillers.py tests/test_layout.py -q`.
  Expected: FAIL, because `speaker.fillers` does not exist.

- [ ] **Step 3: Implement `speaker/fillers.py`** as specified above, with a
  module docstring stating SI units and the frame.

- [ ] **Step 4: Wire it into `layout.py`.**
  - `compute_layout` builds the layout as now, including `coax_chamber_empty`.
  - It calls `validate`, then `compute_fillers`, with
    `volume_to_fill = coax_chamber_empty − coax_chamber_net`.
  - It then sets `fillers`, `filler_volume` and `coax_chamber_as_built` with
    `dataclasses.replace`. This happens alongside `warnings`.
  - In `validate`, the `filler_volume < 0` check becomes
    `coax_chamber_empty < coax_chamber_net`, keeping its message.
  - `cardioid_params` uses `l.coax_chamber_as_built`.
  - Update the `Layout` field comments to match.

- [ ] **Step 5: Change the `[vents]` defaults in `params/default.toml`.**
  - `slots_per_side = 12`.
  - `specific_resistance_rayl = 190`, with the comment
    `# tuned by the model for 10 L / 24,000 mm2 (2026-10-26)`.

- [ ] **Step 6: Run** `uv run pytest -q`.
  - Expected: everything passes except tests that Tasks 2 and 3 fix:
    `test_side_has_six_vent_slots` and the `test_cli` filler/cardioid
    assertions.
  - Update `test_side_has_six_vent_slots` here. Rename it to
    `test_side_has_twelve_vent_slots` and expect 12.
  - Leave `test_cli` to Task 3, and note the failure in the report.

- [ ] **Step 7: Commit.**

```bash
git add src/speaker/fillers.py src/speaker/layout.py params/default.toml tests/test_fillers.py tests/test_layout.py tests/test_panels.py
git commit -m "feat: coax chamber filler geometry, as-built volume, 12 vent slots per side"
```

---

### Task 2: Filler CAD parts and export

**Files:**
- Create: `src/speaker/cad/fillers.py`
- Modify: `src/speaker/cad/panels.py` (`assembly` adds the fillers)
- Modify: `src/speaker/cad/export.py` (`export_all`)
- Test: `tests/test_cad_fillers.py` (new), `tests/test_panels.py` and
  `tests/test_export.py`

**Interfaces:**
- Consumes:
  - `Layout.fillers` (Task 1): its polygons, `corner_y` and the block box.
  - The `_mm` convention.
- Produces: `filler_solids(l: Layout) -> dict[str, cq.Workplane]`, keyed by
  the five names. The fifth is `"filler_rear_block"`. Solids are in global
  millimetre coordinates, already placed, so the assembly adds them with no
  `loc`.

**Notes:**
- Build each corner as a closed polyline of the `(x, z)` vertices in an
  XZ-plane workplane at `y = corner_y[0]`, extruded to `corner_y[1]`.
  CadQuery's `"XZ"` normal is −y, so check the extrusion direction against
  the bbox test, for example by extruding a negative distance or using a
  custom `Plane`.
- The block is a box from `block_min` to `block_max`.

- [ ] **Step 1: Write the failing tests** in `tests/test_cad_fillers.py`.
  Use the default `L`, `P`, `CX` and `B`.

```python
def test_filler_volumes_match_layout():
    s = filler_solids(L)
    assert set(s) == {n for n, _ in L.fillers.corners} | {"filler_rear_block"}
    corners = sum(s[n].val().Volume() for n, _ in L.fillers.corners)
    assert corners == pytest.approx(L.fillers.corner_volume * 1e9, rel=1e-4)
    assert s["filler_rear_block"].val().Volume() == pytest.approx(L.fillers.block_volume * 1e9, rel=1e-6)

def test_fillers_clear_coax_sweep():
    cx, cz = (v * 1000 for v in L.coax_centre)
    sweep = cq.Solid.makeCylinder(CX.outer_d * 500, 296, cq.Vector(cx, 36, cz), cq.Vector(0, 1, 0))
    for name, wp in filler_solids(L).items():
        assert wp.val().intersect(sweep).Volume() < 1.0, name

def test_corner_pieces_stop_before_vents_and_block_clears_sides():
    s = filler_solids(L)
    for n, _ in L.fillers.corners:
        bb = s[n].val().BoundingBox()
        assert bb.ymin == pytest.approx(36, abs=1e-6) and bb.ymax == pytest.approx(205, abs=1e-6)
    bb = s["filler_rear_block"].val().BoundingBox()
    assert bb.xmin >= 18 + 25 - 1e-6 and bb.xmax <= 302 - 25 + 1e-6

def test_default_parts_fit_250mm_printer():
    for name, wp in filler_solids(L).items():
        bb = wp.val().BoundingBox()
        assert max(bb.xlen, bb.ylen, bb.zlen) < 250, name
```

  In `tests/test_panels.py`, `_assert_no_intersections` now expects
  `len(solids) == 14`: 8 panels, the flare insert, 4 corners and the block.

  In `tests/test_export.py`:
  - `test_exports` also asserts that for each of the five names,
    `step/<name>_PROVISIONAL.step` exists and `stl/<name>_PROVISIONAL.stl`
    is non-empty.
  - `test_measured_drops_provisional` asserts
    `stl/filler_rear_block.stl` exists.

- [ ] **Step 2: Run**
  `uv run pytest tests/test_cad_fillers.py tests/test_panels.py tests/test_export.py -q`.
  Expected: FAIL, with an import error or missing files.

- [ ] **Step 3: Implement `filler_solids`** in `src/speaker/cad/fillers.py`.
  Add each part to `assembly()` by name. In `export_all`, write
  `step/<name><suffix>.step` and `stl/<name><suffix>.stl` for each part, and
  append them to the returned lists. Add a comment that fillers depend on the
  provisional coax frame and magnet sizes.

- [ ] **Step 4: Run** `uv run pytest -q`.
  Expected: everything passes except the `test_cli` assertions left for
  Task 3.

- [ ] **Step 5: Commit.**

```bash
git add src/speaker/cad/fillers.py src/speaker/cad/panels.py src/speaker/cad/export.py tests/test_cad_fillers.py tests/test_panels.py tests/test_export.py
git commit -m "feat: printed coax chamber filler parts in the assembly and STEP/STL export"
```

---

### Task 3: As-built cardioid report with hard and aim bands

**Files:**
- Modify: `src/speaker/report.py`
- Test: `tests/test_report.py` and `tests/test_cli.py`

**Interfaces:**
- Consumes:
  - `Layout.coax_chamber_as_built`, `.coax_chamber_empty` and `.fillers`
    (Task 1);
  - `tune_resistance(d, p)`, which tunes over 150–300 Hz by default.
- Produces, in `report.py`:
  ```python
  FB_TARGET_DB = 10.0; FB_BAND = (150.0, 300.0)      # existing
  FB_AIM_DB = 6.0;     FB_AIM_BAND = (300.0, 500.0)  # new
  def cardioid_warnings(fb_min_hard: float, fb_min_aim: float) -> list[str]
  ```

**Behaviour:**
- **Cases.** The two cases become `"as_built"` (with fillers) and `"empty"`
  (no fillers, for comparison). Each is tuned with `tune_resistance` over
  the 150–300 Hz band only. Each also records `fb_min_aim`, the minimum over
  `FB_AIM_BAND`.
- The polar plot uses `as_built`.
- **`front_back.png`.**
  - Labels: "as built N L (fillers)" and "no fillers N L".
  - Keep the 10 dB line and the 150–300 Hz band.
  - Add a 6 dB line and a shaded band for the 300–500 Hz aim.
- **`cardioid_lines`.** Add the line
  `"  - Minimum front/back over 300-500 Hz: X dB"`. When the aim is missed,
  append `"  (BELOW the 6 dB aim)"`.
- **Cardioid section headings:**
  - `"- As built (with printed fillers: net X L = gross less coax displacement less Y L of fillers):"`
  - `"- Without fillers (for comparison: net X L):"`
- **`cardioid_warnings` (as-built case only)** returns:
  - `"as-built min front/back over 150-300 Hz is X dB, below the 10 dB hard requirement"`
    when `fb_min_hard < 10`;
  - `"as-built min front/back over 300-500 Hz is X dB, below the 6 dB aim"`
    when `fb_min_aim < 6`.
- **Summary.**
  - The summary prints `"- Cardioid warnings: none"`, or the heading
    followed by `"  - <warning>"` lines, mirroring the layout warnings.
  - It replaces the "NOT in the CAD" filler line in the Cabinet section with:
    - `"- Coax chamber fillers (printed ABS, in the CAD): 4 corner pieces {corner L} L (y {a}-{b} mm), rear block {w} x {d} x {h} mm = {block L} L; as-built net {as_built L} L (target {net L} L)"`
  - It changes the coax chamber line to report gross, empty (no fillers)
    and as-built.
  - It changes the Confidence line to
    `"- CAD / volumes: geometry high; net volumes medium (displacement estimated)"`.

- [ ] **Step 1: Write the failing tests.**

```python
# tests/test_report.py
def test_cardioid_warnings():
    assert cardioid_warnings(12.0, 7.0) == []
    w = cardioid_warnings(9.4, 5.2)
    assert len(w) == 2 and "9.4 dB" in w[0] and "10 dB hard requirement" in w[0]
    assert "5.2 dB" in w[1] and "6 dB aim" in w[1]
```

  In `tests/test_cli.py::test_build_end_to_end`, replace the cardioid and
  filler assertions with:

```python
    assert "As built (with printed fillers: net 10.0" in s and "Without fillers" in s
    assert s.count("Minimum front/back over 150-300 Hz") == 2
    assert s.count("Minimum front/back over 300-500 Hz") == 2
    assert "- Cardioid warnings: none" in s
    assert "NOT in the CAD" not in s and "rear block" in s
    assert "geometry high; net volumes medium (displacement estimated)" in s
    assert (tmp_path / "stl/filler_rear_block_PROVISIONAL.stl").exists()
```

- [ ] **Step 2: Run** `uv run pytest tests/test_report.py tests/test_cli.py -q`.
  Expected: FAIL.

- [ ] **Step 3: Implement** the behaviour above in `report.py`.

- [ ] **Step 4: Run** `uv run pytest -q`, then `uv run speaker build`.
  Expected: all tests pass, and `out/summary.md` shows the following for
  the as-built case:
  - tuned R about 190 Pa·s/m;
  - min F/B about 22–23 dB over 150–300 Hz, and about 9 dB over 300–500 Hz;
  - "Cardioid warnings: none".

  Paste those lines into the report.

- [ ] **Step 5: Commit.**

```bash
git add src/speaker/report.py tests/test_report.py tests/test_cli.py
git commit -m "feat: report as-built cardioid with fillers against hard and aim bands"
```
