# Foundation Pipeline (Plan 1 of 3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a parametric Python pipeline that, from one parameter file,
does all of the following:
- models the drivers
- simulates the sealed bass, the room/boundary effects and a lumped cardioid
  model of the coax
- builds the cabinet in CadQuery
- exports STEP and CNC-ready DXF files
- writes a summary report

**Architecture:** This is a `uv`-managed Python package called `speaker`, laid
out under `src/`, with one module per concern:
- params
- drivers
- bass
- room
- cardioid
- layout
- CAD panels
- export
- report/CLI

Physics modules are pure NumPy, use SI units, and know nothing about CAD. The
layout module turns parameters and driver data into positions and volumes,
which both the physics and the CAD consume. `speaker build` regenerates
everything into `out/`.

**Tech Stack:** Python 3.12 (needed for the CadQuery wheels), uv, NumPy,
SciPy, matplotlib, CadQuery 2.8, ezdxf 1.4, pytest. The install of all of these
has already been verified on this machine with Python 3.12.

**Spec:** `docs/superpowers/specs/2026-10-22-cardioid-coax-floorstander-design.md`

**Later plans (not this one):**
- Plan 2: BEM directivity of the baffle and flare, spinorama on the 0–30°
  window, VituixCAD FRD/ZMA export, and crossover targets.
- Plan 3: panel-resonance FEM and the bracing study.

This plan produces a working, testable pipeline on its own.

## Global Constraints

- Python `==3.12.*`, managed with `uv`. Commands are run as `uv run ...`.
- **Units:**
  - Code is SI throughout: metres, Hz, Pa, m³, kg.
  - The parameter TOML uses millimetres and litres, with `_mm` / `_l` key
    suffixes. `params.py` is the only place that converts them.
  - CadQuery geometry is in millimetres, and `cad/` is the only place that
    converts.
- Physical constants (in `speaker/constants.py`): `C = 343.0` m/s and
  `RHO = 1.204` kg/m³.
- External width 320 mm. **External depth must be ≤ 350 mm: loading a
  parameter set with depth > 350 raises `ValueError`.**
- Coax acoustic centre at seated ear height, default 950 mm.
- Panels default to 18 mm. The front baffle is two 18 mm layers (36 mm).
- No grille. The cardioid vents are arrays of several small slots on both side
  panels, near the rear.
- Both cabinets are identical (one model).
- **DXF conventions:**
  - One layer per cut depth. Through-cuts go on layer `THROUGH`; pockets go on
    `POCKET_<depth>MM`, with the depth formatted with one decimal, e.g.
    `POCKET_10.0MM`.
  - Notes go on layer `NOTES`.
  - Inside corners need a tool radius of at least 3 mm.
- While `drivers_measured = false` in the parameters, every DXF carries the
  text `PROVISIONAL - drivers not yet measured` on `NOTES`, and the filename
  ends in `_PROVISIONAL.dxf`.
- The coax flare is a separate 3D-printed insert (exported as STEP/STL). The
  baffle only gets simple 2.5D openings.
- Generated outputs go in `out/`, which is git-ignored.

## Review Focus

1. **Parameter edits that break physical fit** must raise a clear `ValueError`
   naming the clash, never produce silently broken CAD. Examples: a coax OD
   wider than the interior, a bass basket deeper than the interior, bass
   drivers overlapping each other or the divider. → Task 6.
2. **Changing ear height or ceiling height** (e.g. 850 / 1050 mm) must move
   the coax, the divider and the cabinet height consistently, and keep the
   volumes valid. → Task 6.
3. **Frequency arrays containing 0 Hz or negative values** must raise
   `ValueError` and never return NaN or inf. → Task 3.
4. **Extreme vent resistance** (0 or `math.inf`) must give a finite dipole or
   sealed result, not a division error. → Task 5.
5. **Unit slips between mm and m** must be caught. Exported STEP bounding
   boxes are checked in mm and layout volumes in m³. → Tasks 6 and 8.

---

## File Structure

```
pyproject.toml                 package metadata, deps, `speaker` console script
.gitignore                     out/, .venv/, __pycache__/
params/default.toml            every tunable dimension, room and drive value
data/drivers/bms_10c262.toml
data/drivers/dayton_rss265ho_44.toml
src/speaker/__init__.py
src/speaker/constants.py       C, RHO
src/speaker/params.py          TOML -> frozen dataclasses (SI), validation
src/speaker/drivers.py         Driver dataclass, derived T/S, loader
src/speaker/bass.py            sealed alignment, SPL, excursion, max SPL
src/speaker/room.py            image-source boundary model, axial modes
src/speaker/cardioid.py        lumped coax+chamber+vent model, polar
src/speaker/layout.py          derived positions/volumes, fit validation
src/speaker/cad/__init__.py
src/speaker/cad/panels.py      CadQuery panels + flare insert + assembly
src/speaker/cad/export.py      STEP / STL / DXF writers
src/speaker/report.py          plots + summary.md
src/speaker/cli.py             `speaker build [--params P] [--out DIR]`
tests/test_params.py ... tests/test_cli.py   (one per module)
```

---

### Task 1: Project scaffold and parameters

**Files:**
- Create: `pyproject.toml`, `.gitignore`, `params/default.toml`, `src/speaker/__init__.py`, `src/speaker/constants.py`, `src/speaker/params.py`
- Test: `tests/test_params.py`

**Interfaces:**
- Produces:
  - `load_params(path: Path | None = None) -> Params`. With `None` it loads
    `params/default.toml`, resolved relative to the repo root.
  - Frozen dataclasses, all lengths in metres and volumes in m³:
    - `Cabinet(width, depth, panel_t, baffle_t, coax_centre_z, top_margin, divider_t, side_edge_radius, flare_mouth_d, bass_centres_z: tuple[float, float], coax_chamber_net_volume)`
    - `Vents(slots_per_side: int, slot_w, slot_h, slot_pitch, centre_from_front, specific_resistance)`.
      `specific_resistance` is in Pa·s/m.
    - `Room(length_x, width_y, height_z, back_gap, speaker_y: tuple[float, float], seat: tuple[float, float, float])`
    - `Drive(amp_v_rms, inductor_dcr, active_v_rms)`
    - `Params(cabinet, vents, room, drive, drivers_measured: bool, coax: str, bass: str)`.
      `coax` and `bass` are driver file stems.

**Default values for `params/default.toml`.** All are from the spec §5, §5.1
and §7, or from the cardioid prototype noted below.

| Key | Value |
|---|---|
| `cabinet.width_mm` | 320 |
| `cabinet.depth_mm` | 350 |
| `cabinet.panel_t_mm` | 18 |
| `cabinet.baffle_t_mm` | 36 |
| `cabinet.coax_centre_z_mm` | 950 |
| `cabinet.top_margin_mm` | 20 (clearance above the coax frame) |
| `cabinet.divider_t_mm` | 18 |
| `cabinet.side_edge_radius_mm` | 0 |
| `cabinet.flare_mouth_d_mm` | 290 |
| `cabinet.bass_centres_z_mm` | [200, 520] |
| `cabinet.coax_chamber_net_volume_l` | 10 |
| `vents.slots_per_side` | 6 |
| `vents.slot_w_mm` | 100 (along the depth, y) |
| `vents.slot_h_mm` | 10 (tall, z) |
| `vents.slot_pitch_mm` | 22 (stacked vertically) |
| `vents.centre_from_front_mm` | 270 |
| `vents.specific_resistance_rayl` | 100 |
| `room.length_x_mm` | 4000 |
| `room.width_y_mm` | 4600 |
| `room.height_z_mm` | 2000 |
| `room.back_gap_mm` | 75 |
| `room.speaker_y_mm` | [960, 3460] |
| `room.seat_mm` | [3425, 2210, 950] |
| `drive.amp_v_rms` | 25.0 (SA550: about 150 W into 4 Ω) |
| `drive.inductor_dcr_ohm` | 0.4 |
| `drive.active_v_rms` | 40.0 |
| `drivers_measured` | false |
| `coax` | "bms_10c262" |
| `bass` | "dayton_rss265ho_44" |

Where the cardioid defaults come from: in a prototype lumped model, a 10 L
chamber with 12 000 mm² of vent area (12 slots of 10 × 100 mm), placed 0.30 m
behind the cone, gave 19 dB minimum front/back from 150 to 300 Hz. A 20 L
chamber gave only 3–6 dB.

- [ ] **Step 1: Write the failing tests**

```python
from pathlib import Path
import pytest
from speaker.params import load_params

def test_defaults_load_in_si():
    p = load_params()
    assert p.cabinet.width == pytest.approx(0.320)
    assert p.cabinet.depth == pytest.approx(0.350)
    assert p.cabinet.coax_centre_z == pytest.approx(0.950)
    assert p.cabinet.coax_chamber_net_volume == pytest.approx(0.010)
    assert p.vents.slots_per_side == 6
    assert p.room.seat == pytest.approx((3.425, 2.210, 0.950))
    assert p.drivers_measured is False

def test_depth_over_350_rejected(tmp_path):
    text = Path("params/default.toml").read_text().replace("depth_mm = 350", "depth_mm = 351")
    f = tmp_path / "p.toml"; f.write_text(text)
    with pytest.raises(ValueError, match="depth"):
        load_params(f)
```

- [ ] **Step 2:** Run `uv run pytest tests/test_params.py -v`. Expected: FAIL (module not found).
- [ ] **Step 3: Create the pyproject, gitignore, constants and default TOML, and implement `load_params`**

  `pyproject.toml`:
  - `requires-python = "==3.12.*"`
  - deps: `numpy`, `scipy`, `matplotlib`, `cadquery>=2.8`, `ezdxf>=1.4`
  - dev group: `pytest`
  - `[project.scripts] speaker = "speaker.cli:main"`
  - src layout

  Parse the TOML with `tomllib`. Convert `_mm` keys to m and `_l` keys to m³.

- [ ] **Step 4:** Run `uv sync && uv run pytest tests/test_params.py -v`. Expected: PASS.
- [ ] **Step 5: Commit** with `git add -A && git commit -m "feat: project scaffold and parameter loading"`

---

### Task 2: Driver models

**Files:**
- Create: `data/drivers/bms_10c262.toml`, `data/drivers/dayton_rss265ho_44.toml`, `src/speaker/drivers.py`
- Test: `tests/test_drivers.py`

**Interfaces:**
- Produces:
  - `Driver`, a frozen dataclass with these fields:
    - `name`, `fs`, `qes`, `qms`, `vas` (m³), `sd` (m²), `re`, `le` (H), `xmax` (m)
    - `outer_d`, `cutout_d`, `depth`, `flange_t`, `basket_d`, `magnet_d` (m)
    - `sens_db` (datasheet, informational)
    - optional `hf_sens_db`
  - Derived properties: `qts`, `cms`, `mms`, `bl`, `rms`. Formulas:
    - `cms = vas / (RHO*C**2*sd**2)`
    - `mms = 1/((2πfs)²·cms)`
    - `bl = sqrt(2πfs·mms·re/qes)`
    - `rms = 2πfs·mms/qms`
  - `load_driver(stem: str) -> Driver`

**Driver data** (spec §6.1). The Dayton is stored **with its coils in
series** (the datasheet values are for parallel coils: Re 1.6 Ω, Le 0.77 mH).
Add a TOML comment explaining the conversion.

| Key | bms_10c262 | dayton_rss265ho_44 (series) |
|---|---|---|
| fs | 65 | 27.6 |
| qes / qms | 0.40 / 4.6 | 0.50 / 4.42 |
| vas_l | 25.5 | 34.8 |
| sd_cm2 | 320 | 352.3 |
| re / le_mh | 5.9 / 0.5 | 6.4 / 3.08 |
| xmax_mm | 5.5 | 13.25 |
| outer_d / cutout_d / depth / flange_t (mm) | 260 / 230 / 146 / 10 | 266 / 235 / 136.5 / 10 |
| basket_d / magnet_d (mm) | 228 / 156 | 235 / 180 |
| sens_db (2.83 V) | 96 (hf_sens_db 112) | 84.5 |

The coax cutout of 230 mm and the Dayton flange and magnet sizes are
provisional estimates. Mark them `# PROVISIONAL` in the TOML.

- [ ] **Step 1: Write the failing tests.** These anchors were computed from
  the formulas above, and they agree with the datasheets: coax Cms
  0.176 mm/N; Dayton Bl 2 × 9.39 = 18.8.

```python
import pytest
from speaker.drivers import load_driver

def test_coax_derived():
    d = load_driver("bms_10c262")
    assert d.qts == pytest.approx(0.367, abs=0.005)
    assert d.cms == pytest.approx(1.758e-4, rel=0.01)
    assert d.mms == pytest.approx(0.0341, rel=0.01)
    assert d.bl == pytest.approx(14.33, rel=0.01)

def test_dayton_series_derived():
    d = load_driver("dayton_rss265ho_44")
    assert d.qts == pytest.approx(0.450, abs=0.005)
    assert d.bl == pytest.approx(19.31, rel=0.01)   # datasheet 2*9.39=18.78, within 3%
    assert d.mms == pytest.approx(0.168, rel=0.01)  # datasheet 161 g
    assert d.re == pytest.approx(6.4)

def test_unknown_driver():
    with pytest.raises(FileNotFoundError):
        load_driver("nope")
```

- [ ] **Step 2:** Run `uv run pytest tests/test_drivers.py -v`. Expected: FAIL.
- [ ] **Step 3:** Create the two TOMLs and implement `Driver` / `load_driver` in `src/speaker/drivers.py`.
- [ ] **Step 4:** Run `uv run pytest tests/test_drivers.py -v`. Expected: PASS.
- [ ] **Step 5: Commit** with message `feat: driver data and derived T/S parameters`.

---

### Task 3: Sealed bass model

**Files:**
- Create: `src/speaker/bass.py`
- Test: `tests/test_bass.py`

**Interfaces:**
- Consumes: `Driver` (Task 2), and `C` and `RHO`.
- Produces:
  - `Alignment(fc: float, qtc: float, f3: float)`
  - `sealed_alignment(d: Driver, vb: float, series_r: float = 0.0) -> Alignment`.
    `series_r` raises Qes to `qes*(re+series_r)/re`.
  - `sealed_spl(d, vb, f, v_rms, r=1.0, n_drivers=1, series_r=0.0) -> np.ndarray`
    - Returns dB SPL in **half-space** (radiation into 2π, a pressure factor
      of 2 over free field).
    - Includes `re + series_r + jωLe`.
    - `n_drivers` drivers share one voltage and sum coherently, each in
      volume `vb`.
  - `sealed_excursion(d, vb, f, v_rms, series_r=0.0) -> np.ndarray`, returning
    peak displacement in metres.
  - `max_spl(d, vb, f, v_rms_limit: float | None, r=1.0, n_drivers=1, series_r=0.0) -> np.ndarray`
    - The lower of the Xmax-limited SPL and, if `v_rms_limit` is not `None`,
      the voltage-limited SPL.
  - All functions raise `ValueError` if any `f <= 0`.

  How to compute it: solve the electro-mechanical series circuit for cone
  velocity:
  - `Zm = jωMms + Rms + 1/(jωCms) + 1/(jωCab)·Sd²`, with `Cab = vb/(RHO C²)`
  - velocity `u = Bl·V / ((Re+Rs+jωLe)·Zm + Bl²)`
  - `p_rms = RHO·ω·Sd·|u| / (2π r)` per driver, using the half-space
    point-source formula.

- [ ] **Step 1: Write the failing tests.** The anchors come from spec §6.1 and
  a hand calculation.

```python
import numpy as np, pytest
from speaker.drivers import load_driver
from speaker.bass import sealed_alignment, sealed_spl, max_spl, sealed_excursion

D = load_driver("dayton_rss265ho_44")

@pytest.mark.parametrize("vb,fc,qtc,f3", [(0.025, 42.7, 0.70, 43.4), (0.030, 40.6, 0.66, 43.6)])
def test_alignment_matches_spec(vb, fc, qtc, f3):
    a = sealed_alignment(D, vb)
    assert a.fc == pytest.approx(fc, abs=0.3)
    assert a.qtc == pytest.approx(qtc, abs=0.01)
    assert a.f3 == pytest.approx(f3, abs=0.5)

def test_inductor_dcr_raises_qtc():
    assert sealed_alignment(D, 0.025, series_r=0.4).qtc > sealed_alignment(D, 0.025).qtc

def test_passband_sensitivity():          # datasheet 84.5 dB @ 2.83 V series-wired
    spl = sealed_spl(D, 0.030, np.array([150.0]), 2.83)
    assert spl[0] == pytest.approx(84.7, abs=1.0)

def test_xmax_limited_spl_40hz():          # half-space 1 m, single driver, hand calc 106.0 dB
    spl = max_spl(D, 0.030, np.array([40.0]), None)
    assert spl[0] == pytest.approx(106.0, abs=0.3)

def test_two_drivers_add_6db():
    f = np.array([60.0])
    assert (sealed_spl(D, 0.03, f, 2.83, n_drivers=2) - sealed_spl(D, 0.03, f, 2.83))[0] == pytest.approx(6.02, abs=0.05)

def test_excursion_below_xmax_at_amp_limit_80hz():
    assert sealed_excursion(D, 0.03, np.array([80.0]), 25.0)[0] < D.xmax

def test_zero_frequency_rejected():
    with pytest.raises(ValueError):
        sealed_spl(D, 0.03, np.array([0.0, 50.0]), 2.83)
```

- [ ] **Step 2:** Run `uv run pytest tests/test_bass.py -v`. Expected: FAIL.
- [ ] **Step 3:** Implement `src/speaker/bass.py` as specified. For `f3`, use a
  closed form or `scipy.optimize.brentq` on the second-order high-pass
  magnitude `|H|² = x⁴/((1−x²)² + x²/Q²) = 0.5`.
- [ ] **Step 4:** Run `uv run pytest tests/test_bass.py -v`. Expected: PASS.
- [ ] **Step 5: Commit** with message `feat: sealed bass alignment, SPL, excursion`.

---

### Task 4: Room and boundary model

**Files:**
- Create: `src/speaker/room.py`
- Test: `tests/test_room.py`

**Interfaces:**
- Consumes: `Room` (Task 1).
- Produces:
  - `Surface = Literal["front", "back", "left", "right", "floor", "ceiling"]`
    - "front" is the wall behind the speakers, at x=0.
    - "back" is the wall behind the listener, at x=length_x.
    - y=0 is "left", y=width_y is "right"; z=0 is "floor", z=height_z is
      "ceiling".
  - `boundary_response(src: tuple[float, float, float], rcv: tuple[float, float, float], f: np.ndarray, room: Room, surfaces: Iterable[Surface] = ALL, reflection: float = 1.0) -> np.ndarray`
    - Returns the complex pressure transfer: the direct sound plus
      first-order image sources, each `e^{-jkr}/r`, normalised so the direct
      sound alone gives `1/r_direct` magnitude.
  - `sbir_null(d: float) -> float` returns `C/(4d)`.
  - `axial_modes(room: Room, n_max: int = 3) -> dict[str, list[float]]` with
    keys `"x"`, `"y"` and `"z"`; each mode is `n·C/(2L)`.

- [ ] **Step 1: Write the failing tests.** The anchors come from spec §5.1.

```python
import numpy as np, pytest
from speaker.params import load_params
from speaker.room import boundary_response, axial_modes, sbir_null

R = load_params().room

def test_front_wall_null_matches_quarter_wave():
    d = 0.6
    f = np.linspace(80, 260, 2000)
    h = boundary_response((d, 2.3, 1.0), (d + 3.0, 2.3, 1.0), f, R, surfaces=["front"])
    assert f[np.argmin(np.abs(h))] == pytest.approx(sbir_null(d), rel=0.03)
    assert sbir_null(0.6) == pytest.approx(142.9, abs=0.1)
    assert sbir_null(0.4) == pytest.approx(214.4, abs=0.1)

def test_floor_bounce_null():               # ear and coax at 0.95 m, 3 m apart -> 311 Hz
    f = np.linspace(200, 450, 2000)
    h = boundary_response((1.0, 2.3, 0.95), (4.0, 2.3, 0.95), f, R, surfaces=["floor"])
    assert f[np.argmin(np.abs(h))] == pytest.approx(311.2, rel=0.02)

def test_axial_modes():
    m = axial_modes(R)
    assert m["x"] == pytest.approx([42.9, 85.8, 128.6], abs=0.1)
    assert m["y"][:2] == pytest.approx([37.3, 74.6], abs=0.1)
    assert m["z"][0] == pytest.approx(85.8, abs=0.1)

def test_no_surfaces_is_direct_only():
    f = np.array([100.0, 1000.0])
    h = boundary_response((0.5, 1, 1), (3.5, 1, 1), f, R, surfaces=[])
    assert np.abs(h) == pytest.approx([1/3.0, 1/3.0])
```

- [ ] **Step 2:** Run `uv run pytest tests/test_room.py -v`. Expected: FAIL.
- [ ] **Step 3:** Implement `src/speaker/room.py`.
- [ ] **Step 4:** Run `uv run pytest tests/test_room.py -v`. Expected: PASS.
- [ ] **Step 5: Commit** with message `feat: image-source boundary model and axial modes`.

---

### Task 5: Lumped cardioid coax model

**Files:**
- Create: `src/speaker/cardioid.py`
- Test: `tests/test_cardioid.py`

**Interfaces:**
- Consumes: `Driver` (Task 2), and `C` and `RHO`.
- Produces:
  - `CardioidParams(chamber_volume: float, vent_area: float, vent_length: float, vent_resistance: float, vent_positions: tuple[tuple[float, float, float], ...], n_slots: int = 12)`
    - `vent_area` is the total over all slots (m²).
    - `vent_resistance` is the specific resistance in Pa·s/m. It may be `0`
      or `math.inf`.
    - Vent positions are relative to the cone centre, with x forward.
      Total vent flow is split equally between the positions.
  - `CardioidResult(f, u_cone: np.ndarray, u_vent: np.ndarray, excursion: np.ndarray)`
    - The `u_*` values are complex volume velocities (m³/s).
    - `excursion` is the peak cone displacement in metres.
  - `simulate(d: Driver, p: CardioidParams, f: np.ndarray, v_rms: float) -> CardioidResult`
  - `polar(res: CardioidResult, p: CardioidParams, angles_deg: np.ndarray, plane: Literal["horizontal", "vertical"] = "horizontal") -> np.ndarray`
    - Free-field far-field magnitude in dB, with shape `(len(f), len(angles))`,
      normalised to 0° at each frequency.
    - Horizontal angles rotate about z; vertical angles rotate about y.
  - `front_back(res, p) -> np.ndarray`, in dB (0° minus 180°).
  - `tune_resistance(d, p, f_lo=150.0, f_hi=300.0) -> float`
    - Returns the specific resistance that maximises the minimum
      `front_back` over `[f_lo, f_hi]`.
    - Use a log-spaced search from 10 to 5000 Pa·s/m, then
      `scipy.optimize.minimize_scalar` bounded.

**Algorithm.** This is not implied by the signatures, so use it as written:

```
Cab = V/(RHO C²);  Rv = r_spec/A;  a = sqrt(A/n_slots/π);  Mv = RHO (L + 1.7a)/A
Yv  = 0 if r_spec is inf else 1/(Rv + jωMv)
Zbox = 1/(jωCab + Yv)
Zm  = jωMms + Rms + 1/(jωCms) + Sd²·Zbox + Bl²/(Re + jωLe)
u   = Bl·V/((Re + jωLe)·Zm)        # cone velocity
Uc  = Sd·u ;  Uv = −Uc·Zbox·Yv     # vent flow (outward); equals −Uc when open, 0 when sealed
p(θ) ∝ | Uc + Σ_i (Uv/N)·exp(j k (r̂(θ)·pos_i)) |
```

This is a **low-frequency estimate only**. It ignores the baffle and the
cabinet, so above about 400 Hz the front/back figure it gives is
pessimistic. Plan 2's BEM replaces it there. Say this in the module
docstring.

- [ ] **Step 1: Write the failing tests.** The anchors come from the prototype
  noted in Task 1.

```python
import math, numpy as np, pytest
from speaker.drivers import load_driver
from speaker.cardioid import CardioidParams, simulate, front_back, tune_resistance, polar

D = load_driver("bms_10c262")
F = np.geomspace(100, 1000, 200)
VENTS = ((-0.30, 0.16, 0.0), (-0.30, -0.16, 0.0))

def P(r):
    return CardioidParams(0.010, 0.012, 0.018, r, VENTS)

def test_sealed_is_omni_and_matches_box():
    res = simulate(D, P(math.inf), F, 2.83)
    assert np.allclose(front_back(res, P(math.inf)), 0.0, atol=1e-6)
    assert np.all(res.u_vent == 0)

def test_open_back_is_symmetric_front_back():
    res = simulate(D, P(0.0), F, 2.83)
    assert np.all(np.isfinite(res.u_vent))
    assert np.allclose(front_back(res, P(0.0)), 0.0, atol=0.5)

def test_tuned_vents_reach_10db_150_to_300():
    r = tune_resistance(D, P(100.0))
    res = simulate(D, P(r), F, 2.83)
    fb = front_back(res, P(r))
    band = (F >= 150) & (F <= 300)
    assert fb[band].min() >= 10.0
    assert 30 < r < 400

def test_excursion_tiny_at_1w():
    res = simulate(D, P(100.0), np.array([150.0]), 2.83)
    assert res.excursion[0] < 0.5e-3

def test_polar_shape():
    res = simulate(D, P(100.0), F, 2.83)
    assert polar(res, P(100.0), np.arange(0, 181, 10)).shape == (200, 19)
```

- [ ] **Step 2:** Run `uv run pytest tests/test_cardioid.py -v`. Expected: FAIL.
- [ ] **Step 3:** Implement `src/speaker/cardioid.py` exactly per the algorithm above.
- [ ] **Step 4:** Run `uv run pytest tests/test_cardioid.py -v`. Expected: PASS.
- [ ] **Step 5: Commit** with message `feat: lumped passive-cardioid coax model`.

---

### Task 6: Cabinet layout and fit validation

**Files:**
- Create: `src/speaker/layout.py`
- Test: `tests/test_layout.py`

**Interfaces:**
- Consumes: `Params` (Task 1), `Driver` (Task 2), `CardioidParams` (Task 5).
- Produces:
  - `Layout`, a frozen dataclass, all in SI.
  - **Coordinate frame:** origin at the bottom-left-front corner of the
    cabinet; x is across the width, y points backwards (depth), z points up.
    The front face is at y=0.
  - Fields:
    - `height` (external)
    - `interior_w`, `interior_d`
    - `coax_centre: tuple[float, float]`, as (x, z)
    - `bass_centres: tuple[tuple[float, float], tuple[float, float]]`
    - `divider_z`, the bottom face of the divider
    - `coax_chamber_gross`, `coax_chamber_net`
    - `filler_volume`, the volume to take out of the coax chamber to reach
      `coax_chamber_net_volume`
    - `bass_gross`, `bass_net` (total for both drivers)
    - `vent_slots: tuple[tuple[float, float, float, float], ...]`, as
      (y_centre, z_centre, w, h) per slot, the same on both sides
  - `compute_layout(p: Params, coax: Driver, bass: Driver) -> Layout`. Rules:
    - `height = coax_centre_z + coax.outer_d/2 + top_margin + panel_t`
    - `interior_w = width − 2·panel_t`
    - `interior_d = depth − baffle_t − panel_t`
    - divider top face `= coax_centre_z − coax.outer_d/2 − 10 mm`
    - Coax chamber spans from the divider top to the underside of the top
      panel.
    - Driver displacement is estimated as a cylinder of
      `basket_d × (depth − flange_t)` per driver.
    - Vents: `slots_per_side` horizontal slots, all centred at
      `y = centre_from_front`, stacked at `slot_pitch` along z and vertically
      centred on the coax chamber.
  - `cardioid_params(p: Params, l: Layout, coax: Driver) -> CardioidParams`
    - Vent positions are relative to the cone, which is taken as a point
      source on the coax axis at `y = baffle_t` (the mounting plane). x is
      forward, so a vent's x is `−(centre_from_front − baffle_t)`.
    - Vent lateral offsets are the side panels at ±`width/2`.
    - `vent_length = panel_t`
    - `n_slots = 2·slots_per_side`
  - `validate(l: Layout, p: Params, coax: Driver, bass: Driver) -> None`.
    Raises `ValueError`, naming the problem, when:
    - `coax.outer_d > interior_w`, or `bass.outer_d > interior_w`
    - `coax.depth − coax.flange_t > interior_d`, and likewise for the bass
      driver
    - the bass frames overlap each other
    - the bass frames overlap the floor panel (`z − outer_d/2 < panel_t`) or
      the divider
    - `coax_chamber_gross < coax_chamber_net_volume`
    - the vent slots run past the back panel
    - `bass_net < 0.050` m³, the spec minimum of 25 L per driver

    `compute_layout` calls `validate`.

- [ ] **Step 1: Write the failing tests**

```python
import dataclasses, pytest
from speaker.params import load_params
from speaker.drivers import load_driver
from speaker.layout import compute_layout

P = load_params(); CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")

def test_default_layout():
    l = compute_layout(P, CX, B)
    assert l.height == pytest.approx(0.950 + 0.130 + 0.020 + 0.018, abs=1e-6)  # 1.118 m
    assert l.interior_w == pytest.approx(0.284)
    assert l.interior_d == pytest.approx(0.296)
    assert l.bass_net >= 0.050                 # ≥25 L per driver
    assert l.coax_chamber_net == pytest.approx(0.010)
    assert len(l.vent_slots) == 6

@pytest.mark.parametrize("ear", [0.85, 1.05])
def test_ear_height_moves_everything_consistently(ear):
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, coax_centre_z=ear))
    l = compute_layout(p, CX, B)
    assert l.coax_centre[1] == pytest.approx(ear)
    assert l.divider_z < ear - CX.outer_d / 2
    assert l.height == pytest.approx(ear + 0.168, abs=1e-6)

def test_too_low_ear_height_rejected():     # divider would hit the upper bass driver
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, coax_centre_z=0.70))
    with pytest.raises(ValueError, match="divider|bass"):
        compute_layout(p, CX, B)

def test_coax_wider_than_interior_rejected():
    big = dataclasses.replace(CX, outer_d=0.300)
    with pytest.raises(ValueError, match="coax"):
        compute_layout(P, big, B)

def test_overlapping_bass_rejected():
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, bass_centres_z=(0.20, 0.30)))
    with pytest.raises(ValueError, match="overlap"):
        compute_layout(p, CX, B)
```

- [ ] **Step 2:** Run `uv run pytest tests/test_layout.py -v`. Expected: FAIL.
- [ ] **Step 3:** Implement `src/speaker/layout.py`.
- [ ] **Step 4:** Run `uv run pytest tests/test_layout.py -v`. Expected: PASS. If
  `bass_net` comes out below 50 L with the defaults, report it rather than
  weakening the test. The fix is a parameter decision for the owner.
- [ ] **Step 5: Commit** with message `feat: cabinet layout, volumes and fit validation`.

---

### Task 7: CadQuery panels, flare insert and assembly

**Files:**
- Create: `src/speaker/cad/__init__.py`, `src/speaker/cad/panels.py`
- Test: `tests/test_panels.py`

**Interfaces:**
- Consumes: `Layout`, `Params`, `Driver`.
- Produces:
  - `Feature(kind: Literal["through", "pocket"], depth_mm: float | None, shape: Literal["circle", "rect", "slot"], cx_mm, cy_mm, w_mm, h_mm, corner_r_mm)`
    - `depth_mm` is `None` for through-cuts.
    - Coordinates are panel-local, in mm, with the origin at the panel's
      lower-left corner as cut.
  - `Panel(name: str, outline_w_mm: float, outline_h_mm: float, thickness_mm: float, features: tuple[Feature, ...], placement: cq.Location)`
  - **Joinery (decided):**
    - `baffle_outer` is full width (320) and full height, so it covers the
      side panels' front edges and gives a clean face.
    - `side_left` and `side_right` are full height and run from behind the
      outer baffle to the rear (332 mm).
    - `baffle_inner`, `top`, `bottom`, `back` and `divider` fit between the
      sides (284 mm wide).
  - `build_panels(l: Layout, p: Params, coax: Driver, bass: Driver) -> list[Panel]`.
    Names are fixed:
    - `baffle_outer`: through-cuts for the flare insert opening
      (Ø `flare_mouth_d`) and for the bass openings (Ø `bass.cutout_d`)
    - `baffle_inner`:
      - a through-cut for the coax at Ø `coax.cutout_d`, with a pocket
        `coax.flange_t` deep and Ø `coax.outer_d + 2 mm` on the rear face
        for rear mounting
      - bass openings at Ø `bass.cutout_d`, with pockets Ø `bass.outer_d + 2`
        and `flange_t` deep
      - bolt holes at Ø 6 mm on each driver's bolt circle. The circle is at
        `(outer_d + cutout_d)/2`; the coax has 8 holes and the bass drivers 6.
    - `side_left` and `side_right`: vent slots as through-cuts with a 3 mm
      corner radius
    - `back`: a through-cut for a 100 × 160 mm recessed terminal cup
    - `top`, `bottom`, `divider`
  - `panel_solid(panel: Panel) -> cq.Workplane`, the panel with its features
    cut.
  - `flare_insert(l: Layout, p: Params, coax: Driver) -> cq.Workplane`
    - A revolved flare: throat Ø `coax.cutout_d` at the rear, mouth Ø
      `flare_mouth_d` at the front face.
    - Axial length `baffle_t`.
    - The profile is a quarter-ellipse.
    - Wall thickness at least 3 mm.
    - The outer surface is a cylinder of Ø `flare_mouth_d − 0.2 mm`, for a
      push fit.
  - `assembly(l, p, coax, bass) -> cq.Assembly`, all panels placed, plus the
    insert.

- [ ] **Step 1: Write the failing tests**

```python
import pytest
from speaker.params import load_params
from speaker.drivers import load_driver
from speaker.layout import compute_layout
from speaker.cad.panels import build_panels, panel_solid, assembly, flare_insert

P = load_params(); CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")
L = compute_layout(P, CX, B)

def test_panel_names():
    names = {pn.name for pn in build_panels(L, P, CX, B)}
    assert names == {"baffle_outer", "baffle_inner", "side_left", "side_right",
                     "back", "top", "bottom", "divider"}

def test_assembly_bbox_mm():
    bb = assembly(L, P, CX, B).toCompound().BoundingBox()
    assert (bb.xlen, bb.ylen, bb.zlen) == pytest.approx((320, 350, L.height * 1000), abs=0.5)

def test_side_has_six_vent_slots():
    side = next(pn for pn in build_panels(L, P, CX, B) if pn.name == "side_left")
    slots = [f for f in side.features if f.shape == "slot"]
    assert len(slots) == 6 and all(f.corner_r_mm >= 3 for f in slots)

def test_inner_corner_radii_at_least_3mm():
    for pn in build_panels(L, P, CX, B):
        for f in pn.features:
            if f.shape in ("rect", "slot"):
                assert f.corner_r_mm >= 3.0

def test_flare_insert_dims():
    bb = flare_insert(L, P, CX).val().BoundingBox()
    assert max(bb.xlen, bb.ylen) == pytest.approx(289.8, abs=0.2)
    assert bb.zlen == pytest.approx(36.0, abs=0.2)
```

- [ ] **Step 2:** Run `uv run pytest tests/test_panels.py -v`. Expected: FAIL.
- [ ] **Step 3:** Implement `src/speaker/cad/panels.py`. Convert m→mm only in this module.
- [ ] **Step 4:** Run `uv run pytest tests/test_panels.py -v`. Expected: PASS.
- [ ] **Step 5: Commit** with message `feat: parametric CadQuery panels, flare insert and assembly`.

---

### Task 8: STEP, STL and DXF export

**Files:**
- Create: `src/speaker/cad/export.py`
- Test: `tests/test_export.py`

**Interfaces:**
- Consumes: `build_panels`, `panel_solid`, `assembly` and `flare_insert`
  (Task 7).
- Produces:
  - `export_all(l, p, coax, bass, out_dir: Path) -> dict[str, list[Path]]`,
    with keys `"step"`, `"stl"` and `"dxf"`. It writes:
    - `step/cabinet_assembly.step`
    - `step/<panel>.step` for each panel
    - `step/flare_insert.step` and `stl/flare_insert.stl`
    - `dxf/<panel>[_PROVISIONAL].dxf` for each panel

  DXF content, in ezdxf with units set to mm:
  - The outline as a closed polyline on `THROUGH`.
  - Each through-feature on `THROUGH`.
  - Each pocket on `POCKET_<depth>MM`, with the depth formatted `f"{d:.1f}"`.
  - If `not p.drivers_measured`: an MTEXT `PROVISIONAL - drivers not yet measured`
    on `NOTES`, and the `_PROVISIONAL` filename suffix.

- [ ] **Step 1: Write the failing tests**

```python
import dataclasses, ezdxf, cadquery as cq, pytest
from speaker.params import load_params
from speaker.drivers import load_driver
from speaker.layout import compute_layout
from speaker.cad.export import export_all

P = load_params(); CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")
L = compute_layout(P, CX, B)

def test_exports(tmp_path):
    files = export_all(L, P, CX, B, tmp_path)
    assert len(files["dxf"]) == 8 and all(f.name.endswith("_PROVISIONAL.dxf") for f in files["dxf"])
    bb = cq.importers.importStep(str(tmp_path / "step/cabinet_assembly.step")).val().BoundingBox()
    assert (bb.xlen, bb.ylen) == pytest.approx((320, 350), abs=0.5)
    assert (tmp_path / "stl/flare_insert.stl").stat().st_size > 0

def test_dxf_layers(tmp_path):
    export_all(L, P, CX, B, tmp_path)
    doc = ezdxf.readfile(tmp_path / "dxf/baffle_inner_PROVISIONAL.dxf")
    layers = {e.dxf.layer for e in doc.modelspace()}
    assert {"THROUGH", "POCKET_10.0MM", "NOTES"} <= layers
    ext = ezdxf.bbox.extents(e for e in doc.modelspace() if e.dxf.layer == "THROUGH")
    assert ext.size.x == pytest.approx(320 - 2 * 18, abs=0.5)   # inner baffle fits between sides

def test_measured_drops_provisional(tmp_path):
    files = export_all(L, dataclasses.replace(P, drivers_measured=True), CX, B, tmp_path)
    assert not any("PROVISIONAL" in f.name for f in files["dxf"])
```

- [ ] **Step 2:** Run `uv run pytest tests/test_export.py -v`. Expected: FAIL.
- [ ] **Step 3:** Implement `src/speaker/cad/export.py`.
- [ ] **Step 4:** Run `uv run pytest tests/test_export.py -v`. Expected: PASS.
- [ ] **Step 5: Commit** with message `feat: STEP/STL/DXF export with CNC layer conventions`.

---

### Task 9: Report and `speaker build` CLI

**Files:**
- Create: `src/speaker/report.py`, `src/speaker/cli.py`
- Modify: `README.md` (create it). Add a short "Usage" section:
  `uv sync`, `uv run speaker build`, and what lands in `out/`.
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `build(params_path: Path | None, out_dir: Path) -> Path`, returning the
    path of `summary.md`.
  - `main(argv: list[str] | None = None) -> int`, which accepts `build`,
    `--params PATH` and `--out DIR` (default `out/`).
  - `report.write_summary(...)` and plotting helpers (the builder chooses the
    signatures).
- What `out_dir` must contain:
  - Everything from Task 8.
  - `plots/bass_spl.png`: the 2-driver sealed SPL at 2.83 V, plus max SPL
    for passive (`amp_v_rms`, `inductor_dcr`) and active (`active_v_rms`),
    at the seat distance of 3 m.
  - `plots/bass_in_room.png`: both bass drivers at their `bass_centres_z`,
    through `boundary_response` with all surfaces, at the seat.
  - `plots/cardioid_polar.png`: horizontal polar at 150, 200, 300, 450 and
    600 Hz, with the tuned resistance.
  - `plots/front_back.png`: front/back ratio from 100 Hz to 1 kHz.
  - `summary.md`, with these headings, each containing the listed numbers:
    - `## Bass`: Fc, Qtc and F3 (with and without the inductor DCR); net
      volume; max SPL at 40 and 80 Hz at the seat, passive and active.
    - `## Cardioid`: tuned specific resistance; minimum front/back over
      150–300 Hz; coax excursion at the passive amp limit at 150 Hz versus
      Xmax.
    - `## Room`: the SBIR null for the current `back_gap`; the seat-wall
      null; axial modes; the floor-bounce null.
    - `## Cabinet`: external H×W×D; chamber volumes; the list of exported
      files.
    - `## Confidence`: one line per result, marked
      high / medium / needs-measurement. The cardioid is "needs-measurement"
      above 400 Hz.

- [ ] **Step 1: Write the failing test**

```python
from speaker.cli import main

def test_build_end_to_end(tmp_path):
    assert main(["build", "--out", str(tmp_path)]) == 0
    s = (tmp_path / "summary.md").read_text()
    for h in ("## Bass", "## Cardioid", "## Room", "## Cabinet", "## Confidence"):
        assert h in s
    for png in ("bass_spl", "bass_in_room", "cardioid_polar", "front_back"):
        assert (tmp_path / "plots" / f"{png}.png").stat().st_size > 0
    assert (tmp_path / "step/cabinet_assembly.step").exists()
```

- [ ] **Step 2:** Run `uv run pytest tests/test_cli.py -v`. Expected: FAIL.
- [ ] **Step 3:** Implement `report.py`, `cli.py` and the README section. Use
  the matplotlib `Agg` backend.
- [ ] **Step 4:** Run `uv run pytest -v` (the full suite). Expected: all PASS.
  Then run `uv run speaker build` and read `out/summary.md` to sanity-check
  the numbers against the spec:
  - Bass F3 is about 43–44 Hz.
  - The SBIR null for a 75 mm gap is about 200 Hz.
  - Front/back is ≥ 10 dB over 150–300 Hz.

  With the real layout, the vents sit only about 0.23 m behind the cone,
  where the prototype used 0.30 m, so front/back may fall short. **If it does,
  report the number in the summary and to the owner. Don't tune the
  parameters to force it.** It is a design finding.
- [ ] **Step 5: Commit** with message `feat: speaker build CLI with plots and summary report`.
