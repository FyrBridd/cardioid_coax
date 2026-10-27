# Directivity Simulation and Scoring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Simulate the off-axis response of the whole speaker (coax LF with its cardioid vents, coax HF, both bass drivers) on the real cabinet, and score round 0 (the current design) against the directivity targets.

**Architecture:**
- A validated 3D boundary-element solver (bempp-cl, Burton–Miller) radiates each source at unit volume velocity from a gmsh surface mesh of the cabinet exterior. Each frequency band gets its own mesh, sized for that band's top frequency.
- A coax-only "puck" model covers the coax above the trusted 3D limit, up to 12 kHz.
- Source strengths come from the existing lumped models (cardioid chamber, sealed bass) and from the traced datasheet curves. The ways are summed through a fixed ideal LR4 crossover, then turned into a re-centred spinorama and the spec §4 scores.
- `speaker directivity` writes a simple HTML report. `speaker directivity --mesh-check` measures the trusted upper limit.

**Tech Stack:**
- Python 3.12 (uv), numpy and scipy.
- bempp-cl 0.4.x, using the OpenCL backend via `pocl-binary-distribution`, with numba as the fallback.
- gmsh 4.15 (Python API, OCC kernel) and matplotlib.
- pymupdf, for the tracing tool only (dev group).

**Spec:** `docs/superpowers/specs/2026-10-26-directivity-optimisation-design.md` (this plan is its §9 item 1). The parent spec is `docs/superpowers/specs/2026-10-22-cardioid-coax-floorstander-design.md`.

## Global Constraints

- Open-source tools only, running on Linux. The machine has 16 cores, 30 GB RAM, no NVIDIA GPU and no Wine.
- SI units inside the code. Millimetres and litres appear only in TOML files, converted by `params.py` and `drivers.py`. Use `C` and `RHO` from `speaker.constants`.
- **Frame:** the Layout frame. The origin is at the bottom-left-front corner, x runs across the width, y backwards (front face at y=0), and z up.
- **Reference point** for all angles: the coax centre on the front face, `(W/2, 0, coax_centre_z)`.
- **Direction convention:**
  - Horizontal angle θ: `(sin θ, −cos θ, 0)`. Positive θ is towards +x, and the listener sits at **+23°**.
  - Vertical angle φ: `(0, −cos φ, sin φ)`. Positive φ is up.
  - Sphere direction (az, el): `(sin az cos el, −cos az cos el, sin el)`.
- **Observation radius 2 m** for scores. Datasheet curves are at 1 m, 2.83 V.
- **Sources**, named exactly: `coax_lf`, `coax_hf`, `vent_left` (the x=0 side), `vent_right` (the x=W side), `bass_lower`, `bass_upper`. The source field is pressure per **unit volume velocity (1 m³/s)**.
- **Cardioid bands** (spec §4): hard ≥ 10 dB at 180° over 150–300 Hz, measured on the coax section alone (coax_lf plus vents, before any crossover filter). Aim ≥ 6 dB over 300–500 Hz. Rear hemisphere 120–180° over 150–500 Hz is reported.
- **Fixed provisional crossover** (spec §9): ideal acoustic LR4 at 170 Hz (bass/coax) and 1.4 kHz (coax LF/HF). No EQ, no delays, one level-matching gain per way.
- **Listening window** (spec §4, my call): power average of horizontal 0°, 10°, 20°, 30°, plus the two directions at az 23° with el ±10°.
- **The trusted 3D limit** comes from `--mesh-check`. When it hasn't been run, the report uses 3000 Hz and marks it "unverified".
- **Do not commit** the BMS datasheet PDF or images extracted from it. Commit only the traced CSVs, with their source URL and the PDF's sha256 in a header comment.
- **Tests:** BEM-running tests carry `@pytest.mark.bem`. The default `pytest` run skips them (`addopts = "-m 'not bem'"`). Each task names the `-m bem` command to run.
- **gmsh needs `libGLU.so.1`.** The owner installs it once with `sudo apt install libglu1-mesa` before Task 3. If that isn't possible, pin the no-X wheel instead: `gmsh @ https://gmsh.info/python-packages-dev-nox/gmsh/gmsh-4.15.0.dev1+nox-py2.py3-none-manylinux_2_24_x86_64.whl`. Record which one was used as a ledger ruling.

## Review Focus

1. **Burton–Miller sign and interior resonances.** A wrong sign of the hypersingular term gives correct answers at most frequencies but fails near interior resonances: the spike was 5 dB off at 3 kHz on a 0.15 m sphere. The test at ka = 4.493 (Task 1) pins it. Box resonances of the cabinet fall right in the band, so a reviewer should confirm that the solver is the validated one.
2. **Unit-volume-velocity normalisation.** A source patch's ∫v_n dA is computed from the same P1 function that drives the solve, not from the nominal area. Otherwise the cone (curved, inside the cavity) and the vents are mis-scaled against the lumped volume velocities, and the cardioid ratio is wrong. Task 4's monopole test pins it.
3. **Sign of the vent volume velocity.** `cardioid.simulate` returns `u_vent` as an outward volume velocity, already negative relative to the cone. Using it unchanged is right; flipping it destroys the cardioid. Task 7's test compares the BEM rear rejection with the lumped `front_back` at 150 Hz, within 3 dB.
4. **Memory blow-up on a fine mesh.** The dense operators need about 64·N² bytes (N = P1 dofs). A refined mesh check could exhaust 30 GB and take down the session. `solve_band` refuses before assembly when the estimate exceeds `BEM_MEMORY_BUDGET` (Task 1 test).
5. **Stale cache after a code or param change.** The cache key must include the params file content, the mesh settings and the code version (git HEAD, or "dirty" plus the diff hash). Otherwise a changed design silently reuses an old field. Task 4 has a test.

---

## File Structure

| File | Responsibility |
|---|---|
| `src/speaker/bem/__init__.py` | empty package marker |
| `src/speaker/bem/solver.py` | `SurfaceMesh`, `solve_band` (Burton–Miller, multi-source, LU), memory guard, backend selection |
| `src/speaker/bem/analytic.py` | pulsating sphere, cap on a rigid sphere, baffled piston on-axis |
| `src/speaker/coax_geometry.py` | flare and coax-cavity profile shared by the CAD and the mesher |
| `src/speaker/bem/mesh.py` | gmsh: `cabinet_mesh`, `puck_mesh`, `box_piston_mesh` (test), mesh checks |
| `src/speaker/bem/run.py` | frequency grid, observation set, banded runs, `SourceField`, npz cache |
| `src/speaker/curves.py` | load traced curves, minimum phase, fractional-octave smoothing |
| `tools/trace_curve.py` | one-off raster tracer (not part of the package) |
| `data/curves/*.csv` | traced BMS 10C262 LF/HF SPL and LF/HF impedance |
| `src/speaker/system.py` | source strengths, de-embedding, 3D/puck blend, `Ways` |
| `src/speaker/crossover.py` | `lr4`, `FixedCrossover`, `apply_crossover` |
| `src/speaker/spinorama.py` | re-centred CTA-2034-style curves and DIs |
| `src/speaker/scores.py` | spec §4 scores |
| `src/speaker/directivity_report.py` | run pipeline, plots, `report.html`, `scores.json`, mesh check |
| `src/speaker/cli.py` | `directivity` subcommand |
| `docs/rounds/round-0.md` | committed summary of round 0 |

---

### Task 1: BEM solver core, validated against analytic solutions

**Files:**
- Modify: `pyproject.toml` (add the dependencies, the dev-group `pymupdf`, the pytest marker and addopts)
- Create: `src/speaker/bem/__init__.py`, `src/speaker/bem/solver.py`, `src/speaker/bem/analytic.py`
- Test: `tests/test_bem_solver.py`

**Interfaces:**
- Produces:
  - `SurfaceMesh(vertices: np.ndarray (nv,3) m, triangles: np.ndarray (nt,3) int, tags: np.ndarray (nt,) int, tag_names: tuple[str, ...])`. Tag 0 is rigid, and tag i ≥ 1 is source `tag_names[i-1]`.
  - `solve_band(mesh: SurfaceMesh, freqs: np.ndarray, points: np.ndarray (3,np), *, memory_budget: float = BEM_MEMORY_BUDGET) -> np.ndarray (nf, nsrc, np) complex`. The result is pressure per unit volume velocity of each source, with all other surfaces rigid.
  - `BEM_MEMORY_BUDGET = 22e9` (bytes); `estimate_bytes(n_dofs: int) -> float` returns `64.0 * n_dofs**2`.
  - `analytic.pulsating_sphere(f, a, r) -> complex`, `analytic.cap_on_sphere(f, a, alpha, r, theta, n_terms=60) -> complex` (both per unit volume velocity), and `analytic.baffled_piston_on_axis(f, a, r) -> complex` (per unit volume velocity).

- [ ] **Step 1: Add the dependencies.** In `[project].dependencies`, add `"bempp-cl>=0.4.2,<0.5"`, `"pyopencl"`, `"pocl-binary-distribution"`, `"gmsh>=4.15"`. Add `"pymupdf"` to the dev group. Add `markers = ["bem: runs the BEM solver (slow)"]` and `addopts = "-m 'not bem'"` under `[tool.pytest.ini_options]`. Then run `uv sync`.

- [ ] **Step 2: Write the failing tests**

```python
import numpy as np, pytest
from speaker.bem import analytic
from speaker.bem.solver import SurfaceMesh, solve_band, estimate_bytes
import bempp_cl.api as bempp

def _sphere(n, a, cap_alpha=None):
    g = bempp.shapes.regular_sphere(n)
    v, t = g.vertices.T * a, g.elements.T
    cz = v[t].mean(axis=1)[:, 1]  # cap axis is -y (forward), matching the direction convention
    tags = np.ones(len(t), int) if cap_alpha is None else \
        np.where(-cz / np.linalg.norm(v[t].mean(axis=1), axis=1) >= np.cos(cap_alpha), 1, 0)
    return SurfaceMesh(v, t, tags, ("src",))

@pytest.mark.bem
@pytest.mark.parametrize("f", [300.0, 343 * 4.493 / (2 * np.pi * 0.15), 3000.0])  # middle: interior resonance ka=4.493
def test_pulsating_sphere_matches_analytic(f):
    m = _sphere(4, 0.15)
    p = solve_band(m, np.array([f]), np.array([[0.0], [-2.0], [0.0]]))[0, 0, 0]
    assert 20 * np.log10(abs(p) / abs(analytic.pulsating_sphere(f, 0.15, 2.0))) == pytest.approx(0, abs=0.2)

@pytest.mark.bem
def test_cap_on_sphere_polar_matches_analytic():
    m = _sphere(4, 0.15, cap_alpha=np.radians(30))
    th = np.radians([0, 30, 60, 90, 150])
    pts = np.stack([2 * np.sin(th), -2 * np.cos(th), 0 * th])
    p = solve_band(m, np.array([1000.0]), pts)[0, 0]
    ref = np.array([analytic.cap_on_sphere(1000.0, 0.15, np.radians(30), 2.0, t) for t in th])
    err = 20 * np.log10(np.abs(p) / np.abs(ref))
    assert np.all(np.abs(err[:4]) < 0.5) and abs(err[4]) < 1.0

def test_memory_guard_refuses_before_assembly():
    v = np.random.rand(3 * 40000, 3); t = np.arange(3 * 40000).reshape(-1, 3)
    m = SurfaceMesh(v, t, np.ones(40000, int), ("src",))
    with pytest.raises(ValueError, match="memory"):
        solve_band(m, np.array([100.0]), np.zeros((3, 1)), memory_budget=1e9)
    assert estimate_bytes(1000) == 64e6

def test_cap_series_reduces_to_pulsating_sphere_for_full_cap():
    assert abs(analytic.cap_on_sphere(500.0, 0.15, np.pi, 2.0, 0.3)) == \
        pytest.approx(abs(analytic.pulsating_sphere(500.0, 0.15, 2.0)), rel=1e-6)
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `uv run pytest tests/test_bem_solver.py -m "bem or not bem" -v`
Expected: FAIL with `ModuleNotFoundError: speaker.bem`.

- [ ] **Step 4: Implement `analytic.py`**

Use the e^{−iωt} convention (outgoing e^{ikr}). Volume velocities: Q = 4πa²v for the sphere, and Q = 2πa²(1−cos α)v for the cap.

- `pulsating_sphere`: `RHO*C*(-1j*k*a)/(1-1j*k*a) * a/r * exp(1j*k*(r-a)) * v`.
- `cap_on_sphere` (Morse): `p = 1j*RHO*C*v * Σ_n c_n h_n(kr)/h_n'(ka) P_n(cos θ)`, with `c_0 = (1−cos α)/2` and `c_n = (P_{n−1}(cos α) − P_{n+1}(cos α))/2`. Use `scipy.special.spherical_jn/yn` for the spherical Hankel functions of the first kind.
- `baffled_piston_on_axis`: `RHO*C*v*(exp(1j*k*r) - exp(1j*k*sqrt(r**2+a**2)))`, with `v = 1/(π a²)`. Check the sign of the phase prefactor against `pulsating_sphere` at low ka. Magnitudes are what the tests pin.

- [ ] **Step 5: Implement `solve_band` in `solver.py`**

- **Backend:** set `bempp_cl.api.DEFAULT_DEVICE_INTERFACE = "opencl"`. If no OpenCL platform is found, fall back to `"numba"` with a one-time `warnings.warn`.
- **Memory guard:** n_dofs = number of distinct vertices referenced by `triangles` (the P1 nodes). Compute it before any bempp `Grid` is created. Raise `ValueError` mentioning "memory", with the estimate and the budget in GB, if `estimate_bytes(n_dofs) > memory_budget`. Do this before any assembly.
- **Neumann data:** per source, a P1 grid function g, with node value = the mean, over adjacent elements, of (1 on source-tagged elements, 0 elsewhere). Scale g by `1j*RHO*ω` (e^{−iωt}: ∂p/∂n = iωρ v_n). Normalise each source by its volume velocity `Q = ∫ v dA` of that same P1 function, computed with the P1 mass matrix.
- **Assembly**, per frequency (the spike verified this sign):

```python
eta = 1j / k
A = 0.5*I - K + eta*W          # K: double_layer, W: hypersingular (bempp sign convention)
B = -V - eta*(Kp + 0.5*I)      # V: single_layer, Kp: adjoint_double_layer; all P1/P1/P1
```

  Accumulate the dense weak forms into one array, deleting each operator as you go. Factor `A` with `scipy.linalg.lu_factor(overwrite_a=True, check_finite=False)`. Solve all sources' right-hand sides `B @ g_i` in one `lu_solve`.
- **Field:** evaluate `p = DLpot(φ) − SLpot(g)` at `points`. Return `(nf, nsrc, np)`.

- [ ] **Step 6: Run the tests to verify they pass**

Run: `uv run pytest tests/test_bem_solver.py -m "bem or not bem" -v`
Expected: 6 passed (the OpenCL kernel compile takes about 15 s on the first test).

- [ ] **Step 7: Commit**

```bash
git add pyproject.toml uv.lock src/speaker/bem tests/test_bem_solver.py
git commit -m "feat: Burton-Miller BEM solver validated against sphere and cap analytics"
```

---

### Task 2: Shared coax geometry (flare and cavity profile)

**Files:**
- Create: `src/speaker/coax_geometry.py`
- Modify: `src/speaker/cad/panels.py` (import `FLARE_FIT_CLEARANCE`, `FLARE_LIP` and the mouth/throat radii from `coax_geometry`, so there is one source of truth)
- Modify: `src/speaker/drivers.py`, `data/drivers/bms_10c262.toml`, `data/drivers/dayton_rss265ho_44.toml`
- Test: `tests/test_coax_geometry.py`

**Interfaces:**
- Produces:
  - `Driver` gains optional fields:
    - `diaphragm_d: float | None` (m);
    - `cone_depth: float | None` (m, surround plane to cone inner edge);
    - `cone_inner_d: float | None` (m);
    - `hf_exit_d: float | None` (m).
  - `coax_geometry.flare_radii(p, coax) -> tuple[float, float, float]` returns `(r_throat, r_mouth, r_out)` in metres, as the CAD uses them now.
  - `coax_geometry.cavity_profile(p, coax, n_arc: int = 16) -> tuple[tuple[str, tuple[tuple[float, float], ...]], ...]` returns labelled (r, y) polylines from the mouth inwards. The labels are, in order: `"flare"` (quarter-ellipse from (r_mouth, 0) to (r_throat, panel_t)), `"bore"` (r_throat, panel_t → baffle_t), `"surround"` (r_throat → diaphragm_d/2 at y = baffle_t), `"cone"` (diaphragm_d/2 at baffle_t → cone_inner_d/2 at baffle_t + cone_depth), `"pole"` (cone_inner_d/2 → hf_exit_d/2) and `"hf_exit"` (hf_exit_d/2 → 0), the last two at y = baffle_t + cone_depth.

- [ ] **Step 1: Write the failing tests**

```python
import pytest
from speaker.coax_geometry import cavity_profile, flare_radii
from speaker.drivers import load_driver
from speaker.params import load_params

P = load_params(); CX = load_driver("bms_10c262")

def test_new_driver_fields():
    assert CX.diaphragm_d == pytest.approx(0.202)          # datasheet "Effect. Diaphragm Diameter"
    assert CX.cone_depth == pytest.approx(0.065)           # PROVISIONAL estimate
    assert CX.cone_inner_d == pytest.approx(0.0636)        # voice coil diameter
    assert CX.hf_exit_d == pytest.approx(0.0254)           # PROVISIONAL
    assert load_driver("dayton_rss265ho_44").diaphragm_d == pytest.approx(0.2118, abs=5e-4)  # from Sd 352.3 cm2

def test_flare_radii_match_cad_constants():
    r_throat, r_mouth, r_out = flare_radii(P, CX)
    assert r_throat == pytest.approx(0.115) and r_out == pytest.approx(0.1449)
    assert r_mouth == pytest.approx(0.1419)

def test_cavity_profile_is_continuous_and_ordered():
    segs = cavity_profile(P, CX)
    assert [n for n, _ in segs] == ["flare", "bore", "surround", "cone", "pole", "hf_exit"]
    pts = [pt for _, s in segs for pt in s]
    for (n0, a), (n1, b) in zip(segs, segs[1:]):
        assert a[-1] == pytest.approx(b[0])
    assert pts[0] == pytest.approx((0.1419, 0.0)) and pts[-1] == pytest.approx((0.0, 0.036 + 0.065))
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_coax_geometry.py -v`
Expected: FAIL with `ModuleNotFoundError: speaker.coax_geometry`.

- [ ] **Step 3: Implement**

- **TOML:** add `diaphragm_d_mm = 202.0`, `cone_depth_mm = 65.0  # PROVISIONAL: estimate, measure with a profile gauge`, `cone_inner_d_mm = 63.6` and `hf_exit_d_mm = 25.4  # PROVISIONAL`.
- **Dayton:** `diaphragm_d` is derived in `load_driver` from Sd when it's absent: `sqrt(4*sd/pi)`.
- **Flare radii:** use the CAD's existing values (`r_out = (flare_mouth_d − 0.2 mm)/2`, `r_mouth = r_out − 3 mm`, `r_throat = cutout_d/2`).
- **Ellipse points:** generate them exactly as the CAD's quarter-ellipse (centre (r_mouth, panel_t), semi-axes (r_mouth − r_throat, panel_t)).
- **CAD:** keep the CAD geometry unchanged; only move the constants.

- [ ] **Step 4: Run all the tests**

Run: `uv run pytest -v`
Expected: all pass, including the existing CAD/export tests.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/coax_geometry.py src/speaker/cad/panels.py src/speaker/drivers.py data/drivers tests/test_coax_geometry.py
git commit -m "feat: shared coax flare/cavity profile and provisional cone geometry"
```

---

### Task 3: Cabinet and puck surface meshes (gmsh)

**Files:**
- Create: `src/speaker/bem/mesh.py`
- Test: `tests/test_bem_mesh.py`

**Interfaces:**
- Consumes: `compute_layout`, `cavity_profile` (Task 2), `SurfaceMesh` (Task 1).
- Produces:
  - `cabinet_mesh(p, coax, bass, *, h: float, h_coax: float = 0.010, h_vent: float = 0.020) -> SurfaceMesh` with tag_names `("coax_lf", "coax_hf", "vent_left", "vent_right", "bass_lower", "bass_upper")`.
  - `puck_mesh(p, coax, *, h: float, h_coax: float) -> SurfaceMesh` with tag_names `("coax_lf", "coax_hf")`.
  - `box_piston_mesh(size: tuple[float, float, float], piston_r: float, *, h: float) -> SurfaceMesh`, tag `("piston",)`. The box spans `[0, sx] × [0, sy] × [0, sz]` and the flush piston is centred at `(sx/2, 0, sz/2)` on the front face. Used by Task 4's tests.
  - `edge_check(m: SurfaceMesh) -> None` raises `ValueError` unless every edge is shared by exactly two triangles and the signed volume is positive (outward normals).

**Geometry** *(the decisions)*:
- **Exterior:** the cabinet box (width × depth × layout height) less the coax cavity. The cavity is a solid of revolution of `cavity_profile` about the coax axis (along +y from the reference point), closed along the axis.
- **Sources:**
  - `coax_lf` is the `"cone"` segment surface.
  - `coax_hf` is the `"hf_exit"` disc.
  - The `"surround"`, `"pole"`, `"flare"` and `"bore"` surfaces are rigid.
- **Bass drivers:** flush discs of `diaphragm_d` on the front face at `layout.bass_centres` *(my call: rear-mount recess ignored below 250 Hz)*.
- **Vents:** one flush rectangle per side panel, covering the envelope of that side's slot stack: y from `centre_from_front − slot_w/2` to `+ slot_w/2`, z from the lowest slot's bottom to the highest slot's top *(my call: the slot pitch is ≤ λ/4 at 4 kHz)*.
- **Imprinting:** source surfaces are imprinted with `gmsh.model.occ.fragment` and identified by their centre of mass.
- **Element size:** a `Box`/`Distance` field gives `h_coax` within 20 mm of the cavity, `h_vent` on the vent rectangles, and `h` elsewhere, using `min(h, h_coax)` and `min(h, h_vent)`.
- **Puck:** a cylinder of radius `width/2`, coaxial with the coax, depth `baffle_t + cone_depth + 0.030`, with its front rim filleted at 20 mm. The same cavity is subtracted.

- [ ] **Step 1: Write the failing tests**

```python
import math, numpy as np, pytest
from speaker.bem.mesh import cabinet_mesh, puck_mesh, box_piston_mesh, edge_check
from speaker.drivers import load_driver
from speaker.params import load_params

P = load_params(); CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")

def _area(m, name):
    i = m.tag_names.index(name) + 1
    v = m.vertices[m.triangles[m.tags == i]]
    return 0.5 * np.linalg.norm(np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=1).sum()

def test_cabinet_mesh_closed_and_tagged():
    m = cabinet_mesh(P, CX, B, h=0.05)
    edge_check(m)
    assert m.tag_names == ("coax_lf", "coax_hf", "vent_left", "vent_right", "bass_lower", "bass_upper")
    assert _area(m, "bass_lower") == pytest.approx(math.pi * 0.2118**2 / 4, rel=0.03)
    assert _area(m, "coax_hf") == pytest.approx(math.pi * 0.0254**2 / 4, rel=0.05)
    assert _area(m, "vent_left") == pytest.approx(0.100 * (11 * 0.022 + 0.010), rel=0.03)
    slant = math.hypot(0.065, (0.202 - 0.0636) / 2)
    assert _area(m, "coax_lf") == pytest.approx(math.pi * (0.101 + 0.0318) * slant, rel=0.03)

def test_vents_on_correct_sides_and_coax_region_refined():
    m = cabinet_mesh(P, CX, B, h=0.05)
    c = m.vertices[m.triangles].mean(axis=1)
    assert np.allclose(c[m.tags == 3][:, 0], 0.0, atol=1e-6)     # vent_left at x = 0
    assert np.allclose(c[m.tags == 4][:, 0], 0.320, atol=1e-6)   # vent_right at x = W
    e = np.linalg.norm(m.vertices[m.triangles[m.tags == 1][:, 0]] - m.vertices[m.triangles[m.tags == 1][:, 1]], axis=1)
    assert e.max() < 0.010 * 1.5

def test_puck_and_test_box_closed():
    edge_check(puck_mesh(P, CX, h=0.02, h_coax=0.01))
    edge_check(box_piston_mesh((0.2, 0.15, 0.3), 0.015, h=0.02))

def test_edge_check_rejects_open_mesh():
    from speaker.bem.solver import SurfaceMesh
    m = SurfaceMesh(np.eye(3), np.array([[0, 1, 2]]), np.zeros(1, int), ())
    with pytest.raises(ValueError, match="edge"):
        edge_check(m)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_bem_mesh.py -v`
Expected: FAIL with `ModuleNotFoundError: speaker.bem.mesh`.

- [ ] **Step 3: Implement `mesh.py`** as described under Geometry. Initialise and finalise gmsh inside each function, and set `General.Terminal = 0`. Meshing algorithm: 2D Frontal-Delaunay (`Mesh.Algorithm = 6`). Read the nodes and triangles back with `gmsh.model.mesh.getNodes` and `getElementsByType(2)`, per surface entity, to assign tags.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_bem_mesh.py -v`
Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/bem/mesh.py tests/test_bem_mesh.py
git commit -m "feat: gmsh surface meshes of the cabinet exterior and coax puck with tagged sources"
```

---

### Task 4: Banded cabinet runs, observation set and cache

**Files:**
- Create: `src/speaker/bem/run.py`
- Test: `tests/test_bem_run.py`

**Interfaces:**
- Consumes: `solve_band` (Task 1); `cabinet_mesh`, `puck_mesh`, `box_piston_mesh` (Task 3).
- Produces:
  - `freq_grid(f_lo=50.0, f_hi=20000.0, ppo=12) -> np.ndarray` gives log-spaced points that include f_lo.
  - `Observation`, a dataclass with:
    - `dirs: np.ndarray (nd,3)` (unit vectors);
    - `r: float`;
    - `h_angles: np.ndarray` (−180…175 step 5) and `h_idx: np.ndarray` (indices into dirs);
    - `v_angles` and `v_idx`, the same for the vertical ring;
    - `lw_extra_idx: tuple[int, int]` (az 23°, el +10/−10);
    - `sphere_idx: np.ndarray`, `sphere_w: np.ndarray` (solid-angle weights summing to 1).
  - `observation_set(r: float = 2.0) -> Observation`. The sphere grid is az 0…350 step 10 and el −80…80 step 10, plus both poles, with weights ∝ cos(el); the pole weights are ∝ the polar cap area.
  - `SourceField`, a dataclass with: `f: np.ndarray`, `names: tuple[str, ...]`, `obs: Observation`, `p: np.ndarray (nf, nsrc, nd) complex` (per unit volume velocity), `meta: dict`.
  - `cabinet_field(p, coax, bass, f: np.ndarray, *, params_text: str, per_wavelength: float = 6.0, h_max: float = 0.05, cache_dir: Path | None = Path("out/cache/bem")) -> SourceField`.
  - `puck_field` is added in Task 6, reusing this task's banding and cache.
  - `BAND_EDGES = (250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0, 12000.0)`. Frequencies above the last edge are an error.

**Banding:** each frequency is solved on the mesh of the smallest band edge ≥ f, with `h = min(h_max, C / (edge * per_wavelength))`. Frequencies ≤ 250 Hz use the 250 Hz mesh.

**Cache key:** sha256 over:
- `params_text`;
- the driver TOML text;
- `per_wavelength` and `h_max`;
- the function name (cabinet or puck);
- the code version: `git rev-parse HEAD` plus a hash of `git diff HEAD -- src/` when the tree is dirty.

Store one `.npz` per key. Each band's mesh element count and wall time go in `meta`.

- [ ] **Step 1: Write the failing tests**

```python
import numpy as np, pytest
from speaker.bem import analytic
from speaker.bem.mesh import box_piston_mesh
from speaker.bem.run import freq_grid, observation_set, cabinet_field
from speaker.bem.solver import solve_band
from speaker.constants import C, RHO
from speaker.drivers import load_driver
from speaker.params import load_params, DEFAULT_PATH

def test_observation_set():
    o = observation_set(2.0)
    assert len(o.h_idx) == 72 and len(o.v_idx) == 72
    assert np.allclose(o.dirs[o.h_idx[o.h_angles == 0][0]], [0, -1, 0])
    assert np.allclose(o.dirs[o.h_idx[o.h_angles == 90][0]], [1, 0, 0])
    assert o.sphere_w.sum() == pytest.approx(1.0)
    # Sphere weights integrate cos^2 (dipole power) to 1/3.
    d = o.dirs[o.sphere_idx]
    assert (o.sphere_w * d[:, 1] ** 2).sum() == pytest.approx(1 / 3, abs=0.01)

def test_freq_grid():
    f = freq_grid(50, 20000, 12)
    assert f[0] == 50 and f[-1] <= 20000 and np.allclose(np.diff(np.log2(f)), 1 / 12)

@pytest.mark.bem
def test_cabinet_radiates_as_monopole_at_low_frequency(tmp_path):
    P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
    fld = cabinet_field(P, CX, B, np.array([20.0]), params_text=DEFAULT_PATH.read_text(),
                        h_max=0.08, cache_dir=tmp_path)   # ka ~ 0.2 for the whole cabinet
    mono = abs(-1j * 2 * np.pi * 20 * RHO / (4 * np.pi * 2.0))   # |p| per unit volume velocity at 2 m
    lv = 20 * np.log10(np.abs(fld.p[0]) / mono)
    assert np.all(np.abs(lv) < 0.5)              # every source, every direction

@pytest.mark.bem
def test_small_box_baffle_step():
    m = box_piston_mesh((0.2, 0.15, 0.3), 0.015, h=0.012)
    f = np.array([100.0, 2000.0, 3000.0, 4000.0])
    p = np.abs(solve_band(m, f, np.array([[0.0], [-2.0], [0.0]]) + np.array([[0.1], [0], [0.15]]))[:, 0, 0])
    step = 20 * np.log10(p[1:].mean() / p[0]) - 20 * np.log10(f[1:].mean() / f[0])  # remove the ω rise
    assert step == pytest.approx(6.0, abs=1.5)

@pytest.mark.bem
def test_cache_hit_skips_solver(tmp_path, monkeypatch):
    P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
    a = cabinet_field(P, CX, B, np.array([60.0]), params_text="x", h_max=0.1, cache_dir=tmp_path)
    import speaker.bem.run as run
    monkeypatch.setattr(run, "solve_band", lambda *a, **k: (_ for _ in ()).throw(AssertionError("solver called")))
    b = cabinet_field(P, CX, B, np.array([60.0]), params_text="x", h_max=0.1, cache_dir=tmp_path)
    assert np.array_equal(a.p, b.p)
    with pytest.raises(AssertionError, match="solver called"):
        cabinet_field(P, CX, B, np.array([60.0]), params_text="y", h_max=0.1, cache_dir=tmp_path)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_bem_run.py -m "bem or not bem" -v`
Expected: FAIL with `ModuleNotFoundError: speaker.bem.run`.

- [ ] **Step 3: Implement `run.py`** per the Interfaces, Banding and Cache key sections. Evaluate the points as `reference_point + r * dirs`, using the Layout reference point (puck: its own coax centre).

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_bem_run.py -m "bem or not bem" -v`
Expected: 5 passed.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/bem/run.py tests/test_bem_run.py
git commit -m "feat: banded BEM runs over the observation sphere with a keyed result cache"
```

---

### Task 5: Traced datasheet curves

**Files:**
- Create: `tools/trace_curve.py`, `data/curves/bms_10c262_lf_spl.csv`, `data/curves/bms_10c262_hf_spl.csv`, `data/curves/bms_10c262_lf_imp.csv`, `data/curves/bms_10c262_hf_imp.csv`
- Create: `src/speaker/curves.py`
- Test: `tests/test_curves.py`

**Interfaces:**
- Produces:
  - `Curve(f: np.ndarray, value: np.ndarray, unit: str, source: str)`.
  - `load_curve(stem: str) -> Curve` reads `data/curves/<stem>.csv`: `#` comment header lines, a `# unit: dB SPL` or `# unit: ohm` line, then `f,value` rows.
  - `smooth(f, db, fraction: float) -> np.ndarray` gives fractional-octave power smoothing on the given grid.
  - `minimum_phase(f, db) -> np.ndarray` gives the phase in radians.
  - `interp_db(curve, f) -> np.ndarray` interpolates in log f; outside the traced range it holds the end value.

**Tracing** *(the decisions)*:
- **Source:** the datasheet at `https://bmsspeakers.com/wp-content/uploads/2022/01/bms_10c262_datasheet.pdf`, downloaded into a scratch directory and never committed. The PDF's embedded images are image 0 (SPL) and image 1 (impedance), extracted with pymupdf.
- **What the tool takes:** `tools/trace_curve.py` takes the image path, two x-gridline pixel/frequency pairs (log axis), two y-gridline pixel/value pairs, a seed pixel on the wanted curve and an x range. It follows the dark curve column by column, choosing the dark run nearest the predicted y (previous slope), so the LF and HF curves are kept apart where they cross.
- **What it writes:** a CSV at 1/24-octave points, plus an overlay PNG in the scratch directory, which the implementer views and checks against the original.
- **Calibration:** read it off the gridlines: SPL x from 20 Hz to 20 kHz and y from 70 to 120 dB; impedance x from 10 Hz to 20 kHz and y from 0 to 100 Ω.
- **Ranges:**
  - LF SPL: 48 Hz–10.5 kHz.
  - HF SPL: 480 Hz–19.5 kHz. Stop before the plot's 20 kHz cliff and hold beyond.
  - LF impedance: 10 Hz–20 kHz.
  - HF impedance: 200 Hz–20 kHz.

- [ ] **Step 1: Write the failing tests** (spot values read from the datasheet plots)

```python
import numpy as np, pytest
from speaker.curves import load_curve, interp_db, minimum_phase, smooth

@pytest.mark.parametrize("stem,f,val,tol", [
    ("bms_10c262_lf_spl", 200, 95.0, 1.0), ("bms_10c262_lf_spl", 1000, 93.5, 1.0),
    ("bms_10c262_lf_spl", 1600, 89.5, 1.0), ("bms_10c262_lf_spl", 2700, 96.0, 1.0),
    ("bms_10c262_hf_spl", 2000, 111.0, 1.0), ("bms_10c262_hf_spl", 10000, 103.5, 1.5),
    ("bms_10c262_lf_imp", 70, 95.0, 4.0), ("bms_10c262_lf_imp", 300, 7.5, 1.0),
    ("bms_10c262_hf_imp", 750, 26.0, 2.0), ("bms_10c262_hf_imp", 5000, 8.5, 1.0),
])
def test_traced_spot_values(stem, f, val, tol):
    assert interp_db(load_curve(stem), np.array([f]))[0] == pytest.approx(val, abs=tol)

def test_headers_record_provenance():
    c = load_curve("bms_10c262_hf_spl")
    assert "bms_10c262_datasheet.pdf" in c.source and "sha256" in c.source and c.unit == "dB SPL"

def test_minimum_phase_of_first_order_highpass():
    f = np.geomspace(5, 40000, 2000); fc = 500.0
    h = (1j * f / fc) / (1 + 1j * f / fc)
    ph = minimum_phase(f, 20 * np.log10(np.abs(h)))
    m = (f > 50) & (f < 5000)
    assert np.degrees(np.abs(ph[m] - np.angle(h[m]))).max() < 2.0

def test_smooth_keeps_flat_flat():
    f = np.geomspace(20, 20000, 500)
    assert np.allclose(smooth(f, np.full_like(f, 90.0), 1 / 3), 90.0)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_curves.py -v`
Expected: FAIL with `ModuleNotFoundError: speaker.curves`.

- [ ] **Step 3: Write `tools/trace_curve.py`**, download the PDF to the scratchpad, trace the four curves, and view each overlay. Iterate the seed or range until the overlay follows the printed curve everywhere. Write the CSVs with this header:
  - `# source: <url> sha256=<hex>`
  - `# unit: ...`
  - `# traced: tools/trace_curve.py <args>`

- [ ] **Step 4: Implement `curves.py`.**
  - **Minimum phase:** the Hilbert transform of ln|H| on a uniform FFT grid. Extend beyond the data with the end slopes and hold the values. Interpolate back to f.
  - **Smoothing:** a power average over ±fraction/2 octave around each point.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `uv run pytest tests/test_curves.py -v`
Expected: 13 passed.

- [ ] **Step 6: Commit**

```bash
git add tools/trace_curve.py data/curves src/speaker/curves.py tests/test_curves.py
git commit -m "feat: traced BMS 10C262 SPL and impedance curves with minimum-phase helper"
```

---

### Task 6: Coax puck run above the 3D limit

**Files:**
- Modify: `src/speaker/bem/run.py` (add `puck_field`)
- Test: `tests/test_bem_puck.py`

**Interfaces:**
- Consumes: `puck_mesh` (Task 3) and the banding and cache (Task 4).
- Produces: `puck_field(p, coax, f, *, params_text: str, per_wavelength: float = 6.0, cache_dir: Path | None = Path("out/cache/bem")) -> SourceField`, for f from 500 Hz to 12 kHz.
  - The observation set is centred on the puck's coax centre.
  - The on-axis point at 1 m is evaluated as well and stored as `meta["axis_1m"]`, with shape `(nf, nsrc)`.

- [ ] **Step 1: Write the failing test**

```python
import numpy as np, pytest
from speaker.bem.run import puck_field
from speaker.drivers import load_driver
from speaker.params import load_params, DEFAULT_PATH

@pytest.mark.bem
def test_hf_beam_narrows_with_frequency(tmp_path):
    P, CX = load_params(), load_driver("bms_10c262")
    fld = puck_field(P, CX, np.array([2000.0, 8000.0]), params_text=DEFAULT_PATH.read_text(),
                     per_wavelength=5.0, cache_dir=tmp_path)
    o = fld.obs; hf = fld.names.index("coax_hf")
    ang = o.h_angles[(o.h_angles >= 0) & (o.h_angles <= 90)]
    lv = 20 * np.log10(np.abs(fld.p[:, hf, o.h_idx[(o.h_angles >= 0) & (o.h_angles <= 90)]]))
    lv -= lv[:, :1]
    bw = [ang[np.argmax(row < -6.0)] for row in lv]   # first angle 6 dB down
    assert bw[1] < bw[0]
    assert fld.meta["axis_1m"].shape == (2, 2)
```

- [ ] **Step 2: Run it to verify it fails**

Run: `uv run pytest tests/test_bem_puck.py -m bem -v`
Expected: FAIL with `ImportError: cannot import name 'puck_field'`.

- [ ] **Step 3: Implement `puck_field`** with the same banding as `cabinet_field`. Bands above 4 kHz use `h_coax = h` (the coax region needs no extra refinement), and `h_max = 0.03`.

- [ ] **Step 4: Run it to verify it passes**

Run: `uv run pytest tests/test_bem_puck.py -m bem -v`
Expected: PASS (a few minutes).

- [ ] **Step 5: Commit**

```bash
git add src/speaker/bem/run.py tests/test_bem_puck.py
git commit -m "feat: coax-only puck BEM model to 12 kHz for directivity above the cabinet limit"
```

---

### Task 7: Source strengths and way responses

**Files:**
- Modify: `src/speaker/bass.py` (rename `_velocity` to public `cone_velocity`, and update its callers)
- Create: `src/speaker/system.py`
- Test: `tests/test_system.py` (no BEM: synthetic `SourceField`s), and one `bem` test

**Interfaces:**
- Consumes: `SourceField`, `observation_set` (Task 4); `load_curve`, `interp_db`, `minimum_phase`, `smooth` (Task 5); `cardioid.simulate`, `layout.cardioid_params`, `bass.cone_velocity`; `analytic.baffled_piston_on_axis` (Task 1).
- Produces:
  - `Ways`, a dataclass with: `f`, `obs`, `ways: dict[str, np.ndarray (nf, nd) complex]`, and `meta: dict`. The ways are `"bass"`, `"coax_lf"` (cone plus both vents) and `"coax_hf"`, each pressure in Pa at `obs.r` for 2.83 V into that way, unfiltered. `meta` holds `"lf_join_mismatch_db"` (lumped against datasheet-derived level at 400 Hz) and `"trusted_limit"`.
  - `coax_section(ways) -> np.ndarray` returns `ways["coax_lf"]`, named for the cardioid gate.
  - `way_responses(p, coax, bass, layout, cab: SourceField, puck: SourceField | None, f: np.ndarray, *, trusted_limit: float) -> Ways`.
    - `f` is the output grid. Each field is interpolated onto it in log f, using dB and unwrapped phase per direction.
    - Raise `ValueError` if `cab` doesn't cover `[f[0], min(f[-1], trusted_limit)]`.
    - Raise `ValueError` if a given `puck` doesn't cover `[trusted_limit/2, min(f[-1], 12 kHz)]`.
  - Constants: `LF_JOIN = (300.0, 500.0)`, `PUCK_HOLD_F = 12000.0` and `REF_SMOOTH = 1/3`.

**Source strengths** *(the decisions; spec §6.1)*:

**Coax LF volume velocity U_lf(f):**
- Below 300 Hz: the lumped `cardioid.simulate(coax, cardioid_params(...), f, 2.83).u_cone`. This includes the chamber.
- Above 500 Hz: from the datasheet. `|U| = 10^(LF_SPL/20)·20e-6·√2 / |baffled_piston_on_axis(f, diaphragm_d/2, 1.0)|`, with minimum phase.
- Between 300 and 500 Hz: blend the dB values and the unwrapped phase linearly in log f.
- Vents: `U_vent = u_vent/u_cone · U_lf`, split equally between `vent_left` and `vent_right`.

**Coax HF:** `|U_hf| = 10^(HF_SPL/20)·20e-6·√2 / |smooth(puck.meta["axis_1m"][:, hf], REF_SMOOTH)|`, with minimum phase.
- If no puck is available, use a baffled piston of `hf_exit_d` as the reference and set `meta["hf_reference"] = "piston"`.
- Below the puck's lowest frequency (500 Hz), the piston reference is used too. The HF way is more than 20 dB down there after the crossover.

**Bass:** each driver's volume velocity is `Sd · bass.cone_velocity(bass, layout.bass_net/2, f, 2.83)`. The two drivers are wired in parallel at 2.83 V.

**Fields:**
- **Up to `trusted_limit/2`:** use the cabinet field `cab` only.
- **From `trusted_limit/2` to `trusted_limit`:** blend `cab` into `puck`, for `coax_lf` and `coax_hf` only. Blend the dB values and the unwrapped phase linearly in log f, per direction.
- **Above `trusted_limit`:** use `puck` only.
- **Above `PUCK_HOLD_F`:** hold each direction's level relative to on-axis at its 12 kHz value. The on-axis level comes from the datasheet, and the phase is the minimum phase of the on-axis level.
- **Bass and vents** are set to 0 above `trusted_limit`. After the crossover they are more than 40 dB down there.

- [ ] **Step 1: Write the failing tests**

```python
import numpy as np, pytest
from speaker.bem import analytic
from speaker.bem.run import SourceField, observation_set
from speaker.cardioid import front_back, simulate
from speaker.drivers import load_driver
from speaker.layout import cardioid_params, compute_layout
from speaker.params import load_params
from speaker.system import way_responses, coax_section, LF_JOIN
from speaker.curves import load_curve, interp_db

P, CX, B = load_params(), load_driver("bms_10c262"), load_driver("dayton_rss265ho_44")
L = compute_layout(P, CX, B)
F = np.geomspace(100, 16000, 120)
NAMES = ("coax_lf", "coax_hf", "vent_left", "vent_right", "bass_lower", "bass_upper")

def _piston_field(f, r=2.0, scale=1.0):
    o = observation_set(r)
    pp = np.array([analytic.baffled_piston_on_axis(x, CX.diaphragm_d / 2, r) for x in f]) * scale
    p = np.broadcast_to(pp[:, None, None], (len(f), len(NAMES), len(o.dirs))).copy()
    p[:, 2:4, :] = 0.0   # vents silent, so coax_lf is the cone alone
    return SourceField(f, NAMES, o, p, {"axis_1m": np.stack([np.array([analytic.baffled_piston_on_axis(x, CX.hf_exit_d / 2, 1.0) for x in f])] * 2, 1)})

def test_deembedding_reproduces_datasheet_at_2m():
    fld = _piston_field(F)
    w = way_responses(P, CX, B, L, fld, fld, F, trusted_limit=3000.0)
    m = (F > 600) & (F < 2500)
    got = 20 * np.log10(np.abs(w.ways["coax_lf"][m, w.obs.h_idx[w.obs.h_angles == 0][0]]) / 20e-6 / np.sqrt(2))
    # Vents are silent in this field: the cone must match the datasheet minus 6 dB (1 m -> 2 m) within 1 dB.
    want = interp_db(load_curve("bms_10c262_lf_spl"), F[m]) - 20 * np.log10(2)
    assert np.abs(got - want).max() < 1.0

def test_blend_moves_monotonically_from_cabinet_to_puck():
    cab, puck = _piston_field(F), _piston_field(F, scale=2.0)   # puck 6 dB hotter everywhere
    w0 = way_responses(P, CX, B, L, cab, cab, F, trusted_limit=3000.0)
    w1 = way_responses(P, CX, B, L, cab, puck, F, trusted_limit=3000.0)
    i0 = w0.obs.h_idx[w0.obs.h_angles == 0][0]
    d = 20 * np.log10(np.abs(w1.ways["coax_hf"][:, i0]) / np.abs(w0.ways["coax_hf"][:, i0]))
    assert np.allclose(d[F < 1500], 0, atol=0.01) and np.allclose(d[(F > 3000) & (F < 12000)], 6.02, atol=0.05)
    m = (F >= 1500) & (F <= 3000)
    assert np.all(np.diff(d[m]) >= -1e-9)
    assert np.all(w1.ways["bass"][F > 3000] == 0)

def test_lf_join_mismatch_reported():
    fld = _piston_field(F)
    w = way_responses(P, CX, B, L, fld, fld, F, trusted_limit=3000.0)
    assert abs(w.meta["lf_join_mismatch_db"]) < 6.0

@pytest.mark.bem
def test_bem_rear_rejection_close_to_lumped_at_150hz(tmp_path):
    from speaker.bem.run import cabinet_field
    from speaker.params import DEFAULT_PATH
    f = np.array([150.0, 200.0])
    cab = cabinet_field(P, CX, B, f, params_text=DEFAULT_PATH.read_text(), h_max=0.06, cache_dir=tmp_path)
    w = way_responses(P, CX, B, L, cab, None, f, trusted_limit=3000.0)
    o = w.obs; s = coax_section(w)
    bem_fb = 20 * np.log10(np.abs(s[:, o.h_idx[o.h_angles == 0][0]]) / np.abs(s[:, o.h_idx[o.h_angles == -180][0]]))
    cp = cardioid_params(P, L, CX)
    lump = front_back(simulate(CX, cp, f, 2.83), cp)
    assert np.all(np.abs(bem_fb - lump) < 3.0)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_system.py -m "bem or not bem" -v`
Expected: FAIL with `ModuleNotFoundError: speaker.system`.

- [ ] **Step 3: Implement** the `cone_velocity` rename, then `system.py`, per the Source strengths section.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_system.py -m "bem or not bem" -v && uv run pytest`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/bass.py src/speaker/report.py src/speaker/system.py tests/test_system.py
git commit -m "feat: source strengths from lumped models and traced datasheet, joined 3D/puck way responses"
```

---

### Task 8: Fixed crossover and re-centred spinorama

**Files:**
- Create: `src/speaker/crossover.py`, `src/speaker/spinorama.py`
- Test: `tests/test_crossover.py`, `tests/test_spinorama.py`

**Interfaces:**
- Consumes: `Ways` (Task 7) and `Observation` (Task 4).
- Produces:
  - `lr4(f, fc, kind: Literal["low", "high"]) -> np.ndarray` complex (the squared 2nd-order Butterworth).
  - `FixedCrossover(bass_coax: float = 170.0, lf_hf: float = 1400.0)`.
  - `apply_crossover(w: Ways, xo: FixedCrossover) -> tuple[np.ndarray (nf, nd), dict[str, float]]` returns the total and the per-way gains in dB.
    - Bass gets LR4 low at bass_coax.
    - coax_lf gets LR4 high at bass_coax times LR4 low at lf_hf.
    - coax_hf gets LR4 high at lf_hf.
  - Level matching *(my call)*, at 23° horizontal, using the dB mean over a band:
    - coax_lf is the reference, measured over 300–1000 Hz with gain 0;
    - coax_hf is matched over 2–8 kHz;
    - bass is matched over 80–150 Hz.
  - `Spin`, a dataclass with: `f`, `on_axis`, `at23`, `lw`, `er`, `sp`, `erdi`, `spdi` (all dB, `np.ndarray (nf,)`), and `er_parts: dict[str, np.ndarray]`.
  - `spinorama(total: np.ndarray, obs: Observation) -> Spin`.
    - **LW:** the re-centred window from Global Constraints.
    - **ER parts** (CTA-2034):
      - floor: V −20, −30, −40;
      - ceiling: V 40, 50, 60;
      - front: H 0, ±10, ±20, ±30;
      - side: H ±40…±80, step 10;
      - rear: H ±90, 180.
    - **ER** is the power mean of the five parts.
    - **SP** is the sphere-weighted power mean.
    - **DIs:** `erdi = lw − er`, `spdi = lw − sp`.
    - All levels are dB SPL re 20 µPa (rms).

- [ ] **Step 1: Write the failing tests**

```python
import numpy as np, pytest
from speaker.crossover import lr4, FixedCrossover, apply_crossover
from speaker.spinorama import spinorama
from speaker.bem.run import observation_set
from speaker.system import Ways

def test_lr4_sums_to_allpass():
    f = np.geomspace(20, 20000, 300)
    s = lr4(f, 1400, "low") + lr4(f, 1400, "high")
    assert np.allclose(np.abs(s), 1.0, atol=1e-9)
    assert 20 * np.log10(abs(lr4(np.array([1400.0]), 1400, "low")[0])) == pytest.approx(-6.02, abs=0.01)

def _ways(fn):
    o = observation_set(2.0); f = np.geomspace(50, 16000, 200)
    d = o.dirs
    w = {k: np.broadcast_to(fn(d)[None, :], (len(f), len(d))).astype(complex) for k in ("bass", "coax_lf", "coax_hf")}
    return Ways(f, o, w, {})

def test_level_match_and_omni_spin():
    w = _ways(lambda d: np.ones(len(d)))
    w.ways["coax_hf"] *= 10.0
    total, gains = apply_crossover(w, FixedCrossover())
    assert gains["coax_hf"] == pytest.approx(-20.0, abs=0.01) and gains["coax_lf"] == 0.0
    s = spinorama(total, w.obs)
    assert np.allclose(s.lw, s.on_axis, atol=1e-6) and np.allclose(s.sp, s.on_axis, atol=1e-6)
    assert np.allclose(s.spdi, 0.0, atol=1e-6)

def test_dipole_power_di():
    # |cos| pattern along the forward axis: DI = 10 log10(3) = 4.77 dB on axis.
    o = observation_set(2.0); p = np.abs(-o.dirs[:, 1])
    s = spinorama(np.broadcast_to(p, (5, len(p))).astype(complex), o)
    assert (s.on_axis - s.sp)[0] == pytest.approx(4.77, abs=0.2)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_crossover.py tests/test_spinorama.py -v`
Expected: FAIL with `ModuleNotFoundError`.

- [ ] **Step 3: Implement** both modules per the Interfaces.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_crossover.py tests/test_spinorama.py -v`
Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/crossover.py src/speaker/spinorama.py tests/test_crossover.py tests/test_spinorama.py
git commit -m "feat: fixed LR4 crossover with level matching and re-centred spinorama"
```

---

### Task 9: Spec §4 scores

**Files:**
- Create: `src/speaker/scores.py`
- Test: `tests/test_scores.py`

**Interfaces:**
- Consumes: `Spin` (Task 8), `Observation`, `smooth` (Task 5).
- Produces:
  - `Score(key: str, label: str, value: float, unit: str, status: Literal["pass", "fail", "warn", "info"], band: tuple[float, float], confidence: Literal["trusted", "low"], formula: str, detail: dict)`.
  - `score_all(f, total, coax_sec, obs, spin, *, trusted_limit, xo: FixedCrossover) -> list[Score]`. 
  - Each band-based score is computed separately over `band ∩ [f[0], trusted_limit]` (confidence "trusted") and over `band ∩ [trusted_limit, f[-1]]` (confidence "low").
  - A part is skipped when its intersection is narrower than 1/3 octave.
  - `Score.band` is the intersected band.

**Score keys and formulas** (spec §4; copy the formula text into `Score.formula`):

| key | Formula | Status |
|---|---|---|
| `flat23` | The 23° H response in 1/3-octave band levels over 150 Hz–16 kHz, minus their mean. Value = max \|dev\|; `detail["rms"]` holds the RMS. | `info` |
| `window_spread` | For H 0, 5, …, 30, the 1/6-octave-smoothed level minus the 23° curve. Value = worst \|dev\| over 150 Hz–16 kHz. `detail["bands_within_1db"]` lists the contiguous (lo, hi) ranges where every angle is within ±1 dB. | `pass` if value ≤ 1, else `warn` |
| `presence_lift` | Over the LW curves (H 0–30 step 5, plus the two elevated directions), each 1/6-octave smoothed: value = max over f in 2–5 kHz of max(c(f) − mean₂₋₅ₖ(c23), c(f) − mean₁₋₂ₖ(c)). | `warn` if value > 0.5 dB, else `pass` |
| `di_smooth_bass`, `di_smooth_lfhf` | For SPDI and ERDI over [fc/2, 2fc]: the max residual from a least-squares line in log f. Value = the worse of the two. | `info` |
| `rolloff_monotonic` | For each 1/3-octave band over 150 Hz–16 kHz and H 40, 50, …, 90: count the bands where the level rises by > 0.5 dB with increasing angle. Value = that count; `detail["bands"]` lists them. | `pass` if 0, else `warn` |
| `cardioid_hard` | The coax section: min over 150–300 Hz of L(0°) − L(180°). | `pass` if ≥ 10, else `fail` |
| `cardioid_aim` | The same over 300–500 Hz. | `pass` if ≥ 6, else `warn` |
| `rear_hemisphere` | The coax section: min over 150–500 Hz and over H ±120…180 of L(0°) − L(θ). | `info` |
| `whole_rear` | The total after the crossover: min over 150–300 Hz of L(0°) − L(180°). | `info` |

- [ ] **Step 1: Write the failing tests**, using synthetic inputs built with `observation_set` and simple patterns:

```python
import numpy as np, pytest
from speaker.bem.run import observation_set
from speaker.crossover import FixedCrossover
from speaker.scores import score_all
from speaker.spinorama import spinorama

O = observation_set(2.0); F = np.geomspace(50, 16000, 300)

def _run(pattern, coax_pattern=None, limit=3000.0):
    p = np.stack([pattern(O.dirs, f) for f in F]).astype(complex)
    c = p if coax_pattern is None else np.stack([coax_pattern(O.dirs, f) for f in F]).astype(complex)
    return {s.key + "/" + s.confidence: s for s in score_all(F, p, c, O, spinorama(p, O), trusted_limit=limit, xo=FixedCrossover())}

def test_omni_passes_window_and_fails_cardioid():
    s = _run(lambda d, f: np.ones(len(d)))
    assert s["window_spread/trusted"].value == pytest.approx(0, abs=1e-6) and s["window_spread/trusted"].status == "pass"
    assert s["cardioid_hard/trusted"].status == "fail" and s["cardioid_hard/trusted"].value == pytest.approx(0, abs=1e-6)
    assert s["window_spread/low"].band == (3000.0, 16000.0)

def test_ideal_cardioid_passes_hard():
    card = lambda d, f: 0.5 * (1 - d[:, 1]) + 1e-3          # forward is -y
    s = _run(lambda d, f: np.ones(len(d)), coax_pattern=card)
    assert s["cardioid_hard/trusted"].status == "pass" and s["cardioid_hard/trusted"].value > 40

def test_presence_lift_flags_on_axis_hot_2_to_5k():
    def pat(d, f):
        az = np.degrees(np.arctan2(d[:, 0], -d[:, 1]))
        hot = (2000 <= f <= 5000) * (np.abs(az) < 7) * 10 ** (2 / 20)   # 0 deg 2 dB hot in 2-5 kHz
        return np.where(hot > 0, hot, 1.0)
    s = _run(pat, limit=8000.0)
    assert s["presence_lift/trusted"].status == "warn" and s["presence_lift/trusted"].value == pytest.approx(2.0, abs=0.3)

def test_rolloff_counts_rising_bands():
    def pat(d, f):
        az = np.abs(np.degrees(np.arctan2(d[:, 0], -d[:, 1])))
        return np.where((az >= 60) & (az <= 70) & (f > 1000), 2.0, 1.0)
    assert _run(pat)["rolloff_monotonic/trusted"].value >= 1
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_scores.py -v`
Expected: FAIL with `ModuleNotFoundError: speaker.scores`.

- [ ] **Step 3: Implement `scores.py`** per the table.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_scores.py -v`
Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/scores.py tests/test_scores.py
git commit -m "feat: spec section 4 directivity and cardioid scores with confidence split"
```

---

### Task 10: `speaker directivity` report and mesh check

**Files:**
- Create: `src/speaker/directivity_report.py`
- Modify: `src/speaker/cli.py`
- Test: `tests/test_directivity_cli.py`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `run_directivity(params_path: Path | None, out_dir: Path, *, fast: bool = False, f_max_3d: float = 4000.0, per_wavelength: float = 6.0) -> Path` writes and returns `out_dir/report.html`.
    - Full run: output grid `freq_grid(50, 20000, 12)`. The cabinet runs to `f_max_3d` and the puck from 500 Hz to 12 kHz.
    - With `fast=True`: 1/6-octave points, 50 Hz–1 kHz, `h_max = 0.08`, and no puck.
  - `mesh_check(params_path, out_dir, *, f_max: float = 4000.0, fine: float = 8.0) -> Path` writes `out_dir/mesh_check.json`. Its keys are `trusted_limit_hz`, `band_max_dev_db` (a dict from 1/3-octave centre to dB), `per_wavelength: [6.0, 8.0]` and `elements`.
  - CLI: `speaker directivity [--params P] [--out DIR] [--fast] [--mesh-check]`. The default out is `out/directivity`. Exit code 2 with `error: ...` on `ValueError`, as `build` does.

**Mesh check rule:**
- Run `cabinet_field` at 6 and 8 elements per wavelength, up to f_max. Feed each through `way_responses`, `apply_crossover` and `spinorama`.
- For each 1/3-octave band, take the max |dB difference| across: at23, LW, ER, SP and the H 0…30 curves.
- `trusted_limit_hz` is the upper edge of the highest band such that it and every band below it differ by ≤ 0.5 dB.
- If the fine mesh for a band would exceed `BEM_MEMORY_BUDGET`, the check stops there. `trusted_limit_hz` is then capped at the last completed band, with `"capped_by_memory": true`.
- `run_directivity` reads `mesh_check.json` from `out_dir` if it's present. Otherwise it uses 3000 Hz and writes "unverified" next to the limit.

**Report contents** (spec §7, round-0 subset). PNGs go in `out_dir/plots/`, embedded with `<img>`:
1. **Header:** the params file path and sha256, the code commit, the date, "Round 0 — current design (baseline)", and the trusted limit (with "unverified" if no mesh check).
2. **Plots:**
   - H and V polar maps: normalised to on-axis, −180…180°, log f, with a vertical line at the trusted limit.
   - The spinorama: on-axis, at23, LW, ER, SP, with the DIs on a second axis.
   - The H 0…30 family normalised to 23°, with ±1 dB band shading and the 2–5 kHz band shaded.
   - Coax-section rear rejection, BEM and lumped, with the 10 and 6 dB lines and the 150–300 and 300–500 Hz bands. The plot is titled with "lumped vs BEM differ by > 3 dB" if they do anywhere in 150–500 Hz.
   - The whole-speaker 180° rejection, marked info.
3. **Score table:** one row per `Score`, with key, value, unit, status, band, confidence and formula. Low-confidence rows are styled grey.
4. **Confidence notes:**
   - the trusted limit;
   - the ±30% vent resistance band. Run `cardioid_hard` at R × 0.7 and R × 1.3 by rescaling the lumped vent ratio only. No new BEM run is needed, because the vent source field is resistance-independent.
   - provisional cone geometry;
   - the LF join mismatch.
5. **Crossover:** the fixed points and the per-way gains.

`scores.json` is the list of `Score` as dicts.

- [ ] **Step 1: Write the failing tests**

```python
import json
from pathlib import Path
import pytest
from speaker.cli import main

@pytest.mark.bem
def test_directivity_fast_end_to_end(tmp_path):
    assert main(["directivity", "--fast", "--out", str(tmp_path)]) == 0
    html = (tmp_path / "report.html").read_text()
    for s in ("Round 0", "unverified", "cardioid_hard", "Listening window", "Confidence"):
        assert s in html
    keys = {d["key"] for d in json.loads((tmp_path / "scores.json").read_text())}
    assert {"cardioid_hard", "cardioid_aim", "window_spread", "rear_hemisphere", "whole_rear"} <= keys
    for png in ("polar_h", "polar_v", "spinorama", "window_family", "rear_rejection"):
        assert (tmp_path / "plots" / f"{png}.png").stat().st_size > 0

def test_directivity_rejects_bad_params(tmp_path, capsys):
    text = Path("params/default.toml").read_text().replace("depth_mm = 350", "depth_mm = 351")
    pf = tmp_path / "p.toml"; pf.write_text(text)
    assert main(["directivity", "--params", str(pf), "--out", str(tmp_path / "o")]) == 2
    assert capsys.readouterr().err.startswith("error: ")
```

- [ ] **Step 2: Run them to verify they fail**

Run: `uv run pytest tests/test_directivity_cli.py -m "bem or not bem" -v`
Expected: FAIL (`invalid choice: 'directivity'`).

- [ ] **Step 3: Implement** `directivity_report.py` and the CLI subcommand.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_directivity_cli.py -m "bem or not bem" -v && uv run pytest`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add src/speaker/directivity_report.py src/speaker/cli.py tests/test_directivity_cli.py
git commit -m "feat: speaker directivity report and mesh-convergence check"
```

---

### Task 11: Round 0 run and summary

**Files:**
- Create: `docs/rounds/round-0.md`
- Modify: `docs/superpowers/specs/2026-10-26-directivity-optimisation-design.md` (§8: replace "about 3–4 kHz" with the measured trusted limit and the date)

- [ ] **Step 1: Run the full BEM test suite**

Run: `uv run pytest -m bem -v 2>&1 | tail -20`
Expected: all bem tests pass.

- [ ] **Step 2: Run the mesh check** (expect hours; run it in the background and log to a file)

Run: `uv run speaker directivity --mesh-check --out out/directivity > out/mesh_check.log 2>&1`
Expected: `out/directivity/mesh_check.json` exists with `trusted_limit_hz` > 0.

- [ ] **Step 3: Run round 0 at full resolution**

Run: `uv run speaker directivity --out out/directivity`
Expected: prints `wrote out/directivity/report.html`. The report header shows the measured trusted limit, not "unverified".

- [ ] **Step 4: View the report's plots** (read the PNGs) and check that they are physically plausible:
  - The H polar map narrows with frequency.
  - The rear rejection shows the cardioid band.
  - The spinorama SP sits below LW.

  Anything implausible goes to superpowers:systematic-debugging before Step 5.

- [ ] **Step 5: Write `docs/rounds/round-0.md`:**
  - the params commit;
  - the trusted limit and how it was measured;
  - the score table (trusted and low rows);
  - the run times per band and the element counts (from `meta`);
  - the three most important observations, in plain words for the owner;
  - the provisional assumptions in force (cone depth 65 mm, HF exit 25.4 mm, fixed crossover).

  Update the spec's §8 trusted limit.

- [ ] **Step 6: Commit**

```bash
git add docs/rounds/round-0.md docs/superpowers/specs/2026-10-26-directivity-optimisation-design.md
git commit -m "docs: round 0 directivity baseline and measured trusted 3D limit"
```
