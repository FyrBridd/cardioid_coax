import numpy as np, pytest
from speaker.bem import analytic
from speaker.bem.solver import SurfaceMesh, solve_band, estimate_bytes, _assemble
import bempp_cl.api as bempp
from bempp_cl.api.operators.boundary import sparse

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
    ref = analytic.pulsating_sphere(f, 0.15, 2.0)
    assert 20 * np.log10(abs(p) / abs(ref)) == pytest.approx(0, abs=0.2)
    assert abs(np.degrees(np.angle(p / ref))) < 2  # both e^{-iωt}

@pytest.mark.bem
def test_cap_on_sphere_polar_matches_analytic():
    m = _sphere(4, 0.15, cap_alpha=np.radians(30))
    th = np.radians([0, 30, 60, 90, 150])
    pts = np.stack([2 * np.sin(th), -2 * np.cos(th), 0 * th])
    p = solve_band(m, np.array([1000.0]), pts)[0, 0]
    ref = np.array([analytic.cap_on_sphere(1000.0, 0.15, np.radians(30), 2.0, t) for t in th])
    err = 20 * np.log10(np.abs(p) / np.abs(ref))
    assert np.all(np.abs(err[:4]) < 0.5) and abs(err[4]) < 1.0
    assert np.all(np.abs(np.degrees(np.angle(p[:4] / ref[:4]))) < 2)  # both e^{-iωt}

@pytest.mark.bem
@pytest.mark.filterwarnings("ignore:Non-empty compiler output")
@pytest.mark.parametrize("ka", [np.pi, 4.493])  # zeros of j0 (where plain CBIE fails for monopole data) and of j0'
def test_burton_miller_system_is_well_conditioned_at_interior_resonances(ka):
    m = _sphere(4, 0.15)
    grid = bempp.Grid(np.ascontiguousarray(m.vertices.T), np.ascontiguousarray(m.triangles.T, dtype=np.uint32))
    space = bempp.function_space(grid, "P", 1)
    a, _ = _assemble(space, ka / 0.15, sparse.identity(space, space, space).weak_form().to_sparse())
    cond = np.linalg.cond(a)
    assert cond < 20

def test_memory_guard_refuses_before_assembly():
    v = np.random.rand(3 * 40000, 3); t = np.arange(3 * 40000).reshape(-1, 3)
    m = SurfaceMesh(v, t, np.ones(40000, int), ("src",))
    with pytest.raises(ValueError, match="memory"):
        solve_band(m, np.array([100.0]), np.zeros((3, 1)), memory_budget=1e9)
    assert estimate_bytes(1000) == 64e6

def test_cap_series_reduces_to_pulsating_sphere_for_full_cap():
    assert abs(analytic.cap_on_sphere(500.0, 0.15, np.pi, 2.0, 0.3)) == \
        pytest.approx(abs(analytic.pulsating_sphere(500.0, 0.15, 2.0)), rel=1e-6)

def test_baffled_piston_low_ka_is_twice_pulsating_sphere_with_same_phase():
    # Low ka: both reduce to -i*rho*c*k*Q*e^{ikr}/(4 pi r) times 1 (sphere) or 2 (half-space).
    f, a, r = 20.0, 0.01, 2.0
    assert analytic.baffled_piston_on_axis(f, a, r) == \
        pytest.approx(2 * analytic.pulsating_sphere(f, a, r), rel=1e-3)
