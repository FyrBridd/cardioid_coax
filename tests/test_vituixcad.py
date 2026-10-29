import re

import numpy as np
import pytest

from speaker.bem.run import observation_set
from speaker.cardioid import electrical_impedance
from speaker.curves import interp_db, load_curve
from speaker.evaluate import load_design
from speaker.layout import cardioid_params
from speaker.constants import C
from speaker.params import DEFAULT_PATH
from speaker.spinorama import P_REF_PEAK
from speaker.system import LF_JOIN, Ways
from speaker.vituixcad import ZMA_GRID, export_vituixcad, impedances


def test_impedances_join_the_model_and_the_datasheet():
    d = load_design(DEFAULT_PATH)  # no BEM: load_design builds only the layout
    f = np.geomspace(*ZMA_GRID[:2], 400)
    z = impedances(d, f)
    assert set(z) == {"coax_lf", "coax_hf", "bass_driver", "bass_pair"}
    db = 20 * np.log10(np.abs(z["coax_lf"]))
    lump = 20 * np.log10(np.abs(electrical_impedance(d.coax, cardioid_params(d.p, d.layout, d.coax), f)))
    sheet = 20 * np.log10(interp_db(load_curve("bms_10c262_lf_imp"), f))
    lo, hi = f < LF_JOIN[0], f > LF_JOIN[1]
    assert lo.any() and hi.any()
    assert np.allclose(db[lo], lump[lo], atol=0.01)
    assert np.allclose(db[hi], sheet[hi], atol=0.01)
    assert np.allclose(20 * np.log10(np.abs(z["coax_hf"])), 20 * np.log10(interp_db(load_curve("bms_10c262_hf_imp"), f)), atol=0.01)
    assert np.allclose(z["bass_pair"], z["bass_driver"] / 2)


@pytest.fixture(scope="module")
def d():
    return load_design(DEFAULT_PATH)


@pytest.fixture(scope="module")
def ways():
    obs = observation_set(2.0)
    f = np.geomspace(50, 20000, 120)
    p = P_REF_PEAK * 10 ** (80 / 20) * np.exp(-2j * np.pi * f * obs.r / C)
    full = np.repeat(p[:, None], len(obs.dirs), axis=1)
    bass = full * (f <= 2851.0)[:, None]
    return Ways(f=f, obs=obs, ways={"bass": bass, "coax_lf": full, "coax_hf": full},
                meta={"bass_each": {"bass_lower": bass, "bass_upper": bass}})


def _rows(path):
    return np.array([l.split("\t") for l in path.read_text().splitlines() if not l.startswith("*")], dtype=float)


def test_frd_names_are_unique_and_parseable(tmp_path, ways, d):
    out = export_vituixcad(tmp_path, ways, d, 2851.0, "2851 Hz (test)")
    names = sorted(p.name for p in out.glob("coax_lf *.frd"))
    assert len(names) == 72 + 71 == len(set(names))
    assert "coax_lf hor 180.frd" in names and "coax_lf hor -175.frd" in names
    assert "coax_lf hor -180.frd" not in names and "coax_lf ver 0.frd" not in names
    assert all(re.fullmatch(r"coax_lf (hor|ver) -?\d+\.frd", n) for n in names)


def test_no_non_finite_values_and_bass_stops_at_the_limit(tmp_path, ways, d):
    out = export_vituixcad(tmp_path, ways, d, 2851.0, "2851 Hz (test)")
    assert len(list(out.glob("*.frd"))) == 5 * 143
    for p in out.glob("*.frd"):
        vals = _rows(p)
        assert np.isfinite(vals).all()
        if p.name.startswith("bass"):
            assert vals[-1, 0] <= 2851.0


def test_level_and_phase(tmp_path, ways, d):
    out = export_vituixcad(tmp_path, ways, d, 2851.0, "2851 Hz (test)")
    vals = _rows(out / "coax_hf hor 30.frd")
    assert np.allclose(vals[:, 1], 80.0, atol=1e-3) and np.allclose(vals[:, 2], 0.0, atol=0.01)


def test_every_file_is_marked_simulated(tmp_path, ways, d):
    out = export_vituixcad(tmp_path, ways, d, 2851.0, "2851 Hz (test)")
    files = list(out.glob("*.frd")) + list(out.glob("*.zma"))
    assert sorted(p.name for p in out.glob("*.zma")) == ["bass_driver.zma", "bass_pair.zma", "coax_hf.zma", "coax_lf.zma"]
    assert (out / "README.txt").read_text().strip()
    for p in files:
        first = p.read_text().splitlines()[0]
        assert first.startswith("*") and "SIMULATED" in first, p.name
    frd = (out / "coax_lf hor 0.frd").read_text().splitlines()
    assert all(l.startswith("*") for l in frd[:5]) and not frd[5].startswith("*")
    assert frd[4].count("trusted below 2851 Hz (test)") == 1 and "2851 Hz (2851" not in frd[4]
    for z in out.glob("*.zma"):
        rows = _rows(z)
        assert rows[0, 0] == ZMA_GRID[0] and rows[-1, 0] == ZMA_GRID[1] and np.isfinite(rows).all()


def test_a_zero_in_a_set_is_refused(tmp_path, ways, d):
    bad = ways.ways["coax_lf"].copy()
    bad[7, :] = 0.0
    w = Ways(f=ways.f, obs=ways.obs, ways={**ways.ways, "coax_lf": bad}, meta=ways.meta)
    with pytest.raises(ValueError, match=r"coax_lf .*\.frd: non-finite value at"):
        export_vituixcad(tmp_path, w, d, 2851.0, "2851 Hz (test)")
