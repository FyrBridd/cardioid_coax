import numpy as np

from speaker.cardioid import electrical_impedance
from speaker.curves import interp_db, load_curve
from speaker.evaluate import load_design
from speaker.layout import cardioid_params
from speaker.params import DEFAULT_PATH
from speaker.system import LF_JOIN
from speaker.vituixcad import ZMA_GRID, impedances


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
