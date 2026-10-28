import numpy as np, pytest
from speaker.bem.run import puck_field
from speaker.drivers import load_driver
from speaker.params import load_params, DEFAULT_PATH

@pytest.mark.bem
def test_hf_beam_narrows_with_frequency(tmp_path):
    P, CX = load_params(), load_driver("bms_10c262")
    fld = puck_field(P, CX, np.array([2000.0, 8000.0]),
                     per_wavelength=5.0, cache_dir=tmp_path)
    o = fld.obs; hf = fld.names.index("coax_hf")
    ang = o.h_angles[(o.h_angles >= 0) & (o.h_angles <= 90)]
    lv = 20 * np.log10(np.abs(fld.p[:, hf, o.h_idx[(o.h_angles >= 0) & (o.h_angles <= 90)]]))
    lv -= lv[:, :1]
    bw = [ang[np.argmax(row < -6.0)] for row in lv]   # first angle 6 dB down
    print("beamwidth -6 dB (deg):", dict(zip((2000, 8000), bw)))
    for b in fld.meta["bands"]:
        print(b)
    assert bw[1] < bw[0]
    assert fld.meta["axis_1m"].shape == (2, 2)


@pytest.mark.parametrize("f", [[400.0], [12500.0], [1000.0, 13000.0]])
def test_puck_field_rejects_frequencies_outside_range(f):
    with pytest.raises(ValueError, match="500 Hz"):
        puck_field(load_params(), load_driver("bms_10c262"), np.array(f), cache_dir=None)
