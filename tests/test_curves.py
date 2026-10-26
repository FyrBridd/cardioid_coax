import numpy as np, pytest
from speaker.curves import load_curve, interp_db, minimum_phase, smooth

@pytest.mark.parametrize("stem,f,val,tol", [
    ("bms_10c262_lf_spl", 200, 95.0, 1.0), ("bms_10c262_lf_spl", 1000, 93.5, 1.0),
    ("bms_10c262_lf_spl", 1600, 89.5, 1.0), ("bms_10c262_lf_spl", 2700, 96.0, 1.0),
    ("bms_10c262_hf_spl", 2000, 111.0, 1.0), ("bms_10c262_hf_spl", 10000, 103.5, 1.5),
    ("bms_10c262_lf_imp", 69.7118, 90.49, 2.0),  # datasheet cursor reading; peak ~95 ohm at ~66 Hz
    ("bms_10c262_lf_imp", 300, 7.5, 1.0),
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
