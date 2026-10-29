"""VituixCAD exports: electrical impedance of every way (FRD/ZMA writers follow)."""
from __future__ import annotations

import numpy as np

from speaker.bass import sealed_impedance
from speaker.cardioid import electrical_impedance
from speaker.curves import interp_db, load_curve, minimum_phase
from speaker.evaluate import Design
from speaker.layout import cardioid_params
from speaker.system import COAX_CURVES, LF_JOIN, _log_ramp

ZMA_GRID = (10.0, 20000.0, 48)  # f_lo, f_hi (Hz), points per octave


def _from_magnitude(f: np.ndarray, mag: np.ndarray) -> np.ndarray:
    """Complex impedance of minimum phase with magnitude mag (ohm)."""
    return mag * np.exp(1j * minimum_phase(f, 20 * np.log10(mag)))


def impedances(d: Design, f: np.ndarray) -> dict[str, np.ndarray]:
    """Complex impedance (ohm) of the coax LF and HF and of one and both woofers on the grid f.

    The coax LF is the lumped cone-in-chamber model below LF_JOIN and the datasheet trace above,
    blended in dB across it; the HF and the phase of the joined LF are the datasheet's minimum phase."""
    f = np.asarray(f, dtype=float)
    lump = 20 * np.log10(np.abs(electrical_impedance(d.coax, cardioid_params(d.p, d.layout, d.coax), f)))
    sheet = 20 * np.log10(interp_db(load_curve(f"{COAX_CURVES}_lf_imp"), f))  # the traces are in ohm
    r = _log_ramp(f, *LF_JOIN)
    lf = 10 ** (((1 - r) * lump + r * sheet) / 20)
    hf = interp_db(load_curve(f"{COAX_CURVES}_hf_imp"), f)  # ohm
    driver = sealed_impedance(d.bass, d.layout.bass_net / 2, f)
    return {
        "coax_lf": _from_magnitude(f, lf),
        "coax_hf": _from_magnitude(f, hf),
        "bass_driver": driver,
        "bass_pair": driver / 2,
    }
