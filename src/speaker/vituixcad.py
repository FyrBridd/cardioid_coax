"""VituixCAD exports of a full run: the polar response of every way (FRD) and its electrical impedance (ZMA).

The FRD files hold the simulated pressure at the observation distance, one file per ring angle, with the
common flight time removed; the ZMA files hold the impedance on ZMA_GRID."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from speaker.bass import sealed_impedance
from speaker.cardioid import electrical_impedance
from speaker.curves import interp_db, load_curve, minimum_phase
from speaker.constants import C
from speaker.evaluate import Design
from speaker.layout import cardioid_params
from speaker.spinorama import P_REF_PEAK
from speaker.system import COAX_CURVES, LF_JOIN, Ways, _log_ramp

ZMA_GRID = (10.0, 20000.0, 48)  # f_lo, f_hi (Hz), points per octave
FRD_SETS = ("coax_lf", "coax_hf", "bass_pair", "bass_lower", "bass_upper")
FRD_COMMENT = "*"


def _from_magnitude(f: np.ndarray, mag: np.ndarray) -> np.ndarray:
    """Complex impedance of minimum phase with magnitude mag (ohm)."""
    return mag * np.exp(1j * minimum_phase(f, 20 * np.log10(mag)))


def _lf_models(d: Design, f: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """|Z| (dB re 1 ohm) of the coax LF from the lumped chamber model and from the datasheet trace."""
    lump = 20 * np.log10(np.abs(electrical_impedance(d.coax, cardioid_params(d.p, d.layout, d.coax), f)))
    sheet = 20 * np.log10(interp_db(load_curve(f"{COAX_CURVES}_lf_imp"), f))  # the traces are in ohm
    return lump, sheet


def impedances(d: Design, f: np.ndarray) -> dict[str, np.ndarray]:
    """Complex impedance (ohm) of the coax LF and HF and of one and both woofers on the grid f.

    The coax LF is the lumped cone-in-chamber model below LF_JOIN and the datasheet trace above,
    blended in dB across it; the HF and the phase of the joined LF are the datasheet's minimum phase."""
    f = np.asarray(f, dtype=float)
    lump, sheet = _lf_models(d, f)
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


_ZMA_NOTES = {
    "coax_lf": "the coax LF (BMS 10C262 low-frequency way) in its cardioid chamber. |Z|: lumped chamber model below "
               "{lo:g} Hz, BMS datasheet trace above {hi:g} Hz, blended in dB between, where the two differ by up to "
               "{gap:.1f} dB (mean {mean:.1f} dB). Phase: minimum phase from |Z|.",
    "coax_hf": "the coax HF (BMS 10C262 high-frequency way). |Z|: BMS datasheet trace, held flat beyond the trace's "
               "ends. Phase: minimum phase from |Z|.",
    "bass_driver": "one bass driver in the sealed box (lumped Thiele-Small model). Voice-coil inductance is a pure "
                   "Le with no semi-inductance, so |Z| is overstated above about 2 kHz.",
    "bass_pair": "the two bass drivers in parallel (half one driver's impedance), in the sealed box (lumped "
                 "Thiele-Small model). Voice-coil inductance is a pure Le with no semi-inductance, so |Z| is "
                 "overstated above about 2 kHz.",
}

_README = """\
VituixCAD export of the simulated cardioid coax speaker.
SIMULATED, not measured. Design file {name}, sha256 {sha}.

Importing the polar responses (.frd)
Each set is one driver: coax_lf, coax_hf, bass_pair (both woofers in parallel),
and bass_lower and bass_upper (each woofer alone). Load the files of a set as
one driver's measurements: the keyword "hor" marks the horizontal plane and
"ver" the vertical plane, and the angle multiplier is 1. The file name carries
the angle in degrees, for example "coax_lf hor -30.frd".

Set every driver's X, Y and Z to 0. The positions are already in the phase.

The data are at {r:g} m in free space, so do not add VituixCAD's own baffle or
diffraction simulation. The cabinet's diffraction is already in the data.

Level is SPL in dB re 20 uPa at {r:g} m for 2.83 V RMS into the way, unfiltered.
Phase is in degrees with the common {r:g} m flight time removed (e^+jwt
convention, so a delay shows as falling phase). The rotation centre is the
coax centre on the front face.

Angles
Horizontal angles run from -175 to 180 degrees in 5 degree steps; 0 is straight
ahead on the axis. Positive horizontal angles are toward +x of the cabinet,
which is the listener's right when facing the speaker (the speaker's left).
VituixCAD's own sign convention should be checked against this before use.
Vertical angles are positive above the axis. There is no "ver 0" file; it is
the same direction as "hor 0".

What each set contains
The coax LF set includes both vents. The bass sets stop at the trusted limit
({tl:g} Hz), so a bass driver has no data above it. The bass pair is the two
drivers in parallel at 2.83 V.

Confidence
The simulation is trusted below {tl_text}. Above that it has low
confidence: the cone and flare are estimated there.

Impedances (.zma)
coax_lf.zma, coax_hf.zma, bass_driver.zma and bass_pair.zma give the electrical
impedance (frequency, ohm, degrees) from {z_lo:g} to {z_hi:g} Hz. The header of
each file says how it was made and where it is limited.
"""


def _write_frd(path: Path, head: list[str], f: np.ndarray, p: np.ndarray, r: float) -> None:
    p = p * np.exp(2j * np.pi * f * r / C)  # remove the common flight time
    with np.errstate(divide="ignore"):
        db = 20 * np.log10(np.abs(p) / P_REF_PEAK)
    deg = np.degrees(np.angle(p))
    bad = ~(np.isfinite(db) & np.isfinite(deg))
    if bad.any():
        raise ValueError(f"{path.name}: non-finite value at {f[bad][0]:.4f} Hz")
    rows = "\n".join(f"{fi:.4f}\t{m:.3f}\t{q:.2f}" for fi, m, q in zip(f, db, deg))
    path.write_text("\n".join(head) + "\n" + rows + "\n")


def export_vituixcad(out_dir: Path, ways: Ways, d: Design, tl: float, tl_text: str) -> Path:
    """Write out_dir/vituixcad/: an FRD per set, plane and ring angle, the four ZMA files and README.txt."""
    out = Path(out_dir) / "vituixcad"
    out.mkdir(parents=True, exist_ok=True)
    obs, f = ways.obs, ways.f
    sha = hashlib.sha256(d.path.read_bytes()).hexdigest()[:8]
    each = ways.meta["bass_each"]
    fields = {"coax_lf": ways.ways["coax_lf"], "coax_hf": ways.ways["coax_hf"], "bass_pair": ways.ways["bass"],
              "bass_lower": each["bass_lower"], "bass_upper": each["bass_upper"]}
    planes = (("hor", obs.h_angles, obs.h_idx), ("ver", obs.v_angles, obs.v_idx))
    for name in FRD_SETS:
        keep = f <= tl if name.startswith("bass") else np.ones(len(f), bool)
        head = [
            f"{FRD_COMMENT} SIMULATED, not measured: cardioid_coax, {name} (FRD for VituixCAD)",
            f"{FRD_COMMENT} design file {d.path.name}, sha256 {sha}",
            f"{FRD_COMMENT} SPL dB re 20 uPa at {obs.r:g} m for 2.83 V RMS into the way, unfiltered, free space",
            f"{FRD_COMMENT} phase deg, common {obs.r:g} m flight time removed; rotation centre: the coax centre on the front face",
            f"{FRD_COMMENT} trusted below {tl_text}; low confidence above",
        ]
        for plane, angles, idx in planes:
            for a, i in zip(angles.astype(int), idx):
                if plane == "ver" and a == 0:
                    continue
                _write_frd(out / f"{name} {plane} {180 if a == -180 else a}.frd", head, f[keep], fields[name][keep, i], obs.r)

    fz = np.geomspace(ZMA_GRID[0], ZMA_GRID[1], round(np.log2(ZMA_GRID[1] / ZMA_GRID[0]) * ZMA_GRID[2]) + 1)
    lump, sheet = _lf_models(d, fz)
    band = (fz >= LF_JOIN[0]) & (fz <= LF_JOIN[1])
    diff = np.abs(lump[band] - sheet[band])
    gap, mean = (float(diff.max()), float(diff.mean())) if band.any() else (0.0, 0.0)
    for name, z in impedances(d, fz).items():
        note = _ZMA_NOTES[name].format(lo=LF_JOIN[0], hi=LF_JOIN[1], gap=gap, mean=mean)
        rows = "\n".join(f"{fi:.4f}\t{abs(zi):.3f}\t{np.degrees(np.angle(zi)):.2f}" for fi, zi in zip(fz, z))
        head = [f"{FRD_COMMENT} SIMULATED, not measured: electrical impedance, {note}",
                f"{FRD_COMMENT} design file {d.path.name}, sha256 {sha}; columns: Hz, ohm, degrees"]
        (out / f"{name}.zma").write_text("\n".join(head) + "\n" + rows + "\n")

    (out / "README.txt").write_text(_README.format(
        name=d.path.name, sha=sha, r=obs.r, tl=tl, tl_text=tl_text, z_lo=ZMA_GRID[0], z_hi=ZMA_GRID[1]))
    return out
