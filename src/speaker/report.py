"""Build all outputs (plots, CAD, summary.md) from a parameter set."""
from __future__ import annotations

import dataclasses
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .bass import max_spl, sealed_alignment, sealed_excursion, sealed_spl  # noqa: E402
from .cad.export import export_all  # noqa: E402
from .cardioid import front_back, polar, simulate, tune_resistance  # noqa: E402
from .constants import C  # noqa: E402
from .drivers import load_driver  # noqa: E402
from .layout import cardioid_params, compute_layout  # noqa: E402
from .params import Params, load_params  # noqa: E402
from .room import axial_modes, boundary_response, sbir_null  # noqa: E402

SEAT_R = 3.0  # m, seat distance used for max-SPL numbers
PAIR_GAIN_DB = 20 * math.log10(2)  # two coherent speakers
POLAR_FREQS = (150, 200, 300, 450, 600)


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def _at(f: np.ndarray, y: np.ndarray, f0: float) -> float:
    return float(np.interp(f0, f, y))


def build(params_path: Path | None, out_dir: Path) -> Path:
    p = load_params(params_path)
    coax, bass = load_driver(p.coax), load_driver(p.bass)
    layout = compute_layout(p, coax, bass)
    out_dir = Path(out_dir)
    plots = out_dir / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    drive, room = p.drive, p.room

    # ---- Bass ----
    vb = layout.bass_net / 2
    f_bass = np.geomspace(20, 300, 300)
    al0 = sealed_alignment(bass, vb)
    al_dcr = sealed_alignment(bass, vb, series_r=drive.inductor_dcr)
    spl283 = sealed_spl(bass, vb, f_bass, 2.83, r=1.0, n_drivers=2)
    mx_pass = max_spl(bass, vb, f_bass, drive.amp_v_rms, r=SEAT_R, n_drivers=2, series_r=drive.inductor_dcr)
    mx_act = max_spl(bass, vb, f_bass, drive.active_v_rms, r=SEAT_R, n_drivers=2)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    a1.semilogx(f_bass, spl283, label="2 drivers, 2.83 V, 1 m")
    a1.set_ylabel("SPL (dB)")
    a1.set_title("Sealed bass, one speaker (half-space)")
    a2.semilogx(f_bass, mx_pass, label=f"passive max ({drive.amp_v_rms:g} V, {drive.inductor_dcr:g} ohm DCR)")
    a2.semilogx(f_bass, mx_act, label=f"active max ({drive.active_v_rms:g} V)")
    a2.set_ylabel(f"Max SPL at {SEAT_R:g} m (dB)")
    a2.set_xlabel("Frequency (Hz)")
    for a in (a1, a2):
        a.grid(True, which="both", alpha=0.3)
        a.legend()
    _save(fig, plots / "bass_spl.png")

    # ---- In-room bass ----
    depth_x = room.back_gap + p.cabinet.depth
    f_room = np.linspace(20, 300, 561)
    pw = 10 ** (sealed_spl(bass, vb, f_room, 2.83, r=1.0) / 20)
    tot = np.zeros(len(f_room), dtype=complex)
    tot_free = np.zeros(len(f_room), dtype=complex)
    for y in room.speaker_y:
        for _, z in layout.bass_centres:
            src = (depth_x, y, z)
            tot += pw * boundary_response(src, room.seat, f_room, room)
            tot_free += pw * boundary_response(src, room.seat, f_room, room, surfaces=[])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(f_room, 20 * np.log10(np.abs(tot)), label="in room (all surfaces)")
    ax.plot(f_room, 20 * np.log10(np.abs(tot_free)), "--", label="free field")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Level at seat (dB, relative)")
    ax.set_title("Bass drivers of both speakers at the seat")
    ax.grid(True, alpha=0.3)
    ax.legend()
    _save(fig, plots / "bass_in_room.png")

    # ---- Cardioid ----
    cp0 = cardioid_params(p, layout, coax)
    r_tuned = tune_resistance(coax, cp0)
    cp = dataclasses.replace(cp0, vent_resistance=r_tuned)
    f_pol = np.array(POLAR_FREQS, dtype=float)
    ang = np.linspace(-180, 180, 361)
    pol = polar(simulate(coax, cp, f_pol, drive.amp_v_rms), cp, ang)
    fig = plt.figure(figsize=(6.5, 6.5))
    ax = fig.add_subplot(projection="polar")
    for fi, row in zip(POLAR_FREQS, pol):
        ax.plot(np.radians(ang), np.maximum(row, -40), label=f"{fi} Hz")
    ax.set_theta_zero_location("N")
    ax.set_ylim(-40, 0)
    ax.set_title(f"Horizontal polar (dB, tuned R = {r_tuned:.0f} Pa s/m)")
    ax.legend(loc="lower left", bbox_to_anchor=(1.0, 0.0))
    _save(fig, plots / "cardioid_polar.png")

    f_fb = np.geomspace(100, 1000, 300)
    fb = front_back(simulate(coax, cp, f_fb, drive.amp_v_rms), cp)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.semilogx(f_fb, fb)
    ax.axhline(10, color="grey", ls=":", label="10 dB")
    ax.axvspan(150, 300, color="C1", alpha=0.12, label="150-300 Hz target band")
    ax.axvspan(400, 1000, color="grey", alpha=0.12, label="needs measurement")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Front/back (dB)")
    ax.set_title("Front/back ratio (tuned resistance)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    _save(fig, plots / "front_back.png")

    f_band = np.geomspace(150, 300, 100)
    fb_min = float(np.min(front_back(simulate(coax, cp, f_band, drive.amp_v_rms), cp)))
    exc150 = float(simulate(coax, cp, np.array([150.0]), drive.amp_v_rms).excursion[0])

    # ---- Room numbers ----
    sbir = sbir_null(depth_x)
    wall_d = room.length_x - room.seat[0]
    seat_wall = sbir_null(wall_d)
    cx, cz = layout.coax_centre[0], layout.coax_centre[1]
    horiz = math.hypot(room.seat[0] - depth_x, room.seat[1] - room.speaker_y[0])
    dz_img = cz + room.seat[2]
    dz_dir = room.seat[2] - cz
    delta = math.hypot(horiz, dz_img) - math.hypot(horiz, dz_dir)
    floor_null = C / (2 * delta)
    modes = axial_modes(room)

    # ---- CAD ----
    files = export_all(layout, p, coax, bass, out_dir)

    # ---- Summary ----
    def pair(x):
        return x + PAIR_GAIN_DB

    def spl_line(label, y):
        s = _at(f_bass, y, 40.0), _at(f_bass, y, 80.0)
        return (f"- Max SPL at {SEAT_R:g} m, {label}: {s[0]:.1f} dB @ 40 Hz, {s[1]:.1f} dB @ 80 Hz "
                f"(one speaker); {pair(s[0]):.1f} dB @ 40 Hz, {pair(s[1]):.1f} dB @ 80 Hz (pair)")

    h_mm = layout.height * 1000
    lines = [
        "# Speaker design summary",
        "",
        f"Drivers measured: {'yes' if p.drivers_measured else 'NO - results are provisional'}",
        "",
        "## Bass",
        "",
        f"- Net volume per driver: {vb * 1000:.1f} L ({layout.bass_net * 1000:.1f} L both drivers, "
        f"gross {layout.bass_gross * 1000:.1f} L)",
        f"- Without inductor DCR: Fc {al0.fc:.1f} Hz, Qtc {al0.qtc:.3f}, F3 {al0.f3:.1f} Hz",
        f"- With inductor DCR ({drive.inductor_dcr:g} ohm): Fc {al_dcr.fc:.1f} Hz, "
        f"Qtc {al_dcr.qtc:.3f}, F3 {al_dcr.f3:.1f} Hz",
        spl_line(f"passive ({drive.amp_v_rms:g} V rms, DCR in series)", mx_pass),
        spl_line(f"active ({drive.active_v_rms:g} V rms)", mx_act),
        "- Assumptions: half-space (2 pi) radiation, no room gain; one speaker = 2 bass drivers; "
        f"pair = one speaker + {PAIR_GAIN_DB:.2f} dB (coherent sum); max SPL limited by the voltage limit "
        "or Xmax, whichever is lower.",
        "",
        "## Cardioid",
        "",
        f"- Tuned vent specific resistance: {r_tuned:.0f} Pa s/m (parameter file: "
        f"{p.vents.specific_resistance:g} Pa s/m)",
        f"- Minimum front/back over 150-300 Hz: {fb_min:.1f} dB"
        + ("" if fb_min >= 10 else "  (BELOW the 10 dB target - design finding, parameters not tuned to hide it)"),
        f"- Coax excursion at 150 Hz at {drive.amp_v_rms:g} V rms: {exc150 * 1000:.2f} mm peak "
        f"vs Xmax {coax.xmax * 1000:.2f} mm ({exc150 / coax.xmax * 100:.0f}% of Xmax)",
        "",
        "## Room",
        "",
        f"- SBIR null (d = back gap + cabinet depth = {depth_x:.3f} m): {sbir:.0f} Hz",
        f"- Seat-wall null (d = {wall_d:.3f} m): {seat_wall:.0f} Hz",
        f"- Floor-bounce null (coax at z = {cz:.3f} m, seat z = {room.seat[2]:.3f} m, "
        f"path difference {delta:.3f} m, left speaker to seat): {floor_null:.0f} Hz",
        "- Axial modes (Hz): " + "; ".join(
            f"{k}: " + ", ".join(f"{v:.1f}" for v in vs) for k, vs in modes.items()
        ),
        "",
        "## Cabinet",
        "",
        f"- External H x W x D: {h_mm:.0f} x {p.cabinet.width * 1000:.0f} x {p.cabinet.depth * 1000:.0f} mm",
        f"- Coax chamber: gross {layout.coax_chamber_gross * 1000:.2f} L, net {layout.coax_chamber_net * 1000:.2f} L, "
        f"filler {layout.filler_volume * 1000:.2f} L",
        f"- Bass chamber: gross {layout.bass_gross * 1000:.2f} L, net {layout.bass_net * 1000:.2f} L (both drivers)",
        "- Layout warnings: " + ("none" if not layout.warnings else ""),
    ]
    lines += [f"  - {w}" for w in layout.warnings]
    lines.append("- Exported files:")
    for kind, paths in files.items():
        for fp in paths:
            lines.append(f"  - {os.path.relpath(fp, out_dir)}")
    lines += [
        "",
        "## Confidence",
        "",
        "- Sealed bass response (Fc, Qtc, F3): medium (Thiele-Small from datasheet"
        + ("" if p.drivers_measured else ", drivers not yet measured") + ")",
        "- Bass max SPL: medium (voltage and Xmax limits only, no thermal or compression model)",
        "- Cardioid front/back 150-400 Hz: medium (lumped model, point-source vents)",
        "- Cardioid above 400 Hz: needs-measurement",
        "- Coax excursion at 150 Hz: medium",
        "- SBIR, seat-wall and floor-bounce nulls: high (geometry only, first-order)",
        "- Axial modes: high (rigid-wall rectangular room)",
        "- In-room bass response: needs-measurement (rigid, fully reflective boundaries)",
        "- CAD / volumes: high",
        "",
    ]
    summary = out_dir / "summary.md"
    summary.write_text("\n".join(lines))
    return summary
