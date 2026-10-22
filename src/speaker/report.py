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
from .room import ALL, axial_modes, boundary_response, sbir_null  # noqa: E402

SEAT_R = 3.0  # m, seat distance used for max-SPL numbers
PAIR_GAIN_DB = 20 * math.log10(2)  # two coherent speakers
HALF_SPACE_DB = 20 * math.log10(2)  # half-space (2 pi) over free-field (4 pi) level
N_BASS_PARALLEL = 2  # bass drivers in parallel behind the one shared inductor
POLAR_FREQS = (150, 200, 300, 450, 600)
FB_TARGET_DB = 10.0
FB_BAND = (150.0, 300.0)


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def _at(f: np.ndarray, y: np.ndarray, f0: float) -> float:
    return float(np.interp(f0, f, y))


def fb_crossing(
    f: np.ndarray, fb: np.ndarray, threshold: float = FB_TARGET_DB, f_min: float = FB_BAND[0]
) -> float | None:
    """Lowest frequency >= f_min where F/B first drops below ``threshold``.

    Linear interpolation between grid points; returns ``f_min`` if F/B is
    already below the threshold there, None if it never drops below.
    """
    f = np.asarray(f, dtype=float)
    fb = np.asarray(fb, dtype=float)
    keep = f > f_min
    ff = np.concatenate(([f_min], f[keep]))
    yy = np.concatenate(([_at(f, fb, f_min)], fb[keep]))
    below = np.nonzero(yy < threshold)[0]
    if len(below) == 0:
        return None
    i = int(below[0])
    if i == 0:
        return float(f_min)
    f0, f1, y0, y1 = ff[i - 1], ff[i], yy[i - 1], yy[i]
    return float(f0 + (y0 - threshold) / (y0 - y1) * (f1 - f0))


def in_room_bass(p: Params, layout, bass, vb: float, f: np.ndarray, surfaces=ALL) -> np.ndarray:
    """dB SPL at the seat from the bass drivers of both speakers, 2.83 V per speaker.

    Each driver is a free-field point source (half-space SPL less 6.02 dB); the
    floor and the other surfaces enter only as first-order image sources, so the
    floor is not counted twice. No room modes.
    """
    room = p.room
    x = room.back_gap + p.cabinet.depth
    src = 10 ** ((sealed_spl(bass, vb, f, 2.83, r=1.0) - HALF_SPACE_DB) / 20)
    tot = np.zeros(len(f), dtype=complex)
    for y in room.speaker_y:
        for _, z in layout.bass_centres:
            tot += src * boundary_response((x, y, z), room.seat, f, room, surfaces=surfaces)
    return 20 * np.log10(np.abs(tot))


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
    # One inductor feeds both drivers in parallel: it carries both currents, so
    # each driver sees N x DCR as source resistance.
    r_eff = N_BASS_PARALLEL * drive.inductor_dcr
    dcr_text = (f"{drive.inductor_dcr:g} ohm shared by {N_BASS_PARALLEL} drivers in parallel -> "
                f"{r_eff:g} ohm effective per driver")
    al_dcr = sealed_alignment(bass, vb, series_r=r_eff)
    spl283 = sealed_spl(bass, vb, f_bass, 2.83, r=1.0, n_drivers=2)
    mx_pass = max_spl(bass, vb, f_bass, drive.amp_v_rms, r=SEAT_R, n_drivers=2, series_r=r_eff)
    mx_act = max_spl(bass, vb, f_bass, drive.active_v_rms, r=SEAT_R, n_drivers=2)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    a1.semilogx(f_bass, spl283, label="2 drivers, 2.83 V, 1 m")
    a1.set_ylabel("SPL (dB)")
    a1.set_title("Sealed bass, one speaker (half-space)")
    a2.semilogx(f_bass, mx_pass, label=f"passive max ({drive.amp_v_rms:g} V, {r_eff:g} ohm effective DCR per driver)")
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
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(f_room, in_room_bass(p, layout, bass, vb, f_room), label="in room (all six surfaces)")
    ax.plot(f_room, in_room_bass(p, layout, bass, vb, f_room, surfaces=[]), "--",
            label="free field (no boundaries)")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("dB SPL at seat, 2.83 V per speaker")
    ax.set_title("Bass drivers of both speakers at the seat\n"
                 "first-order reflections only, no room modes")
    ax.grid(True, alpha=0.3)
    ax.legend()
    _save(fig, plots / "bass_in_room.png")

    # ---- Cardioid ----
    # Two chambers: the target net volume (needs the filler, which is NOT in the
    # CAD) and the as-built CAD chamber (gross less coax displacement). Each is
    # tuned separately.
    cp0 = cardioid_params(p, layout, coax)
    f_fb = np.geomspace(100, 1000, 300)
    f_band = np.geomspace(*FB_BAND, 100)
    cases = {}
    for key, vol in (("target", layout.coax_chamber_net), ("as_built", layout.coax_chamber_as_built)):
        q = dataclasses.replace(cp0, chamber_volume=vol)
        r = tune_resistance(coax, q)
        q = dataclasses.replace(q, vent_resistance=r)
        fb_c = front_back(simulate(coax, q, f_fb, drive.amp_v_rms), q)
        cases[key] = {
            "vol": vol,
            "r": r,
            "cp": q,
            "fb": fb_c,
            "fb_min": float(np.min(front_back(simulate(coax, q, f_band, drive.amp_v_rms), q))),
            "cross": fb_crossing(f_fb, fb_c),
            "exc150": float(simulate(coax, q, np.array([150.0]), drive.amp_v_rms).excursion[0]),
        }
    cp, r_tuned = cases["target"]["cp"], cases["target"]["r"]
    f_pol = np.array(POLAR_FREQS, dtype=float)
    ang = np.linspace(-180, 180, 361)
    pol = polar(simulate(coax, cp, f_pol, drive.amp_v_rms), cp, ang)
    fig = plt.figure(figsize=(6.5, 6.5))
    ax = fig.add_subplot(projection="polar")
    for fi, row in zip(POLAR_FREQS, pol):
        ax.plot(np.radians(ang), np.maximum(row, -40), label=f"{fi} Hz")
    ax.set_theta_zero_location("N")
    ax.set_ylim(-40, 0)
    ax.set_title(f"Horizontal polar, target {cases['target']['vol'] * 1000:.1f} L chamber\n"
                 f"(dB, tuned R = {r_tuned:.0f} Pa s/m)")
    ax.legend(loc="lower left", bbox_to_anchor=(1.0, 0.0))
    _save(fig, plots / "cardioid_polar.png")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.semilogx(f_fb, cases["target"]["fb"],
                label=f"target {cases['target']['vol'] * 1000:.1f} L net (R {cases['target']['r']:.0f} Pa s/m)")
    ax.semilogx(f_fb, cases["as_built"]["fb"], "--",
                label=f"as built {cases['as_built']['vol'] * 1000:.1f} L, no filler "
                      f"(R {cases['as_built']['r']:.0f} Pa s/m)")
    ax.axhline(10, color="grey", ls=":", label="10 dB")
    ax.axvspan(150, 300, color="C1", alpha=0.12, label="150-300 Hz target band")
    ax.axvspan(400, 1000, color="grey", alpha=0.12, label="needs measurement")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Front/back (dB)")
    ax.set_title("Front/back ratio (resistance tuned per chamber volume)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    _save(fig, plots / "front_back.png")


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

    def cardioid_lines(cs):
        fb_min, cross, exc150 = cs["fb_min"], cs["cross"], cs["exc150"]
        if cross is None:
            cross_text = f"not within {f_fb[0]:.0f}-{f_fb[-1]:.0f} Hz"
        elif cross <= FB_BAND[0]:
            cross_text = f"{cross:.0f} Hz (already below at the band start)"
        else:
            cross_text = f"{cross:.0f} Hz"
        return [
            f"  - Tuned vent specific resistance: {cs['r']:.0f} Pa s/m (parameter file: "
            f"{p.vents.specific_resistance:g} Pa s/m)",
            f"  - Minimum front/back over 150-300 Hz: {fb_min:.1f} dB"
            + ("" if fb_min >= FB_TARGET_DB
               else "  (BELOW the 10 dB target - design finding, parameters not tuned to hide it)"),
            f"  - Lowest frequency above 150 Hz where F/B first drops below 10 dB: {cross_text}",
            f"  - Coax excursion at 150 Hz at {drive.amp_v_rms:g} V rms: {exc150 * 1000:.2f} mm peak "
            f"vs Xmax {coax.xmax * 1000:.2f} mm ({exc150 / coax.xmax * 100:.0f}% of Xmax)",
        ]

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
        f"- With inductor DCR ({dcr_text}): Fc {al_dcr.fc:.1f} Hz, "
        f"Qtc {al_dcr.qtc:.3f}, F3 {al_dcr.f3:.1f} Hz",
        spl_line(f"passive ({drive.amp_v_rms:g} V rms, {r_eff:g} ohm effective DCR in series)", mx_pass),
        spl_line(f"active ({drive.active_v_rms:g} V rms)", mx_act),
        "- Assumptions: half-space (2 pi) radiation, no room gain; one speaker = 2 bass drivers; "
        f"pair = one speaker + {PAIR_GAIN_DB:.2f} dB (coherent sum); max SPL limited by the voltage limit "
        "or Xmax, whichever is lower.",
        "",
        "## Cardioid",
        "",
        f"- Target chamber (target net volume {layout.coax_chamber_net * 1000:.2f} L; "
        "needs the filler, see Cabinet):",
        *cardioid_lines(cases["target"]),
        f"- As built (CAD chamber with no filler: net {layout.coax_chamber_as_built * 1000:.2f} L = "
        "gross less coax displacement):",
        *cardioid_lines(cases["as_built"]),
        "- The polar plot shows the target chamber; front_back.png shows both.",
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
        f"- Coax chamber: gross {layout.coax_chamber_gross * 1000:.2f} L, as built (gross less coax "
        f"displacement) {layout.coax_chamber_as_built * 1000:.2f} L, target net "
        f"{layout.coax_chamber_net * 1000:.2f} L",
        f"- Coax chamber filler: {layout.filler_volume * 1000:.2f} L is NOT in the CAD. The builder must add "
        "it (e.g. glued-in blocks) to reach the target net volume; without it the chamber is the "
        "as-built volume and the cardioid result is the as-built one.",
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
        "- CAD / volumes: geometry high; net volumes medium (filler not modelled; displacement estimated)",
        "",
    ]
    summary = out_dir / "summary.md"
    summary.write_text("\n".join(lines))
    return summary
