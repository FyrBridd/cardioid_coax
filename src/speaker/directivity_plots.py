"""Plots for the directivity report (``speaker directivity``). Levels are dB; f in Hz, log axis."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .report import FB_AIM_BAND, FB_AIM_DB, FB_BAND, FB_TARGET_DB, _save  # noqa: E402
from .spinorama import Spin  # noqa: E402

MAP_RANGE = (-30.0, 5.0)  # dB re on-axis, colour scale of the polar maps
MAP_CONTOURS = (-20.0, -10.0, -6.0, -3.0)
PRESENCE = (2000.0, 5000.0)  # Hz


def _freq_axis(ax, f: np.ndarray, tl: float, tl_text: str) -> None:
    """Log f axis over the grid, with the trusted 3D limit as a vertical line when it falls inside."""
    ax.set_xscale("log")
    if f[0] < tl < f[-1]:
        ax.axvline(tl, color="k", ls="--", lw=1.2, label=f"trusted 3D limit {tl_text}")
    ax.set_xlim(f[0], f[-1])
    ax.set_xlabel("Frequency (Hz)")
    ax.grid(True, which="both", alpha=0.3)


def polar_map(path: Path, f: np.ndarray, angles: np.ndarray, rel_db: np.ndarray, plane: str,
              tl: float, tl_text: str) -> None:
    """rel_db (nf, na): level re on-axis over angles -180..175; drawn -180..180."""
    order = np.argsort(angles)
    a, z = angles[order], rel_db[:, order]
    a, z = np.append(a, 180.0), np.concatenate([z, z[:, :1]], axis=1)  # -180 is also +180
    fig, ax = plt.subplots(figsize=(9, 5.5))
    mesh = ax.pcolormesh(f, a, np.clip(z.T, *MAP_RANGE), shading="nearest", cmap="viridis",
                         vmin=MAP_RANGE[0], vmax=MAP_RANGE[1])
    cs = ax.contour(f, a, z.T, levels=MAP_CONTOURS, colors="white", linewidths=0.7)
    ax.clabel(cs, fmt="%g", fontsize=7)
    fig.colorbar(mesh, ax=ax, label="dB re on-axis")
    _freq_axis(ax, f, tl, tl_text)
    ax.grid(False, which="both")
    ax.set_yticks(np.arange(-180, 181, 30))
    ax.set_ylabel(f"{plane} angle (deg)")
    ax.set_title(f"{plane} polar map, normalised to on-axis (2 m); trusted to {tl_text}")
    if f[0] < tl < f[-1]:
        ax.legend(loc="lower left")
    _save(fig, path)


def spinorama(path: Path, spin: Spin, f: np.ndarray, tl: float, tl_text: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for y, label, kw in ((spin.on_axis, "on-axis", {}), (spin.at23, "23° (listener)", {}),
                         (spin.lw, "listening window", {}), (spin.er, "early reflections", {}),
                         (spin.sp, "sound power", {})):
        ax.plot(f, y, label=label, **kw)
    top = float(np.max(spin.on_axis))
    ax.set_ylim(top - 45, top + 5)
    ax.set_ylabel("dB SPL at 2 m, 2.83 V")
    _freq_axis(ax, f, tl, tl_text)
    di = ax.twinx()
    di.plot(f, spin.erdi, color="C5", ls="--", label="ERDI")
    di.plot(f, spin.spdi, color="C6", ls="--", label="SPDI")
    di.set_ylim(-5, max(20.0, float(np.max([spin.erdi, spin.spdi])) + 2))
    di.set_ylabel("Directivity index (dB)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = di.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="lower left", fontsize=8, ncol=2)
    ax.set_title("Spinorama, listening window re-centred on 23°")
    _save(fig, path)


def window_family(path: Path, f: np.ndarray, rel: dict[int, np.ndarray], tl: float, tl_text: str) -> None:
    """rel: H angle -> 1/6-octave smoothed level minus the smoothed 23° level."""
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axhspan(-1, 1, color="C2", alpha=0.12, label="±1 dB")
    ax.axvspan(*PRESENCE, color="C1", alpha=0.10, label="2-5 kHz")
    ax.axhline(0, color="k", lw=1.5, label="23° (reference)")
    for i, (a, y) in enumerate(sorted(rel.items())):
        ax.plot(f, y, color=plt.cm.plasma(i / max(len(rel) - 1, 1) * 0.85), label=f"H {a}°")
    _freq_axis(ax, f, tl, tl_text)
    ax.set_ylabel("dB re 23° (1/6-oct smoothed)")
    ax.set_title("Horizontal 0-30° family normalised to the listener direction (23°)")
    ax.legend(loc="lower left", fontsize=8, ncol=3)
    _save(fig, path)


def rear_rejection(path: Path, f: np.ndarray, bem: np.ndarray, lumped: np.ndarray, band: tuple[np.ndarray, np.ndarray],
                   r_text: str, title: str) -> None:
    """Coax-section L(0) - L(180): BEM, the BEM at vent resistance x0.7 and x1.3 (shaded), lumped."""
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.fill_between(f, *band, color="C0", alpha=0.2, label=f"3D (BEM), vent resistance {r_text}")
    ax.plot(f, bem, color="C0", lw=2, label="3D (BEM), nominal vent resistance")
    ax.plot(f, lumped, color="C3", ls="--", label="lumped model (free field, no cabinet)")
    ax.axhline(FB_TARGET_DB, color="grey", ls=":", label=f"{FB_TARGET_DB:g} dB hard")
    ax.axhline(FB_AIM_DB, color="grey", ls="-.", label=f"{FB_AIM_DB:g} dB aim")
    ax.axvspan(*FB_BAND, color="C1", alpha=0.12, label="150-300 Hz hard band")
    ax.axvspan(*FB_AIM_BAND, color="C2", alpha=0.10, label="300-500 Hz aim band")
    _freq_axis(ax, f, np.inf, "")
    ax.set_ylabel("L(0°) - L(180°) (dB)")
    ax.set_title(title)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8)
    _save(fig, path)


def whole_rear(path: Path, f: np.ndarray, rej: np.ndarray, fc: float) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(f, rej, label="whole speaker after the crossover")
    ax.axvline(fc, color="grey", ls=":", label=f"bass-coax crossover {fc:g} Hz")
    ax.axvspan(*FB_BAND, color="C1", alpha=0.12, label="150-300 Hz")
    _freq_axis(ax, f, np.inf, "")
    ax.set_ylabel("L(0°) - L(180°) (dB)")
    ax.set_title("Whole-speaker rear rejection (info only: the bass drivers radiate both ways)")
    ax.legend(loc="upper right", fontsize=8)
    _save(fig, path)
