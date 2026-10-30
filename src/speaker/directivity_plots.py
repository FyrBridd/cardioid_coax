"""Plots for the directivity report (``speaker directivity``). Levels are dB; f in Hz, log axis."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

from .curves import smooth  # noqa: E402
from .report import FB_AIM_BAND, FB_AIM_DB, FB_BAND, FB_TARGET_DB, _save  # noqa: E402
from .spinorama import Spin  # noqa: E402

MAP_RANGE = (-30.0, 5.0)  # dB re on-axis, colour scale of the polar maps
MAP_CONTOURS = (-20.0, -10.0, -6.0, -3.0)
PRESENCE = (2000.0, 5000.0)  # Hz


def _freq_axis(ax, f: np.ndarray, tl: float, tl_text: str, blend: tuple[float, float] | None = None,
               hatch_color: str = "0.35") -> None:
    """Log f axis over the grid, with the trusted 3D limit as a vertical line when it falls inside and
    the cabinet-to-puck blend range (lo, hi) hatched when given."""
    ax.set_xscale("log")
    if blend is not None:
        ax.axvspan(*blend, facecolor="none", edgecolor=hatch_color, hatch="//", lw=0,
                   label=f"cabinet-to-puck blend {blend[0]:.0f}-{blend[1]:.0f} Hz")
    if f[0] < tl < f[-1]:
        ax.axvline(tl, color="k", ls="--", lw=1.2, label=f"trusted 3D limit {tl_text}")
    ax.set_xlim(f[0], f[-1])
    ax.set_xlabel("Frequency (Hz)")
    ax.grid(True, which="both", alpha=0.3)


def polar_map(path: Path, f: np.ndarray, angles: np.ndarray, rel_db: np.ndarray, plane: str,
              tl: float, tl_text: str, blend: tuple[float, float] | None = None) -> None:
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
    _freq_axis(ax, f, tl, tl_text, blend, hatch_color="white")
    ax.grid(False, which="both")
    ax.set_yticks(np.arange(-180, 181, 30))
    ax.set_ylabel(f"{plane} angle (deg)")
    ax.set_title(f"{plane} polar map, normalised to on-axis (2 m); trusted to {tl_text}")
    if f[0] < tl < f[-1] or blend is not None:
        ax.legend(loc="lower left", fontsize=8)
    _save(fig, path)


def spinorama(path: Path, spin: Spin, f: np.ndarray, tl: float, tl_text: str,
              blend: tuple[float, float] | None = None) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for y, label, kw in ((spin.on_axis, "on-axis", {}), (spin.at23, "23° (listener)", {}),
                         (spin.lw, "listening window", {}), (spin.er, "early reflections", {}),
                         (spin.sp, "sound power", {})):
        ax.plot(f, y, label=label, **kw)
    top = float(np.max(spin.on_axis))
    ax.set_ylim(top - 45, top + 5)
    ax.set_ylabel("dB SPL at 2 m, 2.83 V")
    _freq_axis(ax, f, tl, tl_text, blend)
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


def window_family(path: Path, f: np.ndarray, rel: dict[int, np.ndarray], tl: float, tl_text: str,
                  blend: tuple[float, float] | None = None) -> None:
    """rel: H angle -> 1/6-octave smoothed level minus the smoothed 23° level."""
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axhspan(-1, 1, color="C2", alpha=0.12, label="±1 dB")
    ax.axvspan(*PRESENCE, color="C1", alpha=0.10, label="2-5 kHz")
    ax.axhline(0, color="k", lw=1.5, label="23° (reference)")
    for i, (a, y) in enumerate(sorted(rel.items())):
        ax.plot(f, y, color=plt.cm.plasma(i / max(len(rel) - 1, 1) * 0.85), label=f"H {a}°")
    _freq_axis(ax, f, tl, tl_text, blend)
    ax.set_ylabel("dB re 23° (1/6-oct smoothed)")
    ax.set_title("Horizontal 0-30° family normalised to the listener direction (23°)")
    ax.legend(loc="lower left", fontsize=8, ncol=3)
    _save(fig, path)


def rear_rejection(path: Path, f: np.ndarray, bem: np.ndarray, lumped: np.ndarray, band: tuple[np.ndarray, np.ndarray],
                   r_text: str, title: str, far: np.ndarray | None = None) -> None:
    """Coax-section L(0) - L(180): BEM, the BEM at vent resistance x0.7 and x1.3 (shaded), lumped, and
    the BEM at 50 m when given (NaN where not computed)."""
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.fill_between(f, *band, color="C0", alpha=0.2, label=f"3D (BEM), vent resistance {r_text}")
    ax.plot(f, bem, color="C0", lw=2, label="3D (BEM), nominal vent resistance")
    if far is not None:
        ax.plot(f, far, color="C0", ls=":", lw=1.5, label="3D (BEM) at 50 m (far field, reference only)")
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


def lf_reference(path: Path, f: np.ndarray, raw_db: np.ndarray, fs: np.ndarray, applied_db: np.ndarray) -> None:
    """The coax LF recess ratio 20·log10|R| at its BEM points, and the smoothed curve the de-embedding applies."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(f, raw_db, "o-", ms=3, lw=0.8, label="20·log10|R|, BEM (recess on a disc / flat piston on a disc)")
    ax.plot(fs, applied_db, lw=2, label="1/3-octave smoothed, as applied to the LF reference")
    ax.axhline(0.0, color="grey", lw=0.8)
    _freq_axis(ax, fs, np.inf, "")
    ax.set_ylabel("dB")
    ax.set_title("Coax LF reference: the cone's own recess relative to a flat piston, on axis at 1 m")
    ax.legend(loc="best", fontsize=8)
    _save(fig, path)


def crossover_23(path: Path, f: np.ndarray, curves: dict[str, np.ndarray], tl: float, tl_text: str,
                 blend: tuple[float, float] | None = None,
                 title: str = "Each way at 23° after its crossover filter, and the sum") -> None:
    """Each way's level at 23° after its crossover filter, and the sum (the "sum" curve), 1/6-octave smoothed."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    sm = {}
    for name, y in curves.items():  # a way with no data (-inf or NaN) is left out there, not smoothed
        ok = np.isfinite(y)
        sm[name] = np.full(len(f), np.nan)
        sm[name][ok] = smooth(f[ok], y[ok], 1 / 6)
    for name, y in sm.items():
        ax.plot(f, y, **({"color": "k", "lw": 2} if name == "sum" else {"lw": 1.2}), label=name)
    top = float(np.max(sm.get("sum", np.nanmax(list(sm.values()), axis=0))))
    ax.set_ylim(top - 40, top + 5)
    ax.set_ylabel("dB SPL at 2 m, 2.83 V, 23°")
    _freq_axis(ax, f, tl, tl_text, blend)
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title(title)
    _save(fig, path)


def impedance(path: Path, f: np.ndarray, curves: dict[str, np.ndarray], floor: float, floor_text: str,
              dashed: tuple[str, ...] = ()) -> None:
    """|Z| (ohm) of each labelled curve against frequency, with the minimum-impedance line at floor; the curves
    named in dashed are drawn dashed (for reference, not what the amplifier sees in the main mode)."""
    fig, ax = plt.subplots(figsize=(9, 5))
    top = 0.0
    for name, z in curves.items():
        ax.plot(f, z, lw=1.2 if name in dashed else 2, ls="--" if name in dashed else "-", label=name)
        top = max(top, float(np.nanmax(z)))
    ax.axhline(floor, color="C3", lw=1.5, ls=":", label=floor_text)
    ax.set_ylim(0, min(max(top * 1.1, 2 * floor), 60.0))
    _freq_axis(ax, f, np.inf, "")
    ax.set_ylabel("|Z| (Ω)")
    ax.set_title("Impedance the amplifier sees (modelled, nominal part values)")
    ax.legend(loc="upper left", fontsize=8)
    _save(fig, path)


def whole_rear_compare(path: Path, f: np.ndarray, rej: dict[str, np.ndarray], fc: dict[str, float]) -> None:
    """Whole-speaker L(0°) - L(180°) of several crossovers on one axis; fc: each one's bass-to-coax point (Hz)."""
    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.axvspan(*FB_BAND, color="C1", alpha=0.12, label="150-300 Hz")
    for i, (name, y) in enumerate(rej.items()):
        ax.plot(f, y, color=f"C{i}", lw=1.8, label=name)
        if name in fc:
            ax.axvline(fc[name], color=f"C{i}", ls=":", lw=1, label=f"{name}: bass-coax {fc[name]:g} Hz")
    _freq_axis(ax, f, np.inf, "")
    ax.set_ylabel("L(0°) - L(180°) (dB)")
    ax.set_title("Whole-speaker rear rejection by crossover (info only: the bass drivers radiate both ways)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8)
    _save(fig, path)


def tolerance_spread(path: Path, flat23: np.ndarray, flat23_nominal: float, flat23_p95: float, flat_limit_db: float,
                     z_min: np.ndarray, z_limit: float) -> None:
    """The part-tolerance draws: a histogram of flat23 (with the nominal value, the 95th percentile and the limit
    for it), and each draw's minimum |Z| (with the lowest allowed)."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ok = np.isfinite(flat23)
    a1.hist(flat23[ok], bins=20, color="C0", alpha=0.8)
    a1.axvline(flat23_nominal, color="k", lw=1.5, label=f"nominal values {flat23_nominal:.2f} dB")
    a1.axvline(flat23_p95, color="C1", ls="--", lw=1.5, label=f"95th percentile {flat23_p95:.2f} dB")
    a1.axvline(flat23_nominal + flat_limit_db, color="C3", ls=":", lw=1.5,
               label=f"limit for the 95th percentile (+{flat_limit_db:g} dB)")
    a1.set_xlabel("flat23 (dB, lower is flatter)")
    a1.set_ylabel("draws")
    a1.yaxis.set_major_locator(MaxNLocator(integer=True))
    a1.set_title(f"23° flatness over {ok.sum()} draws")
    a1.legend(loc="upper right", fontsize=8)
    a2.plot(np.arange(1, len(z_min) + 1), z_min, "o", ms=3, color="C0")
    a2.axhline(z_limit, color="C3", ls=":", lw=1.5, label=f"lowest allowed {z_limit:g} Ω")
    a2.set_xlabel("draw")
    a2.set_ylabel("minimum |Z| (Ω)")
    a2.set_ylim(min(z_limit, float(np.min(z_min))) - 0.3, float(np.max(z_min)) + 0.3)
    a2.set_title("Minimum impedance of each draw")
    a2.legend(loc="lower right", fontsize=8)
    for a in (a1, a2):
        a.grid(True, alpha=0.3)
    _save(fig, path)
