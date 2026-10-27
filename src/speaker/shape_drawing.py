"""Drawings of cabinet designs side by side, at one scale: front elevation above, plan section below.

Front elevation (looking at the baffle): the outline is solid; the edge treatment's inner lines, the flare
mouth, the bass cutouts and the vent stack (on the side faces, so hidden from the front) are dashed.
Plan section at the coax centre height: front_shape.section, the flare mouth across the front, the
coax centre and the vent slots on the side faces. Lengths are in mm.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Rectangle  # noqa: E402

from .drivers import load_driver  # noqa: E402
from .front_shape import section  # noqa: E402
from .layout import compute_layout  # noqa: E402
from .params import Params  # noqa: E402

MM = 1000.0
PANEL_W_IN = 3.6  # inches of figure width for the widest design
INK, DASH = "#222222", "#1f5fa8"


def _dashed(ax, xs, ys, **kw) -> None:
    ax.plot(xs, ys, linestyle=(0, (4, 3)), color=kw.pop("color", DASH), linewidth=1.1, **kw)


def _edge_lines(ax, p: Params, height: float) -> None:
    c = p.cabinet
    W = c.width * MM
    if c.front_edge in ("roundover", "chamfer"):
        s = c.edge_size * MM
        for x in (s, W - s):
            _dashed(ax, [x, x], [0, height])
    elif c.front_edge == "facets":
        fw, base, ramp = c.facet_width * MM, c.facet_base_z * MM, c.facet_ramp * MM
        for sign in (1, -1):
            x0, x1 = (0, fw) if sign == 1 else (W, W - fw)
            # from the side edge at the base, up the ramp (a step when the ramp is 0), then straight up
            _dashed(ax, [x0, x1, x1], [base, base + ramp, height])


def _elevation(ax, p: Params, title: str) -> None:
    d = compute_layout(p, load_driver(p.coax), load_driver(p.bass))
    bass = load_driver(p.bass)
    c = p.cabinet
    W, H = c.width * MM, d.height * MM
    ax.add_patch(Rectangle((0, 0), W, H, fill=False, edgecolor=INK, linewidth=1.8))
    _edge_lines(ax, p, H)
    cx, cz = d.coax_centre[0] * MM, d.coax_centre[1] * MM
    ax.add_patch(Circle((cx, cz), c.flare_mouth_d * MM / 2, fill=False, edgecolor=DASH, linestyle=(0, (4, 3))))
    ax.plot([cx], [cz], marker="+", color=DASH, markersize=7)
    for bx, bz in d.bass_centres:
        ax.add_patch(Circle((bx * MM, bz * MM), bass.cutout_d * MM / 2, fill=False, edgecolor=DASH,
                            linestyle=(0, (4, 3))))
    if d.vent_slots:  # on the side faces: a dashed strip along each edge
        z0 = min(s[1] - s[3] / 2 for s in d.vent_slots) * MM
        z1 = max(s[1] + s[3] / 2 for s in d.vent_slots) * MM
        for x0 in (0.0, W - 10.0):
            ax.add_patch(Rectangle((x0, z0), 10.0, z1 - z0, fill=False, edgecolor="#b3261e",
                                   linestyle=(0, (4, 3)), linewidth=1.1))
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(f"{W:.0f} wide x {H:.0f} high x {c.depth * MM:.0f} deep", fontsize=8)


def _plan(ax, p: Params) -> None:
    c = p.cabinet
    d = compute_layout(p, load_driver(p.coax), load_driver(p.bass))
    pts = [(x * MM, y * MM) for x, y in section(c, c.coax_centre_z)]
    xs, ys = zip(*pts, pts[0])
    ax.fill(xs, ys, facecolor="#e8eef6", edgecolor=INK, linewidth=1.8)
    cx = d.coax_centre[0] * MM
    half = c.flare_mouth_d * MM / 2
    ax.plot([cx - half, cx + half], [0, 0], color=DASH, linewidth=3, solid_capstyle="butt")
    ax.plot([cx], [c.baffle_t * MM], marker="+", color=DASH, markersize=7)
    if d.vent_slots:
        y, _, w, _ = d.vent_slots[0]
        for x in (0.0, c.width * MM):
            ax.plot([x, x], [(y - w / 2) * MM, (y + w / 2) * MM], color="#b3261e", linewidth=3,
                    solid_capstyle="butt")
    ax.set_xlabel("plan at the coax centre height", fontsize=8)


def draw_designs(path: Path, designs: list[tuple[str, Params]]) -> None:
    """One column per design, all at the same scale: the front elevation above, the plan section below."""
    if not designs:
        raise ValueError("no designs to draw")
    dims = [(p.cabinet.width * MM, compute_layout(p, load_driver(p.coax), load_driver(p.bass)).height * MM,
             p.cabinet.depth * MM) for _, p in designs]
    max_w, max_h, max_d = (max(x[i] for x in dims) for i in range(3))
    pad = 0.04 * max_w
    scale = PANEL_W_IN / (max_w + 2 * pad)
    row_h = [(max_h + 2 * pad) * scale + 0.6, (max_d + 2 * pad) * scale + 0.6]
    fig, axes = plt.subplots(2, len(designs), figsize=(PANEL_W_IN * len(designs) + 0.4, sum(row_h)),
                             gridspec_kw={"height_ratios": row_h}, squeeze=False)
    for col, (title, p) in enumerate(designs):
        top, bot = axes[0][col], axes[1][col]
        _elevation(top, p, title)
        _plan(bot, p)
        for ax, ymax in ((top, max_h), (bot, max_d)):
            ax.set_xlim(-pad, max_w + pad)
            ax.set_ylim(-pad, ymax + pad)
            ax.set_aspect("equal", adjustable="box")
            ax.tick_params(labelsize=7)
            ax.grid(alpha=0.25)
        top.set_ylabel("height (mm)", fontsize=8)
        bot.set_ylabel("back from the front face (mm)", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.text(0.5, 0.008, "dashed blue: edge treatment, flare mouth, bass cutouts; dashed red: vent slots "
             "(on the side faces); solid: outline", ha="center", fontsize=8, color="#555555")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=110)
    plt.close(fig)
