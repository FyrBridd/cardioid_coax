"""Schematics of the passive crossover, for the person who builds the boards: one ladder per sub-network
(series parts along the top, shunt branches down to ground, drivers as speakers, every part labelled with its
name and value in engineering units, inductors with their DCR), and a table of how the three modes use the
terminals and links.

The drawing is derived from the parts' nodes, so any Topology draws: a part with a node "0" is a shunt; the
chain of series parts from the source node to a driver is the spine; parts between the same two nodes (the
notch) are drawn in parallel; a branch off the spine that ends at ground (R_bz + C_bz, R_tr + L_tr + C_tr) is
one leg; where the spine forks to two drivers, the second runs on a row below the first.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import schemdraw  # noqa: E402
import schemdraw.elements as elm  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from schemdraw.backends.svg import config as svg_config  # noqa: E402

from speaker.network import Part  # noqa: E402
from speaker.passive import Topology  # noqa: E402

schemdraw.use("svg")
svg_config.text = "text"  # labels stay <text> elements, not glyph paths

GND = "0"
STEP = 4.4  # drawing units: a series part, and the spacing of legs at one node
LEG = 2.7  # the height of one part in a leg
ROW_GAP = 1.5  # clear space between a row's deepest leg and the next row
PARALLEL_GAP = 2.8  # spacing of parts in parallel


def fmt_value(kind: str, value: float) -> str:
    """Engineering units: "6.8 µF", "2.5 mH", "3.3 Ω" (kind "C", "L" or "R"; "Z" is a driver and has none)."""
    if kind == "C":
        v, unit = value * 1e6, "µF"
    elif kind == "L":
        v, unit = value * 1e3, "mH"
    else:
        v, unit = value, "Ω"
    return f"{v:.3g} {unit}"


def _label(p: Part) -> str:
    if p.kind == "Z":
        return p.name.replace("_", " ")
    s = f"{p.name}\n{fmt_value(p.kind, p.value)}"
    return s + f"\nDCR {p.dcr:.2f} Ω" if p.kind == "L" else s


def _element(p: Part):
    return {"R": elm.Resistor, "L": elm.Inductor2, "C": elm.Capacitor}[p.kind]()


class _Net:
    """The parts' graph: series edges between non-ground nodes, and what hangs on each node."""

    def __init__(self, parts: list[Part]):
        self.series = [p for p in parts if GND not in (p.a, p.b)]
        self.shunt = [p for p in parts if GND in (p.a, p.b)]
        self.kids: dict[str, list[list[Part]]] = {}
        self.nodes = list(dict.fromkeys(n for p in self.series for n in (p.a, p.b)))
        for p in self.shunt:
            n = p.b if p.a == GND else p.a
            if n not in self.nodes:
                self.nodes.append(n)
        self._memo: dict[tuple[str, str | None], bool] = {}

    def _other(self, p: Part, n: str) -> str:
        return p.b if p.a == n else p.a

    def reach(self, n: str, came: str | None) -> bool:
        """Whether a driver hangs on node n or anywhere beyond it, not going back the way it came."""
        key = (n, came)
        if key not in self._memo:
            self._memo[key] = False  # cycle guard
            self._memo[key] = (any(p.kind == "Z" for p in self.shunt if n in (p.a, p.b))
                               or any(self.reach(self._other(p, n), n) for p in self.series
                                      if n in (p.a, p.b) and self._other(p, n) != came))
        return self._memo[key]

    def source(self) -> str:
        """The node no series part ends on; "in" or "cx" when there is a choice."""
        ends = {p.b for p in self.series}
        cands = [n for n in dict.fromkeys(p.a for p in self.series) if n not in ends] or self.nodes
        for pref in ("in", "cx"):
            if pref in cands:
                return pref
        return cands[0]

    def groups(self, n: str, came: str | None) -> list[tuple[str, list[Part]]]:
        """Series edges leaving n, parts between the same two nodes grouped (the parallel notch)."""
        out: dict[str, list[Part]] = {}
        for p in self.series:
            if n in (p.a, p.b) and self._other(p, n) != came:
                out.setdefault(self._other(p, n), []).append(p)
        return list(out.items())

    def legs(self, n: str, came: str | None, seen: frozenset[str] = frozenset()) -> list[list[Part]]:
        """Chains of parts from n down to ground: shunt parts, and series edges to a node with no driver beyond
        it followed by that node's own legs."""
        seen = seen | {n}
        chains = [[p] for p in self.shunt if n in (p.a, p.b)]
        for m, grp in self.groups(n, came):
            if not self.reach(m, n) and m not in seen:
                chains += [grp[:1] + rest for rest in self.legs(m, n, seen)]
        return chains

    def forks(self, n: str, came: str | None) -> list[tuple[str, list[Part]]]:
        return [(m, g) for m, g in self.groups(n, came) if self.reach(m, n)]

    def extent(self, n: str, came: str | None) -> float:
        """How far below its row the part of the drawing from n on reaches."""
        own = max((len(c) * LEG + 1.0 for c in self.legs(n, came)), default=0.0)
        cur = own
        for i, (m, _) in enumerate(self.forks(n, came)):
            e = self.extent(m, n)
            cur = max(cur, e) if i == 0 else (cur + ROW_GAP) + e
        return cur

    def row_gap(self, n: str, came: str | None, first: str) -> float:
        """The depth reached by the first branch at n: where the next row can start."""
        return max(self.extent(first, n), max((len(c) * LEG + 1.0 for c in self.legs(n, came)), default=0.0)) + ROW_GAP


def _draw_leg(d: schemdraw.Drawing, x: float, y: float, chain: list[Part]) -> None:
    for p in chain:
        if p.kind == "Z":
            d += elm.Line().at((x, y)).down(0.6)
            d += elm.Speaker().at((x, y - 0.6))
            d += elm.Label().at((x + 0.9, y - 0.85)).label(_label(p), halign="left")
            y -= 1.1
            d += elm.Line().at((x, y)).down(0.4)
            y -= 0.4
        else:
            e = _element(p).at((x, y)).down().length(LEG).label(_label(p), loc="bottom", ofst=0.2)  # "bottom" is the right-hand side of a downward part
            d += e
            y -= LEG
    d += elm.Ground().at((x, y))


def _draw_group(d: schemdraw.Drawing, x: float, y: float, grp: list[Part]) -> float:
    """Series part(s) from (x, y) rightwards; returns the x after them."""
    for k, p in enumerate(grp):
        yy = y + k * PARALLEL_GAP
        if k:
            d += elm.Line().at((x, y)).to((x, yy))
            d += elm.Line().at((x + STEP, yy)).to((x + STEP, y))
        d += _element(p).at((x, yy)).right().length(STEP).label(_label(p), loc="top", ofst=0.15)
    return x + STEP


def _draw_row(d: schemdraw.Drawing, net: _Net, n: str, came: str | None, x: float, y: float) -> None:
    """The spine from node n rightwards along row y, with its legs and forks."""
    while True:
        legs = net.legs(n, came)
        for j, chain in enumerate(legs):
            if j:
                d += elm.Line().at((x + (j - 1) * STEP, y)).to((x + j * STEP, y))
            d += elm.Dot().at((x + j * STEP, y))
            _draw_leg(d, x + j * STEP, y, chain)
        x += max(len(legs) - 1, 0) * STEP
        forks = net.forks(n, came)
        if not forks:
            return
        first_node, first_grp = forks[0]
        drop = net.row_gap(n, came, first_node)
        if legs and forks[1:]:  # the drop to the next row clears this node's legs
            d += elm.Line().at((x, y)).right(1.2)
            x += 1.2
        for m, grp in forks[1:]:
            d += elm.Dot().at((x, y))
            d += elm.Line().at((x, y)).to((x, y - drop))
            x2 = _draw_group(d, x, y - drop, grp)
            _draw_row(d, net, m, n, x2, y - drop)
            drop += net.extent(m, n) + ROW_GAP
        d += elm.Line().at((x, y)).right(0.4)
        x = _draw_group(d, x + 0.4, y, first_grp)
        n, came = first_node, n


def draw_network(parts: list[Part], title: str, path: Path) -> Path:
    """An SVG ladder of one sub-network; every part labelled with its name and value (inductors with DCR)."""
    net = _Net(parts)
    d = schemdraw.Drawing(show=False)
    src = net.source()
    d += elm.Dot(open=True).at((0, 0)).label(f"{src}", loc="left")
    d += elm.Line().at((0, 0)).right(0.6)
    _draw_row(d, net, src, None, 0.6, 0.0)
    d += elm.Label().at((0.0, 2.0 + PARALLEL_GAP)).label(title, halign="left", fontsize=16)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    d.save(str(path))
    return path


def _links_rows(topo: Topology) -> list[tuple[str, str, str, str]]:
    wiring = "two woofers in series" if topo.bass_wiring == "series" else "two woofers in parallel"
    return [
        ("Passive", "1 terminal pair\nIN (amplifier)",
         f"No links fitted.\nBass low-pass and coax high-pass in circuit.\nBass: {wiring}.",
         "None bypassed.\nEvery network is in circuit."),
        ("Semi-active", "2 terminal pairs\nCOAX, BASS",
         "Coax high-pass links fitted:\nlink across C_h1, L_h1 left out.\nCOAX feeds the coax LF and HF networks.",
         "Bass low-pass and coax high-pass bypassed.\nBASS drives the woofers directly\n(the DSP supplies the split)."),
        ("Fully active", "3 terminal pairs\nHF, LF, BASS",
         "All series parts linked out,\nall shunt parts left out.\nHF, LF and BASS go straight to their drivers.",
         "Every network bypassed\n(bass, coax LF, coax HF)."),
    ]


def draw_links(topo: Topology, path: Path) -> Path:
    """An SVG table of the three modes: terminal pairs, which links are fitted, which networks are bypassed."""
    rows = _links_rows(topo)
    heads = ("Mode", "Terminals", "Links", "Networks")
    xs = (0.02, 0.16, 0.38, 0.70)
    fig = Figure(figsize=(14, 6.2))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(xs[0], 0.95, "Terminals and links by mode", fontsize=16, fontweight="bold", va="center")
    for x, h in zip(xs, heads):
        ax.text(x, 0.86, h, fontsize=12, fontweight="bold", va="center")
    ax.plot([0.02, 0.98], [0.82, 0.82], color="black", lw=1.2)
    for i, row in enumerate(rows):
        y = 0.66 - i * 0.27
        for x, cell in zip(xs, row):
            ax.text(x, y, cell, fontsize=11 if x != xs[0] else 13, fontweight="bold" if x == xs[0] else "normal",
                    va="center", linespacing=1.5)
        ax.plot([0.02, 0.98], [y - 0.13, y - 0.13], color="0.7", lw=0.8)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig.savefig(path, format="svg")
    return path
