"""Schematics of the passive crossover, for the person who builds the boards: one ladder per sub-network
(series parts along the top, shunt branches down to ground, drivers as speakers, every part labelled with its
name and value in engineering units, inductors with their DCR), and a table of how the three modes use the
terminals and links.

The drawing is derived from the parts' nodes, so any Topology draws: a part with a node "0" is a shunt; the
chain of series parts from the source node to a driver is the spine; parts between the same two nodes (the
notch) are drawn in parallel; a branch off the spine that ends at ground (R_bz + C_bz, R_bc + L_bc + C_bc,
R_tr + L_tr + C_tr) is one leg; where the spine forks to two drivers, the second runs on a row below the first.
"""
from __future__ import annotations

import textwrap
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
LABEL_CLEAR = 3.6  # how far the spine runs on past a node's legs before it drops to the next row (clears the labels)
NODE_LABELS = ("cx",)  # nodes named on the drawing where they are not the source


def fmt_value(kind: str, value: float) -> str:
    """Engineering units: "6.8 µF", "2.5 mH", "3.3 Ω" (kind "C", "L" or "R"; "Z" is a driver and has none)."""
    if kind == "C":
        v, unit = value * 1e6, "µF"
    elif kind == "L":
        v, unit = value * 1e3, "mH"
    else:
        v, unit = value, "Ω"
    return f"{v:.4g} {unit}" if v >= 1000 else f"{v:.3g} {unit}"  # "1200 µF", not "1.2e+03 µF"


def _label(p: Part, inverted: frozenset[str] = frozenset()) -> str:
    if p.kind == "Z":
        return p.name.replace("_", " ") + ("\n(inverted)" if p.name in inverted else "")
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


def _draw_leg(d: schemdraw.Drawing, x: float, y: float, chain: list[Part], inverted: frozenset[str]) -> None:
    for p in chain:
        if p.kind == "Z":
            d += elm.Line().at((x, y)).down(0.6)
            d += elm.Speaker().at((x, y - 0.6))
            d += elm.Label().at((x + 0.9, y - 0.85)).label(_label(p, inverted), halign="left")
            # the driver's terminals: + on the network side unless the network wires it inverted
            top, bottom = ("−", "+") if p.name in inverted else ("+", "−")
            d += elm.Label().at((x - 0.3, y - 0.3)).label(top, fontsize=14)
            d += elm.Label().at((x - 0.3, y - 1.75)).label(bottom, fontsize=14)
            y -= 1.1
            d += elm.Line().at((x, y)).down(0.8)
            y -= 0.8
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


def _draw_row(d: schemdraw.Drawing, net: _Net, n: str, came: str | None, x: float, y: float,
              inverted: frozenset[str]) -> None:
    """The spine from node n rightwards along row y, with its legs and forks."""
    while True:
        legs = net.legs(n, came)
        if came is not None and n in NODE_LABELS:
            d += elm.Dot().at((x, y))
            d += elm.Label().at((x - 0.35, y + 0.45)).label(n, fontsize=13, color="blue")
        for j, chain in enumerate(legs):
            if j:
                d += elm.Line().at((x + (j - 1) * STEP, y)).to((x + j * STEP, y))
            d += elm.Dot().at((x + j * STEP, y))
            _draw_leg(d, x + j * STEP, y, chain, inverted)
        x += max(len(legs) - 1, 0) * STEP
        forks = net.forks(n, came)
        if not forks:
            return
        first_node, first_grp = forks[0]
        drop = net.row_gap(n, came, first_node)
        if legs and forks[1:]:  # the drop to the next row clears this node's legs and their labels
            d += elm.Line().at((x, y)).right(LABEL_CLEAR)
            x += LABEL_CLEAR
        for m, grp in forks[1:]:
            d += elm.Dot().at((x, y))
            d += elm.Line().at((x, y)).to((x, y - drop))
            x2 = _draw_group(d, x, y - drop, grp)
            _draw_row(d, net, m, n, x2, y - drop, inverted)
            drop += net.extent(m, n) + ROW_GAP
        d += elm.Line().at((x, y)).right(0.4)
        x = _draw_group(d, x + 0.4, y, first_grp)
        n, came = first_node, n


def draw_network(parts: list[Part], title: str, path: Path, inverted: set[str] | frozenset[str] = frozenset()) -> Path:
    """An SVG ladder of one sub-network; every part labelled with its name and value (inductors with DCR), every
    driver with its + and − terminals, those in inverted (by part name) wired the other way round and labelled so,
    and the nodes of NODE_LABELS named."""
    net = _Net(parts)
    d = schemdraw.Drawing(show=False)
    d.config(margin=0.8)  # room for the driver labels, whose width the bounding box underestimates
    src = net.source()
    d += elm.Dot(open=True).at((0, 0)).label(f"{src}", loc="left")
    d += elm.Line().at((0, 0)).right(0.6)
    _draw_row(d, net, src, None, 0.6, 0.0, frozenset(inverted))
    d += elm.Label().at((0.0, 2.0 + PARALLEL_GAP)).label(title, halign="left", fontsize=16)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    d.save(str(path))
    return path


LINK_POLES = {"K1": 1, "K2": 1, "K3": 2, "K4": 2, "K5": 2}
LINK_TEXT = {
    "K1": "one-pole, between IN + (node in) and C_h1, the input of the coax high-pass",
    "K2": "one-pole, in series with L_h1 (between L_h1 and node 0)",
    "K3": "two-pole, between the bass network's output (its last node and node 0) and the woofers' + and − "
          "terminals: it breaks both the hot and the return (the Zobel and the impedance-peak compensator stay on "
          "the network side)",
    "K4": "two-pole, between the coax LF network's output (node lp and node 0, after R_lp) and the coax LF "
          "driver's + and − terminals: it breaks both the hot and the return",
    "K5": "two-pole, between the coax HF network's output (node tp and node 0, after R_tp and the RLC shunt) and "
          "the coax HF driver's + and − terminals: it breaks both the hot and the return",
}
_DSP_POLARITY = "polarity is set in the DSP, relative to each driver's own + terminal"


def _driver_terminal(name: str, driver: str) -> dict[str, str]:
    return {"name": name, "plus": f"{driver} + terminal", "minus": f"{driver} − terminal"}


def link_states(topo: Topology) -> dict[str, dict]:
    """Per mode: the terminal pairs (each terminal's name, where its + and its − land, and what it feeds), the
    state of links K1-K5, and where the drivers' polarity is set. K3-K5 are two-pole (LINK_POLES), so in the active
    modes each terminal pair reaches its driver alone, with no return shared with another."""
    wiring = f"the bare woofer pair (wired in {topo.bass_wiring}), no network"
    inv = [n for n, i in (("the woofers", topo.bass_inverted), ("the coax HF", topo.hf_inverted)) if i]
    passive_pol = ("the network wires " + " and ".join(inv) + " inverted (driver + to node 0), in the wiring between "
                   "the network and the driver behind " + ", ".join(k for k, i in (("K3", topo.bass_inverted),
                                                                                   ("K5", topo.hf_inverted)) if i)
                   if inv else "every driver is wired normally (driver + to the network output)")
    hf = "inverted" if topo.hf_inverted else "normal"
    return {
        "Passive": {"pairs": "1 terminal pair",
                    "terminals": [{"name": "IN", "plus": "node in", "minus": "node 0",
                                   "feeds": "the whole network (one amplifier channel)"}],
                    "K": dict(K1="closed", K2="closed", K3="closed", K4="closed", K5="closed"),
                    "polarity": passive_pol},
        "Semi-active": {"pairs": "2 terminal pairs",
                        "terminals": [{"name": "COAX", "plus": "node cx",
                                       "minus": "node 0 (the ground of the coax boards)",
                                       "feeds": "the coax LF and HF networks, with the coax high-pass (C_h1, L_h1) "
                                                "out of circuit"},
                                      _driver_terminal("BASS", "the woofers'") | {"feeds": wiring}],
                        "K": dict(K1="open", K2="open", K3="open", K4="closed", K5="closed"),
                        "polarity": f"{_DSP_POLARITY}: BASS + on the woofers' + terminal; the coax boards keep "
                                    f"their own wiring of the drivers through K4 and K5 (coax LF normal, coax HF "
                                    f"{hf}). "
                                    "The passive-mode wiring of the woofers is behind K3, which is open"},
        "Fully active": {"pairs": "3 terminal pairs",
                         "terminals": [_driver_terminal("HF", "the coax HF driver's") | {"feeds": "the coax HF driver"},
                                       _driver_terminal("LF", "the coax LF driver's") | {"feeds": "the coax LF driver"},
                                       _driver_terminal("BASS", "the woofers'") | {"feeds": wiring}],
                         "K": dict(K1="either", K2="either", K3="open", K4="open", K5="open"),
                         "polarity": f"{_DSP_POLARITY}; K3, K4 and K5 open take the networks and their shared "
                                     "ground (node 0) off every driver"},
    }


def links_text(mode: dict) -> str:
    """"K1 open, K2 open, ..." for one mode of link_states."""
    return ", ".join(f"{k} {v}" for k, v in mode["K"].items())


def terminal_text(t: dict) -> str:
    """"COAX + to node cx, COAX − to node 0 (...)" for one terminal of link_states."""
    return f"{t['name']} + to {t['plus']}, {t['name']} − to {t['minus']}"


def draw_links(topo: Topology, path: Path) -> Path:
    """An SVG table of the three modes: the terminal pairs, where each terminal's + and − land, links K1-K5 open or
    closed, and where the drivers' polarity is set; then what each link breaks (K3-K5 two-pole)."""
    modes = link_states(topo)
    ks = list(LINK_TEXT)
    xs = (0.02, 0.14, 0.62, 0.69, 0.76, 0.83, 0.90, 0.97)
    rows = [max(len(m["terminals"]), 1) * 2 + 2 for m in modes.values()]  # text lines per row
    line = 0.024
    top = 0.86
    fig = Figure(figsize=(16, 10))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(xs[0], 0.965, "Terminals and links by mode", fontsize=16, fontweight="bold", va="center")
    for x, h in zip(xs, ("Mode", "Terminals (+ and − of each pair)", *(f"{k}\n{LINK_POLES[k]}-pole" for k in ks))):
        ax.text(x, 0.915, h, fontsize=12, fontweight="bold", va="center", ha="left" if x < 0.6 else "center")
    ax.plot([0.02, 0.98], [top + 0.01, top + 0.01], color="black", lw=1.2)
    y0 = top
    for (mode, m), n in zip(modes.items(), rows):
        y = y0 - n * line / 2
        ax.text(xs[0], y, f"{mode}\n{m['pairs']}", fontsize=12, fontweight="bold", va="center", linespacing=1.5)
        lines = []
        for t in m["terminals"]:
            lines += [terminal_text(t), f"    feeds {t['feeds']}"]
        lines.append(_wrap(m["polarity"][0].upper() + m["polarity"][1:] + ".", 95))
        ax.text(xs[1], y, "\n".join(lines), fontsize=10, va="center", linespacing=1.45)
        for x, k in zip(xs[2:], ks):
            ax.text(x, y, m["K"][k], fontsize=11, va="center", ha="center")
        y0 -= n * line + 0.02
        ax.plot([0.02, 0.98], [y0 + 0.01, y0 + 0.01], color="0.7", lw=0.8)
    ax.text(0.02, y0 - 0.02, "Links (removable):", fontsize=11, fontweight="bold", va="center")
    for j, (k, t) in enumerate(LINK_TEXT.items()):
        ax.text(0.02, y0 - 0.055 - j * 0.032, f"{k}: {t}", fontsize=10, va="center")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with matplotlib.rc_context({"svg.fonttype": "none"}):
        fig.savefig(path, format="svg")
    return path


def _wrap(s: str, width: int) -> str:
    return "\n".join(textwrap.wrap(s, width))
