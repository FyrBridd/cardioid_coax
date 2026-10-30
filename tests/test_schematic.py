"""Passive crossover schematics and the three-mode link diagram."""
from __future__ import annotations

import re

import pytest
from fakes import SAMPLE

from speaker.passive import Topology, build
from speaker.schematic import draw_links, draw_network, fmt_value

FULL = Topology(bass_order=3, zobel=True, lf_order=3, lf_notch=True, hf_shunt=True)


def _text(svg: str) -> str:
    """The SVG's <text> content, one element per line (a label's lines are <tspan>s, joined by spaces)."""
    return "\n".join(re.sub(r"<[^>]+>", " ", t) for t in re.findall(r"<text[^>]*>(.*?)</text>", svg, re.S))


def _flat(s: str) -> str:
    return s.replace("_", "").replace(" ", "")


def test_fmt_value():
    assert fmt_value("C", 6.8e-6) == "6.8 µF"
    assert fmt_value("L", 2.5e-3) == "2.5 mH"
    assert fmt_value("R", 3.3) == "3.3 Ω"
    assert fmt_value("L", 0.3e-3) == "0.3 mH"
    assert fmt_value("C", 100e-6) == "100 µF"


def test_network_svg_names_every_part(tmp_path):
    parts = build(dict(SAMPLE), Topology(hf_shunt=True), "passive")
    hf = [p for p in parts if p.name.startswith(("C_t", "L_t", "R_t")) or p.name == "coax_hf"]
    svg = draw_network(hf, "Coax HF", tmp_path / "hf.svg").read_text()
    text = _flat(_text(svg))
    for p in hf:
        if p.kind != "Z":
            assert _flat(p.name) in text
            assert _flat(fmt_value(p.kind, p.value)) in text
    assert "Coax HF" in _text(svg)


@pytest.mark.parametrize("subset", ["bass", "coax", "all"])
def test_full_network_draws_every_part(tmp_path, subset):
    parts = build(dict(SAMPLE), FULL, "passive")
    if subset == "bass":
        parts = [p for p in parts if p.name.startswith(("L_b", "C_b", "R_bz")) or p.name == "bass"]
    elif subset == "coax":
        parts = [p for p in parts if not p.name.startswith(("L_b", "C_b", "R_bz", "bass"))]
    text = _flat(_text(draw_network(parts, subset, tmp_path / f"{subset}.svg").read_text()))
    for p in parts:
        assert _flat(p.name) in text
        if p.kind == "L":
            assert _flat(f"DCR {p.dcr:.2f} Ω") in text


def test_semi_active_network_draws(tmp_path):
    parts = build(dict(SAMPLE), Topology(), "semi_active")
    text = _flat(_text(draw_network(parts, "Coax", tmp_path / "c.svg").read_text()))
    for p in parts:
        assert _flat(p.name) in text


def test_links_svg_lists_three_modes(tmp_path):
    text = _text(draw_links(Topology(), tmp_path / "links.svg").read_text())
    assert "Passive" in text and "Semi-active" in text and "Fully active" in text
    assert "1 terminal pair" in text and "2 terminal pairs" in text and "3 terminal pairs" in text
