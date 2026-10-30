"""Passive crossover schematics and the three-mode link diagram."""
from __future__ import annotations

import re

import pytest
from fakes import SAMPLE

from speaker.passive import Topology, build
from speaker.schematic import draw_links, draw_network, fmt_value, link_states

FULL = Topology(bass_order=3, zobel=True, bass_comp=True, lf_order=3, lf_notch=True, lf_shelf=True, hf_shunt=True)


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
    assert fmt_value("C", 1.2e-3) == "1200 µF" and fmt_value("C", 3000e-6) == "3000 µF"   # the compensator's range


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
        parts = [p for p in parts if p.name.startswith(("L_b", "C_b", "R_bz", "R_bc")) or p.name == "bass"]
    elif subset == "coax":
        parts = [p for p in parts if not p.name.startswith(("L_b", "C_b", "R_bz", "R_bc", "bass"))]
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
    for k in ("K1", "K2", "K3", "K4", "K5"):
        assert text.count(k) >= 2  # the header and the legend
    assert "linked out" not in text


def test_link_states_pin_each_mode():
    m = link_states(Topology())
    assert list(m) == ["Passive", "Semi-active", "Fully active"]
    assert all(set(v["K"]) == {"K1", "K2", "K3", "K4", "K5"} for v in m.values())
    assert set(m["Passive"]["K"].values()) == {"closed"}
    semi = m["Semi-active"]
    assert semi["K"] == dict(K1="open", K2="open", K3="open", K4="closed", K5="closed")
    assert [t["name"] for t in semi["terminals"]] == ["COAX", "BASS"]
    active = m["Fully active"]
    assert active["K"]["K4"] == active["K"]["K5"] == active["K"]["K3"] == "open"
    assert active["K"]["K1"] == active["K"]["K2"] == "either"
    assert [t["name"] for t in active["terminals"]] == ["HF", "LF", "BASS"]


def test_k3_to_k5_are_two_pole_links():
    from speaker.schematic import LINK_POLES, LINK_TEXT
    assert LINK_POLES == {"K1": 1, "K2": 1, "K3": 2, "K4": 2, "K5": 2}
    for k in ("K3", "K4", "K5"):
        assert LINK_TEXT[k].startswith("two-pole") and "both" in LINK_TEXT[k]
    for k in ("K1", "K2"):
        assert LINK_TEXT[k].startswith("one-pole")


def test_every_active_mode_terminal_names_where_its_plus_and_minus_land():
    m = link_states(Topology())
    term = {(mode, t["name"]): (t["plus"], t["minus"]) for mode, v in m.items() for t in v["terminals"]}
    assert term[("Passive", "IN")] == ("node in", "node 0")
    assert term[("Semi-active", "COAX")] == ("node cx", "node 0 (the ground of the coax boards)")
    assert term[("Semi-active", "BASS")] == ("the woofers' + terminal", "the woofers' − terminal")
    assert term[("Fully active", "HF")] == ("the coax HF driver's + terminal", "the coax HF driver's − terminal")
    assert term[("Fully active", "LF")] == ("the coax LF driver's + terminal", "the coax LF driver's − terminal")
    assert term[("Fully active", "BASS")] == ("the woofers' + terminal", "the woofers' − terminal")
    # in the active modes the DSP sets polarity against each driver's own + terminal
    for mode in ("Semi-active", "Fully active"):
        assert "relative to each driver's own + terminal" in m[mode]["polarity"]
    inv = link_states(Topology(bass_inverted=True, hf_inverted=True))["Passive"]["polarity"]
    assert "the woofers and the coax HF inverted" in inv and "behind K3, K5" in inv


def test_link_and_terminal_sentences():
    from speaker.schematic import links_text, terminal_text
    m = link_states(Topology())
    assert links_text(m["Semi-active"]) == "K1 open, K2 open, K3 open, K4 closed, K5 closed"
    assert terminal_text(m["Semi-active"]["terminals"][0]) == (
        "COAX + to node cx, COAX − to node 0 (the ground of the coax boards)")


def test_links_svg_shows_the_states(tmp_path):
    text = _text(draw_links(Topology(), tmp_path / "links.svg").read_text())
    assert "COAX + to node cx" in text and "open" in text and "either" in text
    assert "two-pole" in text and "own + terminal" in text


def test_speakers_are_marked_plus_and_inverted_ones_say_so(tmp_path):
    parts = build(dict(SAMPLE), Topology(), "passive")
    hf = [p for p in parts if p.name.startswith(("C_t", "L_t", "R_t")) or p.name == "coax_hf"]
    text = _text(draw_network(hf, "HF", tmp_path / "a.svg").read_text())
    assert "+" in text.split() and "(inverted)" not in text
    text = _text(draw_network(hf, "HF", tmp_path / "b.svg", inverted={"coax_hf"}).read_text())
    assert "+" in text.split() and "coax hf (inverted)" in " ".join(text.split())


def test_the_coax_drawing_labels_node_cx(tmp_path):
    parts = build(dict(SAMPLE), Topology(), "passive")
    coax = [p for p in parts if not p.name.startswith(("L_b", "C_b", "R_bz", "R_bc", "bass"))]
    text = _text(draw_network(coax, "Coax", tmp_path / "c.svg").read_text())
    assert "cx" in text.split()


def test_the_compensator_is_one_three_part_leg_across_the_woofers(tmp_path):
    from speaker.schematic import LINK_TEXT, _Net
    t = Topology(zobel=True, bass_comp=True)
    parts = [p for p in build(dict(SAMPLE), t, "passive")
             if p.name.startswith(("L_b", "C_b", "R_b")) or p.name == "bass"]
    legs = [[p.name for p in c] for c in _Net(parts).legs("b1", "in")]
    assert ["R_bc", "L_bc", "C_bc"] in legs and ["R_bz", "C_bz"] in legs
    text = _flat(_text(draw_network(parts, "Bass", tmp_path / "b.svg").read_text()))
    for n in ("R_bc", "L_bc", "C_bc"):
        assert _flat(n) in text
    assert "compensator" in LINK_TEXT["K3"]


def test_the_woofers_hang_after_the_compensator_and_zobel_behind_k3(tmp_path):
    """The network's legs come first at a driver's node, then the driver behind its link (K3 for the woofers):
    the compensator and the Zobel are on the network side, so they drop out with the network in the active modes."""
    from speaker.schematic import DRIVER_LINK, _Net
    t = Topology(zobel=True, bass_comp=True)
    parts = [p for p in build(dict(SAMPLE), t, "passive")
             if p.name.startswith(("L_b", "C_b", "R_b")) or p.name == "bass"]
    legs = [[p.name for p in c] for c in _Net(parts).legs("b1", "in")]
    assert legs[-1] == ["bass"] and legs.index(["R_bz", "C_bz"]) < len(legs) - 1
    assert DRIVER_LINK == {"bass": "K3", "coax_lf": "K4", "coax_hf": "K5"}
    assert "K3" in _text(draw_network(parts, "Bass", tmp_path / "b.svg").read_text())
    coax = [p for p in build(dict(SAMPLE), Topology(hf_shunt=True), "passive") if not p.name.startswith(("L_b", "C_b"))
            and p.name != "bass"]
    hf = [[p.name for p in c] for c in _Net(coax).legs("tp", "t2")]
    assert hf[-1] == ["coax_hf"] and ["R_tr", "L_tr", "C_tr"] in hf   # the HF's RLC shunt is before K5 too
    text = _text(draw_network(coax, "Coax", tmp_path / "c.svg").read_text())
    assert "K4" in text and "K5" in text


def test_links_svg_wraps_every_line_inside_the_canvas(tmp_path):
    from speaker.schematic import LINK_TEXT
    svg = draw_links(Topology(), tmp_path / "l.svg").read_text()
    lines = [re.sub(r"<[^>]+>", "", t).strip() for t in re.findall(r"<text[^>]*>(.*?)</text>", svg, re.S)]
    assert max(len(x) for x in lines) <= 150
    flat = " ".join(lines)
    for t in LINK_TEXT.values():
        assert " ".join(t.split()[-4:]) in flat   # each link's sentence is there to its end
