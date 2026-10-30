"""The passive crossover topology: which parts exist for a set of options, the netlist they form in passive
and semi-active modes, the driver loads, and each way's terminal voltage per input volt.

Passive: a bass low-pass on `in`, and a coax high-pass from `in` to node `cx`, which feeds the coax LF
(low-pass, optional notch, L-pad) and the coax HF (high-pass, L-pad, optional series RLC shunt).
Semi-active: the coax network alone, with `cx` being `in` (the DSP supplies the split, the bass runs active).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from speaker.network import Part, inductor_dcr, solve

Mode = Literal["passive", "semi_active"]


@dataclass(frozen=True)
class Topology:
    bass_order: int = 2  # 2 or 3
    zobel: bool = False
    lf_order: int = 2  # 2 or 3
    lf_notch: bool = False
    hf_shunt: bool = False
    bass_inverted: bool = False
    hf_inverted: bool = False
    bass_wiring: Literal["parallel", "series"] = "parallel"


_PREFIX_BOUNDS = {"L": (0.05e-3, 20e-3), "C": (0.5e-6, 400e-6), "R": (0.1, 50.0)}  # H, F, ohm


def value_names(topo: Topology, mode: Mode) -> list[str]:
    """Names of the fitted values, in a fixed order; semi-active has no bass and no coax high-pass."""
    names: list[str] = []
    if mode == "passive":
        names += ["L_b1", "C_b1"] + (["L_b2"] if topo.bass_order == 3 else [])
        names += ["R_bz", "C_bz"] if topo.zobel else []
        names += ["C_h1", "L_h1"]
    names += ["R_ls", "R_lp", "L_l1", "C_l1"] + (["L_l2"] if topo.lf_order == 3 else [])
    names += ["L_ln", "C_ln", "R_ln"] if topo.lf_notch else []
    names += ["C_t1", "L_t1", "C_t2", "R_ts", "R_tp"]
    names += ["R_tr", "L_tr", "C_tr"] if topo.hf_shunt else []
    return names


_ALL = value_names(Topology(bass_order=3, zobel=True, lf_order=3, lf_notch=True, hf_shunt=True), "passive")
BOUNDS: dict[str, tuple[float, float]] = {
    n: ((0.5, 100.0) if n in ("R_lp", "R_tp") else _PREFIX_BOUNDS[n[0]]) for n in _ALL}


def build(values: dict[str, float], topo: Topology, mode: Mode) -> list[Part]:
    """The netlist: unit-source node `in`, loads `bass`, `coax_lf` and `coax_hf` as "Z" parts."""
    parts: list[Part] = []

    def r(name, a, b):
        parts.append(Part(name, "R", a, b, values[name]))

    def c(name, a, b):
        parts.append(Part(name, "C", a, b, values[name]))

    def ind(name, a, b):
        v = values[name]
        parts.append(Part(name, "L", a, b, v, dcr=inductor_dcr(v, cored=name.startswith("L_b"))))

    if mode == "passive":
        ind("L_b1", "in", "b1")
        c("C_b1", "b1", "0")
        end = "b1"
        if topo.bass_order == 3:
            ind("L_b2", "b1", "b2")
            end = "b2"
        if topo.zobel:
            r("R_bz", end, "bz")
            c("C_bz", "bz", "0")
        parts.append(Part("bass", "Z", end, "0", 1.0))
        c("C_h1", "in", "cx")
        ind("L_h1", "cx", "0")
        cx = "cx"
    else:
        cx = "in"

    # Coax LF: low-pass, notch, then the L-pad next to the driver.
    ind("L_l1", cx, "l1")
    c("C_l1", "l1", "0")
    end = "l1"
    if topo.lf_order == 3:
        ind("L_l2", "l1", "l2")
        end = "l2"
    if topo.lf_notch:
        ind("L_ln", end, "ln")
        c("C_ln", end, "ln")
        r("R_ln", end, "ln")
        end = "ln"
    r("R_ls", end, "lp")
    r("R_lp", "lp", "0")
    parts.append(Part("coax_lf", "Z", "lp", "0", 1.0))

    # Coax HF: high-pass, L-pad, optional series RLC across the driver.
    c("C_t1", cx, "t1")
    ind("L_t1", "t1", "0")
    c("C_t2", "t1", "t2")
    r("R_ts", "t2", "tp")
    r("R_tp", "tp", "0")
    if topo.hf_shunt:
        r("R_tr", "tp", "tr1")
        ind("L_tr", "tr1", "tr2")
        c("C_tr", "tr2", "0")
    parts.append(Part("coax_hf", "Z", "tp", "0", 1.0))
    return parts


def driver_loads(z: dict[str, np.ndarray], topo: Topology) -> dict[str, np.ndarray]:
    """Load impedances by netlist name from `speaker.vituixcad.impedances`; two woofers in series see 4x one pair."""
    bass = z["bass_pair"] * (4.0 if topo.bass_wiring == "series" else 1.0)
    return {"bass": bass, "coax_lf": z["coax_lf"], "coax_hf": z["coax_hf"]}


def bass_scale(topo: Topology) -> float:
    """Series woofers each see half the terminal voltage."""
    return 0.5 if topo.bass_wiring == "series" else 1.0


def transfers(parts: list[Part], f: np.ndarray, loads: dict[str, np.ndarray],
              topo: Topology) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Each way's terminal voltage per input volt (polarity applied), and the impedance the amplifier sees."""
    sol = solve(parts, f, loads)
    flip = {"bass": topo.bass_inverted, "coax_hf": topo.hf_inverted}
    out = {}
    for p in parts:
        if p.kind == "Z":
            v = sol.across(p)
            out[p.name] = -v if flip.get(p.name, False) else v
    return out, sol.z_in
