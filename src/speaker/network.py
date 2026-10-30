"""Passive networks: a nodal solver for ladders of R, L (with series DCR), C and driver loads, the
standard part values to snap to, and the series resistance of air-core and cored inductors.

e^{+jωt} convention (inductor +jωL, capacitor 1/(jωC)), as in speaker.crossover.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

# Standard-value mantissas: E12 for capacitors and inductors' ladders, E24 for resistors.
E12 = (1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2)
E24 = (1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0,
       3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2, 9.1)

AIR_DCR = 0.5  # ohm at 1 mH, air core
CORED_DCR = 0.14  # ohm at 1 mH, iron or ferrite core


@dataclass(frozen=True)
class Part:
    name: str
    kind: Literal["R", "L", "C", "Z"]
    a: str
    b: str
    value: float  # ohm, H or F; ignored for "Z" (a driver load, impedance from `loads[name]`)
    dcr: float = 0.0  # ohm, series loss of an L


@dataclass
class Solution:
    v: dict[str, np.ndarray]  # complex node voltages, "0" included
    z_in: np.ndarray  # complex impedance seen by the 1 V source

    def across(self, p: Part) -> np.ndarray:
        return self.v[p.a] - self.v[p.b]


def _admittance(p: Part, w: np.ndarray, loads: dict[str, np.ndarray]) -> np.ndarray:
    if p.kind == "R":
        return np.full(w.shape, 1.0 / p.value, dtype=complex)
    if p.kind == "L":
        return 1.0 / (p.dcr + 1j * w * p.value)
    if p.kind == "C":
        return 1j * w * p.value
    if p.name not in loads:
        raise ValueError(f"no load impedance given for driver {p.name!r}")
    return 1.0 / np.asarray(loads[p.name], dtype=complex)


def solve(parts: list[Part], f: np.ndarray, loads: dict[str, np.ndarray]) -> Solution:
    """Nodal analysis at every frequency; a 1 V source drives node "in" against ground "0"."""
    f = np.asarray(f, dtype=float)
    w = 2 * np.pi * f
    nodes = sorted(({n for p in parts for n in (p.a, p.b)} | {"in"}) - {"0"})
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    # Unknowns: node voltages, then the source current. The last row pins v_in = 1 V.
    m = np.zeros((len(f), n + 1, n + 1), dtype=complex)
    for p in parts:
        y = _admittance(p, w, loads)
        ia, ib = idx.get(p.a), idx.get(p.b)
        if ia is not None:
            m[:, ia, ia] += y
        if ib is not None:
            m[:, ib, ib] += y
        if ia is not None and ib is not None:
            m[:, ia, ib] -= y
            m[:, ib, ia] -= y
    i_in = idx["in"]
    m[:, i_in, n] = -1.0
    m[:, n, i_in] = 1.0
    rhs = np.zeros((len(f), n + 1, 1), dtype=complex)
    rhs[:, n, 0] = 1.0
    x = np.linalg.solve(m, rhs)[:, :, 0]
    v = {name: x[:, i] for name, i in idx.items()}
    v["0"] = np.zeros(len(f), dtype=complex)
    return Solution(v=v, z_in=1.0 / x[:, n])


def _series(kind: Literal["R", "L", "C"]) -> np.ndarray:
    """Every standard value of a kind, ascending, over all useful decades."""
    if kind == "R":
        return np.array([m * 10.0 ** d for d in range(-1, 6) for m in E24])
    if kind == "C":
        return np.array([m * 10.0 ** d for d in range(-9, -2) for m in E12])
    # 0.1 mH steps below 1 mH, then 0.5 mH steps (to 100 mH)
    return np.array([k * 1e-4 for k in range(1, 10)] + [k * 0.5e-3 for k in range(2, 201)])


def snap(kind: Literal["R", "L", "C"], value: float) -> float:
    """Nearest standard value on a log scale (never zero)."""
    s = _series(kind)
    return float(s[np.argmin(np.abs(np.log(s) - np.log(max(value, 1e-12))))])


def neighbours(kind: Literal["R", "L", "C"], value: float) -> tuple[float, float]:
    """The standard values just below and just above a standard value."""
    s = _series(kind)
    i = int(np.argmin(np.abs(np.log(s) - np.log(value))))
    return float(s[max(i - 1, 0)]), float(s[min(i + 1, len(s) - 1)])


def inductor_dcr(value_h: float, cored: bool) -> float:
    """Series resistance (ohm) of a coil: DCR at 1 mH scaled with the square root of the inductance
    (more turns of the same gauge in a bigger winding)."""
    return (CORED_DCR if cored else AIR_DCR) * float(np.sqrt(value_h / 1e-3))
