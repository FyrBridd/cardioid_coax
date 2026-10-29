"""Driver models. The only place driver TOML units (mm, L, cm^2, mH) are converted to SI."""
from __future__ import annotations

import math
import tomllib
from dataclasses import dataclass
from pathlib import Path

from speaker.constants import C, RHO

DRIVERS_DIR = Path(__file__).resolve().parents[2] / "data" / "drivers"


@dataclass(frozen=True)
class Driver:
    name: str
    fs: float  # Hz
    qes: float
    qms: float
    vas: float  # m^3
    sd: float  # m^2
    re: float  # ohm
    le: float  # H
    xmax: float  # m
    outer_d: float  # m
    cutout_d: float  # m
    depth: float  # m
    flange_t: float  # m
    basket_d: float  # m
    magnet_d: float  # m
    sens_db: float  # datasheet, informational
    hf_sens_db: float | None = None
    diaphragm_d: float | None = None  # m
    cone_depth: float | None = None  # m, surround plane to cone inner edge
    cone_inner_d: float | None = None  # m
    hf_exit_d: float | None = None  # m
    hf_xo_min: float | None = None  # Hz, the datasheet's recommended (lowest) HF crossover

    @property
    def qts(self) -> float:
        return self.qes * self.qms / (self.qes + self.qms)

    @property
    def cms(self) -> float:
        return self.vas / (RHO * C**2 * self.sd**2)

    @property
    def mms(self) -> float:
        return 1.0 / ((2 * math.pi * self.fs) ** 2 * self.cms)

    @property
    def bl(self) -> float:
        return math.sqrt(2 * math.pi * self.fs * self.mms * self.re / self.qes)

    @property
    def rms(self) -> float:
        return 2 * math.pi * self.fs * self.mms / self.qms


def _opt_mm(raw: dict, key: str) -> float | None:
    return raw[key] / 1000.0 if key in raw else None


def load_driver(stem: str) -> Driver:
    with open(DRIVERS_DIR / f"{stem}.toml", "rb") as fh:
        raw = tomllib.load(fh)
    sd = raw["sd_cm2"] / 10000.0
    diaphragm_d = _opt_mm(raw, "diaphragm_d_mm")
    if diaphragm_d is None:
        diaphragm_d = math.sqrt(4 * sd / math.pi)
    return Driver(
        name=raw["name"],
        fs=float(raw["fs"]),
        qes=float(raw["qes"]),
        qms=float(raw["qms"]),
        vas=raw["vas_l"] / 1000.0,
        sd=sd,
        re=float(raw["re"]),
        le=raw["le_mh"] / 1000.0,
        xmax=raw["xmax_mm"] / 1000.0,
        outer_d=raw["outer_d_mm"] / 1000.0,
        cutout_d=raw["cutout_d_mm"] / 1000.0,
        depth=raw["depth_mm"] / 1000.0,
        flange_t=raw["flange_t_mm"] / 1000.0,
        basket_d=raw["basket_d_mm"] / 1000.0,
        magnet_d=raw["magnet_d_mm"] / 1000.0,
        sens_db=float(raw["sens_db"]),
        hf_sens_db=float(raw["hf_sens_db"]) if "hf_sens_db" in raw else None,
        diaphragm_d=diaphragm_d,
        cone_depth=_opt_mm(raw, "cone_depth_mm"),
        cone_inner_d=_opt_mm(raw, "cone_inner_d_mm"),
        hf_exit_d=_opt_mm(raw, "hf_exit_d_mm"),
        hf_xo_min=float(raw["hf_xo_min_hz"]) if "hf_xo_min_hz" in raw else None,
    )
