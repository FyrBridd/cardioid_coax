"""Round files: stages of variants, each a set of dotted-path overrides of a base params TOML."""
from __future__ import annotations

import copy
import itertools
import re
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import tomli_w

from speaker.evaluate import SETTINGS

REPO = Path(__file__).resolve().parents[2]
RESERVED = ("stage-base", "round0")
RANK_BY = ("cardioid", "directivity")
_NAME = re.compile(r"[A-Za-z0-9._=-]+")


@dataclass(frozen=True)
class VariantSpec:
    name: str
    overrides: dict[str, Any]


@dataclass(frozen=True)
class StageSpec:
    name: str
    settings: str
    rank_by: str
    include_round0: bool
    variants: tuple[VariantSpec, ...]


@dataclass(frozen=True)
class Compare:
    label: str
    scores: Path


@dataclass(frozen=True)
class RoundSpec:
    number: int
    title: str
    base: Path
    trusted_limit: float
    trusted_note: str
    finals: int
    compare: tuple[Compare, ...]
    stages: tuple[StageSpec, ...]


def read_raw(path: Path) -> dict:
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def write_params(raw: dict, path: Path) -> None:
    with open(path, "wb") as fh:
        tomli_w.dump(raw, fh)


def _has(raw: dict, dotted: str) -> bool:
    node: Any = raw
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return False
        node = node[part]
    return True


def apply_overrides(raw: dict, overrides: dict) -> dict:
    out = copy.deepcopy(raw)
    for dotted, value in overrides.items():
        *path, key = dotted.split(".")
        node = out
        for part in path:
            node = node[part]
        node[key] = value
    return out


def with_tuned(raw: dict, resistance: float, net_volume: float) -> dict:
    """net_volume is in cubic metres; the file stores litres."""
    out = copy.deepcopy(raw)
    out["vents"]["specific_resistance_rayl"] = resistance
    out["cabinet"]["coax_chamber_net_volume_l"] = net_volume * 1000.0
    return out


def _fmt(value: Any) -> str:
    if isinstance(value, float) and value.is_integer():
        value = int(value)
    return str(value)


def _grid_variants(grid: dict) -> list[VariantSpec]:
    keys = sorted(grid)
    out = []
    for combo in itertools.product(*(grid[k] for k in keys)):
        name = "-".join(f"{k.split('.')[-1].removesuffix('_mm')}={_fmt(v)}" for k, v in zip(keys, combo))
        out.append(VariantSpec(name, dict(zip(keys, combo))))
    return out


def load_round(path: Path) -> RoundSpec:
    with open(path, "rb") as fh:
        d = tomllib.load(fh)
    base = REPO / d["base"]
    raw = read_raw(base)
    finals = int(d.get("finals", 0))
    if finals < 0:
        raise ValueError(f"finals must be >= 0, got {finals}")
    stages = []
    for s in d.get("stage", []):
        sname = s["name"]
        if s["settings"] not in SETTINGS:
            raise ValueError(f"stage {sname!r}: unknown settings {s['settings']!r} (known: {sorted(SETTINGS)})")
        if s["rank_by"] not in RANK_BY:
            raise ValueError(f"stage {sname!r}: unknown rank_by {s['rank_by']!r} (known: {list(RANK_BY)})")
        variants = _grid_variants(s["grid"]) if s.get("grid") else []
        variants += [VariantSpec(v["name"], dict(v.get("set", {}))) for v in s.get("variant", [])]
        if not variants:
            raise ValueError(f"stage {sname!r}: no variants")
        seen: set[str] = set()
        for v in variants:
            if v.name in RESERVED:
                raise ValueError(f"stage {sname!r}, variant {v.name!r}: reserved name")
            if not _NAME.fullmatch(v.name):
                raise ValueError(f"stage {sname!r}, variant {v.name!r}: bad name (allowed: A-Za-z0-9._=-)")
            if v.name in seen:
                raise ValueError(f"duplicate variant name {v.name!r} in stage {sname!r}")
            seen.add(v.name)
            for key in v.overrides:
                if not _has(raw, key):
                    raise ValueError(f"override key {key!r} is not in the base file (stage {sname!r}, variant {v.name!r})")
        stages.append(StageSpec(sname, s["settings"], s["rank_by"], bool(s.get("include_round0", False)),
                                tuple(variants)))
    return RoundSpec(
        number=int(d["number"]), title=d["title"], base=base,
        trusted_limit=float(d["trusted_limit_hz"]), trusted_note=d.get("trusted_note", ""),
        finals=finals,
        compare=tuple(Compare(c["label"], REPO / c["scores"]) for c in d.get("compare", [])),
        stages=tuple(stages))
