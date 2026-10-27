"""Ranking a stage's variant results: by the tuned cardioid margin, or by directivity behind the gate.

directivity (spec §4, trusted scores only): a variant whose trusted cardioid_hard fails is out, and the
rest are ordered ascending on

    (floor((window_spread + presence_lift) / 0.25), di_smooth_bass + di_smooth_lfhf, flat23)

If every ok variant other than round0 (a reference, not a candidate) fails the gate, gate_unmet is
set and all ok variants, round0 included, are ranked by cardioid_hard (descending) and then by that
key, so the round still has an order and a best to carry on.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

DIRECTIVITY_KEYS = ("cardioid_hard", "window_spread", "presence_lift", "di_smooth_bass", "di_smooth_lfhf", "flat23")
BUCKET_DB = 0.25
REFERENCE = "round0"  # ranked for comparison, but never a stage's best


@dataclass
class VariantResult:
    name: str
    stage: str
    status: Literal["ok", "rejected", "error"]
    reason: str
    overrides: dict
    params_sha256: str  # of the pre-tune params.toml
    settings: str  # EvalSettings name
    code_version: str
    commit: str
    scores: list[dict]  # dataclasses.asdict(Score)
    tune: dict | None  # dataclasses.asdict(VentTune)
    tuned: bool
    wall_s: float


@dataclass
class Ranking:
    order: list[str]
    out: dict[str, str]  # name: why it is not in the order
    gate_unmet: bool


def _not_ok(r: VariantResult) -> str:
    lines = [x for x in r.reason.splitlines() if x.strip()]
    if r.status == "error":  # a traceback: its last line names the exception
        return f"error: {lines[-1].strip() if lines else ''}"
    return f"{r.status}: {r.reason}"


def trusted(r: VariantResult, key: str) -> dict | None:
    """The trusted Score dict of key, or None (absent, or a non-finite value stored as null)."""
    for s in r.scores:
        if s["key"] == key and s["confidence"] == "trusted" and s["value"] is not None:
            return s
    return None


def _margin(r: VariantResult) -> float | None:
    if r.tune is None:
        return None
    m = r.tune["margin"] if r.tuned else r.tune["base"]["margin"]
    return None if m is None or not math.isfinite(m) else float(m)


def _directivity_key(v: dict[str, float]) -> tuple[int, float, float]:
    return (math.floor((v["window_spread"] + v["presence_lift"]) / BUCKET_DB),
            v["di_smooth_bass"] + v["di_smooth_lfhf"], v["flat23"])


def rank(results: list[VariantResult], rank_by: str) -> Ranking:
    """The order (best first) and the reason each other variant is out. Ties keep the input order."""
    out: dict[str, str] = {}
    ok = []
    for r in results:
        if r.status != "ok":
            out[r.name] = _not_ok(r)
        else:
            ok.append(r)

    if rank_by == "cardioid":
        margins = []
        for r in ok:
            m = _margin(r)
            if m is None:
                out[r.name] = "no cardioid margin from the vent tune"
            else:
                margins.append((r.name, m))
        margins.sort(key=lambda x: -x[1])
        return Ranking([n for n, _ in margins], out, False)

    if rank_by != "directivity":
        raise ValueError(f"unknown rank_by {rank_by!r}")
    rows = []  # (name, values, passes the gate)
    for r in ok:
        found = {k: trusted(r, k) for k in DIRECTIVITY_KEYS}
        missing = [k for k, s in found.items() if s is None]
        if missing:
            out[r.name] = f"missing trusted score {missing[0]}"
            continue
        rows.append((r.name, {k: float(s["value"]) for k, s in found.items()},
                     found["cardioid_hard"]["status"] != "fail"))
    # round0 is a reference, not a candidate: the candidates alone decide whether the gate is unmet,
    # and the resulting rule then applies to every row, round0 included.
    candidates = [x for x in rows if x[0] != REFERENCE]
    if candidates and not any(passes for _, _, passes in candidates):
        rows.sort(key=lambda x: (-x[1]["cardioid_hard"], *_directivity_key(x[1])))
        return Ranking([n for n, _, _ in rows], out, True)
    for name, v, passes in rows:
        if not passes:
            out[name] = f"fails the hard cardioid gate ({v['cardioid_hard']:.1f} dB)"
    passing = sorted((x for x in rows if x[2]), key=lambda x: _directivity_key(x[1]))
    return Ranking([n for n, _, _ in passing], out, False)
