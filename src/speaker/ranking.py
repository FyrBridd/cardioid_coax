"""Ranking a stage's variant results: by the tuned cardioid margin, or by directivity behind the gate.

directivity (spec §4, trusted scores only): a variant whose trusted cardioid_hard fails is out, and the
rest are put in directivity order (directivity_order): sorted on window_spread + presence_lift into
tiers, each holding every row within 0.25 dB of the tier's best, the first tier starting at the best
row overall and each next tier at the best row left. Inside a tier the rows are ordered ascending on
(di_smooth_bass + di_smooth_lfhf, flat23).

If every ok variant other than round0 (a reference, not a candidate) fails the gate, gate_unmet is
set and all ok variants, round0 included, are still ranked, so the round has an order and a best to
carry on. The best candidate's cardioid_hard (round0 excluded) sets the bar: every row whose
cardioid_hard is no more than GATE_UNMET_TOL_DB below it counts as at the gate, and those rows come
first, ordered by that key alone. The other rows follow, by cardioid_hard (descending) and then by
directivity order. Ranking on raw cardioid_hard alone would let tenths of a dB of rear rejection decide the
shape stage, with directivity only breaking exact float ties.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Literal, Sequence, TypeVar

DIRECTIVITY_KEYS = ("cardioid_hard", "window_spread", "presence_lift", "di_smooth_bass", "di_smooth_lfhf", "flat23")
BUCKET_DB = 0.25
REFERENCE = "round0"  # ranked for comparison, but never a stage's best
GATE_UNMET_TOL_DB = 0.5  # gate unmet: rows this close to the best candidate's cardioid_hard rank on directivity
_EPS = 1e-9  # dB, float slack on the tier edge

T = TypeVar("T")


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
    scoring: int = 1  # speaker.evaluate.SCORING_VERSION it was scored with; last, so older files load as 1


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
    return _trusted(r.scores, key)


def _trusted(scores: list[dict], key: str) -> dict | None:
    for s in scores:
        if s["key"] == key and s["confidence"] == "trusted" and s["value"] is not None:
            return s
    return None


def _margin(r: VariantResult) -> float | None:
    if r.tune is None:
        return None
    m = r.tune["margin"] if r.tuned else r.tune["base"]["margin"]
    return None if m is None or not math.isfinite(m) else float(m)


def directivity_order(items: Sequence[T], values: Callable[[T], dict[str, float]]) -> list[T]:
    """items best first: tiers of window_spread + presence_lift, each within BUCKET_DB of its best, then
    di_smooth_bass + di_smooth_lfhf and flat23 inside a tier. Ties keep the input order."""
    def wp(x: T) -> float:
        v = values(x)
        return v["window_spread"] + v["presence_lift"]

    left, out = sorted(items, key=wp), []
    while left:
        top = wp(left[0]) + BUCKET_DB + _EPS
        n = next((i for i, x in enumerate(left) if wp(x) > top), len(left))
        tier, left = left[:n], left[n:]
        out += sorted(tier, key=lambda x: (values(x)["di_smooth_bass"] + values(x)["di_smooth_lfhf"],
                                           values(x)["flat23"]))
    return out


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
    r = rank_directivity([(x.name, x.scores) for x in ok])
    return Ranking(r.order, out | r.out, r.gate_unmet)


def rank_directivity(rows: Sequence[tuple[str, list[dict]]]) -> Ranking:
    """The directivity ranking (module docstring) of (name, Score dicts) rows. Ties keep the input order."""
    out: dict[str, str] = {}
    vals = []  # (name, values, passes the gate)
    for name, scores in rows:
        found = {k: _trusted(scores, k) for k in DIRECTIVITY_KEYS}
        missing = [k for k, s in found.items() if s is None]
        if missing:
            out[name] = f"missing trusted score {missing[0]}"
            continue
        vals.append((name, {k: float(s["value"]) for k, s in found.items()},
                     found["cardioid_hard"]["status"] != "fail"))
    v = lambda x: x[1]  # noqa: E731
    # round0 is a reference, not a candidate: the candidates alone decide whether the gate is unmet,
    # and the resulting rule then applies to every row, round0 included.
    candidates = [x for x in vals if x[0] != REFERENCE]
    if candidates and not any(passes for _, _, passes in candidates):
        bar = max(x[1]["cardioid_hard"] for x in candidates) - GATE_UNMET_TOL_DB
        near = directivity_order([x for x in vals if x[1]["cardioid_hard"] >= bar], v)
        rest = sorted(directivity_order([x for x in vals if x[1]["cardioid_hard"] < bar], v),
                      key=lambda x: -x[1]["cardioid_hard"])
        return Ranking([n for n, _, _ in near + rest], out, True)
    for name, x, passes in vals:
        if not passes:
            out[name] = f"fails the hard cardioid gate ({x['cardioid_hard']:.1f} dB)"
    passing = directivity_order([x for x in vals if x[2]], v)
    return Ranking([n for n, _, _ in passing], out, False)
