"""``speaker round``: run a round file's stages, tuning, scoring and ranking each variant, then the finals.

Layout of out_dir (default out/rounds/<number>):

    <stage>/<name>/params.toml          the variant as run (pre-tune)
    <stage>/<name>/params_tuned.toml    with the tuned vents (round0: the same as params.toml)
    <stage>/<name>/scores.json          [asdict(Score)]
    <stage>/<name>/result.json          asdict(VariantResult), written last and atomically
    <stage>/ranking.json                asdict(Ranking)
    finals/<name>/                      run_directivity's report, scores and plots
    finals.json                         the finals in rank order
    results.json                        every stage's results and ranking, the finals, the spec and its sha256
    report.html, shape.png, plots/      the round report (round_report.write_report); notes.md feeds it

A result.json is reused when its params_sha256, settings, code_version and scoring version match (an
"error" result is always retried, as its cause may be transient), so a killed round
resumes where it stopped; a partial write only ever leaves a result.json.<pid>.tmp, which is never read.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path

from .bem.run import code_version
from .directivity_report import _jsonable, run_directivity
from .evaluate import SCORING_VERSION, SETTINGS, evaluate
from .ranking import REFERENCE, Ranking, VariantResult, rank
from .round_report import write_report
from .variants import REPO, RoundSpec, StageSpec, apply_overrides, load_round, read_raw, with_tuned, write_params

TRACEBACK_LINES = 20


def git_commit() -> str:
    """HEAD of the repository, with "-dirty" when tracked files have changed; "unknown" if git fails."""
    try:
        run = lambda *a: subprocess.run(["git", *a], cwd=REPO, capture_output=True, text=True,  # noqa: E731
                                        check=True, timeout=30).stdout
        head = run("rev-parse", "HEAD").strip()
        dirty = run("status", "--porcelain", "--untracked-files=no").strip()
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    return head + ("-dirty" if dirty else "") if head else "unknown"


def write_json(path: Path, data) -> None:
    """Write JSON to path atomically: a sibling temp file, fsynced, then os.replace."""
    tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    with open(tmp, "w") as fh:
        fh.write(json.dumps(_jsonable(data), indent=2, allow_nan=False))
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _traceback() -> str:
    """The last TRACEBACK_LINES lines of the exception being handled."""
    return "\n".join(traceback.format_exc().rstrip().splitlines()[-TRACEBACK_LINES:])


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_result(path: Path) -> VariantResult | None:
    """The stored result, or None when it is missing or unreadable (never a partial file)."""
    try:
        return VariantResult(**json.loads(path.read_text()))
    except (OSError, ValueError, TypeError):
        return None


def _run_variant(vdir: Path, stage: StageSpec, name: str, overrides: dict, raw: dict, *, tune: bool,
                 spec: RoundSpec, cv: str, commit: str) -> tuple[VariantResult, bool]:
    """Evaluate one variant in vdir, or reuse its result. Returns (result, reused)."""
    vdir.mkdir(parents=True, exist_ok=True)
    params = vdir / "params.toml"
    write_params(raw, params)
    sha = _sha256(params)
    prev = _load_result(vdir / "result.json")
    if (prev is not None and (prev.params_sha256, prev.settings, prev.code_version) == (sha, stage.settings, cv)
            and prev.scoring == SCORING_VERSION  # scored as this code scores: older results are re-scored
            and prev.status != "error"  # errors may be transient (memory, OpenCL): retry them
            and (prev.status != "ok" or (vdir / "params_tuned.toml").exists())):
        return prev, True

    for stale in ("params_tuned.toml", "scores.json"):
        (vdir / stale).unlink(missing_ok=True)
    base = dict(name=name, stage=stage.name, overrides=overrides, params_sha256=sha, settings=stage.settings,
                code_version=cv, commit=commit, scoring=SCORING_VERSION)
    t0 = time.perf_counter()
    try:
        ev = evaluate(params, SETTINGS[stage.settings], trusted_limit=spec.trusted_limit, tune=tune)
        write_params(with_tuned(raw, ev.tune.resistance, ev.tune.net_volume) if ev.tuned else raw,
                     vdir / "params_tuned.toml")
        scores = [dataclasses.asdict(s) for s in ev.scores]
        write_json(vdir / "scores.json", scores)
        res = VariantResult(**base, status="ok", reason="", scores=scores, tune=dataclasses.asdict(ev.tune),
                            tuned=bool(ev.tuned), wall_s=float(ev.wall_s))
    except ValueError as e:
        res = VariantResult(**base, status="rejected", reason=str(e), scores=[], tune=None, tuned=False,
                            wall_s=time.perf_counter() - t0)
    except Exception:
        res = VariantResult(**base, status="error", reason=_traceback(), scores=[], tune=None, tuned=False,
                            wall_s=time.perf_counter() - t0)
    res = VariantResult(**_jsonable(dataclasses.asdict(res)))  # as it reads back: tuples as lists, NaN as None
    write_json(vdir / "result.json", dataclasses.asdict(res))
    return res, False


def _best(ranking: Ranking) -> str | None:
    return next((n for n in ranking.order if n != REFERENCE), None)


def _final(out_dir: Path, stage: str, k: int, name: str, spec: RoundSpec) -> dict:
    params = out_dir / stage / name / "params_tuned.toml"
    entry = {"rank": k, "name": name, "stage": stage, "params": f"{stage}/{name}/params_tuned.toml",
             "dir": f"finals/{name}", "report": f"finals/{name}/report.html"}
    fdir = out_dir / "finals" / name
    fdir.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    try:
        run_directivity(params, fdir, trusted_limit=(spec.trusted_limit, spec.trusted_note))
        entry.update(status="ok", reason="")
    except Exception:
        entry.update(status="error", reason=_traceback())
    entry["wall_s"] = time.perf_counter() - t0
    return entry


def run_round(spec_path: Path, out_dir: Path | None = None) -> Path:
    """Run every stage of the round file at spec_path, then its finals; writes results.json, then the report and returns its path."""
    spec_path = Path(spec_path)
    spec = load_round(spec_path)
    out_dir = Path("out/rounds") / str(spec.number) if out_dir is None else Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cv, commit = code_version(), git_commit()
    round0_raw = read_raw(spec.base)
    base_raw = round0_raw
    stages = []
    for stage in spec.stages:
        runs = [("stage-base", {}, base_raw, True)]
        if stage.include_round0:
            runs.append(("round0", {}, round0_raw, False))
        runs += [(v.name, v.overrides, apply_overrides(base_raw, v.overrides), True) for v in stage.variants]
        results = []
        for i, (name, overrides, raw, tune) in enumerate(runs, 1):
            res, reused = _run_variant(out_dir / stage.name / name, stage, name, overrides, raw, tune=tune,
                                       spec=spec, cv=cv, commit=commit)
            results.append(res)
            print(f"[{stage.name} {i}/{len(runs)}] {name}: {res.status}{' (reused)' if reused else ''}, "
                  f"{res.wall_s:.1f} s", file=sys.stderr, flush=True)
        ranking = rank(results, stage.rank_by)
        write_json(out_dir / stage.name / "ranking.json", dataclasses.asdict(ranking))
        stages.append({"name": stage.name, "settings": stage.settings, "rank_by": stage.rank_by,
                       "results": [dataclasses.asdict(r) for r in results], "ranking": dataclasses.asdict(ranking)})
        best = _best(ranking)
        if best is None:
            raise ValueError(f"stage {stage.name!r}: no ranked variant other than round0, so there is no best "
                             f"to carry on; see {out_dir / stage.name / 'ranking.json'}")
        base_raw = read_raw(out_dir / stage.name / best / "params_tuned.toml")

    finals = []
    if spec.stages:
        last = spec.stages[-1].name
        top = [n for n in stages[-1]["ranking"]["order"] if n != "round0"][: spec.finals]
        for k, name in enumerate(top, 1):
            finals.append(_final(out_dir, last, k, name, spec))
            print(f"[finals {k}/{len(top)}] {name}: {finals[-1]['status']}, {finals[-1]['wall_s']:.1f} s",
                  file=sys.stderr, flush=True)
    write_json(out_dir / "finals.json", finals)

    results = out_dir / "results.json"
    write_json(results, {"spec": str(spec_path), "spec_sha256": _sha256(spec_path), "number": spec.number,
                      "title": spec.title, "code_version": cv, "commit": commit, "stages": stages,
                      "finals": finals})
    return write_report(out_dir, spec)


def report_only(spec_path: Path, out_dir: Path | None = None) -> Path:
    """Re-render the report from the stored results of the round at spec_path, evaluating nothing."""
    spec = load_round(Path(spec_path))
    out_dir = Path("out/rounds") / str(spec.number) if out_dir is None else Path(out_dir)
    return write_report(out_dir, spec)
