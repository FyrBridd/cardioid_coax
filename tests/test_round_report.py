import json
import struct
import zlib
from pathlib import Path

import pytest

import speaker.round as round_mod
from speaker.cli import main
from speaker.evaluate import SETTINGS
from speaker.round_report import write_report
from speaker.variants import Compare, RoundSpec, StageSpec, VariantSpec, apply_overrides, read_raw, write_params

DEFAULT = Path(__file__).resolve().parents[1] / "params" / "default.toml"
PLOTS = ("polar_h", "polar_v", "rear_rejection", "spinorama", "whole_rear", "window_family")
TOML = '''
number = 7
title = "Test round"
base = "params/default.toml"
trusted_limit_hz = 2851
trusted_note = "inherited from round 0's mesh check"
finals = 2

[[stage]]
name = "vents"
settings = "vents"
rank_by = "cardioid"
include_round0 = true
grid = { "vents.centre_from_front_mm" = [190] }

[[stage]]
name = "shape"
settings = "screen"
rank_by = "directivity"
include_round0 = true
[[stage.variant]]
name = "best"
set = { "cabinet.front_edge" = "roundover", "cabinet.edge_size_mm" = 30, "cabinet.width_mm" = 400 }
'''


def png() -> bytes:
    def chunk(tag, data):
        c = struct.pack(">I", len(data)) + tag + data
        return c + struct.pack(">I", zlib.crc32(tag + data))
    raw = zlib.compress(b"\x00\x00\x00\x00")
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", 1, 1, 8, 0, 0, 0, 0))
            + chunk(b"IDAT", raw) + chunk(b"IEND", b""))


def score(key, value, status="pass", conf="trusted"):
    return {"key": key, "label": key.replace("_", " "), "value": value, "unit": "dB", "status": status,
            "band": [150.0, 300.0], "confidence": conf, "formula": "f", "detail": {}}


def scores(hard=12.0, spread=0.5):
    return [score("cardioid_hard", hard, "pass" if hard >= 10 else "fail"), score("cardioid_aim", 8.0),
            score("rear_hemisphere", 15.0), score("window_spread", spread), score("presence_lift", 0.2, "warn"),
            score("di_smooth_bass", 1.0), score("di_smooth_lfhf", 1.1), score("flat23", 0.9),
            score("window_spread", 3.3, "fail", "low")]


def tune(margin=2.0, base=1.0):
    return {"resistance": 250.0, "net_volume": 0.0085, "hard": 12.0, "aim": 8.0, "margin": margin,
            "spread": {"0.7": [11.0, 7.0], "1.3": [12.5, 8.5]}, "base": {"hard": 11.0, "aim": 7.0, "margin": base}}


def vresult(stage, name, status="ok", reason="", overrides=None, sc=None, tuned=True, **kw):
    return {"name": name, "stage": stage, "status": status, "reason": reason, "overrides": overrides or {},
            "params_sha256": "0" * 64, "settings": "screen", "code_version": "cv", "commit": "c" * 40,
            "scores": sc if sc is not None else (scores() if status == "ok" else []),
            "tune": tune(**kw) if status == "ok" else None, "tuned": tuned and status == "ok", "wall_s": 1.0}


SPEC = RoundSpec(7, "Test round", DEFAULT, 2851.0, "inherited from round 0's mesh check", 2, (),
                 (StageSpec("vents", "vents", "cardioid", True, (VariantSpec("cff=190", {}),)),
                  StageSpec("shape", "screen", "directivity", True, (VariantSpec("best", {}),))))


def put(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def build(out: Path, *, gate_unmet=False, notes=None, compare=None):
    base = read_raw(DEFAULT)
    best_raw = apply_overrides(base, {"cabinet.front_edge": "roundover", "cabinet.edge_size_mm": 30,
                                      "cabinet.width_mm": 400, "vents.specific_resistance_rayl": 250})
    stages = []
    vres = [vresult("vents", "stage-base", sc=scores()[:3]), vresult("vents", "round0", tuned=False, sc=scores()[:3]),
            vresult("vents", "cff=190", overrides={"vents.centre_from_front_mm": 190}, sc=scores()[:3])]
    shape = [vresult("shape", "stage-base"), vresult("shape", "round0", tuned=False),
             vresult("shape", "best", overrides={"cabinet.front_edge": "roundover"}),
             vresult("shape", "bad<one>", status="rejected", reason="flat front width <too> narrow"),
             vresult("shape", "boom", status="error", reason="Traceback\nRuntimeError: out of memory")]
    for stage, rank_by, results, order, out_ in (
            ("vents", "cardioid", vres, ["cff=190", "stage-base", "round0"], {}),
            ("shape", "directivity", shape, ["best", "stage-base", "round0"],
             {"bad<one>": "rejected: flat front width <too> narrow", "boom": "error: RuntimeError: out of memory"})):
        for r in results:
            put(out / stage / r["name"] / "result.json", r)
            write_params(best_raw if r["name"] == "best" else base, out / stage / r["name"] / "params_tuned.toml")
        ranking = {"order": order, "out": out_, "gate_unmet": gate_unmet and stage == "shape"}
        put(out / stage / "ranking.json", ranking)
        stages.append({"name": stage, "settings": results[0]["settings"], "rank_by": rank_by, "results": results,
                       "ranking": ranking})
    finals = []
    for k, name in enumerate(["best", "stage-base"], 1):
        fdir = out / "finals" / name
        (fdir / "plots").mkdir(parents=True)
        for pl in PLOTS:
            (fdir / "plots" / f"{pl}.png").write_bytes(png())
        put(fdir / "scores.json", scores(spread=0.2 if name == "best" else 0.6))
        (fdir / "report.html").write_text("<html></html>")
        put(fdir / "gains.json", {"bass": -1.25, "coax_lf": 0.0, "coax_hf": -3.5})
        finals.append({"rank": k, "name": name, "stage": "shape", "params": f"shape/{name}/params_tuned.toml",
                       "dir": f"finals/{name}", "report": f"finals/{name}/report.html", "status": "ok", "reason": "",
                       "wall_s": 1.0})
    put(out / "finals.json", finals)
    put(out / "results.json", {"spec": "round.toml", "spec_sha256": "a" * 64, "number": 7, "title": "Test round",
                               "code_version": "cv", "commit": "c" * 40, "stages": stages, "finals": finals})
    if notes is not None:
        (out / "notes.md").write_text(notes)
    spec = SPEC
    if compare:
        put(compare, scores(hard=9.0))
        spec = RoundSpec(
            SPEC.number, SPEC.title, SPEC.base, SPEC.trusted_limit, SPEC.trusted_note, SPEC.finals,
            (Compare("Round 0 final", compare),), SPEC.stages)
    return spec


@pytest.fixture
def round_dir(tmp_path):
    out = tmp_path / "round"
    build(out)
    return out


def html_of(out, spec=SPEC):
    return write_report(out, spec).read_text()


def test_report_sections_in_spec_order(round_dir):
    html = html_of(round_dir)
    heads = ["Recommendation", "Curves", "Scores", "Confidence", "Crossover", "What I'd try next", "All variants"]
    idx = [html.index(f">{h}<") for h in heads]
    assert idx == sorted(idx)
    assert (round_dir / "report.html").exists() and (round_dir / "shape.png").stat().st_size > 5000


def test_notes_fill_the_recommendation_and_next_steps(tmp_path):
    out = tmp_path / "r"
    build(out, notes="# Notes\n\n## Recommendation\n\nBuild <best> with a 30 mm roundover.\n\n- first & foremost\n"
                     "- second\n\n## What I'd try next\n\nA wider baffle?\n")
    html = html_of(out)
    assert "Build &lt;best&gt; with a 30 mm roundover." in html and "<li>first &amp; foremost</li>" in html
    assert "A wider baffle?" in html
    assert "Draft: the recommendation is written in the next session." not in html


def test_placeholder_without_notes(round_dir):
    assert "Draft: the recommendation is written in the next session." in html_of(round_dir)


def test_rejected_variants_listed_with_reason(round_dir):
    html = html_of(round_dir)
    assert "bad&lt;one&gt;" in html and "flat front width &lt;too&gt; narrow" in html
    assert "RuntimeError: out of memory" in html


def test_gate_unmet_banner(tmp_path):
    ok, bad = tmp_path / "ok", tmp_path / "bad"
    build(ok)
    build(bad, gate_unmet=True)
    banner = "hard cardioid gate is unmet by every variant"
    assert banner not in html_of(ok) and banner in html_of(bad)
    assert ("ranked on directivity among the variants within 0.5 dB of the best rear rejection, then by "
            "rear rejection") in html_of(bad)


def test_curves_are_copied_and_shown(round_dir):
    html = html_of(round_dir)
    for pl in ("polar_h", "polar_v", "spinorama", "window_family", "rear_rejection"):
        assert (round_dir / "plots" / f"{pl}.png").exists() and f'src="plots/{pl}.png"' in html
    assert 'src="shape.png"' in html


def test_score_table_rows_and_greyed_low_table(tmp_path):
    out = tmp_path / "r"
    spec = build(out, compare=tmp_path / "prev.json")
    html = html_of(out, spec)
    assert "Round 0 final" in html and "cardioid_hard" in html
    assert 'num fail">9.00' in html and 'num warn">0.20' in html and 'num pass">12.00' in html
    grey = html[html.index('<table class="grey">'):]
    assert "3.30" in grey and "3.30" not in html[:html.index('<table class="grey">')]  # low scores only in the grey table


def test_score_table_has_round0_rescored_and_flags_compare_entries(tmp_path):
    out = tmp_path / "r"
    spec = build(out, compare=tmp_path / "prev.json")
    html = html_of(out, spec)
    rescored = "Round 0 rescored (this code, screen settings, untuned)"
    compare = "Round 0 final (scored with earlier code; not like-for-like)"
    table = html[html.index(">Scores<"):html.index('<table class="grey">')]
    assert table.index("2. stage-base (runner-up; screen rank 2)") < table.index(rescored) < table.index(compare)


def test_score_table_omits_round0_rescored_when_not_ok(tmp_path):
    out = tmp_path / "r"
    build(out)
    r0 = out / "shape" / "round0" / "result.json"
    put(r0, {**json.loads(r0.read_text()), "status": "error", "scores": []})
    data = json.loads((out / "results.json").read_text())
    for r in data["stages"][-1]["results"]:
        if r["name"] == "round0":
            r.update(status="error", scores=[])
    put(out / "results.json", data)
    assert "Round 0 rescored" not in html_of(out)


def test_confidence_says_the_screen_density_is_outside_the_mesh_check(round_dir):
    html = html_of(round_dir)
    pw = SETTINGS["screen"].per_wavelength
    assert (f"The mesh check that set the trusted limit compared 6 against 8 elements per wavelength, so the "
            f"screen's {pw:g} per wavelength is not covered by it.") in html
    assert "presence lift especially, are therefore less certain than the finals' scores" in html
    assert "The recommendation should be argued from the finals' 6-per-wavelength scores." in html


def test_confidence_states_the_known_limits(round_dir):
    html = html_of(round_dir)
    assert "2715 Hz" in html and "2851 Hz" in html
    assert "3-7 %" in html and "h_coax 10 mm" in html
    assert "not mesh-checked" in html


def test_crossover_lists_points_and_way_gains(round_dir):
    html = html_of(round_dir)
    assert "170" in html and "1400" in html and "-1.25 dB" in html and "-3.50 dB" in html


def test_appendix_cardioid_stage_shows_only_cardioid_scores(round_dir):
    html = html_of(round_dir)
    app = html[html.index(">All variants<"):]
    vents = app[:app.index("shape")]
    assert "rear_hemisphere" in vents and "window_spread" not in vents
    assert "+1.0 dB" in app  # what tuning bought: margin 2.0 against base 1.0
    assert "250" in app and "8.5" in app  # tuned R and V (litres)


def test_no_results_is_a_clear_error(tmp_path):
    with pytest.raises(ValueError, match="no results"):
        write_report(tmp_path, SPEC)


def test_report_only_does_not_evaluate(tmp_path, monkeypatch, capsys):
    out = tmp_path / "round"
    build(out)
    (out / "report.html").unlink(missing_ok=True)
    toml = tmp_path / "round.toml"
    toml.write_text(TOML)

    def boom(*a, **k):
        raise AssertionError("evaluate must not run")
    monkeypatch.setattr(round_mod, "evaluate", boom)
    monkeypatch.setattr(round_mod, "run_directivity", boom)
    assert main(["round", str(toml), "--out", str(out), "--report-only"]) == 0
    assert (out / "report.html").exists() and f"wrote {out / 'report.html'}" in capsys.readouterr().out


def test_report_only_without_results_says_so(tmp_path, capsys):
    toml = tmp_path / "round.toml"
    toml.write_text(TOML)
    assert main(["round", str(toml), "--out", str(tmp_path / "none"), "--report-only"]) == 2
    assert "no results" in capsys.readouterr().err


def test_missing_gains_warn_visibly(round_dir):
    (round_dir / "finals" / "best" / "gains.json").unlink()
    html = html_of(round_dir)
    assert "Way gains not found (no gains.json) for: best" in html and "n/a" in html
    assert "-1.25 dB" in html  # the other final still shows its gains


def test_errored_first_final_is_listed_and_rank_is_real(round_dir):
    finals = json.loads((round_dir / "finals.json").read_text())
    finals[0].update(status="error", reason="Traceback\nRuntimeError: out of <memory>")
    put(round_dir / "finals.json", finals)
    html = html_of(round_dir)
    assert "<b>stage-base</b> (screen rank 2;" in html
    assert "Final 1, best, did not complete (error): RuntimeError: out of &lt;memory&gt;" in html


def test_finals_are_ranked_on_their_own_scores(round_dir):
    # best is screen rank 1, but its own final scores are worse than stage-base's (window spread 0.9 vs 0.6)
    put(round_dir / "finals" / "best" / "scores.json", scores(spread=0.9))
    html = html_of(round_dir)
    assert "<b>stage-base</b> (screen rank 2; the best final on its own scores)" in html
    table = html[html.index(">Scores<"):]
    assert table.index("1. stage-base (recommendation; screen rank 2)") < table.index(
        "2. best (runner-up; screen rank 1)")


def test_base_label_names_the_reference_design(tmp_path):
    import dataclasses
    out = tmp_path / "r"
    build(out)
    html = html_of(out, dataclasses.replace(SPEC, base_label="the round 1 recommendation"))
    for s in ("The round 1 recommendation rescored (this code, screen settings, untuned)",
              "<th>The round 1 recommendation</th>", "against the round 1 recommendation.",
              "Rows that differ from the round 1 recommendation are shaded.",
              "the mesh check was run on the round 1 recommendation"):
        assert s in html, s
    assert "Round 0" not in html.replace("Round 0 final", "") and "round 0" not in html.replace("round 0&#x27;s mesh", "")


def test_round_report_without_crossover_json_falls_back_to_gains(round_dir):
    html = html_of(round_dir)
    sec = html[html.index(">Crossover<"):html.index(">What I'd try next<")]
    assert "<th>bass</th>" in sec and "-1.25 dB" in sec and "-3.50 dB" in sec
    assert "crossover.json" not in sec and "Traceback" not in sec


def test_round_report_shows_the_fitted_crossover(round_dir):
    import dataclasses
    from speaker.crossover import ActiveCrossover, Peq, WayFilter
    xo = ActiveCrossover(190.0, 1600.0, {"bass": WayFilter(-1.0, 0.0, 1, (Peq(100.0, 2.0, 1.0),)),
                                         "coax_lf": WayFilter(0.0, 0.25e-3, 1, ()),
                                         "coax_hf": WayFilter(-3.0, 0.0, 1, (Peq(3500.0, 7.0, 1.0),))}, 16000.0)
    put(round_dir / "finals" / "best" / "crossover.json", dataclasses.asdict(xo) | {"trusted_limit": 2851.0})
    html = html_of(round_dir)
    sec = html[html.index(">Crossover<"):html.index(">What I'd try next<")]
    assert "190 Hz" in sec and "1600 Hz" in sec and "0.25 ms" in sec
    assert "7.0 dB" in sec
    assert "-3.50 dB" in sec   # the other final has no crossover.json and keeps its gains row
