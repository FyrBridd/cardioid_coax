import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import speaker.round as round_mod
from speaker.cli import main
from speaker.evaluate import load_design
from speaker.round import run_round
from speaker.scores import Score
from speaker.variants import read_raw
from speaker.vent_tune import VentTune

ROUND = '''
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
grid = { "vents.centre_from_front_mm" = [140, 190] }

[[stage]]
name = "shape"
settings = "screen"
rank_by = "directivity"
include_round0 = true
[[stage.variant]]
name = "m270"
set = { "cabinet.flare_mouth_d_mm" = 270 }
[[stage.variant]]
name = "m280"
set = { "cabinet.flare_mouth_d_mm" = 280 }
'''

WINNER = "centre_from_front=190"


def fake_tune_of(cff: float) -> tuple[float, float]:
    """The fake's tuned (resistance, net volume) for a vent position in m."""
    return 100.0 + 1000.0 * cff, 0.009 + 0.01 * cff


@pytest.fixture
def fake_eval(monkeypatch):
    """evaluate and run_directivity faked in speaker.round's namespace.

    The design is still loaded (so buildability rules raise their real ValueError). Vents stage: the
    margin peaks at 190 mm from the front. Shape stage: window_spread grows with the flare mouth, and a
    270 mm mouth fails the hard gate. Variant directories named in `fail` raise a RuntimeError.
    """
    st = SimpleNamespace(calls=[], finals=[], fail=set())

    def evaluate(path, settings, *, trusted_limit, tune, cache_dir=None):
        path = Path(path)
        st.calls.append((path.parent.parent.name, path.parent.name, settings.name, tune))
        if path.parent.name in st.fail:
            raise RuntimeError("solver blew up")
        d = load_design(path)
        cff, mouth = d.p.vents.centre_from_front, d.p.cabinet.flare_mouth_d
        margin = -abs(cff - 0.19) * 100.0
        r, v = fake_tune_of(cff)
        hard = 8.0 if mouth < 0.275 else 12.0
        vals = {"cardioid_hard": hard, "cardioid_aim": 8.0, "window_spread": (mouth - 0.26) * 100.0,
                "presence_lift": 0.0, "di_smooth_bass": 1.0, "di_smooth_lfhf": 1.0, "flat23": 1.0}
        scores = [Score(k, k, x, "dB", ("pass" if x >= 10 else "fail") if k == "cardioid_hard" else "pass",
                        (150.0, 300.0), "trusted", "", {}) for k, x in vals.items()]
        vt = VentTune(r, v, hard, 8.0, margin, {"0.7": (hard, 8.0), "1.3": (hard, 8.0)},
                      {"hard": hard, "aim": 8.0, "margin": margin - 1.0})
        return SimpleNamespace(scores=scores, tune=vt, tuned=tune, design=d, wall_s=0.5)

    def run_directivity(params_path, out_dir, *, trusted_limit=None, **kw):
        st.finals.append((Path(params_path), Path(out_dir), trusted_limit))
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        (Path(out_dir) / "report.html").write_text("<html></html>")
        return Path(out_dir) / "report.html"

    monkeypatch.setattr(round_mod, "evaluate", evaluate)
    monkeypatch.setattr(round_mod, "run_directivity", run_directivity)
    monkeypatch.setattr(round_mod, "code_version", lambda: "cv1")
    return st


def spec(tmp_path, text=ROUND):
    p = tmp_path / "round.toml"
    p.write_text(text)
    return p


def result(out, stage, name):
    return json.loads((out / stage / name / "result.json").read_text())


def test_round_writes_results_rankings_and_finals(tmp_path, fake_eval):
    out = tmp_path / "out"
    path = run_round(spec(tmp_path), out)
    assert path == out / "results.json"
    data = json.loads(path.read_text())
    assert data["spec"] == str(tmp_path / "round.toml") and len(data["spec_sha256"]) == 64
    assert [s["name"] for s in data["stages"]] == ["vents", "shape"]
    vents = data["stages"][0]
    assert [r["name"] for r in vents["results"]] == [
        "stage-base", "round0", "centre_from_front=140", WINNER]
    assert vents["ranking"]["order"][0] == WINNER
    assert json.loads((out / "vents" / "ranking.json").read_text()) == vents["ranking"]
    r = result(out, "vents", WINNER)
    assert r["status"] == "ok" and r["stage"] == "vents" and r["settings"] == "vents"
    assert r["code_version"] == "cv1" and r["commit"] and r["tuned"] is True
    assert r["overrides"] == {"vents.centre_from_front_mm": 190}
    assert {s["key"] for s in r["scores"]} >= {"cardioid_hard", "flat23"}
    assert r["tune"]["spread"]["0.7"] == [12.0, 8.0]
    for f in ("params.toml", "params_tuned.toml", "scores.json"):
        assert (out / "vents" / WINNER / f).exists()
    # round0 is the round's base, scored as it stands
    assert ("vents", "round0", "vents", False) in fake_eval.calls
    assert ("shape", "round0", "screen", False) in fake_eval.calls
    assert read_raw(out / "shape" / "round0" / "params.toml") == read_raw(round_mod.load_round(
        tmp_path / "round.toml").base)


def test_round_carries_the_best_vents_into_the_next_stage(tmp_path, fake_eval):
    # stage 2's stage-base params.toml has stage 1's winner's tuned rayl and litres
    out = tmp_path / "out"
    run_round(spec(tmp_path), out)
    raw = read_raw(out / "shape" / "stage-base" / "params.toml")
    r, v = fake_tune_of(0.190)
    assert raw["vents"]["specific_resistance_rayl"] == pytest.approx(r)
    assert raw["cabinet"]["coax_chamber_net_volume_l"] == pytest.approx(v * 1000.0)
    assert raw["vents"]["centre_from_front_mm"] == 190
    assert read_raw(out / "shape" / "m280" / "params.toml")["vents"]["specific_resistance_rayl"] == pytest.approx(r)
    assert read_raw(out / "shape" / "stage-base" / "params.toml") == read_raw(
        out / "vents" / WINNER / "params_tuned.toml")


def test_round_resumes_and_ignores_partial_results(tmp_path, fake_eval):
    # run once; delete one result.json and leave a stray result.json.<pid>.tmp; run again:
    # evaluate is called exactly once more, for that variant only
    out = tmp_path / "out"
    run_round(spec(tmp_path), out)
    first = json.loads((out / "results.json").read_text())
    (out / "shape" / "m280" / "result.json").unlink()
    (out / "shape" / "m280" / "result.json.12345.tmp").write_text('{"name": "m280", "status": "ok"')
    fake_eval.calls.clear()
    run_round(spec(tmp_path), out)
    assert fake_eval.calls == [("shape", "m280", "screen", True)]
    again = json.loads((out / "results.json").read_text())
    assert [s["ranking"] for s in again["stages"]] == [s["ranking"] for s in first["stages"]]


def test_round_reruns_a_truncated_result(tmp_path, fake_eval):
    out = tmp_path / "out"
    run_round(spec(tmp_path), out)
    p = out / "shape" / "m270" / "result.json"
    p.write_text(p.read_text()[:40])
    fake_eval.calls.clear()
    run_round(spec(tmp_path), out)
    assert fake_eval.calls == [("shape", "m270", "screen", True)]


def test_round_reruns_on_a_new_code_version(tmp_path, fake_eval, monkeypatch):
    out = tmp_path / "out"
    run_round(spec(tmp_path), out)
    monkeypatch.setattr(round_mod, "code_version", lambda: "cv2")
    fake_eval.calls.clear()
    run_round(spec(tmp_path), out)
    assert len(fake_eval.calls) == 8
    assert result(out, "vents", WINNER)["code_version"] == "cv2"


def test_round_reruns_a_variant_whose_params_changed(tmp_path, fake_eval):
    out = tmp_path / "out"
    run_round(spec(tmp_path), out)
    fake_eval.calls.clear()
    run_round(spec(tmp_path, ROUND.replace("= 280 }", "= 285 }")), out)
    assert fake_eval.calls == [("shape", "m280", "screen", True)]


def test_round_records_rejected_variants(tmp_path, fake_eval):
    # an override breaking rule (e) gives status "rejected" with the flare-mouth message; the others are ok
    text = ROUND + '[[stage.variant]]\nname = "m400"\nset = { "cabinet.flare_mouth_d_mm" = 400 }\n'
    out = tmp_path / "out"
    run_round(spec(tmp_path, text), out)
    r = result(out, "shape", "m400")
    assert r["status"] == "rejected" and "flare_mouth_d 400.0 mm" in r["reason"]
    assert r["scores"] == [] and r["tune"] is None
    assert not (out / "shape" / "m400" / "params_tuned.toml").exists()
    others = [result(out, "shape", n)["status"] for n in ("stage-base", "round0", "m270", "m280")]
    assert others == ["ok"] * 4
    ranking = json.loads((out / "shape" / "ranking.json").read_text())
    assert ranking["out"]["m400"].startswith("rejected: flat front width")
    assert "m400" not in ranking["order"]


def test_round_records_errors_and_continues(tmp_path, fake_eval, capsys):
    fake_eval.fail = {"m280"}
    out = tmp_path / "out"
    run_round(spec(tmp_path), out)
    r = result(out, "shape", "m280")
    assert r["status"] == "error"
    assert r["reason"].rstrip().endswith("RuntimeError: solver blew up")
    assert len(r["reason"].splitlines()) <= 20
    ranking = json.loads((out / "shape" / "ranking.json").read_text())
    assert ranking["out"]["m280"] == "error: RuntimeError: solver blew up"
    assert ranking["order"] == ["stage-base", "round0"]
    assert [f["name"] for f in json.loads((out / "finals.json").read_text())] == ["stage-base"]
    err = capsys.readouterr().err
    assert "[shape 4/4] m280: error, " in err and "[vents 1/4] stage-base: ok, 0.5 s" in err


def test_finals_are_the_top_of_the_last_stage(tmp_path, fake_eval):
    out = tmp_path / "out"
    run_round(spec(tmp_path), out)
    ranking = json.loads((out / "shape" / "ranking.json").read_text())
    assert ranking["order"] == ["m280", "stage-base", "round0"] and not ranking["gate_unmet"]
    assert ranking["out"] == {"m270": "fails the hard cardioid gate (8.0 dB)"}
    finals = json.loads((out / "finals.json").read_text())
    assert [(f["rank"], f["name"], f["status"]) for f in finals] == [(1, "m280", "ok"), (2, "stage-base", "ok")]
    assert finals[0]["dir"] == "finals/m280" and finals[0]["params"] == "shape/m280/params_tuned.toml"
    tl = (2851.0, "inherited from round 0's mesh check")
    assert fake_eval.finals == [
        (out / "shape" / "m280" / "params_tuned.toml", out / "finals" / "m280", tl),
        (out / "shape" / "stage-base" / "params_tuned.toml", out / "finals" / "stage-base", tl)]
    assert json.loads((out / "results.json").read_text())["finals"] == finals


def test_a_failed_final_is_recorded_and_the_round_finishes(tmp_path, fake_eval, monkeypatch):
    def boom(params_path, out_dir, *, trusted_limit=None, **kw):
        raise RuntimeError("out of memory")
    monkeypatch.setattr(round_mod, "run_directivity", boom)
    out = tmp_path / "out"
    run_round(spec(tmp_path), out)
    finals = json.loads((out / "finals.json").read_text())
    assert [f["status"] for f in finals] == ["error", "error"]
    assert finals[0]["reason"].rstrip().endswith("RuntimeError: out of memory")
    assert (out / "results.json").exists()


def test_round_stops_when_a_stage_has_no_best(tmp_path, fake_eval):
    fake_eval.fail = {"stage-base", "centre_from_front=140", WINNER}
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="stage 'vents'"):
        run_round(spec(tmp_path), out)
    ranking = json.loads((out / "vents" / "ranking.json").read_text())
    assert ranking["order"] == ["round0"]
    assert not (out / "shape").exists() and not (out / "results.json").exists()


def test_round_out_defaults_to_out_rounds_number(tmp_path, fake_eval, monkeypatch):
    s = spec(tmp_path)
    monkeypatch.chdir(tmp_path)
    assert run_round(s) == Path("out/rounds/7/results.json")
    assert (tmp_path / "out/rounds/7/finals.json").exists()


def test_cli_round(tmp_path, fake_eval, capsys):
    out = tmp_path / "o"
    assert main(["round", str(spec(tmp_path)), "--out", str(out)]) == 0
    assert f"wrote {out / 'results.json'}" in capsys.readouterr().out


def test_commit_is_head_or_unknown():
    c = round_mod.git_commit()
    assert c == "unknown" or len(c.removesuffix("-dirty")) == 40


def test_round_reruns_an_ok_result_without_its_tuned_params(tmp_path, fake_eval):
    out = tmp_path / "out"
    run_round(spec(tmp_path), out)
    (out / "vents" / WINNER / "params_tuned.toml").unlink()
    fake_eval.calls.clear()
    run_round(spec(tmp_path), out)
    assert fake_eval.calls == [("vents", WINNER, "vents", True)]
