from speaker.ranking import Ranking, VariantResult, rank


def score(key, value, status="pass", confidence="trusted"):
    return {"key": key, "label": key, "value": value, "unit": "dB", "status": status,
            "band": [150.0, 300.0], "confidence": confidence, "formula": "", "detail": {}}


def res(name, *, hard=12.0, ws=1.0, pl=1.0, dib=1.0, dil=1.0, flat=1.0, status="ok", reason="",
        margin=1.0, base_margin=0.0, tuned=True, drop=(), low=()):
    vals = {"cardioid_hard": hard, "window_spread": ws, "presence_lift": pl,
            "di_smooth_bass": dib, "di_smooth_lfhf": dil, "flat23": flat}
    scores = []
    for k, v in vals.items():
        st = ("pass" if v >= 10.0 else "fail") if k == "cardioid_hard" else "pass"
        if k in low:
            scores.append(score(k, v, st, "low"))
        elif k not in drop:
            scores.append(score(k, v, st))
            scores.append(score(k, v + 100.0, st, "low"))  # the low part must never be read
    tune = {"resistance": 200.0, "net_volume": 0.008, "hard": hard, "aim": 8.0, "margin": margin,
            "spread": {}, "base": {"hard": hard, "aim": 8.0, "margin": base_margin}}
    return VariantResult(name=name, stage="s", status=status, reason=reason, overrides={},
                         params_sha256="x", settings="screen", code_version="c", commit="g",
                         scores=scores if status == "ok" else [], tune=tune if status == "ok" else None,
                         tuned=tuned, wall_s=1.0)


def test_rank_excludes_gate_failures():
    r = rank([res("a", hard=12.0, ws=3, pl=3), res("b", hard=8.04, ws=0, pl=0)], "directivity")
    assert isinstance(r, Ranking)
    assert r.order == ["a"] and not r.gate_unmet
    assert r.out == {"b": "fails the hard cardioid gate (8.0 dB)"}


def test_rank_when_every_variant_fails_the_gate():
    r = rank([res("a", hard=6.0, ws=2, pl=2), res("b", hard=8.0, ws=3, pl=3)], "directivity")
    assert r.gate_unmet and r.order == ["b", "a"]
    assert r.out == {}


def test_rank_gate_unmet_ties_on_hard_go_to_the_key():
    r = rank([res("a", hard=6.0, ws=2, pl=2), res("b", hard=6.0, ws=1, pl=1)], "directivity")
    assert r.gate_unmet and r.order == ["b", "a"]


def test_rank_bucket_ties_go_to_di_smoothness():
    # window+presence 4.05 vs 4.20 (same 0.25 dB bucket); the lower DI sum wins
    a = res("a", ws=2.0, pl=2.05, dib=1.5, dil=1.5)
    b = res("b", ws=2.0, pl=2.20, dib=1.0, dil=1.0)
    c = res("c", ws=2.0, pl=1.90, dib=5.0, dil=5.0)  # 3.90: a lower bucket beats any DI sum
    r = rank([a, b, c], "directivity")
    assert r.order == ["c", "b", "a"]


def test_rank_flat23_breaks_the_last_tie():
    r = rank([res("a", flat=2.0), res("b", flat=1.0)], "directivity")
    assert r.order == ["b", "a"]


def test_rank_reads_trusted_scores_only():
    r = rank([res("a", flat=2.0), res("b", flat=1.0, low=("flat23",))], "directivity")
    assert r.order == ["a"] and r.out == {"b": "missing trusted score flat23"}


def test_rank_missing_trusted_score_is_out():
    r = rank([res("a"), res("b", drop=("cardioid_hard",))], "directivity")
    assert r.order == ["a"] and r.out == {"b": "missing trusted score cardioid_hard"}


def test_rank_non_ok_is_out_with_its_reason():
    rows = [res("a"), res("b", status="rejected", reason="flare too wide"),
            res("c", status="error", reason="Traceback\n  File x\nRuntimeError: boom\n")]
    for by in ("directivity", "cardioid"):
        r = rank(rows, by)
        assert r.order == ["a"]
        assert r.out == {"b": "rejected: flare too wide", "c": "error: RuntimeError: boom"}


def test_rank_cardioid_by_tuned_or_base_margin():
    rows = [res("a", margin=1.0), res("b", margin=3.0),
            res("r0", margin=9.0, base_margin=2.0, tuned=False)]  # untuned: its base margin counts
    r = rank(rows, "cardioid")
    assert r.order == ["b", "r0", "a"] and r.out == {} and not r.gate_unmet


def test_rank_cardioid_ignores_the_gate():
    r = rank([res("a", hard=5.0, margin=-5.0), res("b", hard=6.0, margin=-4.0)], "cardioid")
    assert r.order == ["b", "a"] and not r.gate_unmet
