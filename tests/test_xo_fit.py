import functools

import numpy as np
import pytest

from speaker.bem.run import observation_set
from speaker.crossover import Peq, peq
from speaker.system import Ways
from speaker.xo_fit import (
    BASS_COAX_GRID, EQ_TOL_DB, GAIN_RANGE, LF_HF_GRID, N_PEQ, Q_RANGE, eq_band, fit_active, fit_at, fit_peq,
)


def test_fit_removes_a_single_bump():
    f = np.geomspace(100, 10000, 160)
    bump = Peq(1500.0, 5.0, 2.0)
    db = 80 + 20 * np.log10(np.abs(peq(f, bump)))
    eq = fit_peq(f, db, 80.0, (500.0, 5000.0), 3)
    after = db + sum(20 * np.log10(np.abs(peq(f, p))) for p in eq)
    m = (f >= 500) & (f <= 5000)
    assert np.max(np.abs(after[m] - 80)) < 0.5
    assert any(abs(p.f0 - 1500) < 150 and p.gain_db == pytest.approx(-5.0, abs=0.7) for p in eq)


def test_boost_and_cut_stay_in_range():
    f = np.geomspace(100, 10000, 160)
    db = np.where((f > 1000) & (f < 2000), 50.0, 80.0)  # a 30 dB hole
    for p in fit_peq(f, db, 80.0, (500.0, 5000.0), 4):
        assert GAIN_RANGE[0] <= p.gain_db <= GAIN_RANGE[1] and Q_RANGE[0] <= p.q <= Q_RANGE[1]


def test_short_band_gets_no_eq():
    f = np.geomspace(100, 2851, 60)
    assert fit_peq(f, np.full(len(f), 70.0), 80.0, (2600.0, 16000.0), 5) == ()


def test_eq_bands_never_reach_below_a_high_pass():
    assert eq_band("coax_lf", 170, 1400) == (170, 2800)
    assert eq_band("coax_hf", 170, 1400) == (1400, 16000.0)
    assert eq_band("bass", 170, 1400) == (80.0, 340)


def test_a_single_bump_needs_one_filter():
    f = np.geomspace(100, 10000, 160)
    db = 80 + 20 * np.log10(np.abs(peq(f, Peq(1500.0, 5.0, 2.0))))
    eq = fit_peq(f, db, 80.0, (500.0, 5000.0), 3)
    assert len(eq) == 1
    assert abs(eq[0].f0 - 1500) < 150 and eq[0].gain_db == pytest.approx(-5.0, abs=0.7)


def omni_ways():
    o = observation_set(2.0)
    f = np.geomspace(50, 16000, 200)
    return Ways(f, o, {k: np.ones((len(f), len(o.dirs)), complex) for k in ("bass", "coax_lf", "coax_hf")}, {})


def test_alignment_recovers_a_known_delay_and_inverted_polarity():
    # the LR4 sections have phase of their own, so pin the delay added on top of the flat fit
    xo0 = fit_at(omni_ways(), 16000.0, 170.0, 1400.0)
    d0 = xo0.ways["coax_lf"].delay_s - xo0.ways["coax_hf"].delay_s
    w = omni_ways()
    w.ways["coax_hf"] = -w.ways["coax_hf"] * np.exp(-2j * np.pi * w.f * 0.3e-3)[:, None]
    xo1 = fit_at(w, 16000.0, 170.0, 1400.0)
    d1 = xo1.ways["coax_lf"].delay_s - xo1.ways["coax_hf"].delay_s
    assert xo0.ways["coax_hf"].polarity == 1 and xo1.ways["coax_hf"].polarity == -1
    assert d1 - d0 == pytest.approx(0.3e-3, abs=0.02e-3)
    assert min(x.delay_s for x in xo1.ways.values()) == 0.0


def test_flat_omni_ways_need_no_eq():
    xo = fit_at(omni_ways(), 16000.0, 170.0, 1400.0)
    assert all(x.peq == () for x in xo.ways.values())


def test_search_keys_include_the_whole_system_rear_rejection():
    from speaker.xo_fit import SEARCH_KEYS
    assert "whole_rear" in SEARCH_KEYS


def test_search_stays_in_the_spec_ranges_and_returns_a_grid_point():
    xo = fit_active(omni_ways(), 2851.0)
    assert xo.bass_coax in BASS_COAX_GRID and xo.lf_hf in LF_HF_GRID


def test_directivity_order_is_public():
    from speaker.ranking import directivity_order
    v = dict(window_spread=1.0, presence_lift=1.0, di_smooth_bass=0.5, di_smooth_lfhf=0.5, flat23=2.0)
    rows = [("x", v | {"presence_lift": 1.2}), ("y", v)]
    assert directivity_order(rows, lambda r: r[1]) == [rows[1], rows[0]]


def test_search_floor_keeps_lf_hf_at_or_above_the_driver_recommendation():
    assert _fitted(TL, "").lf_hf == 1600.0  # unfloored
    assert fit_active(featured_ways(TL), TL, lf_hf_min=1800.0).lf_hf >= 1800.0
    with pytest.raises(ValueError, match="2500 Hz"):
        fit_active(omni_ways(), TL, lf_hf_min=2500.0)


def _cand(name, rear, ws):
    return name, dict(window_spread=ws, presence_lift=1.0, di_smooth_bass=0.5, di_smooth_lfhf=0.5, flat23=1.0,
                      whole_rear=rear)


def test_pick_keeps_candidates_within_a_db_of_the_best_rear_rejection():
    from speaker.xo_fit import WHOLE_REAR_TOL_DB, pick
    assert WHOLE_REAR_TOL_DB == 1.0
    # a has the best directivity but 4.2 dB of rear rejection, 5.7 dB below b; c is within 1 dB of b
    # and beats it on directivity
    cands = [_cand("a", 4.2, 0.5), _cand("b", 9.9, 2.0), _cand("c", 9.0, 1.5), _cand("d", 8.8, 1.0)]
    assert pick(cands) == "c"


def test_passive_warnings():
    from speaker.crossover import ActiveCrossover, WayFilter
    from speaker.xo_fit import passive_warnings
    eq = (Peq(3000.0, 8.0, 1.0),)
    xo = ActiveCrossover(170.0, 1200.0, {"bass": WayFilter(0, 0, 1, ()), "coax_lf": WayFilter(0, 0.4e-3, 1, ()),
                                         "coax_hf": WayFilter(-5, 0.0, 1, eq)}, 20000.0)
    w = " ".join(passive_warnings(xo, 2851.0))
    assert "coax_hf" in w and "8.0 dB" in w      # boost over 6 dB
    assert "750 Hz" in w                          # 1.2 kHz is within an octave of the HF resonance
    assert "0.40 ms" in w                         # LF-HF delay a passive network can't add
    assert "3000 Hz" in w                         # PEQ above the trusted limit
    flat = {k: WayFilter(0.0, 0.0, 1, ()) for k in ("bass", "coax_lf", "coax_hf")}
    assert passive_warnings(ActiveCrossover(170.0, 1600.0, flat, 20000.0), 2851.0) == []


TL = 2851.0


def featured_ways(f_top: float = 16000.0, above: str = "") -> Ways:
    """Omni ways with a response feature in every way; above: "shelf" (+10 dB) or "notch" (-15 dB at 8 kHz)
    on the coax ways above TL only."""
    o = observation_set(2.0)
    f = np.geomspace(50, f_top, int(round(12 * np.log2(f_top / 50))) + 1)

    def shape(*ps):
        return np.prod([peq(f, p) for p in ps], axis=0)
    extra = np.ones(len(f))
    if above == "shelf":
        extra = np.where(f > TL, 10 ** (10 / 20), 1.0)
    elif above == "notch":
        extra = np.where(f > TL, np.abs(peq(f, Peq(8000.0, -15.0, 4.0))), 1.0)
    ways = {"bass": shape(Peq(120.0, 4.0, 2.0)) * 2.0,
            "coax_lf": shape(Peq(600.0, 3.0, 2.0), Peq(1800.0, -4.0, 3.0)) * extra,
            "coax_hf": shape(Peq(2300.0, -5.0, 2.0), Peq(6000.0, 6.0, 1.5)) * 0.5 * extra}
    return Ways(f, o, {k: v[:, None] * np.ones((1, len(o.dirs)), complex) for k, v in ways.items()}, {})


@functools.lru_cache(maxsize=None)
def _fitted(f_top: float, above: str):
    return fit_active(featured_ways(f_top, above), TL)


def test_a_fit_on_data_that_stop_at_the_trusted_limit_is_unchanged():
    # pinned from the scoring version 3 search (fixed DI bands, bass from 120 Hz, whole_rear kept, tiers within 0.25 dB);
    # the screen grid stops at the trusted limit, so this is the fit every screen gets
    xo = _fitted(TL, "")
    assert (xo.bass_coax, xo.lf_hf, xo.fit_top) == (170.0, 1600.0, pytest.approx(TL))
    want = {
        "bass": (-7.405096501533659, 0.15e-3, 1, [(120.54775165837202, -12.988997629249576, 0.8843117411617827),
                                                  (120.64398094464696, 11.999999999999998, 0.5000000000000001)]),
        "coax_lf": (0.0, 0.0, 1, [(226.54738993865544, 1.3398126562815094, 0.5000000000000001),
                                  (590.1342936661375, -2.3503634702712537, 2.539103005796063),
                                  (1813.2805195654096, 4.647246916340044, 1.6366165514562367)]),
        "coax_hf": (10.961149650479925, 0.01e-3, 1, [(1600.0000000000014, -2.3315138960690103, 2.5216065096168148),
                                                     (2271.1659110319065, 1.58066905241805, 3.1611612583276827),
                                                     (2809.3328203508295, -1.511420260168743, 5.168117698023399)]),
    }
    for k, (g, d, pol, eq) in want.items():
        wf = xo.ways[k]
        assert wf.gain_db == pytest.approx(g, abs=1e-9) and wf.delay_s == pytest.approx(d, abs=1e-12)
        assert wf.polarity == pol
        assert [(p.f0, p.gain_db, p.q) for p in wf.peq] == [pytest.approx(e, rel=1e-6) for e in eq]


def _trusted_part(xo):
    return (xo.bass_coax, xo.lf_hf,
            {k: (v.gain_db, v.delay_s, v.polarity, tuple(p for p in v.peq if p.f0 <= TL)) for k, v in xo.ways.items()})


@pytest.mark.parametrize("above", ["shelf", "notch"])
def test_data_above_the_trusted_limit_do_not_move_the_crossover(above):
    """Points, gains, delays, polarities and trusted-band EQ come from the trusted band only."""
    assert _trusted_part(_fitted(16000.0, above)) == _trusted_part(_fitted(16000.0, ""))
    w = featured_ways(16000.0)
    m = w.f <= TL  # the same grid, stopping at TL
    stop = fit_active(Ways(w.f[m], w.obs, {k: v[m] for k, v in w.ways.items()}, {}), TL)
    assert _trusted_part(_fitted(16000.0, "")) == _trusted_part(stop) and stop.fit_top <= TL


@pytest.mark.parametrize("above", ["", "shelf", "notch"])
def test_above_limit_peqs_stay_above_it_and_leave_the_trusted_band(above):
    xo = _fitted(16000.0, above)
    f = featured_ways(16000.0).f
    assert xo.fit_top == pytest.approx(16000.0)
    hi = {k: [p for p in v.peq if p.f0 > TL] for k, v in xo.ways.items()}
    assert hi["coax_hf"]  # the HF way's features above TL get EQ
    assert not hi["bass"]
    for k, ps in hi.items():
        assert len(xo.ways[k].peq) <= N_PEQ[k]
        top = eq_band(k, xo.bass_coax, xo.lf_hf)[1]
        assert all(TL < p.f0 <= top for p in ps)
        if ps:
            moved = sum(20 * np.log10(np.abs(peq(f[f <= TL], p))) for p in ps)
            assert np.max(np.abs(moved)) < EQ_TOL_DB


def _hf_peq_xo(hf_peq=(Peq(12800.0, 12.0, 0.65),)):
    from speaker.crossover import ActiveCrossover, WayFilter
    return ActiveCrossover(170.0, 1400.0, {"bass": WayFilter(4.0, 1.25e-3, -1, ()),
                                           "coax_lf": WayFilter(0.0, 0.0, 1, (Peq(600.0, -15.0, 2.0),)),
                                           "coax_hf": WayFilter(-3.0, 0.0, 1, hf_peq)}, 20000.0)


def test_biquad_deviation_reports_a_cramped_hf_peq_and_recommends_96k():
    from speaker.crossover import biquad_deviation
    from speaker.xo_fit import crossover_text
    xo = _hf_peq_xo()
    d48, d96 = biquad_deviation(xo, 48000.0, TL), biquad_deviation(xo, 96000.0, TL)
    assert set(d48) == {"bass", "coax_lf", "coax_hf"}
    assert set(d48["coax_hf"]) == {"trusted", "2-5 kHz", "above 5 kHz"}
    assert d48["coax_hf"]["2-5 kHz"] > 0.25
    for band in ("trusted", "2-5 kHz", "above 5 kHz"):
        assert d96["coax_hf"][band] < d48["coax_hf"][band]
    assert np.isnan(d48["bass"]["above 5 kHz"])  # the bass is more than 20 dB down there
    txt = crossover_text(xo, TL, design="d.toml", sha="0123abcd")
    assert "use 96 kHz" in txt
    assert "# coax_hf\t" + "\t".join(f"{d48['coax_hf'][b]:.2f}" for b in ("trusted", "2-5 kHz", "above 5 kHz")) in txt
    # without the HF PEQ both rates follow the analog ways closely, so no recommendation
    assert "use 96 kHz" not in crossover_text(_hf_peq_xo(()), TL, design="d.toml", sha="0123abcd")


def test_crossover_text_header_and_limit_flags():
    from speaker.xo_fit import crossover_text
    txt = crossover_text(_hf_peq_xo(), TL, design="round-3.toml", sha="0123abcd")
    head = "\n".join(l for l in txt.splitlines() if l.startswith("#"))
    assert "round-3.toml" in head and "sha256 0123abcd" in head and "trusted limit 2851 Hz" in head
    assert "normalise" in head and "+9.00 dB" in head  # coax_hf: -3 dB gain plus the +12 dB PEQ
    assert "-a1" in head and "-a2" in head
    rows = [l.split("\t") for l in txt.splitlines() if "\tpeq\t" in l]
    hf = next(r for r in rows if r[0] == "coax_hf")
    lf = next(r for r in rows if r[0] == "coax_lf")
    assert hf[5:] == ["low", "at the fitter's limit"]   # +12 dB, above TL
    assert lf[5:] == ["at the fitter's limit"]          # -15 dB, trusted
