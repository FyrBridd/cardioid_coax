import numpy as np
import pytest
from fakes import monopole_field, monopole_puck

import speaker.evaluate as ev_mod
from speaker.bem.run import freq_grid
from speaker.evaluate import SCREEN, SETTINGS, VENTS, EvalSettings, evaluate
from speaker.params import DEFAULT_PATH, load_params


@pytest.fixture(autouse=True)
def fake_fields(monkeypatch):
    """The fake fields in evaluate.py's namespace; returns the puck calls."""
    puck_calls = []

    def cab(p, coax, bass, f, *, per_wavelength, h_max, cache_dir, **kw):
        return monopole_field(p, coax, bass, f)

    def puck(p, coax, f, *, per_wavelength, cache_dir):
        puck_calls.append(np.asarray(f, float))
        return monopole_puck(p, coax, f)

    def recess(coax, f=None, *, per_wavelength=6.0, cache_dir=None):
        f = freq_grid(300.0, 4000.0, 6) if f is None else np.asarray(f, float)
        return f, np.ones(len(f), complex)

    monkeypatch.setattr(ev_mod, "cabinet_field", cab)
    monkeypatch.setattr(ev_mod, "puck_field", puck)
    monkeypatch.setattr(ev_mod, "recess_ratio", recess)
    return puck_calls


def test_presets():
    assert SETTINGS == {"vents": VENTS, "screen": SCREEN}
    assert VENTS == EvalSettings("vents", 120.0, 600.0, 12, 6.0, 0.05, None, True)
    assert SCREEN == EvalSettings("screen", 120.0, None, 12, 5.0, 0.05, 4000.0, True)


def test_screen_scores_are_all_trusted(tmp_path, fake_fields):
    ev = evaluate(DEFAULT_PATH, SCREEN, trusted_limit=2851.0, tune=False, cache_dir=tmp_path)
    assert 2851.0 / 2 ** (1 / 12) < ev.f[-1] <= 2851.0
    assert not ev.tuned and ev.design.p == load_params(DEFAULT_PATH)
    assert {s.confidence for s in ev.scores} == {"trusted"}
    assert any(s.key == "presence_lift" for s in ev.scores)
    # The puck runs 500 Hz to puck_top, for the HF reference only: no cabinet-to-puck blend.
    fp = np.concatenate(fake_fields)
    assert fp[0] == 500.0 and fp[-1] == 4000.0
    assert ev.ways.meta["puck_blend"] is None and ev.ways.meta["hf_reference"] == "puck"
    assert ev.ways.meta["lf_reference"] == "recess"
    assert ev.trusted_limit == 2851.0 and ev.wall_s > 0
    assert {b["field"] for b in ev.bands} == {"cabinet", "puck"}
    # Untuned, the tune is still run, and its base margin is the params' own R and V.
    assert ev.tune.base["margin"] <= ev.tune.margin + 1e-9


def test_tuned_design_carries_the_tuned_vents(tmp_path):
    ev = evaluate(DEFAULT_PATH, VENTS, trusted_limit=2851.0, tune=True, cache_dir=tmp_path)
    assert ev.tuned
    assert ev.design.p.vents.specific_resistance == ev.tune.resistance
    assert ev.design.p.cabinet.coax_chamber_net_volume == ev.tune.net_volume
    hard = next(s.value for s in ev.scores if s.key == "cardioid_hard")
    assert hard == pytest.approx(ev.tune.hard, abs=1e-6)


def test_vents_preset_has_no_puck_and_stops_at_600(tmp_path, fake_fields):
    ev = evaluate(DEFAULT_PATH, VENTS, trusted_limit=2851.0, tune=False, cache_dir=tmp_path)
    assert fake_fields == []
    assert ev.f[0] == 120.0 and 600.0 / 2 ** (1 / 12) < ev.f[-1] <= 600.0
    assert {150.0, 300.0, 500.0} <= set(ev.f)
    assert ev.ways.meta["hf_reference"] != "puck"
    assert {s.confidence for s in ev.scores} == {"trusted"}
