import numpy as np
import pytest
from fakes import monopole_field

from speaker.crossover import FixedCrossover
from speaker.layout import compute_layout
from speaker.scores import score_all
from speaker.spinorama import spinorama
from speaker.system import coax_section, way_responses
from speaker.vent_tune import feasible_volumes, tune_vents, with_vents

F_120_600 = np.geomspace(120.0, 600.0, 24)


def test_tuned_margin_not_below_the_base(design):
    cab = monopole_field(*design, F_120_600)
    t = tune_vents(*design, cab)
    assert t.margin >= t.base["margin"] - 1e-9


def test_tuned_hard_matches_the_score(design):
    cab = monopole_field(*design, F_120_600)
    t = tune_vents(*design, cab)
    q = with_vents(design.p, t.resistance, t.net_volume)
    ways = way_responses(q, design.coax, design.bass, compute_layout(q, design.coax, design.bass), cab, None,
                         F_120_600, trusted_limit=600.0)
    sec = coax_section(ways)  # also passed as "total": only the cardioid rows are compared
    sc = score_all(F_120_600, sec, sec, cab.obs, spinorama(sec, cab.obs, f=F_120_600),
                   trusted_limit=600.0, xo=FixedCrossover())
    assert next(s.value for s in sc if s.key == "cardioid_hard") == pytest.approx(t.hard, abs=1e-6)
    assert next(s.value for s in sc if s.key == "cardioid_aim") == pytest.approx(t.aim, abs=1e-6)


def test_feasible_volumes_all_lay_out(design):
    for v in feasible_volumes(*design):
        compute_layout(with_vents(design.p, design.p.vents.specific_resistance, v), design.coax, design.bass)


def test_tune_needs_the_cardioid_band(design):
    with pytest.raises(ValueError, match="150"):
        tune_vents(*design, monopole_field(*design, np.geomspace(200, 600, 10)))
