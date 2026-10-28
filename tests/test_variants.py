import pytest

from speaker.params import DEFAULT_PATH, load_params
from speaker.variants import RESERVED, apply_overrides, load_round, read_raw, with_tuned, write_params

ROUND = '''
number = 1
title = "Cardioid first, then the front shape"
base = "params/default.toml"
trusted_limit_hz = 2851
trusted_note = "inherited from round 0's mesh check (2026-10-27)"
finals = 2

[[compare]]
label = "Round 0 (round 0b run)"
scores = "out/round0b/scores.json"

[[stage]]
name = "vents"
settings = "vents"
rank_by = "cardioid"
include_round0 = true
grid = { "vents.centre_from_front_mm" = [140, 190, 240, 280], "vents.slots_per_side" = [6, 9, 12] }

[[stage]]
name = "shape"
settings = "screen"
rank_by = "directivity"
include_round0 = true
[[stage.variant]]
name = "w330-ro30-m260"
set = { "cabinet.width_mm" = 330, "cabinet.front_edge" = "roundover", "cabinet.edge_size_mm" = 30, "cabinet.flare_mouth_d_mm" = 260 }
'''


def write(tmp_path, text):
    p = tmp_path / "round.toml"
    p.write_text(text)
    return p


def test_loads_fields(tmp_path):
    r = load_round(write(tmp_path, ROUND))
    assert r.number == 1 and r.finals == 2 and r.trusted_limit == 2851.0
    assert r.base.name == "default.toml" and r.base.is_absolute()
    assert r.compare[0].label.startswith("Round 0")
    assert [s.name for s in r.stages] == ["vents", "shape"]
    assert r.stages[1].variants[0].overrides["cabinet.width_mm"] == 330
    assert r.stages[1].include_round0 is True


def test_unknown_override_key_names_it(tmp_path):
    spec = ROUND.replace('"cabinet.width_mm" = 330', '"cabinet.widht_mm" = 330')
    with pytest.raises(ValueError, match=r"cabinet\.widht_mm.*w330-ro30-m260"):
        load_round(write(tmp_path, spec))


def test_grid_expands_in_sorted_key_order(tmp_path):
    r = load_round(write(tmp_path, ROUND))
    names = [v.name for v in r.stages[0].variants]
    assert len(names) == 12 and names[0] == "centre_from_front=140-slots_per_side=6"
    assert names[1] == "centre_from_front=140-slots_per_side=9"


def test_explicit_variants_follow_grid(tmp_path):
    spec = ROUND.replace(
        '[[stage]]\nname = "shape"',
        '[[stage.variant]]\nname = "x"\nset = { "vents.slots_per_side" = 3 }\n\n[[stage]]\nname = "shape"')
    r = load_round(write(tmp_path, spec))
    names = [v.name for v in r.stages[0].variants]
    assert len(names) == 13 and names[-1] == "x"


@pytest.mark.parametrize("old,new,msg", [
    ('settings = "vents"', 'settings = "nope"', "nope"),
    ('rank_by = "cardioid"', 'rank_by = "loud"', "loud"),
    ('finals = 2', 'finals = -1', "finals"),
    ('name = "w330-ro30-m260"', 'name = "round0"', "round0"),
    ('name = "w330-ro30-m260"', 'name = "stage-base"', "stage-base"),
    ('name = "w330-ro30-m260"', 'name = "bad name!"', "bad name"),
])
def test_rejections(tmp_path, old, new, msg):
    with pytest.raises(ValueError, match=msg):
        load_round(write(tmp_path, ROUND.replace(old, new, 1)))


def test_duplicate_name_rejected(tmp_path):
    spec = ROUND + '[[stage.variant]]\nname = "w330-ro30-m260"\nset = {}\n'
    with pytest.raises(ValueError, match="duplicate.*w330-ro30-m260"):
        load_round(write(tmp_path, spec))


def test_stage_without_variants_rejected(tmp_path):
    spec = ROUND.replace('grid = { "vents.centre_from_front_mm" = [140, 190, 240, 280], "vents.slots_per_side" = [6, 9, 12] }', '')
    with pytest.raises(ValueError, match="vents.*no variants"):
        load_round(write(tmp_path, spec))


def test_reserved():
    assert RESERVED == ("stage-base", "round0")


def test_written_variant_loads_with_the_override(tmp_path):
    raw = apply_overrides(read_raw(DEFAULT_PATH), {"cabinet.width_mm": 330, "cabinet.front_edge": "chamfer"})
    write_params(raw, tmp_path / "v.toml")
    p = load_params(tmp_path / "v.toml")
    assert p.cabinet.width == pytest.approx(0.330) and p.cabinet.front_edge == "chamfer"


def test_apply_overrides_deep_copies():
    raw = read_raw(DEFAULT_PATH)
    w = raw["cabinet"]["width_mm"]
    apply_overrides(raw, {"cabinet.width_mm": w + 1})
    assert raw["cabinet"]["width_mm"] == w


def test_with_tuned_writes_rayl_and_litres(tmp_path):
    raw = with_tuned(read_raw(DEFAULT_PATH), 123.0, 0.0042)
    assert raw["vents"]["specific_resistance_rayl"] == 123.0
    assert raw["cabinet"]["coax_chamber_net_volume_l"] == pytest.approx(4.2)
    write_params(raw, tmp_path / "t.toml")
    p = load_params(tmp_path / "t.toml")
    assert p.vents.specific_resistance == pytest.approx(123.0)
    assert p.cabinet.coax_chamber_net_volume == pytest.approx(0.0042)
