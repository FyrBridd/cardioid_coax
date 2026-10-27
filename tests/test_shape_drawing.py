from pathlib import Path

import pytest

from speaker.params import load_params
from speaker.shape_drawing import draw_designs
from speaker.variants import apply_overrides, read_raw, write_params

DEFAULT = Path(__file__).resolve().parents[1] / "params" / "default.toml"

FAMILIES = {
    "square": {},
    "roundover": {"cabinet.front_edge": "roundover", "cabinet.edge_size_mm": 30},
    "chamfer": {"cabinet.front_edge": "chamfer", "cabinet.edge_size_mm": 25},
    "facets": {"cabinet.front_edge": "facets", "cabinet.facet_angle_deg": 30, "cabinet.facet_width_mm": 40,
               "cabinet.facet_base_z_mm": 660, "cabinet.facet_ramp_mm": 100},
}


def params_for(tmp_path, name):
    p = tmp_path / f"{name}.toml"
    write_params(apply_overrides(read_raw(DEFAULT), {"cabinet.width_mm": 400, **FAMILIES[name]}), p)
    return load_params(p)


def test_draws_each_family(tmp_path):
    designs = [(name, params_for(tmp_path, name)) for name in FAMILIES]
    out = tmp_path / "shape.png"
    draw_designs(out, designs)
    assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n" and out.stat().st_size > 5000


@pytest.mark.parametrize("name", list(FAMILIES))
def test_draws_one_family_against_round0(tmp_path, name):
    out = tmp_path / "s.png"
    draw_designs(out, [(name, params_for(tmp_path, name)), ("round 0", load_params(DEFAULT))])
    assert out.stat().st_size > 5000
