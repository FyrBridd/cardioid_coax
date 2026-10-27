import math
from pathlib import Path
import pytest
from speaker.params import load_params

def test_defaults_load_in_si():
    p = load_params()
    assert p.cabinet.width == pytest.approx(0.320)
    assert p.cabinet.depth == pytest.approx(0.350)
    assert p.cabinet.coax_centre_z == pytest.approx(0.950)
    assert p.cabinet.coax_chamber_net_volume == pytest.approx(0.010)
    assert p.vents.slots_per_side == 12
    assert p.vents.specific_resistance == pytest.approx(190.0)
    assert p.room.seat == pytest.approx((3.425, 2.210, 0.950))
    assert p.drivers_measured is False

def test_depth_over_350_rejected(tmp_path):
    text = Path("params/default.toml").read_text().replace("depth_mm = 350", "depth_mm = 351")
    f = tmp_path / "p.toml"; f.write_text(text)
    with pytest.raises(ValueError, match="depth"):
        load_params(f)

def test_facet_angle_converted_to_radians(tmp_path):
    text = Path("params/default.toml").read_text().replace("facet_angle_deg = 0", "facet_angle_deg = 20")
    f = tmp_path / "p.toml"; f.write_text(text)
    c = load_params(f).cabinet
    assert c.facet_angle == pytest.approx(math.radians(20))

def test_front_shape_defaults():
    c = load_params().cabinet
    assert c.outer_baffle_t == pytest.approx(0.018) and c.baffle_t == pytest.approx(0.036)
    assert c.front_edge == "square" and c.flare_profile == "ellipse"
    assert (c.edge_size, c.facet_angle, c.facet_width, c.facet_base_z, c.facet_ramp, c.coax_offset_x) == (0,) * 6
    assert not hasattr(c, "side_edge_radius")
