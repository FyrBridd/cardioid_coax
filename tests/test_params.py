from pathlib import Path
import pytest
from speaker.params import load_params

def test_defaults_load_in_si():
    p = load_params()
    assert p.cabinet.width == pytest.approx(0.320)
    assert p.cabinet.depth == pytest.approx(0.350)
    assert p.cabinet.coax_centre_z == pytest.approx(0.950)
    assert p.cabinet.coax_chamber_net_volume == pytest.approx(0.010)
    assert p.vents.slots_per_side == 6
    assert p.room.seat == pytest.approx((3.425, 2.210, 0.950))
    assert p.drivers_measured is False

def test_depth_over_350_rejected(tmp_path):
    text = Path("params/default.toml").read_text().replace("depth_mm = 350", "depth_mm = 351")
    f = tmp_path / "p.toml"; f.write_text(text)
    with pytest.raises(ValueError, match="depth"):
        load_params(f)

def test_side_edge_radius_not_implemented(tmp_path):
    text = Path("params/default.toml").read_text().replace("side_edge_radius_mm = 0", "side_edge_radius_mm = 6")
    f = tmp_path / "p.toml"; f.write_text(text)
    with pytest.raises(ValueError, match="side_edge_radius is not implemented yet"):
        load_params(f)
