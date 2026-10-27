import math
from pathlib import Path

import numpy as np
import pytest

from speaker.bass import sealed_spl
from speaker.drivers import load_driver
from speaker.layout import compute_layout
from speaker.params import load_params
from speaker.report import build, cardioid_warnings, fb_crossing, in_room_bass

P = load_params(); CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")
L = compute_layout(P, CX, B)


def test_fb_crossing_interpolates_first_drop():
    f = np.array([100.0, 150.0, 200.0, 300.0, 400.0, 500.0])
    fb = np.array([5.0, 20.0, 15.0, 5.0, 12.0, 3.0])
    assert fb_crossing(f, fb) == pytest.approx(250.0)


def test_fb_crossing_already_below_at_band_start():
    f = np.array([100.0, 200.0, 300.0])
    fb = np.array([12.0, 8.0, 6.0])  # 10 dB at 150 Hz is the boundary; below just after
    assert fb_crossing(f, fb) == pytest.approx(150.0, abs=1.0)
    assert fb_crossing(f, np.array([9.0, 8.0, 6.0])) == pytest.approx(150.0)


def test_fb_crossing_none_when_never_below():
    f = np.array([100.0, 200.0, 300.0])
    assert fb_crossing(f, np.array([12.0, 15.0, 11.0])) is None


def _half_space_coherent(f0):
    """Both speakers' bass drivers, half-space source level, coherent at low f."""
    room = P.room
    x = room.back_gap + P.cabinet.depth
    lvl = 10 ** (sealed_spl(B, L.bass_net / 2, np.array([f0]), 2.83, r=1.0)[0] / 20)
    tot = sum(lvl / math.dist((x, y, z), room.seat) for y in room.speaker_y for _, z in L.bass_centres)
    return 20 * math.log10(tot)


def test_in_room_floor_image_gives_half_space_level():
    # At 20 Hz the floor image is coherent with the source, so free-field source
    # level + floor image must reproduce the half-space level (not +6 dB on top).
    f = np.array([20.0])
    floor_only = in_room_bass(P, L, B, L.bass_net / 2, f, surfaces=["floor"])[0]
    assert floor_only == pytest.approx(_half_space_coherent(20.0), abs=0.3)


def test_in_room_free_field_is_6db_below_half_space():
    f = np.array([20.0])
    free = in_room_bass(P, L, B, L.bass_net / 2, f, surfaces=[])[0]
    assert free == pytest.approx(_half_space_coherent(20.0) - 6.02, abs=0.1)


def test_cardioid_warnings():
    assert cardioid_warnings(12.0, 7.0) == []
    w = cardioid_warnings(9.4, 5.2)
    assert len(w) == 2 and "9.4 dB" in w[0] and "10 dB hard requirement" in w[0]
    assert "5.2 dB" in w[1] and "6 dB aim" in w[1]


@pytest.mark.parametrize("old, new, name", [
    ('front_edge = "square"', 'front_edge = "roundover"', "front_edge 'roundover'"),
    ("coax_offset_x_mm = 0", "coax_offset_x_mm = 3", "coax_offset_x 3.0 mm"),
    ("outer_baffle_t_mm = 18", "outer_baffle_t_mm = 24", "outer_baffle_t 24.0 mm"),
    ('flare_profile = "ellipse"', 'flare_profile = "cone"', "flare_profile 'cone'"),
])
def test_build_refuses_bem_only_geometry(tmp_path, old, new, name):
    text = Path("params/default.toml").read_text()
    assert old in text
    text = text.replace(old, new)
    if "roundover" in new:
        text = text.replace("edge_size_mm = 0", "edge_size_mm = 5")
    pf = tmp_path / "p.toml"; pf.write_text(text)
    out = tmp_path / "out"
    with pytest.raises(ValueError, match="the CAD does not model these yet") as e:
        build(pf, out)
    assert name in str(e.value) and "`speaker directivity` and `speaker round` accept them" in str(e.value)
    assert not out.exists()
