import dataclasses, pytest
from speaker.params import load_params
from speaker.drivers import load_driver
from speaker.layout import compute_layout

P = load_params(); CX = load_driver("bms_10c262"); B = load_driver("dayton_rss265ho_44")

def test_default_layout():
    l = compute_layout(P, CX, B)
    assert l.height == pytest.approx(0.950 + 0.130 + 0.020 + 0.018, abs=1e-6)  # 1.118 m
    assert l.interior_w == pytest.approx(0.284)
    assert l.interior_d == pytest.approx(0.296)
    assert l.bass_net >= 0.050                 # >=25 L per driver
    assert l.coax_chamber_net == pytest.approx(0.010)
    assert len(l.vent_slots) == 6
    assert l.warnings == ()

@pytest.mark.parametrize("ear", [0.85, 1.05])
def test_ear_height_moves_everything_consistently(ear):
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, coax_centre_z=ear))
    l = compute_layout(p, CX, B)
    assert l.coax_centre[1] == pytest.approx(ear)
    assert l.divider_z < ear - CX.outer_d / 2
    assert l.height == pytest.approx(ear + 0.168, abs=1e-6)

def test_too_low_ear_height_rejected():     # divider would hit the upper bass driver
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, coax_centre_z=0.70))
    with pytest.raises(ValueError, match="divider|bass"):
        compute_layout(p, CX, B)

def test_coax_wider_than_interior_rejected():
    big = dataclasses.replace(CX, outer_d=0.300)
    with pytest.raises(ValueError, match="coax"):
        compute_layout(P, big, B)

def test_overlapping_bass_rejected():
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, bass_centres_z=(0.20, 0.30)))
    with pytest.raises(ValueError, match="overlap"):
        compute_layout(p, CX, B)

def test_low_ear_height_warns_on_bass_volume():
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, coax_centre_z=0.85))
    l = compute_layout(p, CX, B)
    assert l.bass_net < 0.050
    assert len(l.warnings) == 1 and "bass" in l.warnings[0]

def test_bass_net_under_30l_rejected():
    fat = dataclasses.replace(B, basket_d=0.284, depth=0.300, flange_t=0.010)
    with pytest.raises(ValueError, match="bass"):
        compute_layout(P, CX, fat)
