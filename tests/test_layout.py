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
    p = dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, coax_centre_z=0.85),
                            room=dataclasses.replace(P.room, seat=(3.425, 2.210, 0.85)))
    l = compute_layout(p, CX, B)
    assert l.bass_net < 0.050
    assert len(l.warnings) == 1 and "bass" in l.warnings[0]

def test_bass_net_under_30l_rejected():
    fat = dataclasses.replace(B, basket_d=0.284, depth=0.300, flange_t=0.010)
    with pytest.raises(ValueError, match="bass"):
        compute_layout(P, CX, fat)

def _cab(**kw):
    return dataclasses.replace(P, cabinet=dataclasses.replace(P.cabinet, **kw))

def _vents(**kw):
    return dataclasses.replace(P, vents=dataclasses.replace(P.vents, **kw))

def test_vent_stack_taller_than_coax_chamber_rejected():
    # 20 slots at 22 mm pitch: stack runs below the divider top and above the top panel.
    with pytest.raises(ValueError, match=r"vent slot.*z=.*mm.*divider top face at 810\.0 mm"):
        compute_layout(_vents(slots_per_side=20), CX, B)

def test_vent_slot_crossing_baffle_rejected():
    # centre 60 mm from front, 100 mm wide: front edge at y=10 mm, inside the 36 mm baffle.
    with pytest.raises(ValueError, match=r"vent slot.*y=10\.0 mm.*baffle rear face at 36\.0 mm"):
        compute_layout(_vents(centre_from_front=0.060), CX, B)

def test_vent_slot_into_back_panel_rejected():
    with pytest.raises(ValueError, match=r"vent slot.*y=340\.0 mm.*back panel front face at 332\.0 mm"):
        compute_layout(_vents(centre_from_front=0.290), CX, B)

def test_zero_slots_rejected():
    with pytest.raises(ValueError, match="slots_per_side"):
        compute_layout(_vents(slots_per_side=0), CX, B)

def test_negative_top_margin_rejected():
    with pytest.raises(ValueError, match=r"top_margin -50\.0 mm.*coax frame top at 1080\.0 mm"):
        compute_layout(_cab(top_margin=-0.050), CX, B)

def test_baffle_must_be_two_panel_layers_in_layout():
    with pytest.raises(ValueError, match=r"baffle_t 40\.0 mm.*panel_t 18\.0 mm"):
        compute_layout(_cab(baffle_t=0.040), CX, B)

def test_as_built_coax_chamber_excludes_filler():
    l = compute_layout(P, CX, B)
    assert l.coax_chamber_as_built == pytest.approx(l.coax_chamber_net + l.filler_volume)
    assert l.coax_chamber_as_built == pytest.approx(0.01883, abs=2e-4)

def test_ear_height_mismatch_with_seat_warns():
    p = dataclasses.replace(P, room=dataclasses.replace(P.room, seat=(3.425, 2.210, 0.980)))
    l = compute_layout(p, CX, B)
    assert any("ear height" in w and "950.0 mm" in w and "980.0 mm" in w for w in l.warnings)
    p = dataclasses.replace(P, room=dataclasses.replace(P.room, seat=(3.425, 2.210, 0.965)))
    assert compute_layout(p, CX, B).warnings == ()
