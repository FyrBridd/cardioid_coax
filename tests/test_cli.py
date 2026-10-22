from speaker.cli import main


def test_build_end_to_end(tmp_path):
    assert main(["build", "--out", str(tmp_path)]) == 0
    s = (tmp_path / "summary.md").read_text()
    for h in ("## Bass", "## Cardioid", "## Room", "## Cabinet", "## Confidence"):
        assert h in s
    for png in ("bass_spl", "bass_in_room", "cardioid_polar", "front_back"):
        assert (tmp_path / "plots" / f"{png}.png").stat().st_size > 0
    assert (tmp_path / "step/cabinet_assembly.step").exists()
    # Cardioid: target-volume and as-built (no filler) results, with the <10 dB crossing.
    assert "target net volume 10.00 L" in s and "as built" in s
    assert s.count("Minimum front/back over 150-300 Hz") == 2
    assert s.count("first drops below 10 dB") == 2
    # Cabinet: filler is not in the CAD.
    assert "filler" in s and "NOT in the CAD" in s
    assert "geometry high; net volumes medium (filler not modelled; displacement estimated)" in s
    # Bass: shared inductor DCR doubled per driver.
    assert "0.4 ohm shared by 2 drivers in parallel -> 0.8 ohm effective per driver" in s
