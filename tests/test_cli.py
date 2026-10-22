from speaker.cli import main


def test_build_end_to_end(tmp_path):
    assert main(["build", "--out", str(tmp_path)]) == 0
    s = (tmp_path / "summary.md").read_text()
    for h in ("## Bass", "## Cardioid", "## Room", "## Cabinet", "## Confidence"):
        assert h in s
    for png in ("bass_spl", "bass_in_room", "cardioid_polar", "front_back"):
        assert (tmp_path / "plots" / f"{png}.png").stat().st_size > 0
    assert (tmp_path / "step/cabinet_assembly.step").exists()
