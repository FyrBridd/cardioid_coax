from pathlib import Path

from speaker.cli import main


def test_build_end_to_end(tmp_path):
    assert main(["build", "--out", str(tmp_path)]) == 0
    s = (tmp_path / "summary.md").read_text()
    for h in ("## Bass", "## Cardioid", "## Room", "## Cabinet", "## Confidence"):
        assert h in s
    for png in ("bass_spl", "bass_in_room", "cardioid_polar", "front_back"):
        assert (tmp_path / "plots" / f"{png}.png").stat().st_size > 0
    assert (tmp_path / "step/cabinet_assembly.step").exists()
    # Cardioid: as built with fillers, and without fillers for comparison.
    assert "As built (with printed fillers: net 10.0" in s and "Without fillers" in s
    assert s.count("Minimum front/back over 150-300 Hz") == 2
    assert s.count("Minimum front/back over 300-500 Hz") == 2
    assert s.count("first drops below 10 dB") == 2
    assert "- Cardioid warnings: none" in s
    # Cabinet: the fillers are in the CAD.
    assert "NOT in the CAD" not in s and "rear block" in s
    assert "geometry high; net volumes medium (displacement estimated)" in s
    assert (tmp_path / "stl/filler_rear_block_PROVISIONAL.stl").exists()
    # Bass: shared inductor DCR doubled per driver.
    assert "0.4 ohm shared by 2 drivers in parallel -> 0.8 ohm effective per driver" in s


def test_build_rejects_bad_params_with_exit_code_2(tmp_path, capsys):
    text = Path("params/default.toml").read_text().replace("depth_mm = 350", "depth_mm = 351")
    pf = tmp_path / "p.toml"
    pf.write_text(text)
    out = tmp_path / "out"
    assert main(["build", "--params", str(pf), "--out", str(out)]) == 2
    err = capsys.readouterr().err
    assert err.startswith("error: ") and "depth 351 mm" in err
    assert not out.exists()  # nothing written before validation


def test_cli_passive_subcommand_parses(monkeypatch, tmp_path):
    import speaker.passive_report as pr
    calls = []

    def fake(params_path, out_dir, *, trusted_limit):
        calls.append((params_path, out_dir, trusted_limit))
        return tmp_path / "report.html"
    monkeypatch.setattr(pr, "run_passive", fake)
    assert main(["passive", "--trusted-limit", "2851"]) == 0
    assert calls == [(None, Path("out/passive"), (2851.0, "given on the command line"))]
    assert main(["passive", "--params", "p.toml", "--out", str(tmp_path)]) == 0
    assert calls[-1] == (Path("p.toml"), tmp_path, None)
