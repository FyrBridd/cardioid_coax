# Cardioid coax floorstander

Parametric design pipeline: acoustics model, CadQuery cabinet and CNC export.

## Usage

```
uv sync
uv run speaker build
```

`uv run speaker build` reads `params/default.toml` (override with `--params PATH`)
and writes to `out/` (override with `--out DIR`):

- `summary.md`: bass, cardioid, room and cabinet numbers, plus a confidence table
- `plots/`: `bass_spl.png`, `bass_in_room.png`, `cardioid_polar.png`, `front_back.png`
- `step/`, `stl/`, `dxf/`: cabinet assembly, panels, flare insert and CNC drawings

`out/` is git-ignored. Run the tests with `uv run pytest`.
