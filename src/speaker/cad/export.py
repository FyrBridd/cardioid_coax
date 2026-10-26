"""STEP / STL / DXF export. All CAD units are millimetres."""
from __future__ import annotations

import math
from pathlib import Path

import cadquery as cq
import ezdxf
import ezdxf.units

from ..drivers import Driver
from ..layout import Layout
from ..params import Params
from .fillers import filler_solids
from .panels import Feature, Panel, assembly, build_panels, flare_insert, panel_solid

PROVISIONAL_TEXT = "PROVISIONAL - drivers not yet measured"
_ARC90 = math.tan(math.pi / 8)  # bulge of a 90 degree arc


def _rounded_rect_points(f: Feature) -> list[tuple[float, float, float]]:
    """CCW closed polyline points (x, y, bulge) of a rounded rectangle."""
    x0, x1 = f.cx_mm - f.w_mm / 2, f.cx_mm + f.w_mm / 2
    y0, y1 = f.cy_mm - f.h_mm / 2, f.cy_mm + f.h_mm / 2
    r = min(f.corner_r_mm, f.w_mm / 2, f.h_mm / 2)
    if r <= 0:
        return [(x0, y0, 0), (x1, y0, 0), (x1, y1, 0), (x0, y1, 0)]
    b = _ARC90
    return [
        (x0 + r, y0, 0), (x1 - r, y0, b), (x1, y0 + r, 0), (x1, y1 - r, b),
        (x1 - r, y1, 0), (x0 + r, y1, b), (x0, y1 - r, 0), (x0, y0 + r, b),
    ]


def _write_dxf(panel: Panel, provisional: bool, path: Path) -> None:
    doc = ezdxf.new("R2010", setup=True)
    doc.units = ezdxf.units.MM
    doc.header["$INSUNITS"] = ezdxf.units.MM
    msp = doc.modelspace()

    def layer(name: str) -> str:
        if name not in doc.layers:
            doc.layers.add(name)
        return name

    msp.add_lwpolyline(
        [(0, 0), (panel.outline_w_mm, 0), (panel.outline_w_mm, panel.outline_h_mm),
         (0, panel.outline_h_mm)],
        close=True, dxfattribs={"layer": layer("THROUGH")},
    )
    for f in panel.features:
        name = "THROUGH" if f.kind == "through" else f"POCKET_{f.depth_mm:.1f}MM"
        attrs = {"layer": layer(name)}
        if f.shape == "circle":
            msp.add_circle((f.cx_mm, f.cy_mm), f.w_mm / 2, dxfattribs=attrs)
        else:
            msp.add_lwpolyline(_rounded_rect_points(f), format="xyb", close=True, dxfattribs=attrs)

    face = "rear" if panel.name == "baffle_inner" else "drawing side (as viewed)"
    note_layer = layer("NOTES")
    lines = [f"{panel.name} - {panel.thickness_mm:.1f} mm stock - machining face: {face}"]
    if provisional:
        lines.append(PROVISIONAL_TEXT)
    for i, text in enumerate(lines):
        msp.add_mtext(
            text, dxfattribs={"layer": note_layer, "char_height": 5.0,
                              "insert": (0, -10.0 - 10.0 * i)},
        )
    doc.saveas(path)


def export_all(
    l: Layout, p: Params, coax: Driver, bass: Driver, out_dir: Path
) -> dict[str, list[Path]]:
    out_dir = Path(out_dir)
    for sub in ("step", "stl", "dxf"):
        d = out_dir / sub
        d.mkdir(parents=True, exist_ok=True)
        # These dirs are ours: drop stale files (e.g. _PROVISIONAL variants
        # left over after drivers_measured flips). Only files, no recursion.
        for old in d.iterdir():
            if old.is_file():
                old.unlink()
    step: list[Path] = []
    stl: list[Path] = []
    dxf: list[Path] = []

    asm_path = out_dir / "step/cabinet_assembly.step"
    assembly(l, p, coax, bass).export(str(asm_path))
    step.append(asm_path)

    provisional = not p.drivers_measured
    suffix = "_PROVISIONAL" if provisional else ""
    for panel in build_panels(l, p, coax, bass):
        sp = out_dir / "step" / f"{panel.name}.step"
        cq.exporters.export(panel_solid(panel), str(sp))
        step.append(sp)
        dp = out_dir / "dxf" / f"{panel.name}{suffix}.dxf"
        _write_dxf(panel, provisional, dp)
        dxf.append(dp)

    # The flare throat is the coax cutout_d, which is provisional until measured.
    flare = flare_insert(l, p, coax)
    fs = out_dir / "step" / f"flare_insert{suffix}.step"
    cq.exporters.export(flare, str(fs))
    step.append(fs)
    fl = out_dir / "stl" / f"flare_insert{suffix}.stl"
    cq.exporters.export(flare, str(fl))
    stl.append(fl)

    # Filler profiles follow the coax frame and magnet sizes: provisional too.
    for name, solid in filler_solids(l).items():
        fs = out_dir / "step" / f"{name}{suffix}.step"
        cq.exporters.export(solid, str(fs))
        step.append(fs)
        fl = out_dir / "stl" / f"{name}{suffix}.stl"
        cq.exporters.export(solid, str(fl))
        stl.append(fl)
    return {"step": step, "stl": stl, "dxf": dxf}
