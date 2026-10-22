import dataclasses

import cadquery as cq
import ezdxf
import ezdxf.bbox
import pytest

from speaker.cad.export import export_all
from speaker.drivers import load_driver
from speaker.layout import compute_layout
from speaker.params import load_params

P = load_params()
CX = load_driver("bms_10c262")
B = load_driver("dayton_rss265ho_44")
L = compute_layout(P, CX, B)


def test_exports(tmp_path):
    files = export_all(L, P, CX, B, tmp_path)
    assert len(files["dxf"]) == 8 and all(f.name.endswith("_PROVISIONAL.dxf") for f in files["dxf"])
    assert all(f.exists() for fs in files.values() for f in fs)
    wp = cq.importers.importStep(str(tmp_path / "step/cabinet_assembly.step"))
    bb = cq.Compound.makeCompound(wp.vals()).BoundingBox()
    assert (bb.xlen, bb.ylen) == pytest.approx((320, 350), abs=0.5)
    assert (tmp_path / "stl/flare_insert.stl").stat().st_size > 0
    assert (tmp_path / "step/flare_insert.step").exists()
    assert (tmp_path / "step/baffle_inner.step").exists()


def test_dxf_layers(tmp_path):
    export_all(L, P, CX, B, tmp_path)
    doc = ezdxf.readfile(tmp_path / "dxf/baffle_inner_PROVISIONAL.dxf")
    assert doc.units == ezdxf.units.MM
    layers = {e.dxf.layer for e in doc.modelspace()}
    assert {"THROUGH", "POCKET_10.0MM", "NOTES"} <= layers
    ext = ezdxf.bbox.extents(e for e in doc.modelspace() if e.dxf.layer == "THROUGH")
    assert ext.size.x == pytest.approx(320 - 2 * 18, abs=0.5)  # inner baffle fits between sides
    notes = " ".join(e.text for e in doc.modelspace() if e.dxf.layer == "NOTES")
    assert "PROVISIONAL - drivers not yet measured" in notes
    assert "18.0 mm stock" in notes and "rear" in notes


def test_rounded_rect_has_arcs(tmp_path):
    export_all(L, P, CX, B, tmp_path)
    doc = ezdxf.readfile(tmp_path / "dxf/back_PROVISIONAL.dxf")
    polys = [e for e in doc.modelspace() if e.dxf.layer == "THROUGH" and e.dxftype() == "LWPOLYLINE"]
    assert len(polys) == 2  # outline + terminal cutout
    assert any(any(pt[4] != 0 for pt in e.get_points("xyseb")) for e in polys)


def test_measured_drops_provisional(tmp_path):
    files = export_all(L, dataclasses.replace(P, drivers_measured=True), CX, B, tmp_path)
    assert not any("PROVISIONAL" in f.name for f in files["dxf"])
    doc = ezdxf.readfile(files["dxf"][0])
    assert not any("PROVISIONAL" in e.text for e in doc.modelspace() if e.dxftype() == "MTEXT")
