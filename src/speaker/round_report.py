"""The round report: ``out_dir/report.html``, rendered from a round's stored results (spec §7).

Sections, in order: Recommendation (notes.md, key dimensions against round 0, a drawing), Curves (the
recommendation's plots), Scores, Confidence, Crossover, What I'd try next (notes.md), and an appendix of
every variant per stage. Nothing is evaluated: everything comes from the files speaker.round wrote.
"""
from __future__ import annotations

import datetime
import html
import json
import math
import re
import shutil
from pathlib import Path

from .crossover import ActiveCrossover, FixedCrossover, Peq, WayFilter
from .xo_fit import WAYS, passive_warnings, peak_boost
from .directivity_report import _CSS, COARSE, FINE, VENT_SPREAD, _row
from .evaluate import SETTINGS, grid
from .params import Params, load_params
from .ranking import DIRECTIVITY_KEYS, GATE_UNMET_TOL_DB, REFERENCE
from .shape_drawing import draw_designs
from .variants import RoundSpec

e = html.escape
PLACEHOLDER = "Draft: the recommendation is written in the next session."
NEXT_PLACEHOLDER = "Draft: what to try next is written in the next session."
COMPARE_NOTE = "(scored with earlier code; not like-for-like)"
CARDIOID_KEYS = ("cardioid_hard", "cardioid_aim", "rear_hemisphere")  # a cardioid stage's grid stops at 600 Hz
CURVES = (
    ("polar_h", "Horizontal polar map: level re on-axis at 2 m. The dashed line is the trusted 3D limit."),
    ("polar_v", "Vertical polar map (positive angles up), normalised to on-axis."),
    ("spinorama", "Spinorama at 2 m after the crossover. Listening window re-centred on 23°. Both DIs on the right axis."),
    ("window_family", "H 0-30° in 5° steps minus the 23° listener curve; ±1 dB band and 2-5 kHz shaded."),
    ("rear_rejection", "Coax section alone: front minus back, with the hard and aim lines. Shaded: the vent "
                       "resistance ×0.7 to ×1.3."),
)
_EXTRA_CSS = """
.banner { border: 2px solid #c62828; padding: 6px 10px; color: #c62828; font-weight: bold; }
table.grey { color: #888; } table.grey .pass, table.grey .fail, table.grey .warn { color: #888; }
td.out { color: #555; font-size: 0.85em; max-width: 380px; }
"""


# ---------------------------------------------------------------- reading

def _read_json(path: Path, default=None):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return default


def _notes(out_dir: Path) -> dict[str, str]:
    """The text under each `## ` heading of notes.md."""
    try:
        text = (out_dir / "notes.md").read_text()
    except OSError:
        return {}
    sections: dict[str, list[str]] = {}
    cur = None
    for line in text.splitlines():
        if line.startswith("## "):
            cur = line[3:].strip()
            sections[cur] = []
        elif cur is not None:
            sections[cur].append(line)
    return {k: "\n".join(v).strip() for k, v in sections.items()}


def _text_html(text: str) -> str:
    """Paragraphs and `- ` lists, escaped."""
    out = []
    for block in re.split(r"\n\s*\n", text.strip()):
        lines = [x for x in block.splitlines() if x.strip()]
        if not lines:
            continue
        if all(x.lstrip().startswith("- ") for x in lines):
            out.append("<ul>" + "".join(f"<li>{e(x.lstrip()[2:].strip())}</li>" for x in lines) + "</ul>")
        else:
            out.append(f"<p>{e(' '.join(x.strip() for x in lines))}</p>")
    return "".join(out)


def _section_text(notes: dict[str, str], heading: str, placeholder: str) -> str:
    return _text_html(notes[heading]) if notes.get(heading) else f"<p><i>{e(placeholder)}</i></p>"


# ---------------------------------------------------------------- scores

def _fmt(v, nd=2) -> str:
    return "n/a" if v is None else f"{v:.{nd}f}"


def _by_key(scores: list[dict], conf: str) -> dict[str, dict]:
    return {s["key"]: s for s in scores if s.get("confidence") == conf}


def _cell(s: dict | None) -> str:
    if s is None:
        return "<td>-</td>"
    return f'<td class="num {e(s["status"])}">{_fmt(s["value"])}</td>'


def _score_table(rows: list[tuple[str, list[dict] | None]], conf: str, grey: bool) -> str:
    """One row per (label, scores); one column per score key of this confidence."""
    labels: dict[str, str] = {}
    for _, sc in rows:
        for s in _by_key(sc or [], conf).values():
            labels.setdefault(s["key"], s["label"])
    first = [k for k in ("cardioid_hard", "cardioid_aim") if k in labels]
    keys = first + [k for k in labels if k not in first]
    if not keys:
        return f"<p><i>No {e(conf)}-confidence scores.</i></p>"
    head = "<tr><th>Design</th>" + "".join(f"<th title='{e(labels[k])}'>{e(k)}</th>" for k in keys) + "</tr>"
    body = []
    for label, sc in rows:
        if sc is None:
            body.append(_row([f"<td>{e(label)}</td>", f'<td colspan="{len(keys)}"><i>scores not found</i></td>']))
            continue
        by = _by_key(sc, conf)
        body.append(_row([f"<td>{e(label)}</td>", *(_cell(by.get(k)) for k in keys)]))
    return f"<table{' class="grey"' if grey else ''}>{head}{''.join(body)}</table>"


# ---------------------------------------------------------------- the data of a round

class _Round:
    def __init__(self, out_dir: Path, spec: RoundSpec, data: dict):
        self.out_dir, self.spec, self.data = out_dir, spec, data
        self.stages = data["stages"]
        self.finals = [f for f in _read_json(out_dir / "finals.json", data.get("finals", [])) or []]
        self.last = self.stages[-1] if self.stages else None

    def result(self, stage: dict, name: str) -> dict | None:
        return next((r for r in stage["results"] if r["name"] == name), None)

    def recommendation(self) -> tuple[str, str, int] | None:
        """(name, params path relative to out_dir, rank) of the first final that ran, else the stage's best."""
        for f in self.finals:
            if f.get("status") == "ok":
                return f["name"], f["params"], int(f["rank"])
        if self.last:
            best = next((n for n in self.last["ranking"]["order"] if n != REFERENCE), None)
            if best:
                return best, f"{self.last['name']}/{best}/params_tuned.toml", 1
        return None

    def final_scores(self, name: str) -> list[dict] | None:
        sc = _read_json(self.out_dir / "finals" / name / "scores.json")
        if sc is None and self.last and (r := self.result(self.last, name)):
            sc = r["scores"] or None
        return sc

    def round0_rescored(self) -> tuple[str, list[dict]] | None:
        """(label, scores) of the base design as the last stage ran it with this code, when that result is ok."""
        if not self.last:
            return None
        r = (_read_json(self.out_dir / self.last["name"] / REFERENCE / "result.json")
             or self.result(self.last, REFERENCE))
        if not isinstance(r, dict) or r.get("status") != "ok" or not r.get("scores"):
            return None
        return (f"{_cap(self.spec.base_label)} rescored (this code, {r.get('settings', self.last['settings'])} settings, "
                f"{'tuned' if r.get('tuned') else 'untuned'})", r["scores"])

    def round0_params(self) -> Path:
        if self.last:
            p = self.out_dir / self.last["name"] / REFERENCE / "params_tuned.toml"
            if p.exists():
                return p
        return self.spec.base


# ---------------------------------------------------------------- sections

def _cap(s: str) -> str:
    return s[:1].upper() + s[1:]


def _dimensions(p1: Params, p0: Params, base: str) -> str:
    def shape(p: Params) -> str:
        c = p.cabinet
        mm = lambda x: f"{x * 1000:.1f}".rstrip("0").rstrip(".")  # noqa: E731
        if c.front_edge in ("roundover", "chamfer"):
            return f"{c.front_edge}, {mm(c.edge_size)} mm"
        if c.front_edge == "facets":
            return (f"facets, {math.degrees(c.facet_angle):.0f}° back, {mm(c.facet_width)} mm wide, from "
                    f"{mm(c.facet_base_z)} mm up (ramp {mm(c.facet_ramp)} mm)")
        return "square"

    mm = lambda x: f"{x * 1000:.1f}".rstrip("0").rstrip(".") + " mm"  # noqa: E731
    rows = [
        ("Width", lambda p: mm(p.cabinet.width)),
        ("Front edge and size", shape),
        ("Flare mouth", lambda p: mm(p.cabinet.flare_mouth_d)),
        ("Flare profile", lambda p: p.cabinet.flare_profile),
        ("Outer baffle layer", lambda p: mm(p.cabinet.outer_baffle_t)),
        ("Coax offset", lambda p: mm(p.cabinet.coax_offset_x)),
        ("Top margin", lambda p: mm(p.cabinet.top_margin)),
        ("Vent position (from the front)", lambda p: mm(p.vents.centre_from_front)),
        ("Vent slots (per side)", lambda p: f"{p.vents.slots_per_side} of {p.vents.slot_w * 1000:.0f} x "
                                            f"{p.vents.slot_h * 1000:.0f} mm, pitch {p.vents.slot_pitch * 1000:.0f} mm"),
        ("Vent resistance R", lambda p: f"{p.vents.specific_resistance:.0f} rayl"),
        ("Coax chamber net volume V", lambda p: f"{p.cabinet.coax_chamber_net_volume * 1000:.1f} L"),
        ("Bass driver heights", lambda p: " and ".join(mm(z) for z in p.cabinet.bass_centres_z)),
    ]
    body = []
    for label, fn in rows:
        a, b = fn(p1), fn(p0)
        body.append(_row([f"<td>{e(label)}</td>", f"<td>{e(a)}</td>", f"<td>{e(b)}</td>"],
                         "" if a == b else "chg"))
    return (f"<table><tr><th>Dimension</th><th>Recommendation</th><th>{e(_cap(base))}</th></tr>" + "".join(body)
            + f"</table><div class='cap'>Rows that differ from {e(base)} are shaded.</div>")


def _recommendation(rd: _Round, notes: dict[str, str]) -> str:
    parts = ["<h2>Recommendation</h2>", _section_text(notes, "Recommendation", PLACEHOLDER)]
    rec = rd.recommendation()
    if rec is None:
        parts.append("<p class='fail'>No variant ran to completion, so there is no recommendation.</p>")
        return "".join(parts)
    name, rel, rank = rec
    try:
        p1, p0 = load_params(rd.out_dir / rel), load_params(rd.round0_params())
    except (OSError, ValueError, KeyError) as exc:
        parts.append(f"<p class='fail'>The parameter files could not be read: {e(str(exc))}</p>")
        return "".join(parts)
    label = rd.spec.base_label
    parts.append(f"<p>Recommended: <b>{e(name)}</b> (rank {rank}), against {e(label)}.</p>")
    parts.append(_dimensions(p1, p0, label))
    draw_designs(rd.out_dir / "shape.png", [(f"Recommendation: {name}", p1), (_cap(label), p0)])
    parts.append(f'<img src="shape.png" alt="front elevation and plan section: the recommendation and {e(label)}">'
                 '<div class="cap">Front elevation above, plan section at the coax centre height below, both '
                 'at the same scale.</div>')
    return "".join(parts)


def _curves(rd: _Round) -> str:
    parts = ["<h2>Curves</h2>"]
    rec = rd.recommendation()
    if rec is None:
        return parts[0] + "<p>No recommendation, so no curves.</p>"
    src = rd.out_dir / "finals" / rec[0] / "plots"
    plots = rd.out_dir / "plots"
    plots.mkdir(exist_ok=True)
    for stem, caption in CURVES:
        if (src / f"{stem}.png").exists():
            shutil.copyfile(src / f"{stem}.png", plots / f"{stem}.png")
            parts.append(f'<img src="plots/{stem}.png" alt="{stem}"><div class="cap">{e(caption)}</div>')
        else:
            parts.append(f"<p class='warn'>{e(stem)}.png was not found in {e(str(src))}.</p>")
    return "".join(parts)


def _scores(rd: _Round) -> str:
    parts = ["<h2>Scores</h2>"]
    gated = [s["name"] for s in rd.stages if s["ranking"].get("gate_unmet")]
    if gated:
        parts.append(f"<p class='banner'>The hard cardioid gate is unmet by every variant (stage "
                     f"{e(', '.join(gated))}). They are ranked on directivity among the variants within "
                     f"{GATE_UNMET_TOL_DB:g} dB of the best rear rejection, then by rear rejection.</p>")
    rows: list[tuple[str, list[dict] | None]] = []
    rec = rd.recommendation()
    done = [f for f in rd.finals if f.get("status") == "ok"]
    for f in done:
        tag = "recommendation" if rec and f["name"] == rec[0] else f"runner-up {f['rank']}"
        rows.append((f"{f['rank']}. {f['name']} ({tag})", rd.final_scores(f["name"])))
    if r0 := rd.round0_rescored():
        rows.append(r0)
    for c in rd.spec.compare:
        rows.append((f"{c.label} {COMPARE_NOTE}", _read_json(Path(c.scores))))
    for f in rd.finals:
        if f.get("status") != "ok":
            last = [x for x in str(f.get("reason", "")).splitlines() if x.strip()]
            parts.append(f"<p class='fail'>Final {f.get('rank')}, {e(str(f['name']))}, did not complete "
                         f"({e(str(f.get('status')))}): {e(last[-1].strip() if last else 'no reason recorded')}</p>")
    parts.append(_score_table(rows, "trusted", False))
    parts.append("<p class='cap'>Trusted scores only (below the trusted limit); they decide the ranking. "
                 "Green passes, orange warns, red fails.</p>")
    parts.append("<h3>Low confidence (above the trusted limit; they do not rank designs)</h3>")
    parts.append(_score_table(rows, "low", True))
    return "".join(parts)


def _confidence(rd: _Round) -> str:
    spec = rd.spec
    scr = SETTINGS["screen"]
    top = float(grid(scr.f_lo, spec.trusted_limit, scr.ppo)[-1])
    tune_spread = ""
    rec = rd.recommendation()
    if rec and rd.last and (r := rd.result(rd.last, rec[0])) and r.get("tune"):
        sp = r["tune"]["spread"]
        tune_spread = (" For the recommendation, the cardioid hard / aim rejection is "
                       + "; ".join(f"{sp[k][0]:.1f} / {sp[k][1]:.1f} dB at R × {k}" for k in ("0.7", "1.3"))
                       + f", and {r['tune']['hard']:.1f} / {r['tune']['aim']:.1f} dB at R × 1.")
    items = [
        f"<b>Trusted limit</b>: {spec.trusted_limit:.0f} Hz, {e(spec.trusted_note)}. Scores are split at the "
        "limit; the ones above it are low confidence and do not rank designs. The variants are not mesh-checked: "
        f"the mesh check was run on {e(spec.base_label)}, and the limit is carried over to them.",
        f"<b>Screen against finals</b>: the stages that rank by directivity ran with the screen settings, "
        f"{scr.per_wavelength:g} elements per wavelength up to the trusted limit. The screen grid stops one "
        f"1/{scr.ppo}-octave step below the trusted limit, at {top:.0f} Hz against {spec.trusted_limit:.0f} Hz. "
        f"The finals were rerun at {COARSE:g} elements per wavelength, with the coax puck. The mesh check that "
        f"set the trusted limit compared {COARSE:g} against {FINE:g} elements per wavelength, so the screen's "
        f"{scr.per_wavelength:g} per wavelength is not covered by it. Screen rankings in the "
        f"1–{spec.trusted_limit / 1000:.2f} kHz bands, presence lift especially, are therefore less certain than "
        f"the finals' scores. The recommendation should be argued from the finals' {COARSE:g}-per-wavelength "
        "scores.",
        "<b>Flare mesh</b>: the mesh understates the ellipse-against-cone flare difference by about 3-7 % at "
        "h_coax 10 mm.",
        "<b>Coax LF de-embedding</b>: above the LF join the coax's low-frequency volume velocity is the datasheet "
        "on-axis response divided by that of an infinite-baffle piston of the diaphragm diameter, times the recess "
        "ratio R (the BEM on-axis pressure of the cone in its own frame recess over that of a flat piston, on a "
        "1 m disc baffle). It assumes the datasheet was measured front-mounted on a large baffle.",
        "<b>Felt spread</b>: the vent resistance as built is unknown, so it is taken as ±30 % "
        f"(×{VENT_SPREAD[0]:g} to ×{VENT_SPREAD[1]:g}). Each variant is tuned for its worst case over that spread "
        "(the cardioid margin)." + tune_spread,
        "<b>Estimated</b>: the coax displacement and filler volumes (a basket cylinder until the driver is "
        "measured). Above the trusted limit the coax's horn behaviour rests on an estimated cone profile. "
        "A geometry decision should not rest on low-confidence detail alone.",
    ]
    return "<h2>Confidence</h2><ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>"


def _gains(final_dir: Path) -> dict[str, float] | None:
    """The way gains (dB) run_directivity wrote to gains.json, or None when missing or unreadable."""
    g = _read_json(final_dir / "gains.json")
    return g if isinstance(g, dict) else None


def _fitted(final_dir: Path) -> tuple[ActiveCrossover, float] | None:
    """The ActiveCrossover and the trusted limit (Hz) in crossover.json, or None when missing or unreadable."""
    d = _read_json(final_dir / "crossover.json")
    try:
        ways = {k: WayFilter(float(w["gain_db"]), float(w["delay_s"]), int(w["polarity"]),
                             tuple(Peq(float(p["f0"]), float(p["gain_db"]), float(p["q"])) for p in w["peq"]))
                for k, w in d["ways"].items()}
        xo = ActiveCrossover(float(d["bass_coax"]), float(d["lf_hf"]), ways, float(d["fit_top"]))
        return xo, float(d["trusted_limit"])
    except (KeyError, TypeError, ValueError, AttributeError):
        return None


def _crossover(rd: _Round) -> str:
    fixed = FixedCrossover()
    fitted, plain, missing = [], [], []
    for f in rd.finals:
        if f.get("status") != "ok":
            continue
        fdir = rd.out_dir / f["dir"]
        got = _fitted(fdir)
        if got is not None:
            xo, tl = got
            delays = ", ".join(f"{k} {xo.ways[k].delay_s * 1e3:.2f} ms" for k in WAYS)
            npeq = ", ".join(f"{k} {len(xo.ways[k].peq)}" for k in WAYS)
            boost = max(peak_boost(xo, k)[0] for k in WAYS)
            fitted.append(_row([f"<td>{e(f['name'])}</td>", f'<td class="num">{xo.bass_coax:g} Hz</td>',
                                f'<td class="num">{xo.lf_hf:g} Hz</td>', f"<td>{delays}</td>", f"<td>{npeq}</td>",
                                f'<td class="num">{boost:.1f} dB</td>',
                                f'<td class="num">{len(passive_warnings(xo, tl))}</td>']))
            continue
        g = _gains(fdir)
        if g is None:
            missing.append(f["name"])
        cell = lambda k: "n/a" if g is None or g.get(k) is None else f"{g[k]:+.2f} dB"  # noqa: E731
        plain.append(_row([f"<td>{e(f['name'])}</td>"]
                          + [f'<td class="num">{cell(k)}</td>' for k in ("bass", "coax_lf", "coax_hf")]))
    parts = ["<h2>Crossover</h2>"]
    if fitted:
        parts.append("<p>Active crossover fitted per final: acoustic LR4 at the two points, with per-way parametric "
                     "EQ towards a flat 23° response, delays (relative) and polarity. The warnings count where a "
                     "passive network could not follow (more than 6 dB of summed EQ boost, a coax LF to HF point "
                     "within an octave of the 750 Hz HF resonance, a coax HF-LF delay over 0.1 ms, or EQ above the "
                     "trusted limit); each final's report lists them.</p>"
                     "<table><tr><th>Final</th><th>Bass to coax</th><th>Coax LF to HF</th><th>Delays</th>"
                     "<th>PEQ filters</th><th>Largest boost</th><th>Warnings</th></tr>" + "".join(fitted) + "</table>")
    if plain:
        parts.append(f"<p>Fixed and provisional: acoustic LR4 at {fixed.bass_coax:g} Hz (bass to coax) and "
                     f"{fixed.lf_hf:g} Hz (coax LF to HF). No EQ and no delays. Each way's gain matches its level "
                     "at 23° to the coax LF way.</p>")
    if missing:
        parts.append(f"<p class='warn'>Way gains not found (no gains.json) for: {e(', '.join(missing))}.</p>")
    if plain:
        parts.append("<table><tr><th>Final</th><th>bass</th><th>coax_lf</th><th>coax_hf</th></tr>"
                     + "".join(plain) + "</table>")
    return "".join(parts)


def _tuning_bought(r: dict) -> str:
    t = r.get("tune")
    if not r.get("tuned") or not t or t.get("margin") is None or t["base"].get("margin") is None:
        return "-"
    return f"{t['margin'] - t['base']['margin']:+.1f} dB"


def _appendix(rd: _Round) -> str:
    parts = ["<h2>All variants</h2>"]
    for st in rd.stages:
        keys = CARDIOID_KEYS if st["rank_by"] == "cardioid" else DIRECTIVITY_KEYS
        ranking = st["ranking"]
        by_name = {r["name"]: r for r in st["results"]}
        order = [n for n in ranking["order"] if n in by_name]
        rest = [r["name"] for r in st["results"] if r["name"] not in order]
        why = ("ranked by the tuned cardioid margin; the grid stops at 600 Hz, so the directivity scores are "
               "partial and are not shown" if st["rank_by"] == "cardioid"
               else "ranked by directivity behind the hard cardioid gate")
        parts.append(f"<h3>{e(st['name'])} ({e(st['settings'])}, {why})</h3>")
        head = ("<tr><th>Rank</th><th>Variant</th><th>Status</th><th>Why out</th><th>Changes</th>"
                + "".join(f"<th>{e(k)}</th>" for k in keys)
                + "<th>Tuned R (rayl)</th><th>Tuned V (L)</th><th>Tuning bought</th></tr>")
        body = []
        for rank, name in [(str(i), n) for i, n in enumerate(order, 1)] + [("-", n) for n in rest]:
            r = by_name[name]
            sc = _by_key(r["scores"], "trusted")
            tuned = r.get("tuned") and r.get("tune")
            changes = "; ".join(f"{k} = {v}" for k, v in r["overrides"].items()) or "-"
            body.append(_row([
                f'<td class="num">{rank}</td>', f"<td>{e(name)}</td>", f"<td>{e(r['status'])}</td>",
                f"<td class='out'>{e(ranking['out'].get(name, ''))}</td>", f"<td class='out'>{e(changes)}</td>",
                *(_cell(sc.get(k)) for k in keys),
                f'<td class="num">{r["tune"]["resistance"]:.0f}</td>' if tuned else "<td>-</td>",
                f'<td class="num">{r["tune"]["net_volume"] * 1000:.2f}</td>' if tuned else "<td>-</td>",
                f'<td class="num">{_tuning_bought(r)}</td>']))
        parts.append(f"<table>{head}{''.join(body)}</table>")
    return "".join(parts)


# ---------------------------------------------------------------- the report

def write_report(out_dir: Path, spec: RoundSpec) -> Path:
    """Write out_dir/report.html from the results in out_dir; returns its path. Raises ValueError without results."""
    out_dir = Path(out_dir)
    data = _read_json(out_dir / "results.json")
    if not isinstance(data, dict) or "stages" not in data:
        raise ValueError(f"no results to render: {out_dir / 'results.json'} is missing or unreadable "
                         "(a round that stopped at a stage writes none; see that stage's ranking.json)")
    rd = _Round(out_dir, spec, data)
    notes = _notes(out_dir)
    head = [("Round file", e(str(data.get("spec", "")))), ("Round file sha256", f"<code>{e(str(data.get('spec_sha256', '')))}</code>"),
            ("Code commit", e(str(data.get("commit", "")))), ("Code version", e(str(data.get("code_version", "")))),
            ("Date", e(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))]
    page = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        f"<title>Round {spec.number}</title><style>{_CSS}{_EXTRA_CSS} tr.chg td {{ background: #fff6dd; }}</style></head><body>",
        f"<h1>Round {spec.number}: {e(spec.title)}</h1>",
        "<table class='meta'>" + "".join(_row([f"<td>{k}</td>", f"<td>{v}</td>"]) for k, v in head) + "</table>",
        _recommendation(rd, notes),
        _curves(rd),
        _scores(rd),
        _confidence(rd),
        _crossover(rd),
        "<h2>What I'd try next</h2>" + _section_text(notes, "What I'd try next", NEXT_PLACEHOLDER),
        _appendix(rd),
        "</body></html>",
    ]
    path = out_dir / "report.html"
    path.write_text("\n".join(page))
    return path
