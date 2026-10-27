"""Trace one dark curve off a raster datasheet plot (log-f x axis, linear y axis).

Calibrate from two x gridlines (pixel, Hz) and two y gridlines (pixel, value). From
a seed pixel on the wanted curve the tracer walks column by column in both
directions, picking the dark run nearest the y predicted from the recent slope, so
two curves that cross stay apart. Where two near-tangent curves merge into one
stroke and part again, the slope cannot tell them apart: a --via pixel, in the first
column where they part (walking away from the seed), picks the branch. It writes a
CSV at 1/24-octave points and an overlay PNG for checking the trace by eye.

    uv run python tools/trace_curve.py IMAGE --x PX HZ PX HZ --y PX VAL PX VAL \
        --seed PX PY [--via PX PY ...] --range FMIN FMAX --unit "dB SPL" \
        --pdf PDF --url URL --out data/curves/NAME.csv --overlay OVERLAY.png
"""
from __future__ import annotations

import argparse
import hashlib
import shlex
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

DARK = 0.3  # grey level below which a pixel is curve ink; gridlines are >= 0.36
FRAME = 0.8  # fraction of a row/column that is ink for it to count as the plot frame
THICK = 4.0  # px, the longest run a near-flat stroke makes in one column
SLOPE_N = 4  # columns averaged for the slope prediction


def _grey(path: Path) -> np.ndarray:
    a = mpimg.imread(path)
    if a.ndim == 3:
        a = a[..., :3].mean(axis=2)
    return a.astype(float) / (255.0 if a.dtype == np.uint8 else 1.0)


def _interior(ink: np.ndarray) -> tuple[int, int, int, int]:
    """Rows and columns strictly inside the plot frame."""
    rows = np.where(ink.mean(axis=1) > FRAME)[0]
    cols = np.where(ink.mean(axis=0) > FRAME)[0]
    h, w = ink.shape
    r0 = rows[rows < h // 2].max() + 2 if np.any(rows < h // 2) else 0
    r1 = rows[rows > h // 2].min() - 2 if np.any(rows > h // 2) else h - 1
    c0 = cols[cols < w // 2].max() + 2 if np.any(cols < w // 2) else 0
    c1 = cols[cols > w // 2].min() - 2 if np.any(cols > w // 2) else w - 1
    return r0, r1, c0, c1


def _runs(grey_col: np.ndarray, r0: int, r1: int) -> list[tuple[int, int, float]]:
    """Dark runs (top, bottom, darkness-weighted centroid) in rows r0..r1."""
    ink = grey_col[r0 : r1 + 1] < DARK
    out = []
    i = 0
    while i < ink.size:
        if not ink[i]:
            i += 1
            continue
        j = i
        while j + 1 < ink.size and ink[j + 1]:
            j += 1
        lo, hi = max(i - 1, 0), min(j + 1, ink.size - 1)  # include anti-aliased edges
        w = 1.0 - grey_col[r0 + lo : r0 + hi + 1]
        y = float(np.sum(w * np.arange(r0 + lo, r0 + hi + 1)) / np.sum(w))
        out.append((r0 + i, r0 + j, y))
        i = j + 1
    return out


def _overlaps(run, runs) -> int:
    """How many of `runs` (two or more pixels long) touch `run` row-wise."""
    top, bot = run[0], run[1]
    return sum(1 for t, b, _ in runs if b > t and t <= bot + 1 and top <= b + 1)


def _walk(grey, cols, y_seed, r0, r1, via: dict[int, float]) -> dict[int, float]:
    """Follow the curve across `cols` (ordered away from the seed)."""
    ys: list[float] = [y_seed]
    found: dict[int, float] = {}
    lost = merged = False
    prev_runs: list[tuple[int, int, float]] = []
    prev2_runs: list[tuple[int, int, float]] = []
    prev_best = None
    for c in cols:
        recent = ys[-(SLOPE_N + 1) :]
        slope = (recent[-1] - recent[0]) / (len(recent) - 1) if len(recent) > 1 else 0.0
        pred = via.get(c, ys[-1] + slope)
        runs = _runs(grey[:, c], r0, r1)
        best, best_d = None, np.inf
        for top, bot, cy in runs:
            d = 0.0 if top <= pred <= bot else min(abs(pred - top), abs(pred - bot))
            if d < best_d:
                best, best_d = (top, bot, cy), d
        if best is None or best_d > 6.0 + 2.0 * abs(slope):
            # Lost (e.g. the curve runs off the plot floor): coast on the prediction.
            ys.append(float(np.clip(pred, r0, r1)))
            lost, prev2_runs, prev_runs, prev_best = True, prev_runs, runs, None
            continue
        top, bot, cy = best
        if lost or c in via:
            # Reacquired, or branch forced: the old slope is stale, so restart here.
            ys, lost, merged = [], False, False
        elif prev_best is not None:
            if _overlaps(best, prev_runs) >= 2 and _overlaps(best, prev2_runs) >= 2:
                # This run joins two strokes seen in each of the last two columns (a
                # cursor's dots or a gridline crossing do not persist like that).
                merged = True
            elif merged and _overlaps(prev_best, runs) >= 2:
                merged = False  # the merged stroke has split again
        if merged and bot - top + 1 > THICK + 2.0:
            # Another curve crosses inside this run: keep to the prediction.
            cy = float(np.clip(pred, top, bot))
        ys.append(cy)
        found[c] = cy
        prev2_runs, prev_runs, prev_best = prev_runs, runs, best
    return found


def trace(grey, seed, c_lo, c_hi, via=()) -> tuple[np.ndarray, np.ndarray]:
    """Pixel (x, y) of the curve through `seed` (and `via`) for columns c_lo..c_hi."""
    ink = grey < DARK
    r0, r1, f0, f1 = _interior(ink)
    c_lo, c_hi = max(c_lo, f0), min(c_hi, f1)
    sx, sy = seed
    runs = _runs(grey[:, sx], r0, r1)
    if not runs:
        raise SystemExit(f"no ink in seed column {sx}")
    y0 = min(runs, key=lambda r: abs(r[2] - sy))[2]
    vias = {int(x): float(y) for x, y in via}
    found = {sx: y0}
    found |= _walk(grey, range(sx + 1, c_hi + 1), y0, r0, r1, vias)
    found |= _walk(grey, range(sx - 1, c_lo - 1, -1), y0, r0, r1, vias)
    xs = np.array(sorted(found), dtype=float)
    return xs, np.array([found[int(x)] for x in xs])


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("image", type=Path)
    ap.add_argument("--x", nargs=4, type=float, required=True, metavar=("PX1", "HZ1", "PX2", "HZ2"))
    ap.add_argument("--y", nargs=4, type=float, required=True, metavar=("PY1", "V1", "PY2", "V2"))
    ap.add_argument("--seed", nargs=2, type=int, required=True, metavar=("PX", "PY"))
    ap.add_argument("--via", nargs=2, type=int, action="append", default=[], metavar=("PX", "PY"),
                    help="pixel on the wanted branch where it parts from a merged stroke")
    ap.add_argument("--range", nargs=2, type=float, required=True, metavar=("FMIN", "FMAX"))
    ap.add_argument("--unit", required=True)
    ap.add_argument("--pdf", type=Path, required=True, help="source PDF, hashed into the header")
    ap.add_argument("--url", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--overlay", type=Path, required=True)
    a = ap.parse_args(argv)

    px1, hz1, px2, hz2 = a.x
    py1, v1, py2, v2 = a.y
    to_f = lambda x: hz1 * (hz2 / hz1) ** ((x - px1) / (px2 - px1))  # noqa: E731
    to_v = lambda y: v1 + (y - py1) * (v2 - v1) / (py2 - py1)  # noqa: E731
    to_x = lambda f: px1 + (px2 - px1) * np.log(f / hz1) / np.log(hz2 / hz1)  # noqa: E731
    to_y = lambda v: py1 + (v - v1) * (py2 - py1) / (v2 - v1)  # noqa: E731

    grey = _grey(a.image)
    fmin, fmax = a.range
    xs, ys = trace(grey, tuple(a.seed), int(np.floor(to_x(fmin))), int(np.ceil(to_x(fmax))), a.via)
    k = np.arange(np.ceil(24 * np.log2(fmin / 1000)), np.floor(24 * np.log2(fmax / 1000)) + 1)
    f = 1000.0 * 2.0 ** (k / 24)
    val = np.interp(np.log(f), np.log(to_f(xs)), to_v(ys))

    sha = hashlib.sha256(a.pdf.read_bytes()).hexdigest()
    cmd = list(argv) if argv is not None else sys.argv[1:]
    for p in (a.image, a.pdf, a.overlay):  # scratch paths: keep only the file name
        cmd = [p.name if s == str(p) else s for s in cmd]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w") as fh:
        fh.write(f"# source: {a.url} sha256={sha}\n")
        fh.write(f"# unit: {a.unit}\n")
        fh.write(f"# traced: {shlex.join(['tools/trace_curve.py', *cmd])}\n")
        fh.write("f,value\n")
        for fi, vi in zip(f, val):
            fh.write(f"{fi:.2f},{vi:.2f}\n")

    h, w = grey.shape
    fig = plt.figure(figsize=(2 * w / 100, 2 * h / 100), dpi=100)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.imshow(grey, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax.plot(xs, ys, color="red", lw=0.8, alpha=0.8)
    ax.plot(to_x(f), to_y(val), "o", ms=2.5, mfc="none", mec="dodgerblue", mew=0.8)
    ax.plot(*a.seed, "+", color="lime", ms=12, mew=2)
    for v in a.via:
        ax.plot(*v, "x", color="orange", ms=10, mew=2)
    ax.set_axis_off()
    a.overlay.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.overlay)
    plt.close(fig)
    print(f"{a.out}: {f.size} points {f[0]:.1f}-{f[-1]:.1f} Hz; columns {int(xs[0])}-{int(xs[-1])}")


if __name__ == "__main__":
    main()
