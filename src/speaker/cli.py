"""Command line entry point: ``speaker build``."""
from __future__ import annotations

import argparse
from pathlib import Path

from .report import build


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="speaker")
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build", help="compute, plot, export and summarise the design")
    b.add_argument("--params", type=Path, default=None, help="parameter TOML (default: params/default.toml)")
    b.add_argument("--out", type=Path, default=Path("out"), help="output directory (default: out/)")
    args = parser.parse_args(argv)
    summary = build(args.params, args.out)
    print(f"wrote {summary}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
