"""Command line entry point: ``speaker build``, ``speaker directivity`` and ``speaker round``."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .report import build


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="speaker")
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build", help="compute, plot, export and summarise the design")
    b.add_argument("--params", type=Path, default=None, help="parameter TOML (default: params/default.toml)")
    b.add_argument("--out", type=Path, default=Path("out"), help="output directory (default: out/)")
    d = sub.add_parser("directivity", help="3D directivity simulation, scores and HTML report")
    d.add_argument("--params", type=Path, default=None, help="parameter TOML (default: params/default.toml)")
    d.add_argument("--out", type=Path, default=Path("out/directivity"),
                   help="output directory (default: out/directivity/)")
    d.add_argument("--fast", action="store_true", help="coarse mesh, 50 Hz-1 kHz, no puck (minutes, not hours)")
    d.add_argument("--mesh-check", action="store_true",
                   help="only run the mesh-convergence check, writing mesh_check.json (hours)")
    d.add_argument("--trusted-limit", type=float, default=None, metavar="HZ",
                   help="trusted 3D limit in Hz, in place of mesh_check.json's")
    r = sub.add_parser("round", help="run a round file's stages and finals (resumes finished variants)")
    r.add_argument("round_toml", type=Path, metavar="ROUND_TOML", help="round file, e.g. rounds/round-1.toml")
    r.add_argument("--out", type=Path, default=None, help="output directory (default: out/rounds/<number>/)")
    args = parser.parse_args(argv)
    try:
        if args.command == "build":
            out = build(args.params, args.out)
        elif args.command == "round":
            from .round import run_round
            out = run_round(args.round_toml, args.out)
        elif args.mesh_check:
            from .directivity_report import mesh_check
            out = mesh_check(args.params, args.out)
        else:
            from .directivity_report import run_directivity
            given = None if args.trusted_limit is None else (args.trusted_limit, "given on the command line")
            out = run_directivity(args.params, args.out, fast=args.fast, trusted_limit=given)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
