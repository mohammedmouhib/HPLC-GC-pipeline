"""
Command-line interface.

    hplc process  <experiment_dir> [--config FILE]   # Stage 1
    hplc analyze  <experiment_dir> [--config FILE]   # Stage 2
    hplc run      <experiment_dir> [--config FILE]   # both stages

The experiment directory is the only required argument; all tunable values come
from its ``hplc_config.yaml`` (or the file given with --config). Outputs are
written back into the experiment directory.
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

from .config import load_config
from . import pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hplc",
        description="HPLC (MOCCA2) data-processing pipeline. Runs against any "
                    "experiment folder; configuration lives in that folder's "
                    "hplc_config.yaml.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    for name, help_text in [
        ("process", "Stage 1: raw .D folders -> peak table + chromatogram plots"),
        ("analyze", "Stage 2: peak table -> statistics + Bokeh dashboard"),
        ("run", "Run both stages in sequence"),
    ]:
        sp = sub.add_parser(name, help=help_text)
        sp.add_argument("experiment_dir", type=Path, help="Path to the experiment folder")
        sp.add_argument("--config", type=Path, default=None,
                        help="Config YAML (default: <experiment_dir>/hplc_config.yaml)")
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        cfg = load_config(args.experiment_dir, args.config)
        if args.command == "process":
            pipeline.run_processing(cfg)
        elif args.command == "analyze":
            pipeline.run_analysis(cfg)
        elif args.command == "run":
            pipeline.run_all(cfg)
    except Exception as exc:
        print(f"\n{'=' * 70}\nERROR\n{'=' * 70}\n{exc}\n", file=sys.stderr)
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
