"""
Command-line interface.

    hplc process  <experiment_dir> [--config FILE] [--modality M]   # Stage 1
    hplc analyze  <experiment_dir> [--config FILE] [--modality M]   # Stage 2
    hplc run      <experiment_dir> [--config FILE] [--modality M]   # both stages

The experiment directory is the only required argument; all tunable values come
from its ``hplc_config.yaml`` (or the file given with --config). Outputs are
written back into the experiment directory. ``--modality`` restricts a run to
one instrument (hplc / gc); by default every modality the config defines runs.
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

from .config import load_config
from . import pipeline

# Errors caused by user input (bad paths, missing or invalid config). These get
# a clean one-line message; anything else is treated as an unexpected bug and
# shown with a full traceback (or via --debug).
USER_ERRORS = (FileNotFoundError, NotADirectoryError, ValueError)


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
        sp.add_argument("--modality", choices=["hplc", "gc", "both"], default=None,
                        help="Restrict to one modality (default: every modality the "
                             "config defines).")
        sp.add_argument("--debug", action="store_true",
                        help="Show the full traceback on error")

    ip = sub.add_parser(
        "init",
        help="Auto-discover a GC-MS experiment folder and write starter config files",
    )
    ip.add_argument("experiment_dir", type=Path, help="Path to the experiment folder")
    ip.add_argument(
        "--rt-min", type=float, default=1.5, metavar="MIN",
        help="Start of RT scan window in minutes (default: 1.5; skips solvent front)",
    )
    ip.add_argument(
        "--rt-max", type=float, default=30.0, metavar="MAX",
        help="End of RT scan window in minutes (default: 30.0)",
    )
    ip.add_argument(
        "--rt-margin", type=float, default=0.08, metavar="MIN",
        help="Margin added to each side of a detected peak RT window (default: 0.08 min)",
    )
    ip.add_argument(
        "--force", action="store_true",
        help="Overwrite existing hplc_config.yaml / gc_compounds.csv without prompting",
    )
    ip.add_argument("--debug", action="store_true", help="Show full traceback on error")

    gp = sub.add_parser("gui", help="Launch the web GUI (needs the 'gui' extra)")
    gp.add_argument("experiment_dir", type=Path, nargs="?", default=None,
                    help="Experiment folder to preload (optional; can be set in the GUI)")
    gp.add_argument("--port", type=int, default=8501, help="Local port (default 8501)")
    return parser


def _run_init(args) -> None:
    """Handle the `hplc init` subcommand."""
    from .gc_init import run_gc_init, write_gc_compounds, write_hplc_config, print_summary
    from .gc_calibration import GCCalibrationConfig

    exp_dir = Path(args.experiment_dir).expanduser().resolve()
    if not exp_dir.is_dir():
        raise NotADirectoryError(f"Experiment directory not found: {exp_dir}")

    config_path = exp_dir / "hplc_config.yaml"
    compounds_path = exp_dir / "gc_compounds.csv"

    existing = [p for p in (config_path, compounds_path) if p.exists()]
    if existing and not args.force:
        names = " and ".join(p.name for p in existing)
        raise FileExistsError(
            f"{names} already exist in {exp_dir}. "
            f"Use --force to overwrite."
        )

    std_pattern = GCCalibrationConfig.standard_pattern

    result = run_gc_init(
        exp_dir,
        std_pattern=std_pattern,
        rt_scan_min=args.rt_min,
        rt_scan_max=args.rt_max,
        rt_margin=args.rt_margin,
    )

    write_gc_compounds(result, compounds_path)
    write_hplc_config(result, config_path, std_pattern)
    print_summary(result, config_path, compounds_path)


def _apply_modality_filter(cfg, modality) -> None:
    """Narrow cfg.modalities to the one requested on the CLI, if any."""
    if modality is None or modality == "both":
        return
    if modality not in cfg.modalities:
        raise ValueError(
            f"--modality {modality} requested, but the config defines no "
            f"'{modality}' block (active: {cfg.modalities})."
        )
    cfg.modalities = [modality]


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "gui":
        from .gui.launch import launch_gui
        return launch_gui(args.experiment_dir, port=args.port)

    try:
        if args.command == "init":
            _run_init(args)
            return 0
        cfg = load_config(args.experiment_dir, args.config)
        _apply_modality_filter(cfg, args.modality)
        if args.command == "process":
            pipeline.run_processing(cfg)
        elif args.command == "analyze":
            pipeline.run_analysis(cfg)
        elif args.command == "run":
            pipeline.run_all(cfg)
    except USER_ERRORS as exc:
        if args.debug:
            traceback.print_exc()
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        # Unexpected failure: show the traceback so it can be reported/fixed.
        print(f"\n{'=' * 70}\nUNEXPECTED ERROR\n{'=' * 70}\n{exc}\n", file=sys.stderr)
        traceback.print_exc()
        print("\nRe-run with --debug for the full traceback, or report this.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
