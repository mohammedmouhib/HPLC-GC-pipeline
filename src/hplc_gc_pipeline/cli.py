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
USER_ERRORS = (FileNotFoundError, NotADirectoryError, ValueError, FileExistsError)


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
        help="Auto-discover an experiment folder and write starter config files",
    )
    ip.add_argument("experiment_dir", type=Path, help="Path to the experiment folder")
    ip.add_argument(
        "--mode", choices=["auto", "hplc", "gc", "both"], default="auto",
        help="Which modalities to configure: auto=detect from data (default), "
             "hplc=HPLC only, gc=GC only, both=force both",
    )
    ip.add_argument(
        "--rt-min", type=float, default=1.5, metavar="MIN",
        help="GC peak detection: RT scan start in minutes (default: 1.5)",
    )
    ip.add_argument(
        "--rt-max", type=float, default=30.0, metavar="MAX",
        help="GC peak detection: RT scan end in minutes (default: 30.0)",
    )
    ip.add_argument(
        "--rt-margin", type=float, default=0.08, metavar="MIN",
        help="GC peak detection: margin added to each peak RT edge (default: 0.08 min)",
    )
    ip.add_argument(
        "--force", action="store_true",
        help="Overwrite existing config/CSV files without prompting",
    )
    ip.add_argument("--debug", action="store_true", help="Show full traceback on error")

    gp = sub.add_parser("gui", help="Launch the web GUI (needs the 'gui' extra)")
    gp.add_argument("experiment_dir", type=Path, nargs="?", default=None,
                    help="Experiment folder to preload (optional; can be set in the GUI)")
    gp.add_argument("--port", type=int, default=8501, help="Local port (default 8501)")
    return parser


def _run_init(args) -> None:
    """Handle the `hplc init` subcommand."""
    import shutil
    from .gc_init import (
        run_gc_init, write_gc_compounds, write_hplc_config, print_summary,
        detect_hplc_data_dir, find_gc_injections,
    )
    from .gc_calibration import GCCalibrationConfig

    exp_dir = Path(args.experiment_dir).expanduser().resolve()
    if not exp_dir.is_dir():
        raise NotADirectoryError(f"Experiment directory not found: {exp_dir}")

    config_path = exp_dir / "hplc_config.yaml"
    gc_compounds_path = exp_dir / "gc_compounds.csv"

    existing = [p for p in (config_path, gc_compounds_path) if p.exists()]
    if existing and not args.force:
        names = " and ".join(p.name for p in existing)
        raise FileExistsError(
            f"{names} already exist in {exp_dir}. "
            f"Use --force to overwrite."
        )

    std_pattern = GCCalibrationConfig.standard_pattern

    # Determine which modalities to configure
    mode = getattr(args, "mode", "auto")
    if mode == "auto":
        has_hplc = detect_hplc_data_dir(exp_dir) is not None or any(
            p.is_dir() and p.name.endswith(".D") and not (p / "data.ms").exists()
            for p in exp_dir.iterdir()
        )
        has_gc = bool(find_gc_injections(exp_dir))
        include_hplc = has_hplc or not has_gc  # default to HPLC when ambiguous
        include_gc = has_gc
    elif mode == "hplc":
        include_hplc, include_gc = True, False
    elif mode == "gc":
        include_hplc, include_gc = False, True
    else:  # "both"
        include_hplc, include_gc = True, True

    result = run_gc_init(
        exp_dir,
        std_pattern=std_pattern,
        rt_scan_min=args.rt_min,
        rt_scan_max=args.rt_max,
        rt_margin=args.rt_margin,
        include_gc=include_gc,
    )

    extra_files: list[Path] = []

    if include_gc:
        write_gc_compounds(result, gc_compounds_path)
    elif gc_compounds_path.exists():
        pass  # leave existing file

    write_hplc_config(result, config_path, std_pattern,
                      include_hplc=include_hplc, include_gc=include_gc)

    # Seed HPLC CSV files from package defaults when HPLC is active.
    if include_hplc:
        _gui_dir = Path(__file__).parent / "gui"
        hplc_compounds_path = exp_dir / "compounds.csv"
        standard_path = exp_dir / "standard.csv"
        for src, dst in (
            (_gui_dir / "compounds_default.csv", hplc_compounds_path),
            (_gui_dir / "standard_default.csv", standard_path),
        ):
            if not dst.exists() and src.exists():
                shutil.copy2(src, dst)
                extra_files.append(dst)

    print_summary(
        result,
        config_path,
        gc_compounds_path if include_gc else None,
        extra_files=extra_files,
        include_hplc=include_hplc,
        include_gc=include_gc,
    )


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
