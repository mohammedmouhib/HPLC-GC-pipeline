"""
GC pipeline orchestration.

Mirrors the HPLC stages but for GC-MS:

* Stage 1 (:func:`run_gc_processing`) -- raw GC .D folders -> gc_results/peak_results.csv
  (targeted integration on the quantifier-ion EIC; see :mod:`gc_processing`).
* Stage 2 (:func:`run_gc_analysis`)   -- peak table -> calibration + concentrations,
  replicate statistics, GC dashboard.

The peak table is already one row per (injection x compound) with undetected
compounds zero-filled, so GC skips the HPLC consolidation step: after sample-name
parsing (which drops the standard injections) the table already has the shape the
shared stats/dashboard code expects.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import Config
from . import gc_processing, gc_calibration, sample_names, gc_plotting
from .gc_agilent import find_gc_injections
from . import stats as stats_mod
from . import plotting_dashboard as dash


def run_gc_processing(cfg: Config) -> None:
    """Stage 1 (GC): integrate every compound in every GC injection."""
    results_dir = cfg.gc_results_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    injections = find_gc_injections(cfg.gc_data_root)
    if not injections:
        raise ValueError(f"No GC-MS .D folders (with data.ms) found in {cfg.gc_data_root}")
    print(f"Found {len(injections)} GC-MS injections in {cfg.gc_data_root}")
    compounds = gc_processing.load_gc_compounds(cfg.gc_compounds_path)
    peak_df = gc_processing.build_gc_peak_table(injections, compounds, cfg.gc.processing)

    print("\n" + "=" * 70)
    print("EXPORTING GC STAGE 1 RESULTS")
    print("=" * 70)
    peak_df.to_csv(cfg.gc_peak_results_csv, index=False)
    print(f"  GC peak table: {cfg.gc_peak_results_csv.name} ({len(peak_df)} rows)")

    gc_plotting.render_gc_plots(injections, compounds, cfg.gc.processing, results_dir)
    print("\nGC Stage 1 complete. Results in:", results_dir)


def _to_consolidated(df: pd.DataFrame) -> pd.DataFrame:
    """Reshape the parsed GC sample rows into the shared consolidated schema.

    The GC peak table is already consolidated (one row per injection x compound,
    undetected zero-filled), so this only maps GC column names onto the ones the
    shared stats/dashboard code expects (``Height_mAU``, ``Area_Percent``) and
    keeps the concentration column produced by calibration.
    """
    df = df.copy()
    # GC "Height" (ion counts) reuses the shared column name; the unit differs
    # from HPLC mAU but the downstream code only reads the column, not the unit.
    if "Height_mAU" not in df.columns:
        df["Height_mAU"] = df["Height"] if "Height" in df.columns else 0.0
    totals = df.groupby("Folder")["Area_Integral"].transform("sum")
    df["Area_Percent"] = np.where(totals > 0, (df["Area_Integral"] / totals * 100).round(2), 0.0)
    if "Concentration_uM" not in df.columns:
        df["Concentration_uM"] = np.nan
    return df


def run_gc_analysis(cfg: Config) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stage 2 (GC): calibrate, quantify, and summarise into stats + a dashboard.

    Returns ``(consolidated_df, stats_df)`` so a caller can build a combined
    HPLC+GC dashboard.
    """
    if not cfg.gc_peak_results_csv.exists():
        raise FileNotFoundError(
            f"{cfg.gc_peak_results_csv} not found. Run GC Stage 1 (process) first."
        )
    gc = cfg.gc
    peak_df = pd.read_csv(cfg.gc_peak_results_csv)

    # Calibration (selectable source) + concentrations (dilution applied to samples).
    calibration = gc_calibration.build_calibration(peak_df, gc.calibration, cfg.experiment_dir)
    quant_df = gc_calibration.add_concentrations(peak_df, calibration, gc.calibration)

    analysis_dir = cfg.gc_analysis_dir
    analysis_dir.mkdir(parents=True, exist_ok=True)
    quantified_csv = analysis_dir / "quantified_peaks.csv"
    quant_df.to_csv(quantified_csv, index=False)
    print(f"\n  Quantified peaks (all injections): {quantified_csv.name} ({len(quant_df)} rows)")

    # Parse sample names -> drops the standard injections (they don't match the
    # sample pattern), then pool replicate letters across runs.
    df = sample_names.load_peak_data(quantified_csv, gc.sample_name.pattern)
    df = sample_names.pool_replicates_across_experiments(df)

    consolidated = _to_consolidated(df)
    stats_df = stats_mod.compute_replicate_stats(consolidated, cv_warning_threshold=gc.cv_warning_threshold)

    dash.export_analysis_results(df, consolidated, stats_df, analysis_dir)
    dash.generate_plots(
        consolidated, stats_df, analysis_dir,
        palette=cfg.plots.palette,
        strain_order=cfg.plots.strain_order,
        plot_strains=cfg.plots.plot_strains,
        n_cols=cfg.plots.n_cols,
    )
    print("\nGC Stage 2 complete. Analysis in:", analysis_dir)
    return consolidated, stats_df
