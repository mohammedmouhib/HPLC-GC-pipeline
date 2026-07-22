"""
Pipeline orchestration.

An experiment can include HPLC data, GC-MS data, or both (``cfg.modalities``).
:func:`run_processing` and :func:`run_analysis` dispatch to whichever modalities
are active; when both are present, :func:`run_analysis` also builds a combined
dashboard showing HPLC and GC compounds together.

Per modality, two stages write back into the experiment folder:

* Stage 1 (processing) -- raw .D folders -> peak table + plots.
* Stage 2 (analysis)   -- peak table -> stats + Bokeh dashboard.
"""

from __future__ import annotations

import pandas as pd

from .config import Config
from . import processing
from . import plotting_chromatograms as pc
from . import sample_names, compounds as compounds_mod, calibration, consolidate, stats
from . import plotting_dashboard as pd_plots
from . import gc_pipeline


# ---------------------------------------------------------------------------
# Top-level dispatch
# ---------------------------------------------------------------------------

def run_processing(cfg: Config) -> None:
    """Stage 1 for every active modality."""
    if "hplc" in cfg.modalities:
        _run_hplc_processing(cfg)
    if "gc" in cfg.modalities:
        gc_pipeline.run_gc_processing(cfg)


def run_analysis(cfg: Config) -> None:
    """Stage 2 for every active modality (+ combined dashboard when both)."""
    results: dict[str, tuple[pd.DataFrame, pd.DataFrame]] = {}
    if "hplc" in cfg.modalities:
        results["hplc"] = _run_hplc_analysis(cfg)
    if "gc" in cfg.modalities:
        results["gc"] = gc_pipeline.run_gc_analysis(cfg)

    if "hplc" in results and "gc" in results:
        _combined_dashboard(cfg, results["hplc"], results["gc"])


def run_all(cfg: Config) -> None:
    run_processing(cfg)
    run_analysis(cfg)


# ---------------------------------------------------------------------------
# HPLC stages
# ---------------------------------------------------------------------------

def _run_hplc_processing(cfg: Config) -> None:
    """Stage 1 (HPLC): process raw .D folders into a peak table + plots."""
    results_dir = cfg.results_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    dataset = processing.load_dataset(cfg.data_root, cfg.processing)
    dataset = processing.process_dataset(dataset, cfg.processing)

    print("\n" + "=" * 70)
    print("EXPORTING STAGE 1 RESULTS")
    print("=" * 70)

    peaks_df = processing.build_peak_table(dataset)
    if not peaks_df.empty:
        peaks_df.to_csv(cfg.peak_results_csv, index=False)
        print(f"  Peak table: {cfg.peak_results_csv.name} ({len(peaks_df)} peaks)")

    summary_df = processing.build_sample_summary(dataset)
    summary_df.to_csv(results_dir / "sample_summary.csv", index=False)
    print(f"  Sample summary: sample_summary.csv ({len(summary_df)} injections)")

    # Per-injection PNGs
    plots_dir = results_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    print(f"\n  Rendering {len(dataset.chromatograms)} per-injection plots...")
    for chrom in dataset.chromatograms:
        safe_name = chrom.sample_name.replace(" ", "_").replace(",", "").replace("/", "-")
        safe_exp = getattr(chrom, "experiment_name", "Main").replace(" ", "_").replace("/", "-")
        plot_path = plots_dir / f"{safe_exp}_{chrom.folder_name.replace('.D', '')}_{safe_name}.png"
        pc.plot_chromatogram(chrom, dataset.wavelength, plot_path)

    pc.export_chromatogram_traces_csv(dataset, results_dir)
    pc.export_html_gallery(plots_dir, results_dir)
    pc.plot_overlay(dataset, results_dir)
    pc.plot_overlay_interactive_bokeh(dataset, results_dir, palette=cfg.plots.palette)

    print("\nStage 1 complete. Results in:", results_dir)


def _run_hplc_analysis(cfg: Config) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stage 2 (HPLC): analyse the peak table into stats + plots.

    Returns ``(consolidated_df, stats_df)`` for the combined dashboard.
    """
    if not cfg.peak_results_csv.exists():
        raise FileNotFoundError(
            f"{cfg.peak_results_csv} not found. Run Stage 1 (process) first."
        )

    acfg = cfg.analysis

    df = sample_names.load_peak_data(cfg.peak_results_csv, acfg.sample_name.pattern)
    df = sample_names.pool_replicates_across_experiments(df)

    compounds_df = compounds_mod.load_compounds(cfg.compounds_path)
    standards_df = calibration.load_standards(
        cfg.standard_path,
        peak_results_csv=cfg.peak_results_csv,
        compounds_df=compounds_df,
        cfg=acfg.calibration,
    )
    cal = calibration.build_calibration_functions(standards_df, compounds_df, cfg=acfg.calibration)

    assigned_df = compounds_mod.assign_compounds(df, compounds_df, exclude_unknown=acfg.exclude_unknown)
    consolidated_df = consolidate.consolidate_compound_areas(assigned_df, compounds_df, all_samples_df=df)
    consolidated_df = calibration.add_concentration_column(consolidated_df, cal, cfg=acfg.calibration)

    stats_df = stats.compute_replicate_stats(consolidated_df, cv_warning_threshold=acfg.cv_warning_threshold)

    analysis_dir = cfg.analysis_dir
    pd_plots.export_analysis_results(assigned_df, consolidated_df, stats_df, analysis_dir)

    stats_for_bars = stats_df
    if cfg.plots.bar_time_points is not None:
        stats_for_bars = stats_df[stats_df["Time_Point_h"].isin(cfg.plots.bar_time_points)]

    pd_plots.generate_plots(
        consolidated_df, stats_for_bars, analysis_dir,
        palette=cfg.plots.palette,
        strain_order=cfg.plots.strain_order,
        plot_strains=cfg.plots.plot_strains,
        n_cols=cfg.plots.n_cols,
    )

    print("\nStage 2 complete. Analysis in:", analysis_dir)
    return consolidated_df, stats_df


# ---------------------------------------------------------------------------
# Combined dashboard (HPLC + GC)
# ---------------------------------------------------------------------------

def _combined_dashboard(cfg: Config,
                        hplc_res: tuple[pd.DataFrame, pd.DataFrame],
                        gc_res: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    """Union HPLC and GC results into one modality-tagged dashboard.

    Compounds are disjoint across instruments (each compound is quantified by
    one technique), so the union simply shows every compound. Concentration
    (uM) is the cross-comparable unit; peak-area bar charts mix instrument
    units and are kept only as a per-compound reference.
    """
    print("\n" + "=" * 70)
    print("BUILDING COMBINED HPLC + GC DASHBOARD")
    print("=" * 70)

    # Compounds are normally quantified by a single instrument, but if the same
    # name appears in both tables (e.g. a compound run on both), suffix it with
    # the modality so the two don't merge into one ambiguous plot.
    overlap = set(hplc_res[0]["Compound"]) & set(gc_res[0]["Compound"])
    if overlap:
        print(f"  Note: {sorted(overlap)} quantified by both instruments; "
              f"suffixing compound names with modality in the combined view.")

    def _tag(frames, label):
        cons, sts = frames
        cons = cons.copy(); cons["Modality"] = label
        sts = sts.copy(); sts["Modality"] = label
        if overlap:
            for f in (cons, sts):
                m = f["Compound"].isin(overlap)
                f.loc[m, "Compound"] = f.loc[m, "Compound"].astype(str) + f" ({label})"
        return cons, sts

    h_cons, h_stats = _tag(hplc_res, "HPLC")
    g_cons, g_stats = _tag(gc_res, "GC")
    consolidated = pd.concat([h_cons, g_cons], ignore_index=True)
    stats_df = pd.concat([h_stats, g_stats], ignore_index=True)

    combined_dir = cfg.experiment_dir / "combined_analysis"
    combined_dir.mkdir(parents=True, exist_ok=True)
    consolidated.to_csv(combined_dir / "consolidated_peaks_combined.csv", index=False)
    stats_df.to_csv(combined_dir / "replicate_summary_combined.csv", index=False)
    print(f"  Combined tables: {combined_dir.name}/consolidated_peaks_combined.csv "
          f"({len(consolidated)} rows), replicate_summary_combined.csv ({len(stats_df)} rows)")

    stats_for_bars = stats_df
    if cfg.plots.bar_time_points is not None:
        stats_for_bars = stats_df[stats_df["Time_Point_h"].isin(cfg.plots.bar_time_points)]

    pd_plots.generate_plots(
        consolidated, stats_for_bars, combined_dir,
        palette=cfg.plots.palette,
        strain_order=cfg.plots.strain_order,
        plot_strains=cfg.plots.plot_strains,
        n_cols=cfg.plots.n_cols,
    )
    print("\nCombined dashboard in:", combined_dir)
