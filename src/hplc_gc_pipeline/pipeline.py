"""
Pipeline orchestration.

Two stages, each writing back into the experiment folder:

* Stage 1 (:func:`run_processing`) -- raw .D folders -> results/peak_results.csv,
  chromatogram plots and overlays.
* Stage 2 (:func:`run_analysis`)   -- peak_results.csv + compounds/standards ->
  analysis/ CSVs and the Bokeh dashboard.

:func:`run_all` runs both in sequence.
"""

from __future__ import annotations

from .config import Config
from . import processing
from . import plotting_chromatograms as pc
from . import sample_names, compounds as compounds_mod, calibration, consolidate, stats
from . import plotting_dashboard as pd_plots


def run_processing(cfg: Config) -> None:
    """Stage 1: process raw .D folders into a peak table + plots."""
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


def run_analysis(cfg: Config) -> None:
    """Stage 2: analyse the peak table into stats + plots."""
    if not cfg.peak_results_csv.exists():
        raise FileNotFoundError(
            f"{cfg.peak_results_csv} not found. Run Stage 1 (process) first."
        )

    acfg = cfg.analysis

    df = sample_names.load_peak_data(cfg.peak_results_csv, acfg.sample_name.pattern)
    df = sample_names.pool_replicates_across_experiments(df)

    compounds_df = compounds_mod.load_compounds(cfg.compounds_path)
    standards_df = calibration.load_standards(cfg.standard_path)
    cal = calibration.build_calibration_functions(standards_df, compounds_df)

    assigned_df = compounds_mod.assign_compounds(df, compounds_df, exclude_unknown=acfg.exclude_unknown)
    consolidated_df = consolidate.consolidate_compound_areas(assigned_df, compounds_df, all_samples_df=df)
    consolidated_df = calibration.add_concentration_column(consolidated_df, cal, acfg.calibration)

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


def run_all(cfg: Config) -> None:
    run_processing(cfg)
    run_analysis(cfg)
