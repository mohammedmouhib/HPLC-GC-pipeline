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

import copy
from pathlib import Path
from typing import Optional

import numpy as np
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

    dataset = processing.load_dataset(cfg.data_roots, cfg.processing)
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


# ---------------------------------------------------------------------------
# Targeted re-integration of a single injection
# ---------------------------------------------------------------------------

def reintegrate_injection(
    cfg: Config,
    folder_name: str,
    *,
    min_r2: Optional[float] = None,
    max_comps: Optional[int] = None,
    min_height: Optional[float] = None,
    min_prominence: Optional[float] = None,
    no_deconvolution: bool = False,
    manual_bounds: Optional[dict] = None,
    skip_analysis: bool = False,
) -> None:
    """Re-process one injection and patch the existing Stage 1 outputs in-place.

    Finds ``folder_name`` in the configured data roots, re-runs MOCCA2 (with
    optional parameter overrides that apply only to this injection), then:

    * replaces that folder's rows in ``peak_results.csv``
    * replaces that folder's rows in ``chromatogram_traces.csv``
    * re-renders the per-injection PNG and refreshes the HTML gallery
    * re-runs Stage 2 analysis (unless ``skip_analysis=True``)

    The rest of the dataset is untouched.
    """
    from .agilent import find_injection_folders, load_chromatogram

    if not cfg.peak_results_csv.exists():
        raise FileNotFoundError(
            f"{cfg.peak_results_csv} not found — run Stage 1 (process) first."
        )

    # ------------------------------------------------------------------
    # 1. Locate the target injection
    # ------------------------------------------------------------------
    target_inj = None
    for data_root in cfg.data_roots:
        for inj in find_injection_folders(data_root):
            if inj.folder_name == folder_name:
                target_inj = inj
                break
        if target_inj is not None:
            break

    if target_inj is None:
        raise ValueError(
            f"Injection folder '{folder_name}' not found in {cfg.data_roots}.\n"
            f"Check the folder name (include the .D suffix)."
        )

    print(f"\n{'=' * 70}")
    print(f"RE-INTEGRATING  {folder_name}  '{target_inj.sample_name}'")
    print(f"{'=' * 70}")

    # ------------------------------------------------------------------
    # 2. Build processing config (override only what was explicitly given)
    # ------------------------------------------------------------------
    proc_cfg = copy.deepcopy(cfg.processing)
    if min_r2 is not None:
        proc_cfg.deconvolution.min_r2 = min_r2
        print(f"  Override: min_r2         = {min_r2}")
    if max_comps is not None:
        proc_cfg.deconvolution.max_comps = max_comps
        print(f"  Override: max_comps      = {max_comps}")
    if min_height is not None:
        proc_cfg.peak_detection.min_height = min_height
        print(f"  Override: min_height     = {min_height}")
    if min_prominence is not None:
        proc_cfg.peak_detection.min_prominence = min_prominence
        print(f"  Override: min_prominence = {min_prominence}")
    if no_deconvolution:
        proc_cfg.deconvolution.enabled = False
        print("  Override: deconvolution  = disabled (manual trapezoid)")

    # ------------------------------------------------------------------
    # 3. Load and process the single chromatogram
    # ------------------------------------------------------------------
    chrom = load_chromatogram(target_inj.folder_path, proc_cfg.wavelength_nm)
    chrom.sample_name = target_inj.sample_name
    chrom.folder_name = target_inj.folder_name
    chrom.experiment_name = target_inj.experiment_name

    n_peaks = processing._process_single(chrom, proc_cfg.peak_detection, proc_cfg.deconvolution)
    print(f"  Detected {n_peaks} peaks")

    # ------------------------------------------------------------------
    # 4. Build peak rows for this injection
    # ------------------------------------------------------------------
    components = processing._safe_components(chrom)
    time_data = processing._get_time(chrom)
    experiment_name = getattr(chrom, "experiment_name", "Main")

    if components:
        new_rows = processing._rows_from_components(
            chrom, components, time_data, experiment_name, proc_cfg.wavelength_nm
        )
        method = "MOCCA2_Component"
    elif getattr(chrom, "peaks", None):
        new_rows = processing._rows_from_peaks(chrom, time_data, experiment_name, proc_cfg.wavelength_nm)
        method = "Manual_Trapezoid"
    else:
        new_rows = []
        method = "(none)"
    print(f"  Integration: {method}  ({len(new_rows)} rows)")

    # ------------------------------------------------------------------
    # 5. Patch peak_results.csv
    # ------------------------------------------------------------------
    peaks_df = pd.read_csv(cfg.peak_results_csv)
    n_old = int((peaks_df["Folder"] == folder_name).sum())
    peaks_df = peaks_df[peaks_df["Folder"] != folder_name].copy()

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        peaks_df = pd.concat([peaks_df, new_df], ignore_index=True)

    # Recompute Area_Percent across the full (patched) dataset
    totals = peaks_df.groupby(["Sample", "Folder"])["Area_Integral"].transform("sum")
    with np.errstate(invalid="ignore", divide="ignore"):
        peaks_df["Area_Percent"] = np.where(
            (totals > 0) & peaks_df["Area_Integral"].notna(),
            peaks_df["Area_Integral"] / totals * 100,
            np.nan,
        )
        peaks_df["Area_Percent"] = peaks_df["Area_Percent"].round(2)

    peaks_df.to_csv(cfg.peak_results_csv, index=False)
    print(f"  peak_results.csv: replaced {n_old} old rows with {len(new_rows)} new rows")

    # ------------------------------------------------------------------
    # 5b. Apply manual RT-bound integrations (optional)
    # ------------------------------------------------------------------
    _manual_plot_segs: list = []   # accumulated for preview plot in step 7
    if manual_bounds:
        time_data = processing._get_time(chrom)
        # _get_signal returns chrom.data (single wavelength after extract_wavelength),
        # which is already baseline-corrected in-place by chrom.correct_baseline() in
        # _process_single — no additional baseline subtraction needed.
        signal_raw = processing._get_signal(chrom)
        if time_data is None or signal_raw is None:
            print("  Warning: cannot apply manual bounds — chromatogram signal unavailable")
        else:
            peaks_df = pd.read_csv(cfg.peak_results_csv)
            for compound_name, (left_rt, right_rt) in manual_bounds.items():
                mask = (time_data >= left_rt) & (time_data <= right_rt)
                if not mask.any():
                    print(f"  Manual bounds {compound_name}: no data in {left_rt:.3f}–{right_rt:.3f} min — skipped")
                    continue
                seg_t = time_data[mask]
                seg_s = signal_raw[mask]  # already baseline-corrected
                # np.sum matches MOCCA's component.integral = np.sum(concentration) convention;
                # np.trapezoid gives mAU·min which is ~50–100× smaller and breaks calibration.
                area = float(np.sum(seg_s))
                apex_idx = int(np.argmax(seg_s))
                apex_rt = float(seg_t[apex_idx])
                height = float(seg_s[apex_idx])

                # remove any existing rows that fall within this RT window
                in_window = (
                    (peaks_df["Folder"] == folder_name) &
                    (peaks_df["Retention_Time_min"] >= left_rt) &
                    (peaks_df["Retention_Time_min"] <= right_rt)
                )
                n_removed = int(in_window.sum())
                peaks_df = peaks_df[~in_window].copy()

                new_row = {
                    "Experiment":         experiment_name,
                    "Folder":             folder_name,
                    "Sample":             target_inj.sample_name,
                    "Peak_Number":        999,
                    "Retention_Time_min": round(apex_rt, 3),
                    "Area_Integral":      round(area, 2),
                    "Height_mAU":         round(height, 2),
                    "Width_min":          round(float(right_rt - left_rt), 3),
                    "Wavelength_nm":      proc_cfg.wavelength_nm,
                    "Integration_Method": "Manual_Bounds",
                    "Area_Percent":       float("nan"),
                }
                peaks_df = pd.concat([peaks_df, pd.DataFrame([new_row])], ignore_index=True)
                _manual_plot_segs.append(
                    (compound_name, left_rt, right_rt, seg_t, seg_s, area, apex_rt, height)
                )
                print(
                    f"  Manual bounds {compound_name}: {left_rt:.3f}–{right_rt:.3f} min "
                    f"→ area {area:.2f}, apex {apex_rt:.3f} min  (removed {n_removed} row(s))"
                )

            # recompute Area_Percent across the full patched dataset
            totals = peaks_df.groupby(["Sample", "Folder"])["Area_Integral"].transform("sum")
            with np.errstate(invalid="ignore", divide="ignore"):
                peaks_df["Area_Percent"] = np.where(
                    (totals > 0) & peaks_df["Area_Integral"].notna(),
                    peaks_df["Area_Integral"] / totals * 100,
                    np.nan,
                )
                peaks_df["Area_Percent"] = peaks_df["Area_Percent"].round(2)
            peaks_df.to_csv(cfg.peak_results_csv, index=False)
            print(f"  peak_results.csv: updated with manual bounds")

    # ------------------------------------------------------------------
    # 6. Patch chromatogram_traces.csv (if it exists)
    # ------------------------------------------------------------------
    traces_csv = cfg.results_dir / "chromatogram_traces.csv"
    if traces_csv.exists():
        traces_df = pd.read_csv(traces_csv)
        traces_df = traces_df[traces_df["Folder"] != folder_name].copy()
        if time_data is not None:
            corrected = pc._corrected_signal(chrom)
            new_trace_rows = [
                {
                    "Folder": folder_name,
                    "Sample": target_inj.sample_name,
                    "Time_min": round(float(t), 5),
                    "Signal_mAU": round(float(s), 4),
                }
                for t, s in zip(time_data, corrected)
            ]
            traces_df = pd.concat([traces_df, pd.DataFrame(new_trace_rows)], ignore_index=True)
        traces_df.to_csv(traces_csv, index=False)
        print(f"  chromatogram_traces.csv: patched")

    # ------------------------------------------------------------------
    # 7. Re-render the per-injection PNG and refresh the HTML gallery
    # ------------------------------------------------------------------
    plots_dir = cfg.results_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    safe_name = chrom.sample_name.replace(" ", "_").replace(",", "").replace("/", "-")
    safe_exp = experiment_name.replace(" ", "_").replace("/", "-")
    plot_path = plots_dir / f"{safe_exp}_{folder_name.replace('.D', '')}_{safe_name}.png"
    pc.plot_chromatogram(chrom, proc_cfg.wavelength_nm, plot_path)
    print(f"  Plot updated: {plot_path.name}")

    if _manual_plot_segs:
        _preview_path = plot_path.with_name(plot_path.stem + "_manual_preview.png")
        pc.plot_manual_bounds_preview(
            time_data=time_data,
            signal_raw=signal_raw,
            segments=_manual_plot_segs,
            sample_name=chrom.sample_name,
            folder_name=folder_name,
            wavelength=proc_cfg.wavelength_nm,
            save_path=_preview_path,
        )
        print(f"  Manual bounds preview: {_preview_path.name}")

    pc.export_html_gallery(plots_dir, cfg.results_dir)

    # ------------------------------------------------------------------
    # 8. Re-run Stage 2 analysis
    # ------------------------------------------------------------------
    if not skip_analysis:
        print()
        _run_hplc_analysis(cfg)
    else:
        print("\n  Skipping Stage 2 (--skip-analysis). Run 'hplc analyze' when ready.")


