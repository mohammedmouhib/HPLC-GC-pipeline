"""
Stage 1 core: load .D folders, run MOCCA2 processing, extract a peak table.

This is a faithful refactor of the original ``HPLC_MOCCA_script.py`` processing
path. Two substantive changes:

* The duplicate ``process_single_chrom`` is gone. The original file defined it
  twice and Python silently used the second (deconvolution) version; only that
  behaviour is kept here.
* Peak-detection and deconvolution parameters come from the config instead of
  being hardcoded/commented out.

The integration logic (MOCCA2 components first, manual trapezoid fallback) is
unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .agilent import Injection, find_injection_folders, load_chromatogram
from .config import DeconvolutionConfig, PeakDetectionConfig, ProcessingConfig


@dataclass
class Dataset:
    """Loaded + processed chromatograms plus run-level metadata."""
    chromatograms: list = field(default_factory=list)
    blank: object = None
    blank_name: Optional[str] = None
    wavelength: int = 0


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_dataset(experiment_dir: "Path | list[Path]", cfg: ProcessingConfig) -> Dataset:
    """Discover and load every chromatogram under ``experiment_dir``.

    ``experiment_dir`` may be a single :class:`Path` or a list of paths
    (e.g. when HPLC data spans multiple subfolders).  Injections from all
    directories are combined into a single :class:`Dataset`.
    """
    dirs: list[Path] = (
        [experiment_dir] if isinstance(experiment_dir, Path) else list(experiment_dir)
    )

    print("=" * 70)
    print("LOADING DATA")
    print("=" * 70)
    if len(dirs) == 1:
        print(f"Experiment directory : {dirs[0]}")
    else:
        for i, d in enumerate(dirs, 1):
            print(f"Experiment directory {i}: {d}")
    print(f"Blank folder         : {cfg.blank_folder_name or '(none)'}")
    print(f"Wavelength           : {cfg.wavelength_nm} nm\n")

    injections = []
    for d in dirs:
        injections.extend(find_injection_folders(d))
    if not injections:
        raise ValueError(
            f"No valid .D folders with sample.xml found in "
            f"{', '.join(str(d) for d in dirs)}"
        )

    print(f"Found {len(injections)} .D folders. Loading chromatograms:")

    dataset = Dataset(wavelength=cfg.wavelength_nm, blank_name=cfg.blank_folder_name)
    for i, inj in enumerate(injections, 1):
        try:
            chrom = load_chromatogram(inj.folder_path, cfg.wavelength_nm)
            chrom.sample_name = inj.sample_name
            chrom.folder_name = inj.folder_name
            chrom.experiment_name = inj.experiment_name

            dataset.chromatograms.append(chrom)
            if cfg.blank_folder_name and inj.folder_name == cfg.blank_folder_name:
                dataset.blank = chrom

            print(f"  [{i}/{len(injections)}] OK  {inj.folder_name}  '{inj.sample_name}'")
        except Exception as exc:
            print(f"  [{i}/{len(injections)}] ERR {inj.folder_name}: {exc}")

    if not dataset.chromatograms:
        raise ValueError("No chromatograms were successfully loaded")

    if cfg.blank_folder_name and dataset.blank is None:
        print(f"\n  Warning: blank folder '{cfg.blank_folder_name}' not found among injections")

    print(f"\nLoaded {len(dataset.chromatograms)} chromatograms")
    return dataset


# ---------------------------------------------------------------------------
# Processing
# ---------------------------------------------------------------------------

def process_dataset(dataset: Dataset, cfg: ProcessingConfig) -> Dataset:
    """Baseline-correct, peak-pick, and (optionally) deconvolve every trace."""
    print("\n" + "=" * 70)
    print("PROCESSING CHROMATOGRAMS")
    print("=" * 70)

    if dataset.blank is not None:
        print(f"Blank: {dataset.blank_name}", end="")
        try:
            n = _process_single(dataset.blank, cfg.peak_detection, cfg.deconvolution)
            print(f" -> {n} peaks")
        except Exception as exc:
            print(f" -> error: {exc}")

    print("Samples:")
    for i, chrom in enumerate(dataset.chromatograms, 1):
        print(f"  [{i}/{len(dataset.chromatograms)}] {chrom.folder_name} '{chrom.sample_name}'", end="")
        try:
            n = _process_single(chrom, cfg.peak_detection, cfg.deconvolution)
            print(f" -> {n} peaks")
        except Exception as exc:
            print(f" -> error: {exc}")

    print("\nProcessing complete")
    return dataset


def _process_single(chrom, peak_cfg: PeakDetectionConfig, decon_cfg: DeconvolutionConfig) -> int:
    """Process one chromatogram; return the number of detected peaks."""
    chrom.correct_baseline()

    # Only pass through parameters the user actually set; otherwise MOCCA2
    # picks its own defaults (matching the original commented-out behaviour).
    find_kwargs = {}
    if peak_cfg.min_height is not None:
        find_kwargs["min_height"] = peak_cfg.min_height
    if peak_cfg.min_prominence is not None:
        find_kwargs["min_prominence"] = peak_cfg.min_prominence
    if peak_cfg.min_width is not None:
        find_kwargs["min_width"] = peak_cfg.min_width
    if peak_cfg.distance is not None:
        find_kwargs["distance"] = peak_cfg.distance
    chrom.find_peaks(**find_kwargs)

    if not (hasattr(chrom, "peaks") and chrom.peaks is not None):
        return 0
    n_peaks = len(chrom.peaks)

    if n_peaks > 0 and decon_cfg.enabled:
        _deconvolve(chrom, decon_cfg)
    return n_peaks


def _deconvolve(chrom, decon_cfg: DeconvolutionConfig) -> None:
    """Try MOCCA2's automatic deconvolution, then each configured model.

    Deconvolution failing is non-fatal: detected peaks still get integrated
    manually downstream.
    """
    try:
        chrom.deconvolve_peaks()
        return
    except Exception:
        pass

    for model in (decon_cfg.models or []):
        try:
            chrom.deconvolve_peaks(
                model=model,
                min_r2=decon_cfg.min_r2,
                relaxe_concs=False,
                max_comps=decon_cfg.max_comps,
            )
            return
        except Exception:
            continue
    print(" (deconvolution failed)", end="")


# ---------------------------------------------------------------------------
# Peak-table extraction
# ---------------------------------------------------------------------------

def build_peak_table(dataset: Dataset) -> pd.DataFrame:
    """Build the per-peak table (one row per detected peak/component).

    Uses MOCCA2's deconvolved component integrals when available, otherwise
    falls back to manual trapezoidal integration over the baseline-corrected
    signal. The ``Integration_Method`` column records which path was taken.
    """
    all_peaks: list[dict] = []
    n_with_components = 0
    n_without_components = 0

    for chrom in dataset.chromatograms:
        components = _safe_components(chrom)
        time_data = _get_time(chrom)
        experiment_name = getattr(chrom, "experiment_name", "Main")

        if components:
            n_with_components += 1
            all_peaks.extend(
                _rows_from_components(chrom, components, time_data, experiment_name, dataset.wavelength)
            )
        elif getattr(chrom, "peaks", None):
            n_without_components += 1
            all_peaks.extend(
                _rows_from_peaks(chrom, time_data, experiment_name, dataset.wavelength)
            )

    print("\nIntegration methods:")
    print(f"  MOCCA2 components (deconvolved): {n_with_components} chromatograms")
    print(f"  Manual trapezoid              : {n_without_components} chromatograms")

    if not all_peaks:
        print("\n  Warning: no peaks detected in any chromatogram")
        return pd.DataFrame()

    df = pd.DataFrame(all_peaks)
    df = _add_area_percent(df)
    print(f"  Total peaks detected          : {len(df)}")
    return df


def _safe_components(chrom):
    if hasattr(chrom, "all_components") and callable(chrom.all_components):
        try:
            comps = list(chrom.all_components())
            return comps if comps else None
        except Exception:
            return None
    return None


def _get_time(chrom):
    for attr in ["time", "times", "rt", "retention_time", "x"]:
        if hasattr(chrom, attr):
            return np.asarray(getattr(chrom, attr)).flatten()
    return None


def _get_signal(chrom):
    for attr in ["intensity", "signal", "y", "absorbance", "data", "raw_data"]:
        if hasattr(chrom, attr):
            return np.asarray(getattr(chrom, attr)).flatten()
    return None


def _rows_from_components(chrom, components, time_data, experiment_name, wavelength):
    rows = []
    for idx, component in enumerate(components, 1):
        rt = None
        if hasattr(component, "elution_time") and time_data is not None:
            try:
                elution_idx = int(component.elution_time)
                if 0 <= elution_idx < len(time_data):
                    rt = time_data[elution_idx]
            except Exception:
                pass

        area = getattr(component, "integral", None)

        height = None
        conc = getattr(component, "concentration", None)
        if conc is not None:
            height = np.max(conc) if hasattr(conc, "__len__") else conc

        width = None
        if area is not None and height is not None and height > 0:
            width = area / (0.6 * height)  # rough FWHM estimate

        rows.append(_peak_row(
            experiment_name, chrom, idx, rt, area, height, width, wavelength, "MOCCA2_Component",
            width_round=4,
        ))
    return rows


def _rows_from_peaks(chrom, time_data, experiment_name, wavelength):
    rows = []
    signal_raw = _get_signal(chrom)
    baseline = np.asarray(chrom.baseline).flatten() if getattr(chrom, "baseline", None) is not None else None

    for idx, peak in enumerate(chrom.peaks, 1):
        rt = None
        if hasattr(peak, "maximum") and time_data is not None:
            try:
                max_idx = int(peak.maximum)
                if 0 <= max_idx < len(time_data):
                    rt = time_data[max_idx]
            except Exception:
                pass

        height = getattr(peak, "height", None)

        area = None
        width = None
        if hasattr(peak, "left") and hasattr(peak, "right") and time_data is not None:
            try:
                left_idx, right_idx = int(peak.left), int(peak.right)
                width = time_data[right_idx] - time_data[left_idx]
                if signal_raw is not None:
                    peak_times = time_data[left_idx:right_idx + 1]
                    peak_signal = signal_raw[left_idx:right_idx + 1]
                    if baseline is not None:
                        peak_baseline = baseline[left_idx:right_idx + 1]
                    else:
                        peak_baseline = np.linspace(
                            signal_raw[left_idx], signal_raw[right_idx], len(peak_times)
                        )
                    area = np.trapezoid(peak_signal - peak_baseline, peak_times)
            except Exception:
                pass

        rows.append(_peak_row(
            experiment_name, chrom, idx, rt, area, height, width, wavelength, "Manual_Trapezoid",
            width_round=3,
        ))
    return rows


def _peak_row(experiment_name, chrom, idx, rt, area, height, width, wavelength, method, width_round):
    return {
        "Experiment": experiment_name,
        "Folder": chrom.folder_name,
        "Sample": chrom.sample_name,
        "Peak_Number": idx,
        "Retention_Time_min": round(float(rt), 3) if rt is not None else None,
        "Area_Integral": round(float(area), 2) if area is not None else None,
        "Height_mAU": round(float(height), 2) if height is not None else None,
        "Width_min": round(float(width), width_round) if width is not None else None,
        "Wavelength_nm": wavelength,
        "Integration_Method": method,
    }


def _add_area_percent(df: pd.DataFrame) -> pd.DataFrame:
    """Area as a percentage of the injection's total integrated area."""
    totals = df.groupby(["Sample", "Folder"])["Area_Integral"].transform("sum")
    with np.errstate(invalid="ignore", divide="ignore"):
        pct = np.where(
            (totals > 0) & df["Area_Integral"].notna(),
            df["Area_Integral"] / totals * 100,
            np.nan,
        )
    df["Area_Percent"] = np.round(pct, 2)
    return df


def build_sample_summary(dataset: Dataset) -> pd.DataFrame:
    """One row per injection: peak count and blank flag."""
    rows = []
    for chrom in dataset.chromatograms:
        n_peaks = len(chrom.peaks) if getattr(chrom, "peaks", None) is not None else 0
        rows.append({
            "Experiment": getattr(chrom, "experiment_name", "Main"),
            "Folder": chrom.folder_name,
            "Sample": chrom.sample_name,
            "N_Peaks": n_peaks,
            "Is_Blank": "Yes" if (dataset.blank is not None and chrom is dataset.blank) else "No",
        })
    return pd.DataFrame(rows)
