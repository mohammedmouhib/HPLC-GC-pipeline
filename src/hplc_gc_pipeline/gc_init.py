"""
hplc init: auto-discover an experiment folder and write starter config files.

What is auto-detected:
  - Whether HPLC and/or GC-MS data are present, and in which subfolders
  - Monitored m/z values (SIM header of any GC injection)
  - Which injections are GC standards (names match standard_pattern)
  - Peak RT windows (from averaged EIC across standard injections)
  - Quantifier m/z per peak (highest-signal ion at that RT)
  - Qualifier ion ratios (area ratios averaged across standards)
  - GC sample name convention (comma-separated vs underscore)

What must be filled in manually after init (flagged in the written files):
  - GC compound names (written as Peak_1, Peak_2, ...)
  - calibrate_as (cross-compound calibration; requires chemistry knowledge)
  - dilution_factor (lab protocol, not encoded in data files)
  - HPLC compound RT windows and calibration data
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from .gc_agilent import GCInjection, find_gc_injections, load_gc_traces, GCTraces
from .gc_calibration import parse_standard_name
from .gc_processing import als_baseline, GCBaselineConfig, GCPeakDetectionConfig, integrate_peak

# ---------------------------------------------------------------------------
# HPLC folder detection
# ---------------------------------------------------------------------------

def _is_hplc_d_folder(path: Path) -> bool:
    """A .D folder without ``data.ms`` is an HPLC (DAD) injection, not GC-MS."""
    return (path.is_dir() and path.name.endswith(".D")
            and not (path / "data.ms").exists())


def detect_hplc_data_dir(experiment_dir: Path) -> Optional[str]:
    """Return the overarching HPLC data subfolder, or None when at the experiment root.

    Walks the experiment tree recursively to find all directories that directly
    contain HPLC ``.D`` folders, then returns their lowest common ancestor
    relative to ``experiment_dir``.  Returns None when the LCA is the experiment
    root itself (i.e. ``.D`` files live directly there).
    """
    experiment_dir = Path(experiment_dir)

    hplc_parents: list[Path] = []

    def walk(directory: Path) -> None:
        try:
            entries = sorted(directory.iterdir())
        except PermissionError:
            return
        has_hplc = False
        for item in entries:
            if not item.is_dir():
                continue
            if _is_hplc_d_folder(item):
                has_hplc = True
            elif not item.name.endswith(".D") and not item.name.startswith("."):
                walk(item)
        if has_hplc:
            hplc_parents.append(directory)

    walk(experiment_dir)

    if not hplc_parents:
        return None
    common = Path(os.path.commonpath([str(p) for p in hplc_parents]))
    if common == experiment_dir:
        return None
    return str(common.relative_to(experiment_dir))


# ---------------------------------------------------------------------------
# Sample-name convention detection
# ---------------------------------------------------------------------------

_COMMA_PATTERN = (
    r"^\s*(?P<time>\d+)\s*,\s*(?P<strain>[^,]+?)\s*,\s*(?P<replicate>[A-Za-z][A-Za-z0-9]*)\s*$"
)
_UNDERSCORE_PATTERN = (
    r"^\s*(?P<time>\d+)\s*_\s*(?P<strain>[^_]+?)\s*_\s*(?P<replicate>[A-Za-z][A-Za-z0-9]*)\s*$"
)
_DEFAULT_STD_PATTERN = (
    r"^(?P<conc>\d+(?:\.\d+)?)\s*(?P<unit>[a-zA-Zµ]*M)_(?P<compound>[^_]+)"
)


def _count_matches(names: list[str], pattern: str) -> int:
    rx = re.compile(pattern)
    return sum(1 for n in names if rx.match(n))


def detect_sample_name_convention(sample_names: list[str]) -> tuple[str, str]:
    """Return (convention_label, regex_pattern) for the better-matching convention."""
    n_comma = _count_matches(sample_names, _COMMA_PATTERN)
    n_under = _count_matches(sample_names, _UNDERSCORE_PATTERN)
    if n_comma >= n_under:
        return "comma", _COMMA_PATTERN
    return "underscore", _UNDERSCORE_PATTERN


# ---------------------------------------------------------------------------
# Peak discovery from standard injections
# ---------------------------------------------------------------------------

@dataclass
class _RawPeak:
    rt_apex: float
    rt_low: float          # integration left edge (before margin)
    rt_high: float         # integration right edge (before margin)
    mz: float
    area: float
    apex_height: float


@dataclass
class DiscoveredPeak:
    """One detected chromatographic peak, ready to become a gc_compounds.csv row."""
    index: int             # 1-based counter for auto-naming
    rt_low: float
    rt_high: float
    rt_apex: float
    quantifier_mz: float
    qualifiers: list[tuple[float, float]]   # (mz, ratio)
    n_standards_detected: int
    total_standards: int


def _find_peaks_in_eic(
    time_min: np.ndarray,
    eic: np.ndarray,
    bl_cfg: GCBaselineConfig,
    pk_cfg: GCPeakDetectionConfig,
    rt_scan_min: float,
    rt_scan_max: float,
) -> list[_RawPeak]:
    """Find all peaks in a single EIC trace that clear the SNR floor."""
    from scipy.signal import find_peaks as sp_find_peaks

    bl = als_baseline(eic, bl_cfg.lam, bl_cfg.p, bl_cfg.niter)
    corrected = eic - bl

    # Restrict scan to sensible RT range (skip solvent front + column dead time)
    mask = (time_min >= rt_scan_min) & (time_min <= rt_scan_max)
    if not mask.any():
        return []

    sig = corrected.copy()
    sig[~mask] = 0.0

    noise = float(1.4826 * np.median(np.abs(corrected - np.median(corrected)))) or float(np.std(corrected)) or 1.0
    min_height = pk_cfg.min_snr * noise

    peaks_idx, props = sp_find_peaks(sig, height=min_height, prominence=min_height * 0.5)

    results = []
    for idx in peaks_idx:
        apex_val = float(sig[idx])
        threshold = max(pk_cfg.edge_fraction * apex_val, noise)
        # Walk left
        left = idx
        while left > 0 and sig[left - 1] >= threshold and sig[left - 1] <= sig[left]:
            left -= 1
        # Walk right
        right = idx
        n = len(sig)
        while right < n - 1 and sig[right + 1] >= threshold and sig[right + 1] <= sig[right]:
            right += 1
        seg_t = time_min[left:right + 1]
        seg_y = np.clip(sig[left:right + 1], 0, None)
        area = float(np.trapezoid(seg_y, seg_t)) if right > left else apex_val * float(time_min[1] - time_min[0])
        results.append(_RawPeak(
            rt_apex=float(time_min[idx]),
            rt_low=float(time_min[left]),
            rt_high=float(time_min[right]),
            mz=0.0,   # filled by caller
            area=area,
            apex_height=apex_val,
        ))
    return results


def _cluster_peaks(all_peaks: list[_RawPeak], rt_tol: float = 0.15) -> list[list[_RawPeak]]:
    """Group _RawPeak objects from different m/z channels that share the same RT apex."""
    if not all_peaks:
        return []
    sorted_peaks = sorted(all_peaks, key=lambda p: p.rt_apex)
    clusters: list[list[_RawPeak]] = []
    current = [sorted_peaks[0]]
    for pk in sorted_peaks[1:]:
        if pk.rt_apex - current[-1].rt_apex <= rt_tol:
            current.append(pk)
        else:
            clusters.append(current)
            current = [pk]
    clusters.append(current)
    return clusters


def discover_peaks(
    standard_injections: list[GCInjection],
    standard_concentrations: Optional[list[float]] = None,
    rt_scan_min: float = 1.5,
    rt_scan_max: float = 30.0,
    rt_margin: float = 0.08,
    qualifier_threshold: float = 0.05,
    min_standards_fraction: float = 0.6,
    min_conc_correlation: float = 0.7,
    bl_cfg: Optional[GCBaselineConfig] = None,
    pk_cfg: Optional[GCPeakDetectionConfig] = None,
) -> list[DiscoveredPeak]:
    """
    Load standard injections and discover chromatographic peaks.

    For each standard, find peaks in every EIC channel. Cluster across
    standards/channels by RT, pick the highest-signal m/z as quantifier,
    compute qualifier ratios.

    Two filters are applied:

    1. Reproducibility: the quantifier m/z must appear in at least
       ``min_standards_fraction`` of standard injections.  Eliminates
       single-injection noise.

    2. Concentration scaling: when ``standard_concentrations`` is provided
       (list of uM concentrations aligned with ``standard_injections``),
       clusters whose area does NOT scale with concentration are discarded.
       The Pearson correlation between area and concentration must exceed
       ``min_conc_correlation``.  This is the most reliable way to separate
       real compound peaks from constant-area matrix/solvent peaks.
    """
    bl_cfg = bl_cfg or GCBaselineConfig()
    pk_cfg = pk_cfg or GCPeakDetectionConfig()

    n_standards = len(standard_injections)

    # ── Pass 1: collect all raw peaks, tagged by which injection they came from ──
    # raw_by_inj[i] = list of _RawPeak from standard_injections[i]
    raw_by_inj: list[list[_RawPeak]] = []
    all_raw: list[_RawPeak] = []

    for inj in standard_injections:
        inj_peaks: list[_RawPeak] = []
        try:
            traces = load_gc_traces(inj)
        except Exception as exc:
            print(f"  [init] Warning: could not load {inj.folder_name}: {exc}")
            raw_by_inj.append([])
            continue
        for mz_idx, mz in enumerate(traces.mz):
            eic = traces.eic[:, mz_idx]
            peaks = _find_peaks_in_eic(
                traces.time_min, eic, bl_cfg, pk_cfg, rt_scan_min, rt_scan_max
            )
            for pk in peaks:
                pk.mz = float(mz)
                inj_peaks.append(pk)
                all_raw.append(pk)
        raw_by_inj.append(inj_peaks)

    if not all_raw:
        return []

    clusters = _cluster_peaks(all_raw, rt_tol=0.15)

    # ── Pass 2: build DiscoveredPeak for each cluster, then filter ────────────
    discovered: list[DiscoveredPeak] = []
    peak_idx = 0
    for cluster in clusters:
        # Per-mz: sum area across the cluster
        mz_area: dict[float, float] = {}
        for pk in cluster:
            mz_area[pk.mz] = mz_area.get(pk.mz, 0.0) + pk.area

        quant_mz = max(mz_area, key=lambda m: mz_area[m])
        quant_area = mz_area[quant_mz]

        # Count how many individual standard injections have a peak at this RT
        # for the selected quantifier m/z.
        rt_apex = float(np.median([pk.rt_apex for pk in cluster]))
        n_detected = sum(
            1 for inj_peaks in raw_by_inj
            if any(abs(pk.rt_apex - rt_apex) < 0.2 and pk.mz == quant_mz
                   for pk in inj_peaks)
        )

        # Filter 1: reproducibility across standards.
        min_required = max(1, round(min_standards_fraction * n_standards)) if n_standards > 1 else 1
        if n_detected < min_required:
            continue

        # Filter 2: area must scale with concentration (real compound vs matrix).
        # Applied only when concentrations are known and at least 3 distinct
        # levels are detectable (fewer points make correlation unreliable).
        if standard_concentrations is not None and n_standards > 2:
            areas_for_corr: list[float] = []
            concs_for_corr: list[float] = []
            for inj_peaks, conc in zip(raw_by_inj, standard_concentrations):
                matching = [
                    pk for pk in inj_peaks
                    if pk.mz == quant_mz and abs(pk.rt_apex - rt_apex) < 0.2
                ]
                if matching:
                    areas_for_corr.append(max(pk.area for pk in matching))
                    concs_for_corr.append(conc)
            if len(areas_for_corr) >= 3:
                a = np.array(areas_for_corr)
                c = np.array(concs_for_corr)
                # Pearson r between concentration and area
                if a.std() > 0 and c.std() > 0:
                    corr = float(np.corrcoef(c, a)[0, 1])
                else:
                    corr = 0.0
                if corr < min_conc_correlation:
                    continue

        qualifiers = []
        for mz, area in sorted(mz_area.items()):
            if mz == quant_mz:
                continue
            ratio = area / quant_area
            if ratio >= qualifier_threshold:
                qualifiers.append((mz, round(ratio, 2)))

        rt_low = min(pk.rt_low for pk in cluster) - rt_margin
        rt_high = max(pk.rt_high for pk in cluster) + rt_margin

        peak_idx += 1
        discovered.append(DiscoveredPeak(
            index=peak_idx,
            rt_low=round(max(rt_low, rt_scan_min), 2),
            rt_high=round(min(rt_high, rt_scan_max), 2),
            rt_apex=round(rt_apex, 3),
            quantifier_mz=quant_mz,
            qualifiers=qualifiers,
            n_standards_detected=n_detected,
            total_standards=n_standards,
        ))

    return discovered


# ---------------------------------------------------------------------------
# Main init routine
# ---------------------------------------------------------------------------

@dataclass
class InitResult:
    experiment_dir: Path
    all_injections: list[GCInjection]
    standard_injections: list[GCInjection]
    sample_injections: list[GCInjection]
    sample_convention: str          # "comma" | "underscore"
    sample_pattern: str
    monitored_mz: list[float]
    peaks: list[DiscoveredPeak]
    hplc_data_dir: Optional[str] = None   # overarching HPLC subfolder (LCA)
    gc_data_dir: Optional[str] = None    # overarching GC subfolder (LCA)
    warnings: list[str] = field(default_factory=list)


def run_gc_init(
    experiment_dir: Path,
    std_pattern: str = _DEFAULT_STD_PATTERN,
    rt_scan_min: float = 1.5,
    rt_scan_max: float = 30.0,
    rt_margin: float = 0.08,
    include_gc: bool = True,
) -> InitResult:
    """Scan ``experiment_dir`` and build an :class:`InitResult`.

    When ``include_gc`` is False, skips GC injection discovery and peak
    detection (for HPLC-only experiments). ``hplc_data_dir`` is always
    detected regardless of ``include_gc``.
    """
    experiment_dir = Path(experiment_dir).expanduser().resolve()

    hplc_data_dir = detect_hplc_data_dir(experiment_dir)

    warnings: list[str] = []
    all_injections: list[GCInjection] = []
    standard_injections: list[GCInjection] = []
    sample_injections: list[GCInjection] = []
    convention = "comma"
    pattern = _COMMA_PATTERN
    monitored_mz: list[float] = []
    peaks: list[DiscoveredPeak] = []

    if not include_gc:
        return InitResult(
            experiment_dir=experiment_dir,
            all_injections=[],
            standard_injections=[],
            sample_injections=[],
            sample_convention=convention,
            sample_pattern=pattern,
            monitored_mz=[],
            peaks=[],
            hplc_data_dir=hplc_data_dir,
            gc_data_dir=None,
            warnings=warnings,
        )

    print(f"\nScanning {experiment_dir} for GC-MS .D folders...")
    all_injections = find_gc_injections(experiment_dir)
    if not all_injections:
        warnings.append(
            f"No GC-MS .D folders (containing data.ms) found in {experiment_dir}. "
            "gc_compounds.csv will be empty."
        )
        return InitResult(
            experiment_dir=experiment_dir,
            all_injections=[],
            standard_injections=[],
            sample_injections=[],
            sample_convention=convention,
            sample_pattern=pattern,
            monitored_mz=[],
            peaks=[],
            hplc_data_dir=hplc_data_dir,
            gc_data_dir=None,
            warnings=warnings,
        )
    print(f"  Found {len(all_injections)} injection(s)")

    # Partition into standards vs samples
    standard_injections = [
        inj for inj in all_injections
        if parse_standard_name(inj.sample_name, std_pattern) is not None
    ]
    sample_injections = [
        inj for inj in all_injections
        if inj not in standard_injections
    ]
    print(f"  Standards (matching standard_pattern): {len(standard_injections)}")
    print(f"  Samples (remaining):                   {len(sample_injections)}")

    if not standard_injections:
        warnings.append(
            "No injections matched the standard_pattern. Peaks will be detected from ALL "
            "injections. RT windows may be wider. calibration.source must be set to 'csv'."
        )

    # Detect sample name convention
    sample_names_list = [inj.sample_name for inj in sample_injections]
    convention, pattern = detect_sample_name_convention(sample_names_list)
    n_match = _count_matches(sample_names_list, pattern)
    if n_match < len(sample_injections):
        warnings.append(
            f"{len(sample_injections) - n_match} sample(s) did not match the detected "
            f"'{convention}' name pattern. Check sample_name.pattern in the config."
        )

    # Read monitored m/z from any injection
    probe_inj = (standard_injections or sample_injections)[0]
    try:
        probe_traces = load_gc_traces(probe_inj)
        monitored_mz = [float(m) for m in probe_traces.mz]
        print(f"  Monitored m/z: {[int(m) for m in monitored_mz]}")
    except Exception as exc:
        warnings.append(f"Could not read m/z list from {probe_inj.folder_name}: {exc}")

    # Parse standard concentrations (used for concentration-scaling filter)
    std_concentrations: Optional[list[float]] = None
    if standard_injections:
        parsed = [
            parse_standard_name(inj.sample_name, std_pattern)
            for inj in standard_injections
        ]
        if all(p is not None for p in parsed):
            std_concentrations = [p[1] for p in parsed]  # type: ignore[index]

    # Discover peaks from standards (fall back to all injections if no standards)
    discovery_injections = standard_injections if standard_injections else all_injections
    print(f"\nDetecting peaks from {len(discovery_injections)} standard injection(s)...")
    peaks = discover_peaks(
        discovery_injections,
        standard_concentrations=std_concentrations,
        rt_scan_min=rt_scan_min,
        rt_scan_max=rt_scan_max,
        rt_margin=rt_margin,
    )
    print(f"  Detected {len(peaks)} peak cluster(s)")

    if not peaks:
        warnings.append(
            "No peaks were detected above the SNR threshold. The compounds file will be "
            "empty. Try adjusting rt_scan_min/rt_scan_max or check that standards have "
            "detectable signal."
        )

    gc_parents = {inj.folder_path.parent for inj in all_injections}
    gc_common = Path(os.path.commonpath([str(p) for p in gc_parents]))
    gc_data_dir: Optional[str] = (
        None if gc_common == experiment_dir
        else str(gc_common.relative_to(experiment_dir))
    )

    return InitResult(
        experiment_dir=experiment_dir,
        all_injections=all_injections,
        standard_injections=standard_injections,
        sample_injections=sample_injections,
        sample_convention=convention,
        sample_pattern=pattern,
        monitored_mz=monitored_mz,
        peaks=peaks,
        hplc_data_dir=hplc_data_dir,
        gc_data_dir=gc_data_dir,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# File writers
# ---------------------------------------------------------------------------

def write_gc_compounds(result: InitResult, output_path: Path) -> None:
    """Write a draft gc_compounds.csv from discovered peaks."""
    lines = ["Compound,RT_low,RT_high,quantifier_mz,qualifier_mz,ion_ratio_tol,calibrate_as,Notes"]
    for pk in result.peaks:
        name = f"Peak_{pk.index}"
        qual_str = ";".join(f"{int(mz)}={ratio:.2f}" for mz, ratio in pk.qualifiers)
        lines.append(
            f"{name},{pk.rt_low},{pk.rt_high},{int(pk.quantifier_mz)},"
            f"{qual_str},0.30,,"
            f"auto-detected RT={pk.rt_apex:.3f} min quant=m/z{int(pk.quantifier_mz)}"
        )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_hplc_config(
    result: InitResult,
    output_path: Path,
    std_pattern: str,
    *,
    include_hplc: bool = True,
    include_gc: bool = True,
) -> None:
    """Write a complete hplc_config.yaml annotated with auto-detected values.

    When ``include_hplc`` is True a ``processing:`` / ``analysis:`` block
    (HPLC modality) is written with the auto-detected ``data_dir``.
    When ``include_gc`` is True a ``gc:`` block is written.
    """
    has_standards = bool(result.standard_injections)

    lines = [
        "# hplc_config.yaml  --  auto-generated by `hplc init`",
        "# Lines marked  # <-- REVIEW  should be confirmed before running.",
        "# Lines marked  # auto-detected  were inferred from the data.",
        "",
    ]

    if include_hplc:
        dd = result.hplc_data_dir
        if dd is None:
            data_dir_line = "  data_dir: null  # .D files are directly in the experiment folder"
        else:
            data_dir_line = f"  data_dir: {dd}  # auto-detected"
        lines += [
            "# -- HPLC modality -------------------------------------------------------",
            "processing:",
            "  wavelength_nm: 262  # <-- REVIEW: set to your DAD detection wavelength (nm)",
            data_dir_line,
            "  blank_folder_name: null  # <-- REVIEW: blank .D folder name, or null",
            "",
            "analysis:",
            "  compounds_file: compounds.csv    # edit compound names and RT windows",
            "  standard_file: standard.csv      # edit calibration points",
            "  exclude_unknown: true",
            "  calibration:",
            "    source: csv  # <-- REVIEW: 'csv' reads standard.csv; 'injections' uses name-encoded standards",
            "    dilution_factor: 1.0  # <-- REVIEW: sample prep dilution (e.g. 10.0 for 1:10 dilution)",
            "  sample_name:",
            "    # Default expects: \"48, s11, A\"  ->  time=48, strain=s11, replicate=A",
            r"    pattern: '^\s*(?P<time>\d+)\s*,\s*(?P<strain>[^,]+?)\s*,\s*(?P<replicate>[A-Za-z])\s*$'",
            "",
        ]

    if include_gc:
        convention_note = (
            "comma-separated  e.g. \"48, s11, A\""
            if result.sample_convention == "comma"
            else "underscore-separated  e.g. \"48_s11_A\""
        )
        n_sample = len(result.sample_injections)
        n_match = _count_matches(
            [i.sample_name for i in result.sample_injections], result.sample_pattern
        )
        gc_dir = result.gc_data_dir
        gc_data_dir_line = (
            f"  data_dir: {gc_dir}  # auto-detected"
            if gc_dir else
            "  data_dir: null  # .D files are directly in the experiment folder"
        )
        lines += [
            "# -- GC-MS modality ------------------------------------------------------",
            "gc:",
            gc_data_dir_line,
            "  compounds_file: gc_compounds.csv  # rename Peak_N to real compound names",
            "",
            "  processing:",
            "    quant_channel: eic  # auto-detected: SIM data present",
            "",
            "  calibration:",
            f"    source: {'injections' if has_standards else 'csv'}  # auto-detected",
        ]

        if has_standards:
            lines += [
                f"    standard_pattern: '{std_pattern}'",
                "    force_through_origin: true",
                "    clamp_negative_to_zero: true",
                "    require_ion_ratio_pass: false",
                "    dilution_factor: 1.0  # <-- REVIEW: sample prep dilution (e.g. 10.0 for 1:10)",
            ]
        else:
            lines += [
                "    standard_file: gc_standard.csv  # <-- REVIEW: provide a standard-curve CSV",
                "    force_through_origin: true",
                "    clamp_negative_to_zero: true",
                "    dilution_factor: 1.0  # <-- REVIEW",
            ]

        lines += [
            "",
            f"  # Sample names auto-detected as {convention_note}  # auto-detected",
            f"  # {n_match}/{n_sample} sample(s) matched.",
            "  sample_name:",
            f"    pattern: '{result.sample_pattern}'  # auto-detected",
            "",
            "  cv_warning_threshold: 15.0",
            "  results_dirname: gc_results",
            "  analysis_dirname: gc_analysis",
            "",
        ]

    lines += [
        "# -- Shared plotting settings ---------------------------------------------",
        "plots:",
        "  n_cols: 3",
    ]

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------

def print_summary(
    result: InitResult,
    config_path: Path,
    compounds_path: Optional[Path],
    extra_files: Optional[list[Path]] = None,
    include_hplc: bool = True,
    include_gc: bool = True,
) -> None:
    sep = "=" * 70
    print(f"\n{sep}")
    print("INIT SUMMARY")
    print(sep)
    if result.hplc_data_dir:
        print(f"  HPLC data dir    : {result.hplc_data_dir}/")
    if result.gc_data_dir:
        print(f"  GC data dir      : {result.gc_data_dir}/")
    if result.all_injections:
        print(f"  GC injections    : {len(result.all_injections)}")
        print(f"    Standards      : {len(result.standard_injections)}")
        print(f"    Samples        : {len(result.sample_injections)}")
        print(f"  Monitored m/z    : {[int(m) for m in result.monitored_mz]}")
        print(f"  GC sample naming : {result.sample_convention}-separated")
        print(f"  GC peaks found   : {len(result.peaks)}")
    if result.peaks:
        print()
        print(f"  {'Peak':<12} {'RT window':>14}  {'Quant m/z':>10}  {'Qualifiers'}")
        for pk in result.peaks:
            qual_str = "  ".join(
                f"m/z{int(mz)}={ratio:.2f}" for mz, ratio in pk.qualifiers
            ) or "(none)"
            print(f"  Peak_{pk.index:<7} {pk.rt_low:.2f}-{pk.rt_high:.2f} min  "
                  f"m/z {int(pk.quantifier_mz):>6}      {qual_str}")

    print(f"\n  Written:")
    print(f"    {config_path}")
    if compounds_path is not None:
        print(f"    {compounds_path}")
    for p in (extra_files or []):
        print(f"    {p}")

    if result.warnings:
        print(f"\n  Warnings:")
        for w in result.warnings:
            print(f"    ! {w}")

    print(f"\n  Next steps:")
    steps = []
    if result.peaks:
        steps.append("Rename Peak_N entries in gc_compounds.csv to real compound names")
        steps.append("Set calibrate_as for any compound that borrows another's standards")
    if include_hplc:
        steps.append("Edit compounds.csv with your HPLC compound RT windows")
    if include_hplc:
        steps.append("Set dilution_factor in hplc_config.yaml under analysis.calibration (current: 1.0)")
    if include_gc:
        steps.append("Set dilution_factor in hplc_config.yaml under gc.calibration (current: 1.0)")
    steps.append(f"Run:  hplc run {result.experiment_dir}")
    for i, step in enumerate(steps, 1):
        print(f"    {i}. {step}")
    print(sep)
