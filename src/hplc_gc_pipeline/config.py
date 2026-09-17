"""
Configuration loading and validation.

The pipeline is driven by a single YAML file (the "input interface"). This
module loads that file, fills in defaults for any omitted field, and exposes
the result as nested dataclasses so the rest of the code reads
``cfg.processing.wavelength_nm`` instead of digging through raw dicts.

Resolution order for the config file (see :func:`resolve_config_path`):
1. an explicit path passed on the CLI (``--config``)
2. ``<experiment_dir>/hplc_config.yaml``

Nothing here hardcodes an experiment path. Every path in the returned config
is resolved relative to the experiment directory that was passed in.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any, Optional

import yaml

from .gc_processing import GCProcessingConfig, GCPeakDetectionConfig, GCBaselineConfig
from .gc_calibration import GCCalibrationConfig

# Name of the config file looked up inside an experiment folder when no
# explicit --config path is given.
DEFAULT_CONFIG_FILENAME = "hplc_config.yaml"


# ---------------------------------------------------------------------------
# Dataclasses mirroring the YAML schema. Defaults here ARE the documented
# defaults; config.example.yaml keeps them in sync for humans.
# ---------------------------------------------------------------------------

@dataclass
class PeakDetectionConfig:
    # All None -> use MOCCA2's find_peaks() defaults. Set a number to override.
    min_height: Optional[float] = None
    min_prominence: Optional[float] = None
    min_width: Optional[float] = None
    distance: Optional[float] = None


@dataclass
class DeconvolutionConfig:
    enabled: bool = True
    # Deconvolution models tried in order; MOCCA2's automatic pass is tried
    # first regardless, these are the explicit fallbacks.
    models: list[str] = field(default_factory=lambda: ["FraserSuzuki", "Gaussian"])
    min_r2: float = 0.95
    max_comps: int = 5


@dataclass
class ProcessingConfig:
    wavelength_nm: int = 262
    # .D folder name of the blank injection (with .D suffix); null = no blank.
    blank_folder_name: Optional[str] = None
    # Subfolder (relative to the experiment dir) that holds the raw .D data.
    # null = the experiment dir itself. Use this when raw data lives in a
    # nested folder (e.g. "Data") while config/CSVs/outputs stay at the root.
    data_dir: Optional[str] = None
    peak_detection: PeakDetectionConfig = field(default_factory=PeakDetectionConfig)
    deconvolution: DeconvolutionConfig = field(default_factory=DeconvolutionConfig)


@dataclass
class CalibrationConfig:
    # Clamp negative predicted concentrations to 0 (a linear fit can dip below
    # zero near the origin).
    clamp_negative_to_zero: bool = True
    # When True, a compound with zero integrated area is reported as 0 uM
    # rather than the calibration intercept b. Fixes the "pCA ~15 uM at zero
    # area" artefact. See README "Known fixes".
    zero_area_zero_conc: bool = True
    # Calibration source: "csv" (read standard.csv), "injections" (auto-detect
    # standard .D folders in peak_results.csv by name pattern), or "both"
    # (pool both sources for a single regression).
    source: str = "csv"
    # Regex (with named groups conc, unit, compound) identifying standard
    # injections by their sample name — e.g. "250uM_pCA_hexane" → 250 µM pCA.
    # Used when source is "injections" or "both".
    standard_pattern: str = (
        r"^(?P<conc>\d+(?:\.\d+)?)\s*(?P<unit>[a-zA-Zµ]*M)_(?P<compound>[^_]+)"
    )
    # Multiply *sample* concentrations by this factor to recover the undiluted
    # value. Standards are never scaled.  Example: 1:10 dilution → 10.
    dilution_factor: float = 1.0
    # Fit conc = slope*area (no intercept). Removes the "non-zero concentration
    # at zero area" artefact from a negative intercept.
    force_through_origin: bool = False


@dataclass
class SampleNameConfig:
    # Sample names are encoded as "<time_h>, <strain>, <replicate>".
    # Injections that don't match (blanks, standards) are skipped with a warning.
    pattern: str = (
        r"^\s*(?P<time>\d+)\s*,\s*(?P<strain>[^,]+?)\s*,\s*(?P<replicate>[A-Za-z])\s*$"
    )


@dataclass
class AnalysisConfig:
    # Reference tables, resolved relative to the experiment folder.
    compounds_file: str = "compounds.csv"
    standard_file: str = "standard.csv"
    # Drop peaks that fall outside every compound RT window.
    exclude_unknown: bool = True
    # CV% above which a warning is printed for a replicate group.
    cv_warning_threshold: float = 15.0
    calibration: CalibrationConfig = field(default_factory=CalibrationConfig)
    sample_name: SampleNameConfig = field(default_factory=SampleNameConfig)


@dataclass
class PlotsConfig:
    # Explicit strain plotting order (controls colour assignment).
    # null -> alphabetical.
    strain_order: Optional[list[str]] = None
    # Subset of strains to include in plots. null -> all strains.
    plot_strains: Optional[list[str]] = None
    # Restrict bar charts to these time points (h). null -> all.
    bar_time_points: Optional[list[int]] = None
    # Columns in the plot grid layout.
    n_cols: int = 2
    # Ordered colour palette (hex); strains/replicates cycle through it.
    palette: list[str] = field(
        default_factory=lambda: [
            "#377eb8", "#e41a1c", "#984ea3", "#4daf4a",
            "#ff7f00", "#000000", "#a65628", "#f781bf",
        ]
    )


@dataclass
class OutputConfig:
    # Subfolder names created inside the experiment directory.
    results_dirname: str = "results"
    analysis_dirname: str = "analysis"


@dataclass
class GCSampleNameConfig:
    # GC sample names use either underscore ("50_s11_A") or comma ("50, s11, A")
    # separators. The default pattern accepts both. Standards ("250uM_34DMS_hexane")
    # don't match and are skipped automatically.
    pattern: str = (
        r"^\s*(?P<time>\d+)\s*[_,]\s*(?P<strain>[^_,]+?)\s*[_,]\s*(?P<replicate>[A-Za-z])\s*$"
    )


@dataclass
class GCConfig:
    """Settings for the GC-MS modality (present only when the config has a
    ``gc:`` block). Mirrors the HPLC blocks but with GC-specific processing,
    calibration, and an underscore sample-name convention."""
    # Subfolder holding the raw GC .D data (relative to the experiment dir, or
    # absolute). null = the experiment dir itself.
    data_dir: Optional[str] = None
    # GC compound table (adds quantifier_mz / qualifier_mz columns).
    compounds_file: str = "gc_compounds.csv"
    processing: GCProcessingConfig = field(default_factory=GCProcessingConfig)
    calibration: GCCalibrationConfig = field(default_factory=GCCalibrationConfig)
    sample_name: GCSampleNameConfig = field(default_factory=GCSampleNameConfig)
    cv_warning_threshold: float = 15.0
    results_dirname: str = "gc_results"
    analysis_dirname: str = "gc_analysis"


@dataclass
class Config:
    processing: ProcessingConfig = field(default_factory=ProcessingConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    plots: PlotsConfig = field(default_factory=PlotsConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    # GC modality settings; None when the config has no ``gc:`` block.
    gc: Optional[GCConfig] = None
    # Active modalities for this experiment, derived from which blocks the YAML
    # contains. Set by load_config(); not read directly from YAML.
    modalities: list[str] = field(default_factory=lambda: ["hplc"])

    # Absolute path to the experiment directory this config was loaded for.
    # Populated by load_config(); not read from YAML.
    experiment_dir: Optional[Path] = None

    # --- Convenience path accessors (all under experiment_dir) -------------
    @property
    def data_roots(self) -> list[Path]:
        """All folders scanned for raw HPLC .D injections.

        ``data_dir`` may be a single subfolder name (str), a list of subfolder
        names, or None (meaning the experiment root itself).
        """
        dd = self.processing.data_dir
        if dd is None:
            return [self.experiment_dir]
        if isinstance(dd, list):
            result = []
            for d in dd:
                p = Path(str(d)).expanduser()
                result.append(p if p.is_absolute() else self.experiment_dir / str(d))
            return result or [self.experiment_dir]
        p = Path(str(dd)).expanduser()
        return [p if p.is_absolute() else self.experiment_dir / str(dd)]

    @property
    def data_root(self) -> Path:
        """First data root (single-dir backward-compat accessor)."""
        return self.data_roots[0]

    @property
    def results_dir(self) -> Path:
        return self.experiment_dir / self.output.results_dirname

    @property
    def analysis_dir(self) -> Path:
        return self.experiment_dir / self.output.analysis_dirname

    @property
    def peak_results_csv(self) -> Path:
        return self.results_dir / "peak_results.csv"

    @property
    def compounds_path(self) -> Path:
        return self.experiment_dir / self.analysis.compounds_file

    @property
    def standard_path(self) -> Path:
        return self.experiment_dir / self.analysis.standard_file

    # --- GC path accessors (only valid when self.gc is set) ----------------
    @property
    def gc_data_root(self) -> Path:
        """Folder scanned for raw GC .D injections."""
        if self.gc and self.gc.data_dir:
            p = Path(self.gc.data_dir).expanduser()
            return p if p.is_absolute() else self.experiment_dir / self.gc.data_dir
        return self.experiment_dir

    @property
    def gc_results_dir(self) -> Path:
        return self.experiment_dir / self.gc.results_dirname

    @property
    def gc_analysis_dir(self) -> Path:
        return self.experiment_dir / self.gc.analysis_dirname

    @property
    def gc_peak_results_csv(self) -> Path:
        return self.gc_results_dir / "peak_results.csv"

    @property
    def gc_compounds_path(self) -> Path:
        return self.experiment_dir / self.gc.compounds_file


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def _from_dict(cls, data: Any):
    """Recursively build a (possibly nested) dataclass from a plain dict.

    Unknown keys raise, so typos in the YAML are caught early rather than
    silently ignored.
    """
    if data is None:
        return cls()
    if not isinstance(data, dict):
        raise TypeError(f"Expected mapping for {cls.__name__}, got {type(data).__name__}")

    kwargs = {}
    field_map = {f.name: f for f in fields(cls)}
    unknown = set(data) - set(field_map)
    if unknown:
        raise ValueError(
            f"Unknown config key(s) for '{cls.__name__}': {sorted(unknown)}. "
            f"Allowed: {sorted(field_map)}"
        )

    for name, f in field_map.items():
        if name not in data or data[name] is None or data[name] == "":
            continue
        value = data[name]
        # Nested dataclass? recurse.
        if is_dataclass(f.type) or (isinstance(f.type, type) and is_dataclass(f.type)):
            kwargs[name] = _from_dict(f.type, value)
        else:
            kwargs[name] = value
    return cls(**kwargs)


def resolve_config_path(experiment_dir: Path, explicit: Optional[Path] = None) -> Path:
    """Determine which YAML file to load."""
    if explicit is not None:
        p = Path(explicit).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"--config file not found: {p}")
        return p
    candidate = experiment_dir / DEFAULT_CONFIG_FILENAME
    if not candidate.exists():
        raise FileNotFoundError(
            f"No config file found. Expected '{DEFAULT_CONFIG_FILENAME}' in "
            f"{experiment_dir}, or pass --config /path/to/config.yaml.\n"
            f"Copy config/config.example.yaml as a starting point."
        )
    return candidate


def load_config(experiment_dir: Path, explicit_config: Optional[Path] = None) -> Config:
    """Load and validate the config for a given experiment directory."""
    experiment_dir = Path(experiment_dir).expanduser().resolve()
    if not experiment_dir.is_dir():
        raise NotADirectoryError(f"Experiment directory not found: {experiment_dir}")

    config_path = resolve_config_path(experiment_dir, explicit_config)
    with open(config_path, "r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}

    # An experiment may include HPLC data, GC data, or both. Which modalities
    # are active is inferred from which blocks the YAML contains:
    #   * HPLC: a top-level `processing`/`analysis` (legacy flat form) or an
    #     explicit `hplc:` block.
    #   * GC:   a `gc:` block.
    # A config with none of these defaults to HPLC (backward compatible).
    hplc_block = raw.get("hplc") if isinstance(raw.get("hplc"), dict) else raw
    modalities = []
    if any(k in raw for k in ("hplc", "processing")):
        modalities.append("hplc")
    if "gc" in raw:
        modalities.append("gc")
    if not modalities:
        modalities = ["hplc"]

    # ProcessingConfig / AnalysisConfig need their own nested handling because
    # they contain further dataclasses; build them explicitly.
    cfg = Config(
        processing=_build_processing(hplc_block.get("processing")),
        analysis=_build_analysis(hplc_block.get("analysis")),
        plots=_from_dict(PlotsConfig, raw.get("plots")),
        output=_from_dict(OutputConfig, raw.get("output")),
        gc=_build_gc(raw.get("gc")) if "gc" in raw else None,
    )
    cfg.modalities = modalities
    cfg.experiment_dir = experiment_dir

    # Old config templates (pre-fix) wrote dilution_factor: 0.0 as a placeholder.
    # Treat 0.0 as "not set" and fall back to the default of 1.0 (no dilution).
    if cfg.analysis.calibration.dilution_factor <= 0:
        cfg.analysis.calibration.dilution_factor = 1.0
    if cfg.gc is not None and cfg.gc.calibration.dilution_factor <= 0:
        cfg.gc.calibration.dilution_factor = 1.0

    _validate(cfg)
    return cfg


def _build_processing(data: Any) -> ProcessingConfig:
    data = data or {}
    _check_keys(data, ProcessingConfig, "processing")
    return ProcessingConfig(
        wavelength_nm=data.get("wavelength_nm", ProcessingConfig.wavelength_nm),
        blank_folder_name=data.get("blank_folder_name"),
        data_dir=data.get("data_dir"),
        peak_detection=_from_dict(PeakDetectionConfig, data.get("peak_detection")),
        deconvolution=_from_dict(DeconvolutionConfig, data.get("deconvolution")),
    )


def _build_analysis(data: Any) -> AnalysisConfig:
    data = data or {}
    _check_keys(data, AnalysisConfig, "analysis")
    defaults = AnalysisConfig()
    return AnalysisConfig(
        compounds_file=data.get("compounds_file", defaults.compounds_file),
        standard_file=data.get("standard_file", defaults.standard_file),
        exclude_unknown=data.get("exclude_unknown", defaults.exclude_unknown),
        cv_warning_threshold=data.get("cv_warning_threshold", defaults.cv_warning_threshold),
        calibration=_from_dict(CalibrationConfig, data.get("calibration")),
        sample_name=_from_dict(SampleNameConfig, data.get("sample_name")),
    )


def _build_gc(data: Any) -> GCConfig:
    data = data or {}
    _check_keys(data, GCConfig, "gc")
    defaults = GCConfig()
    return GCConfig(
        data_dir=data.get("data_dir"),
        compounds_file=data.get("compounds_file", defaults.compounds_file),
        processing=_build_gc_processing(data.get("processing")),
        calibration=_from_dict(GCCalibrationConfig, data.get("calibration")),
        sample_name=_from_dict(GCSampleNameConfig, data.get("sample_name")),
        cv_warning_threshold=data.get("cv_warning_threshold", defaults.cv_warning_threshold),
        results_dirname=data.get("results_dirname", defaults.results_dirname),
        analysis_dirname=data.get("analysis_dirname", defaults.analysis_dirname),
    )


def _build_gc_processing(data: Any) -> GCProcessingConfig:
    data = data or {}
    _check_keys(data, GCProcessingConfig, "gc.processing")
    defaults = GCProcessingConfig()
    _qc = data.get("quant_channel")
    return GCProcessingConfig(
        quant_channel=_qc if _qc is not None else defaults.quant_channel,
        peak_detection=_from_dict(GCPeakDetectionConfig, data.get("peak_detection")),
        baseline=_from_dict(GCBaselineConfig, data.get("baseline")),
    )


def _check_keys(data: dict, cls, name: str) -> None:
    allowed = {f.name for f in fields(cls)}
    unknown = set(data) - allowed
    if unknown:
        raise ValueError(
            f"Unknown config key(s) for '{name}': {sorted(unknown)}. "
            f"Allowed: {sorted(allowed)}"
        )


def _validate(cfg: Config) -> None:
    """Cheap sanity checks that fail fast with a clear message."""
    if cfg.processing.wavelength_nm <= 0:
        raise ValueError("processing.wavelength_nm must be positive")
    if cfg.plots.n_cols < 1:
        raise ValueError("plots.n_cols must be >= 1")
    if cfg.analysis.cv_warning_threshold < 0:
        raise ValueError("analysis.cv_warning_threshold must be >= 0")
    if cfg.analysis.calibration.source not in ("csv", "injections", "both"):
        raise ValueError("analysis.calibration.source must be 'csv', 'injections', or 'both'")
    if cfg.analysis.calibration.dilution_factor <= 0:
        raise ValueError("analysis.calibration.dilution_factor must be positive")
    if cfg.gc is not None:
        if cfg.gc.processing.quant_channel not in ("eic", "tic", "fid"):
            raise ValueError("gc.processing.quant_channel must be 'eic', 'tic', or 'fid'")
        if cfg.gc.calibration.source not in ("injections", "csv"):
            raise ValueError("gc.calibration.source must be 'injections' or 'csv'")
        if cfg.gc.calibration.dilution_factor <= 0:
            raise ValueError("gc.calibration.dilution_factor must be positive")
