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
class Config:
    processing: ProcessingConfig = field(default_factory=ProcessingConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    plots: PlotsConfig = field(default_factory=PlotsConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    # Absolute path to the experiment directory this config was loaded for.
    # Populated by load_config(); not read from YAML.
    experiment_dir: Optional[Path] = None

    # --- Convenience path accessors (all under experiment_dir) -------------
    @property
    def data_root(self) -> Path:
        """Folder scanned for raw .D injections (experiment_dir/<data_dir>)."""
        if self.processing.data_dir:
            return self.experiment_dir / self.processing.data_dir
        return self.experiment_dir

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
        if name not in data:
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

    # ProcessingConfig / AnalysisConfig need their own nested handling because
    # they contain further dataclasses; _from_dict recurses one level using the
    # annotated field types, so build the top level explicitly.
    cfg = Config(
        processing=_build_processing(raw.get("processing")),
        analysis=_build_analysis(raw.get("analysis")),
        plots=_from_dict(PlotsConfig, raw.get("plots")),
        output=_from_dict(OutputConfig, raw.get("output")),
    )
    cfg.experiment_dir = experiment_dir

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
