"""
Lightweight unit tests that don't require MOCCA2 or real .D data.

Covers config loading/validation and the calibration behaviour that differs
from the original scripts (the zero-area fix and negative clamping).

Run with:  pytest   (or)   python -m pytest tests/
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from hplc_gc_pipeline.config import load_config, CalibrationConfig
from hplc_gc_pipeline.calibration import area_to_concentration


def _write_config(dir_path: Path, body: str) -> None:
    (dir_path / "hplc_config.yaml").write_text(textwrap.dedent(body))


def test_defaults_apply_for_empty_config(tmp_path):
    _write_config(tmp_path, "")
    cfg = load_config(tmp_path)
    assert cfg.processing.wavelength_nm == 262
    assert cfg.analysis.calibration.zero_area_zero_conc is True
    assert cfg.plots.n_cols == 2
    # Paths resolve under the experiment folder.
    assert cfg.results_dir == tmp_path / "results"
    assert cfg.peak_results_csv == tmp_path / "results" / "peak_results.csv"


def test_overrides_are_read(tmp_path):
    _write_config(tmp_path, """
        processing:
          wavelength_nm: 210
          blank_folder_name: "091-0202.D"
        plots:
          n_cols: 3
    """)
    cfg = load_config(tmp_path)
    assert cfg.processing.wavelength_nm == 210
    assert cfg.processing.blank_folder_name == "091-0202.D"
    assert cfg.plots.n_cols == 3


def test_unknown_key_raises(tmp_path):
    _write_config(tmp_path, "processing:\n  wavelenght_nm: 262\n")  # typo
    with pytest.raises(ValueError, match="Unknown config key"):
        load_config(tmp_path)


def test_invalid_value_raises(tmp_path):
    _write_config(tmp_path, "plots:\n  n_cols: 0\n")
    with pytest.raises(ValueError, match="n_cols"):
        load_config(tmp_path)


def test_missing_config_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path)


# --- Calibration behaviour --------------------------------------------------

CAL = {"FA": (0.18, 5.0)}  # conc = 0.18*area + 5.0


def test_zero_area_gives_zero_when_enabled():
    cfg = CalibrationConfig(zero_area_zero_conc=True)
    assert area_to_concentration(0.0, "FA", CAL, cfg) == 0.0


def test_zero_area_gives_intercept_when_disabled():
    cfg = CalibrationConfig(zero_area_zero_conc=False)
    # original behaviour: 0.18*0 + 5.0 = 5.0
    assert area_to_concentration(0.0, "FA", CAL, cfg) == 5.0


def test_negative_prediction_clamped():
    cfg = CalibrationConfig(zero_area_zero_conc=False, clamp_negative_to_zero=True)
    # small positive area whose fit dips below zero -> clamped
    assert area_to_concentration(1.0, "FA", {"FA": (0.18, -100.0)}, cfg) == 0.0


def test_uncalibrated_compound_is_nan():
    import math
    cfg = CalibrationConfig()
    assert math.isnan(area_to_concentration(100.0, "Unknown", CAL, cfg))


def test_normal_conversion():
    cfg = CalibrationConfig()
    assert area_to_concentration(100.0, "FA", CAL, cfg) == pytest.approx(0.18 * 100 + 5.0)
