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


def test_data_dir_resolves_data_root(tmp_path):
    # When data_dir is set, raw data is discovered under <exp>/<data_dir> while
    # outputs/CSVs stay at the experiment root.
    _write_config(tmp_path, "processing:\n  data_dir: Data\n")
    cfg = load_config(tmp_path)
    assert cfg.processing.data_dir == "Data"
    assert cfg.data_root == tmp_path / "Data"


def test_data_root_defaults_to_experiment_dir(tmp_path):
    _write_config(tmp_path, "")
    cfg = load_config(tmp_path)
    assert cfg.processing.data_dir is None
    assert cfg.data_root == tmp_path


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


# --- GC: config / modality detection ---------------------------------------

def test_no_gc_block_is_hplc_only(tmp_path):
    _write_config(tmp_path, "processing:\n  wavelength_nm: 262\n")
    cfg = load_config(tmp_path)
    assert cfg.modalities == ["hplc"]
    assert cfg.gc is None


def test_gc_only_experiment(tmp_path):
    _write_config(tmp_path, """
        gc:
          data_dir: Data
          calibration:
            source: injections
            dilution_factor: 10.0
    """)
    cfg = load_config(tmp_path)
    assert cfg.modalities == ["gc"]
    assert cfg.gc.calibration.dilution_factor == 10.0
    assert cfg.gc.processing.quant_channel == "eic"
    assert cfg.gc_data_root == tmp_path / "Data"


def test_both_modalities(tmp_path):
    _write_config(tmp_path, """
        processing:
          wavelength_nm: 262
        gc:
          calibration:
            source: csv
    """)
    cfg = load_config(tmp_path)
    assert cfg.modalities == ["hplc", "gc"]


def test_gc_bad_quant_channel_raises(tmp_path):
    _write_config(tmp_path, "gc:\n  processing:\n    quant_channel: mass\n")
    with pytest.raises(ValueError, match="quant_channel"):
        load_config(tmp_path)


def test_gc_unknown_key_raises(tmp_path):
    _write_config(tmp_path, "gc:\n  dilution_factor: 10\n")  # belongs under calibration
    with pytest.raises(ValueError, match="Unknown config key"):
        load_config(tmp_path)


# --- GC: calibration --------------------------------------------------------

def test_parse_standard_name_units():
    from hplc_gc_pipeline.gc_calibration import parse_standard_name, GCCalibrationConfig
    pat = GCCalibrationConfig().standard_pattern
    assert parse_standard_name("250uM_34DMS_hexane", pat) == ("34DMS", 250.0)
    assert parse_standard_name("1mM_34DMS_5etac_95hexane", pat) == ("34DMS", 1000.0)
    # A sample name is not a standard.
    assert parse_standard_name("50_s11_A", pat) is None


def test_gc_calibration_fit_and_predict():
    from hplc_gc_pipeline.gc_calibration import _fit
    # area = 60*conc exactly.
    cal = _fit([(50, 3000), (100, 6000), (250, 15000)])
    assert cal.slope == pytest.approx(60.0)
    assert cal.r2 == pytest.approx(1.0)
    assert cal.predict(6000) == pytest.approx((6000 - cal.intercept) / 60.0)


def test_gc_force_through_origin_removes_floor():
    from hplc_gc_pipeline.gc_calibration import _fit
    # Points implying a negative intercept; forcing origin must give 0 at area 0.
    pts = [(50, 3000), (100, 6100), (250, 14900)]
    free = _fit(pts, force_origin=False)
    forced = _fit(pts, force_origin=True)
    assert forced.intercept == 0.0
    assert forced.predict(0.0) == 0.0
    # Free fit's non-zero intercept would predict a non-zero conc at area 0.
    assert free.predict(0.0, clamp_negative=False) != 0.0
