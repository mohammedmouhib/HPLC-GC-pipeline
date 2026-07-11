"""
Calibration: area -> concentration (uM), Stage 2.

A linear fit (degree 1) per compound from ``standard.csv``. Two config-gated
behaviours applied when converting areas:

* ``clamp_negative_to_zero`` -- a linear fit can predict a negative
  concentration near the origin; clamp it to 0.
* ``zero_area_zero_conc`` -- a zero-filled (undetected) compound has area 0,
  which the fit maps to the intercept ``b`` rather than 0. When enabled, area 0
  is reported as concentration 0. This is the fix for the original "pCA ~15 uM
  at zero area" artefact.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import CalibrationConfig


def load_standards(standard_csv: Path) -> pd.DataFrame:
    """Load calibration points (Compound, Area_Integral, concentration)."""
    standard_csv = Path(standard_csv)
    if not standard_csv.exists():
        raise FileNotFoundError(f"standard.csv not found: {standard_csv}")
    df = pd.read_csv(standard_csv)
    return df.dropna(subset=["Compound", "Area_Integral", "concentration"])


def build_calibration_functions(standards_df: pd.DataFrame, compounds_df: pd.DataFrame | None = None) -> dict:
    """Build per-compound linear calibration: concentration = m*area + b.

    If ``compounds_df`` is given, warn about any compound that has an RT window
    but no calibration data (it will produce NaN concentrations).
    """
    print("\n" + "=" * 70)
    print("BUILDING CALIBRATION CURVES")
    print("=" * 70)

    cal = {}
    for compound, sub in standards_df.groupby("Compound"):
        x = sub["Area_Integral"].values
        y = sub["concentration"].values
        if len(x) < 2:
            print(f"  Warning: {compound} has <2 calibration points; skipped")
            continue
        m, b = np.polyfit(x, y, 1)
        cal[compound] = (float(m), float(b))
        print(f"  {compound:<10}: conc = {m:.4g} * area + {b:.4g}  (uM, n={len(x)})")

    if compounds_df is not None:
        uncalibrated = [c for c in compounds_df["Compound"] if c not in cal]
        if uncalibrated:
            print(f"\n  Warning: no calibration data for {len(uncalibrated)} compound(s); "
                  f"they will report NaN concentration:")
            for c in uncalibrated:
                print(f"    - {c}")
    return cal


def area_to_concentration(area: float, compound: str, calibration: dict, cfg: CalibrationConfig) -> float:
    """Convert one area value to a concentration in uM."""
    if compound not in calibration:
        return np.nan
    if cfg.zero_area_zero_conc and (pd.isna(area) or area == 0):
        return 0.0
    m, b = calibration[compound]
    conc = m * area + b
    if cfg.clamp_negative_to_zero:
        conc = max(conc, 0.0)
    return conc


def add_concentration_column(df: pd.DataFrame, calibration: dict, cfg: CalibrationConfig) -> pd.DataFrame:
    """Add a ``Concentration_uM`` column to a consolidated peak frame."""
    df = df.copy()
    df["Concentration_uM"] = df.apply(
        lambda row: area_to_concentration(row["Area_Integral"], row["Compound"], calibration, cfg),
        axis=1,
    )
    return df
