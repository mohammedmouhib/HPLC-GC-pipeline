"""
Calibration: area -> concentration (uM), Stage 2.

A linear fit per compound from calibration data that may come from:

* ``standard.csv``         -- pre-existing (Compound, Area_Integral, concentration) rows
* standard .D injections   -- sample names matching a regex (e.g. "250uM_pCA_hexane");
                              their RT-matched peaks are extracted from peak_results.csv
* both                     -- pooled from both sources, fitted together

Additional capabilities (matching the GC calibration module):

* ``dilution_factor``     -- multiply *sample* concentrations only; standards unchanged
* ``force_through_origin``-- fit conc = slope*area (b forced to 0)
* ``calibrate_as``        -- a column in compounds.csv; compound X borrows compound Y's
                             calibration curve (shared chromophore / response factor)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .config import CalibrationConfig


# ---------------------------------------------------------------------------
# Unit handling (shared with gc_calibration)
# ---------------------------------------------------------------------------

_UNIT_TO_UM = {"m": 1000.0, "u": 1.0, "µ": 1.0, "n": 0.001, "": 1_000_000.0}


def _unit_to_uM(unit: str) -> float:
    unit = unit.strip()
    if not unit.endswith("M"):
        raise ValueError(f"unrecognised concentration unit: {unit!r}")
    prefix = unit[:-1].lower()
    if prefix not in _UNIT_TO_UM:
        raise ValueError(f"unrecognised concentration unit: {unit!r}")
    return _UNIT_TO_UM[prefix]


def _parse_standard_name(name: str, pattern: str) -> Optional[tuple[str, float]]:
    """Return (compound, conc_uM) if *name* matches *pattern*, else None."""
    m = re.match(pattern, str(name).strip())
    if not m:
        return None
    conc = float(m.group("conc")) * _unit_to_uM(m.group("unit"))
    return m.group("compound"), conc


# ---------------------------------------------------------------------------
# Loading calibration data
# ---------------------------------------------------------------------------

def _load_from_csv(csv_path: Path) -> pd.DataFrame:
    """Load (Compound, Area_Integral, concentration) rows from *csv_path*."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return pd.DataFrame(columns=["Compound", "Area_Integral", "concentration"])
    df = pd.read_csv(csv_path)
    df = df.dropna(subset=["Compound", "Area_Integral", "concentration"])
    df["_source"] = "csv"
    return df[["Compound", "Area_Integral", "concentration", "_source"]]


def _load_from_injections(
    peak_results_csv: Path,
    compounds_df: pd.DataFrame,
    pattern: str,
) -> pd.DataFrame:
    """Detect standard injections in *peak_results_csv* and return calibration rows.

    A row contributes a calibration point when:
    1. Its ``Sample`` name matches *pattern* (yields compound + conc_uM).
    2. The peak's ``Retention_Time_min`` falls inside that compound's RT window
       in *compounds_df* — this selects the right peak from an injection that
       may contain multiple peaks.
    """
    peak_results_csv = Path(peak_results_csv)
    if not peak_results_csv.exists():
        return pd.DataFrame(columns=["Compound", "Area_Integral", "concentration"])

    peak_df = pd.read_csv(peak_results_csv)
    rows: list[dict] = []

    for _, row in peak_df.iterrows():
        parsed = _parse_standard_name(str(row.get("Sample", "")), pattern)
        if parsed is None:
            continue
        std_compound, conc_uM = parsed

        rt = row.get("Retention_Time_min")
        if pd.isna(rt):
            continue
        area = row.get("Area_Integral")
        if pd.isna(area):
            continue

        # Match peak RT to the target compound's RT window.
        matched = None
        for _, c in compounds_df.iterrows():
            if c["RT_low"] <= float(rt) <= c["RT_high"]:
                matched = c["Compound"]
                break
        if matched != std_compound:
            continue

        rows.append({
            "Compound": std_compound,
            "Area_Integral": float(area),
            "concentration": conc_uM,
            "_source": "injection",
        })

    if not rows:
        return pd.DataFrame(columns=["Compound", "Area_Integral", "concentration", "_source"])
    return pd.DataFrame(rows)


def load_standards(
    standard_csv: Path,
    peak_results_csv: Optional[Path] = None,
    compounds_df: Optional[pd.DataFrame] = None,
    cfg: Optional[CalibrationConfig] = None,
) -> pd.DataFrame:
    """Load calibration points, pooling sources according to *cfg.source*.

    * ``source: csv``         — read *standard_csv* only (default, backward-compatible).
    * ``source: injections``  — detect standard .D injections in *peak_results_csv*.
    * ``source: both``        — pool both; regression uses all points together.

    When *cfg* is None, behaves as ``source: csv`` (original behaviour).
    """
    source = getattr(cfg, "source", "csv")
    pattern = getattr(cfg, "standard_pattern",
                      r"^(?P<conc>\d+(?:\.\d+)?)\s*(?P<unit>[a-zA-Zµ]*M)_(?P<compound>[^_]+)")

    frames: list[pd.DataFrame] = []

    if source in ("csv", "both"):
        df_csv = _load_from_csv(standard_csv)
        if not df_csv.empty:
            frames.append(df_csv)

    if source in ("injections", "both"):
        if peak_results_csv is None or compounds_df is None:
            print("  Warning: calibration source includes 'injections' but "
                  "peak_results_csv or compounds_df was not provided; skipping injection detection.")
        else:
            df_inj = _load_from_injections(peak_results_csv, compounds_df, pattern)
            if not df_inj.empty:
                frames.append(df_inj)

    if not frames:
        src_desc = f"source={source!r}"
        if source in ("csv", "both"):
            src_desc += f", file={standard_csv}"
        raise FileNotFoundError(
            f"No calibration data found ({src_desc}). "
            "Provide a standard.csv or standard .D injections matching the standard_pattern."
        )

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.dropna(subset=["Compound", "Area_Integral", "concentration"])

    n_csv = int((combined.get("_source", pd.Series()) == "csv").sum()) if "_source" in combined.columns else len(combined)
    n_inj = int((combined.get("_source", pd.Series()) == "injection").sum()) if "_source" in combined.columns else 0
    if source == "both":
        print(f"  Calibration data: {n_csv} row(s) from CSV, {n_inj} point(s) from injections (pooled).")

    return combined


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------

def _fit_one(x: np.ndarray, y: np.ndarray,
             force_origin: bool = False) -> tuple[float, float]:
    """Fit y = m*x + b (or y = m*x if force_origin). Returns (m, b)."""
    if force_origin:
        m = float((x @ y) / (x @ x)) if (x @ x) != 0 else 0.0
        return m, 0.0
    m, b = np.polyfit(x, y, 1)
    return float(m), float(b)


def build_calibration_functions(
    standards_df: pd.DataFrame,
    compounds_df: Optional[pd.DataFrame] = None,
    cfg: Optional[CalibrationConfig] = None,
) -> dict:
    """Build per-compound linear calibration: concentration = m*area + b.

    When *compounds_df* has a ``calibrate_as`` column, a compound with a
    non-null entry borrows the named compound's calibration curve instead of
    fitting its own.  This handles shared chromophores / response factors.
    """
    force_origin = getattr(cfg, "force_through_origin", False)

    print("\n" + "=" * 70)
    print("BUILDING CALIBRATION CURVES")
    print("=" * 70)

    cal: dict[str, tuple[float, float]] = {}
    for compound, sub in standards_df.groupby("Compound"):
        x = sub["Area_Integral"].values.astype(float)
        y = sub["concentration"].values.astype(float)
        if len(x) < 2:
            print(f"  Warning: {compound} has <2 calibration points; skipped")
            continue
        m, b = _fit_one(x, y, force_origin)
        cal[compound] = (m, b)
        print(f"  {compound:<12}: conc = {m:.4g} * area + {b:.4g}  (µM, n={len(x)})")

    # Apply calibrate_as cross-calibration
    if compounds_df is not None and "calibrate_as" in compounds_df.columns:
        for _, row in compounds_df.iterrows():
            src = row.get("calibrate_as")
            if pd.isna(src) or not str(src).strip():
                continue
            compound = row["Compound"]
            src = str(src).strip()
            if src in cal:
                cal[compound] = cal[src]
                print(f"  {compound:<12}: borrows calibration from {src!r}")
            else:
                print(f"  Warning: {compound} calibrate_as={src!r} but {src!r} "
                      f"has no calibration curve; {compound} will report NaN")

    if compounds_df is not None:
        uncalibrated = [c for c in compounds_df["Compound"] if c not in cal]
        if uncalibrated:
            print(f"\n  Warning: no calibration data for {len(uncalibrated)} compound(s); "
                  "they will report NaN concentration:")
            for c in uncalibrated:
                print(f"    - {c}")

    return cal


# ---------------------------------------------------------------------------
# Applying calibration to samples
# ---------------------------------------------------------------------------

def area_to_concentration(
    area: float,
    compound: str,
    calibration: dict,
    cfg: CalibrationConfig,
) -> float:
    """Convert one area value to a measured concentration (uM, before dilution)."""
    if compound not in calibration:
        return np.nan
    if cfg.zero_area_zero_conc and (pd.isna(area) or area == 0):
        return 0.0
    m, b = calibration[compound]
    conc = m * area + b
    if cfg.clamp_negative_to_zero:
        conc = max(conc, 0.0)
    return conc


def add_concentration_column(
    df: pd.DataFrame,
    calibration: dict,
    cfg: CalibrationConfig,
) -> pd.DataFrame:
    """Add concentration columns to a consolidated peak frame.

    Adds:
    * ``Concentration_uM_measured`` — concentration in the injected solution
    * ``Dilution_Factor``           — the applied factor (cfg.dilution_factor)
    * ``Concentration_uM``          — original-sample concentration (measured × factor)
    """
    df = df.copy()
    df["Concentration_uM_measured"] = df.apply(
        lambda row: area_to_concentration(
            row["Area_Integral"], row["Compound"], calibration, cfg
        ),
        axis=1,
    )
    dilution = getattr(cfg, "dilution_factor", 1.0)
    df["Dilution_Factor"] = dilution
    df["Concentration_uM"] = df["Concentration_uM_measured"] * dilution
    return df
