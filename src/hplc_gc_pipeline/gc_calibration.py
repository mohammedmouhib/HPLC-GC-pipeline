"""
GC calibration: area -> concentration (uM), Stage 2 for the GC path.

Two things distinguish this from the HPLC calibration module:

1. **Selectable source.** GC calibration is not repeated every experiment, so
   the standard curve can come from either:
     * ``injections`` -- dedicated standard ``.D`` folders whose concentration is
       encoded in the sample name (e.g. ``250uM_34DMS_hexane`` -> 250 uM). They
       are integrated by the same Stage-1 path as the samples, then fitted.
     * ``csv`` -- a table of (Compound, concentration, Area_Integral) rows.

2. **Fit direction.** The curve is fitted as ``area = slope * conc + intercept``
   (concentration, the known/error-free quantity, on the x-axis) and inverted for
   prediction: ``conc = (area - intercept) / slope``. This is the conventional
   calibration direction and yields a meaningful R^2. (The HPLC module fits the
   inverse, ``conc = m * area + b``; the GC form is the deliberate improvement.)

A per-sample **dilution factor** is applied to *sample* concentrations only --
standards carry their true prepared concentration and are never scaled.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Config (standalone for now; merged into the main Config in a later stage)
# ---------------------------------------------------------------------------

@dataclass
class GCCalibrationConfig:
    # Where the standard curve comes from.
    source: str = "injections"                 # "injections" | "csv"
    # Regex with named groups (conc, unit, compound) identifying + decoding
    # standard injections by sample name. Non-matching injections are samples.
    standard_pattern: str = r"^(?P<conc>\d+(?:\.\d+)?)\s*(?P<unit>[a-zA-Zµ]*M)_(?P<compound>[^_]+)"
    # CSV of (Compound, concentration, Area_Integral) used when source == "csv".
    standard_file: str = "gc_standard.csv"
    # Multiply *sample* concentrations by this to recover the original,
    # undiluted sample. GC samples are commonly diluted 1:10 -> factor 10.
    dilution_factor: float = 1.0
    # Fit area = slope*conc (intercept forced to 0). Common for trace GC-MS and
    # removes the "non-zero concentration at zero area" artefact a negative
    # intercept would otherwise produce.
    force_through_origin: bool = False
    # Clamp negative predicted concentrations (possible below the intercept) to 0.
    clamp_negative_to_zero: bool = True
    # Report samples whose qualifier ion-ratio failed as 0 uM (unconfirmed).
    require_ion_ratio_pass: bool = False
    # Saturation correction: when a compound has a ``saturation_mz`` backup ion
    # defined, high-concentration standards where the quantifier ion is
    # flat-topped (area > expected from backup × linear ratio) are corrected.
    # This value is the fractional excess above the linear-range quant/backup
    # ratio that flags a point as saturated (e.g. 0.03 = 3% excess triggers correction).
    saturation_ratio_tolerance: float = 0.03


# ---------------------------------------------------------------------------
# Concentration-unit handling
# ---------------------------------------------------------------------------

_UNIT_TO_UM = {"m": 1000.0, "u": 1.0, "µ": 1.0, "n": 0.001, "": 1_000_000.0}


def unit_to_uM(unit: str) -> float:
    """Factor converting a concentration in ``unit`` to uM (e.g. 'mM' -> 1000)."""
    unit = unit.strip()
    if not unit.endswith("M"):
        raise ValueError(f"unrecognised concentration unit: {unit!r}")
    prefix = unit[:-1].lower()
    if prefix not in _UNIT_TO_UM:
        raise ValueError(f"unrecognised concentration unit: {unit!r}")
    return _UNIT_TO_UM[prefix]


def parse_standard_name(name: str, pattern: str) -> Optional[tuple[str, float]]:
    """Decode a standard injection name into (compound, concentration_uM).

    Returns None if ``name`` doesn't match ``pattern`` (i.e. it's a sample),
    or if the pattern is empty / missing the required named groups.
    """
    if not pattern:
        return None
    m = re.match(pattern, str(name).strip())
    if not m:
        return None
    try:
        conc = float(m.group("conc")) * unit_to_uM(m.group("unit"))
        return m.group("compound"), conc
    except (IndexError, KeyError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Calibration model
# ---------------------------------------------------------------------------

@dataclass
class GCCalibration:
    """Linear calibration for one compound: ``area = slope * conc + intercept``."""
    slope: float
    intercept: float
    r2: float
    n: int
    conc_min: float
    conc_max: float
    points: list[tuple[float, float]] = field(default_factory=list)  # (conc_uM, area)

    def predict(self, area: float, clamp_negative: bool = True) -> float:
        if self.slope == 0 or pd.isna(area):
            return np.nan
        conc = (area - self.intercept) / self.slope
        return max(conc, 0.0) if clamp_negative else conc


def _fit(points: list[tuple[float, float]], force_origin: bool = False) -> Optional[GCCalibration]:
    """Fit area = slope*conc (+ intercept) from (conc, area) points."""
    pts = sorted(points)
    x = np.array([c for c, _ in pts], dtype=float)
    y = np.array([a for _, a in pts], dtype=float)
    if len(x) < 2 or np.ptp(x) == 0:
        return None
    if force_origin:
        slope = float((x @ y) / (x @ x))
        intercept = 0.0
    else:
        slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    ss_res = float(((y - yhat) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    return GCCalibration(float(slope), float(intercept), r2, len(x),
                         float(x.min()), float(x.max()), pts)


# ---------------------------------------------------------------------------
# Saturation correction
# ---------------------------------------------------------------------------

def _correct_saturation(
    triples: list[tuple[float, float, float]],
    compound_name: str,
    tolerance: float,
) -> list[tuple[float, float]]:
    """Detect and correct detector saturation in calibration standard points.

    When a quantifier ion saturates (detector clips the signal), the recorded
    peak area is artificially HIGH because the flat-topped peak occupies more
    time at the clipped value.  This function detects that excess by comparing
    the quant/backup ratio at each concentration against the median ratio of the
    lowest-concentration (unsaturated) standards.  Points where the ratio
    exceeds the reference by more than ``tolerance`` are replaced with
    ``backup_area × reference_ratio``.

    Parameters
    ----------
    triples : [(conc_uM, quant_area, backup_area), ...] sorted by concentration
    compound_name : str  -- for diagnostic printing
    tolerance : float   -- fractional excess (e.g. 0.03 = 3%) above reference
                           ratio that flags a point as saturated
    """
    pts = sorted(triples)
    ratios: list[Optional[float]] = [
        q / b if b > 0 else None for _, q, b in pts
    ]

    # Reference ratio from the lowest-concentration half (assumed unsaturated)
    n_ref = max(2, len(pts) // 2)
    ref_values = [r for r in ratios[:n_ref] if r is not None]
    if not ref_values:
        return [(c, q) for c, q, _ in pts]
    ref_ratio = float(np.median(ref_values))

    corrected: list[tuple[float, float]] = []
    n_fixed = 0
    for i, (conc, quant, backup) in enumerate(pts):
        r = ratios[i]
        if r is not None and r > ref_ratio * (1.0 + tolerance) and backup > 0:
            fixed = backup * ref_ratio
            corrected.append((conc, fixed))
            n_fixed += 1
            print(f"    {compound_name} {conc:.0f} uM: area {quant:.0f} → {fixed:.0f} "
                  f"(saturation corr., ratio {r:.3f} vs ref {ref_ratio:.3f})")
        else:
            corrected.append((conc, quant))

    if n_fixed:
        print(f"    → {compound_name}: {n_fixed} standard(s) corrected for quant-ion saturation")

    return corrected


# ---------------------------------------------------------------------------
# Building the calibration (two sources)
# ---------------------------------------------------------------------------

def build_from_injections(peak_df: pd.DataFrame, pattern: str,
                          force_origin: bool = False,
                          sat_tol: float = 0.03) -> dict[str, GCCalibration]:
    """Build calibrations from the standard rows of an integrated peak table.

    ``peak_df`` is the full Stage-1 GC peak table (standards + samples). Rows
    whose ``Sample`` matches ``pattern`` are treated as standards; their decoded
    concentration + integrated ``Area_Integral`` become calibration points.

    When ``Saturation_Backup_Area`` is present, saturation correction is applied
    per compound before fitting: points where the quant/backup ratio exceeds the
    linear-range reference by more than ``sat_tol`` are corrected.
    """
    has_role = "Role" in peak_df.columns
    has_backup = "Saturation_Backup_Area" in peak_df.columns

    # points: compound -> list of (conc, area) for calib_probe and direct rows
    # without backup data.
    # backup_triples: compound -> list of (conc, area, backup_area) for direct
    # standard rows that have backup data (saturation correction candidates).
    points: dict[str, list[tuple[float, float]]] = {}
    backup_triples: dict[str, list[tuple[float, float, float]]] = {}

    for _, row in peak_df.iterrows():
        parsed = parse_standard_name(row["Sample"], pattern)
        if parsed is None:
            continue
        std_compound, conc_uM = parsed
        role = str(row["Role"]) if has_role else "sample"
        area = float(row["Area_Integral"])

        if role == "calib_probe":
            # Surrogate curve: this compound's quantifier ion integrated over the
            # calibrant's peak. Use it only in the calibrant's own standards.
            if str(row.get("Calibrant")) != std_compound:
                continue
            points.setdefault(str(row["Compound"]), []).append((conc_uM, area))
            continue

        # A standard injection is for one compound; only its own row contributes.
        if str(row["Compound"]) != std_compound:
            continue

        compound_name = str(row["Compound"])
        backup = float(row["Saturation_Backup_Area"]) if (
            has_backup and pd.notna(row.get("Saturation_Backup_Area"))
        ) else None

        if backup is not None:
            backup_triples.setdefault(compound_name, []).append((conc_uM, area, backup))
        else:
            points.setdefault(compound_name, []).append((conc_uM, area))

    # Apply saturation correction for compounds that have backup areas
    for compound_name, triples in backup_triples.items():
        corrected = _correct_saturation(triples, compound_name, sat_tol)
        points.setdefault(compound_name, []).extend(corrected)

    return _finalise(points, force_origin)


def build_from_csv(csv_path: Path, force_origin: bool = False) -> dict[str, GCCalibration]:
    """Build calibrations from a (Compound, concentration, Area_Integral) CSV."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"GC standard CSV not found: {csv_path}")
    df = pd.read_csv(csv_path)
    required = {"Compound", "concentration", "Area_Integral"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"GC standard CSV missing columns: {sorted(missing)}")
    df = df.dropna(subset=["Compound", "concentration", "Area_Integral"])
    points: dict[str, list[tuple[float, float]]] = {}
    for _, row in df.iterrows():
        points.setdefault(str(row["Compound"]), []).append(
            (float(row["concentration"]), float(row["Area_Integral"]))
        )
    return _finalise(points, force_origin)


def _finalise(points: dict[str, list[tuple[float, float]]],
              force_origin: bool = False) -> dict[str, GCCalibration]:
    print("\n" + "=" * 70)
    print("GC CALIBRATION CURVES")
    print("=" * 70)
    cal: dict[str, GCCalibration] = {}
    for compound, pts in sorted(points.items()):
        fit = _fit(pts, force_origin)
        if fit is None:
            print(f"  Warning: {compound} has <2 usable calibration points; skipped")
            continue
        cal[compound] = fit
        flag = "" if fit.r2 >= 0.99 else "  <- R2 below 0.99, check standards"
        print(f"  {compound:<12}: area = {fit.slope:.4g} * uM + {fit.intercept:.4g}  "
              f"R2={fit.r2:.4f}  (n={fit.n}, {fit.conc_min:g}-{fit.conc_max:g} uM){flag}")
    return cal


def build_calibration(peak_df: pd.DataFrame, cfg: GCCalibrationConfig,
                      experiment_dir: Optional[Path] = None) -> dict[str, GCCalibration]:
    """Dispatch to the configured calibration source."""
    if cfg.source == "csv":
        csv_path = Path(cfg.standard_file)
        if experiment_dir is not None and not csv_path.is_absolute():
            csv_path = experiment_dir / cfg.standard_file
        return build_from_csv(csv_path, cfg.force_through_origin)
    if cfg.source == "injections":
        return build_from_injections(peak_df, cfg.standard_pattern, cfg.force_through_origin,
                                     sat_tol=cfg.saturation_ratio_tolerance)
    raise ValueError(f"unknown GC calibration source: {cfg.source!r} (use 'injections' or 'csv')")


# ---------------------------------------------------------------------------
# Applying the calibration to samples
# ---------------------------------------------------------------------------

def add_concentrations(peak_df: pd.DataFrame, calibration: dict[str, GCCalibration],
                       cfg: GCCalibrationConfig) -> pd.DataFrame:
    """Add concentration columns to the *sample* rows of a GC peak table.

    Standard injections (name matches ``cfg.standard_pattern``) are labelled and
    left unquantified. Sample rows gain:

    * ``Concentration_uM_measured`` -- concentration in the injected solution,
    * ``Dilution_Factor``           -- the applied factor,
    * ``Concentration_uM``          -- original-sample concentration (measured x factor),
    * ``Below_Cal_Range``           -- True if measured conc is below the lowest
      standard (quantified by extrapolation, less reliable).
    """
    df = peak_df.copy()

    is_standard = df["Sample"].apply(
        lambda s: parse_standard_name(s, cfg.standard_pattern) is not None
    )
    df["Is_Standard"] = is_standard

    def measured(row) -> float:
        # Each compound uses its own curve. For a surrogate (calibrate_as set)
        # that curve was built from the calibrant's standards via probe rows, so
        # it is already keyed by this compound's name.
        cal = calibration.get(row["Compound"])
        if cal is None:
            return np.nan
        if cfg.require_ion_ratio_pass and row.get("Ion_Ratio_Pass") is False:
            return 0.0
        return cal.predict(row["Area_Integral"], clamp_negative=cfg.clamp_negative_to_zero)

    def below_range(row) -> bool:
        cal = calibration.get(row["Compound"])
        m = row["Concentration_uM_measured"]
        if cal is None or pd.isna(m):
            return False
        return 0 < m < cal.conc_min

    df["Concentration_uM_measured"] = df.apply(
        lambda r: np.nan if r["Is_Standard"] else measured(r), axis=1
    )
    df["Below_Cal_Range"] = df.apply(lambda r: False if r["Is_Standard"] else below_range(r), axis=1)
    df["Dilution_Factor"] = np.where(is_standard, np.nan, cfg.dilution_factor)
    df["Concentration_uM"] = df["Concentration_uM_measured"] * df["Dilution_Factor"]

    n_below = int(df["Below_Cal_Range"].sum())
    if n_below:
        print(f"\n  Note: {n_below} sample(s) quantified below the lowest standard "
              f"(extrapolation) -- see Below_Cal_Range flag.")
    return df
