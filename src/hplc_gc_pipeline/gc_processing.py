"""
Stage 1 core for GC-MS: integrate targeted compounds from ``.D`` injections.

Unlike the HPLC path (MOCCA2 deconvolution over a DAD wavelength), GC-MS here is
SIM (Selected Ion Monitoring): each injection records a handful of m/z channels.
Quantification is *targeted* -- for every compound in the GC compound table we
extract its quantifier-ion chromatogram (EIC), correct the baseline, find the
peak inside the compound's retention-time window, integrate it, and confirm
identity from the qualifier-ion ratios.

Design mirror of ``processing.py``: this module reads config + loaded traces and
returns a peak table with the same core columns (``Sample``, ``Folder``,
``Retention_Time_min``, ``Area_Integral`` ...) plus GC-specific columns.

The quantification channel is configurable:
* ``eic`` -- quantifier-ion extracted chromatogram (default, most selective)
* ``tic`` -- total ion current (one trace for all compounds)
* ``fid`` -- flame-ionization detector trace
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .gc_agilent import GCInjection, GCTraces, find_gc_injections, load_gc_traces


# ---------------------------------------------------------------------------
# Config (standalone for now; merged into the main Config in a later stage)
# ---------------------------------------------------------------------------

@dataclass
class GCPeakDetectionConfig:
    # Signal-to-noise a window's apex must clear to count as detected. Noise is
    # a robust (MAD-based) estimate of the baseline-corrected trace.
    min_snr: float = 5.0
    # Fraction of apex height (and multiple of noise) marking the peak edges
    # when walking outward from the apex to set integration bounds.
    edge_fraction: float = 0.01


@dataclass
class GCBaselineConfig:
    # Asymmetric least squares (Eilers) baseline. lam: smoothness; p: asymmetry
    # (small -> baseline hugs the lower envelope). niter: reweighting passes.
    lam: float = 1e5
    p: float = 0.01
    niter: int = 10


@dataclass
class GCProcessingConfig:
    quant_channel: str = "eic"          # "eic" | "tic" | "fid"
    peak_detection: GCPeakDetectionConfig = field(default_factory=GCPeakDetectionConfig)
    baseline: GCBaselineConfig = field(default_factory=GCBaselineConfig)


# ---------------------------------------------------------------------------
# GC compound table
# ---------------------------------------------------------------------------

@dataclass
class GCCompound:
    name: str
    rt_low: float
    rt_high: float
    quantifier_mz: float
    # (m/z, expected_ratio-or-None) for each qualifier ion.
    qualifiers: list[tuple[float, Optional[float]]] = field(default_factory=list)
    ion_ratio_tol: float = 0.30         # relative tolerance on qualifier ratios


def _parse_qualifiers(raw) -> list[tuple[float, Optional[float]]]:
    """Parse a ``qualifier_mz`` cell like ``"149=0.40; 91=0.36"`` or ``"149;91"``."""
    if raw is None or (isinstance(raw, float) and pd.isna(raw)):
        return []
    out: list[tuple[float, Optional[float]]] = []
    for token in str(raw).replace(",", ";").split(";"):
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            mz_s, ratio_s = token.split("=", 1)
            out.append((float(mz_s.strip()), float(ratio_s.strip())))
        else:
            out.append((float(token), None))
    return out


def load_gc_compounds(compounds_csv: Path) -> list[GCCompound]:
    """Load the GC compound table.

    Required columns: Compound, RT_low, RT_high, quantifier_mz.
    Optional columns: qualifier_mz (e.g. "149=0.40;91=0.36"), ion_ratio_tol, Notes.
    """
    compounds_csv = Path(compounds_csv)
    if not compounds_csv.exists():
        raise FileNotFoundError(
            f"GC compounds file not found: {compounds_csv}\n"
            f"Required columns: Compound, RT_low, RT_high, quantifier_mz"
        )
    df = pd.read_csv(compounds_csv)
    required = {"Compound", "RT_low", "RT_high", "quantifier_mz"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"GC compounds file missing required columns: {sorted(missing)}")

    compounds = []
    print("\n" + "=" * 70)
    print("GC COMPOUND DEFINITIONS")
    print("=" * 70)
    for _, r in df.iterrows():
        quals = _parse_qualifiers(r.get("qualifier_mz"))
        tol = float(r["ion_ratio_tol"]) if "ion_ratio_tol" in df.columns and pd.notna(r.get("ion_ratio_tol")) else 0.30
        c = GCCompound(
            name=str(r["Compound"]),
            rt_low=float(r["RT_low"]),
            rt_high=float(r["RT_high"]),
            quantifier_mz=float(r["quantifier_mz"]),
            qualifiers=quals,
            ion_ratio_tol=tol,
        )
        compounds.append(c)
        qual_str = ", ".join(
            f"{int(m)}" + (f"={ratio:.2f}" if ratio is not None else "") for m, ratio in quals
        ) or "(none)"
        print(f"  {c.name:<16} RT {c.rt_low:.2f}-{c.rt_high:.2f} min | "
              f"quant m/z {int(c.quantifier_mz)} | qualifiers {qual_str}")
    return compounds


# ---------------------------------------------------------------------------
# Signal processing
# ---------------------------------------------------------------------------

def als_baseline(y: np.ndarray, lam: float, p: float, niter: int) -> np.ndarray:
    """Asymmetric least squares baseline (Eilers & Boelens 2005).

    Points above the current baseline get weight ``p``; points below get
    ``1-p``, so the fit is pulled toward the lower envelope while staying smooth
    (``lam`` penalises curvature).
    """
    from scipy import sparse
    from scipy.sparse.linalg import spsolve

    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return np.zeros_like(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(n, n - 2))
    DTD = lam * (D @ D.transpose())
    w = np.ones(n)
    z = y.copy()
    for _ in range(max(1, niter)):
        W = sparse.spdiags(w, 0, n, n)
        z = spsolve((W + DTD).tocsc(), w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def _noise_estimate(sig: np.ndarray) -> float:
    """Robust noise level: MAD scaled to a Gaussian sigma equivalent."""
    med = np.median(sig)
    mad = np.median(np.abs(sig - med))
    return float(1.4826 * mad) or float(np.std(sig)) or 1.0


def _integration_bounds(sig: np.ndarray, apex: int, threshold: float) -> tuple[int, int]:
    """Walk outward from ``apex`` until the signal drops below ``threshold`` or
    starts rising again (a valley), giving the peak's integration limits."""
    n = len(sig)
    left = apex
    while left > 0 and sig[left - 1] >= threshold and sig[left - 1] <= sig[left]:
        left -= 1
    right = apex
    while right < n - 1 and sig[right + 1] >= threshold and sig[right + 1] <= sig[right]:
        right += 1
    return left, right


@dataclass
class GCPeakResult:
    detected: bool
    rt_min: Optional[float]
    area: float
    height: float
    width_min: Optional[float]
    left_idx: Optional[int]
    right_idx: Optional[int]
    noise: float


def integrate_peak(time_min: np.ndarray, sig_corrected: np.ndarray,
                   rt_low: float, rt_high: float,
                   pk_cfg: GCPeakDetectionConfig) -> GCPeakResult:
    """Find and integrate the largest peak inside ``[rt_low, rt_high]``.

    ``sig_corrected`` must be baseline-corrected (baseline ~0 between peaks).
    Returns a zero-area, not-detected result when nothing clears the S/N floor.
    """
    noise = _noise_estimate(sig_corrected)
    window = (time_min >= rt_low) & (time_min <= rt_high)
    if not window.any():
        return GCPeakResult(False, None, 0.0, 0.0, None, None, None, noise)

    win_idx = np.where(window)[0]
    apex = win_idx[int(np.argmax(sig_corrected[win_idx]))]
    apex_val = float(sig_corrected[apex])

    if apex_val < pk_cfg.min_snr * noise:
        return GCPeakResult(False, None, 0.0, apex_val, None, None, None, noise)

    threshold = max(pk_cfg.edge_fraction * apex_val, pk_cfg.min_snr * noise * 0.5, noise)
    left, right = _integration_bounds(sig_corrected, apex, threshold)
    seg_t = time_min[left:right + 1]
    seg_y = np.clip(sig_corrected[left:right + 1], 0, None)
    area = float(np.trapezoid(seg_y, seg_t)) if right > left else 0.0
    width = float(seg_t[-1] - seg_t[0]) if right > left else None
    return GCPeakResult(True, float(time_min[apex]), area, apex_val, width, left, right, noise)


def _channel_signal(traces: GCTraces, compound: GCCompound, quant_channel: str):
    """Return (time, raw_signal) for the configured quantification channel."""
    if quant_channel == "tic":
        return traces.time_min, traces.tic
    if quant_channel == "fid":
        if traces.fid_signal is None:
            raise ValueError(f"{traces.sample_name}: no FID channel available")
        return traces.fid_time_min, traces.fid_signal
    # default: quantifier-ion EIC
    return traces.time_min, traces.eic_for(compound.quantifier_mz)


def _qualifier_report(traces: GCTraces, compound: GCCompound,
                      result: GCPeakResult) -> tuple[str, Optional[bool]]:
    """Observed qualifier/quantifier area ratios and a pass/fail verdict.

    Ratios are area ratios over the integrated bounds. The verdict is True/False
    only when at least one qualifier has an expected ratio; otherwise None (the
    ratios are still reported for the record).
    """
    if not compound.qualifiers or result.left_idx is None:
        return "", None
    lo, hi = result.left_idx, result.right_idx + 1
    quant_area = float(np.clip(traces.eic_for(compound.quantifier_mz)[lo:hi], 0, None).sum())
    if quant_area <= 0:
        return "", None

    parts, verdicts = [], []
    for mz, expected in compound.qualifiers:
        try:
            qual_area = float(np.clip(traces.eic_for(mz)[lo:hi], 0, None).sum())
        except ValueError:
            continue
        ratio = qual_area / quant_area
        tag = f"{int(mz)}={ratio:.2f}"
        if expected is not None and expected > 0:
            ok = abs(ratio - expected) / expected <= compound.ion_ratio_tol
            verdicts.append(ok)
            tag += "ok" if ok else "X"
        parts.append(tag)
    verdict = all(verdicts) if verdicts else None
    return "; ".join(parts), verdict


# ---------------------------------------------------------------------------
# Peak-table extraction
# ---------------------------------------------------------------------------

def build_gc_peak_table(injections: list[GCInjection], compounds: list[GCCompound],
                        cfg: GCProcessingConfig) -> pd.DataFrame:
    """Integrate every compound in every injection into a peak table.

    One row per (injection x compound). Undetected compounds get a zero-area row
    (analogous to the HPLC zero-fill), so downstream statistics see every
    compound for every injection.
    """
    print("\n" + "=" * 70)
    print(f"INTEGRATING GC COMPOUNDS (channel: {cfg.quant_channel})")
    print("=" * 70)

    rows: list[dict] = []
    for i, inj in enumerate(injections, 1):
        try:
            traces = load_gc_traces(inj)
        except Exception as exc:
            print(f"  [{i}/{len(injections)}] ERR {inj.folder_name}: {exc}")
            continue

        n_detected = 0
        for compound in compounds:
            try:
                time_axis, raw = _channel_signal(traces, compound, cfg.quant_channel)
                baseline = als_baseline(raw, cfg.baseline.lam, cfg.baseline.p, cfg.baseline.niter)
                corrected = raw - baseline
                res = integrate_peak(time_axis, corrected, compound.rt_low, compound.rt_high,
                                     cfg.peak_detection)
                qual_str, qual_pass = ("", None)
                if res.detected and cfg.quant_channel == "eic":
                    qual_str, qual_pass = _qualifier_report(traces, compound, res)
                if res.detected:
                    n_detected += 1
            except Exception as exc:
                print(f"      {inj.sample_name} / {compound.name}: {exc}")
                res = GCPeakResult(False, None, 0.0, 0.0, None, None, None, 0.0)
                qual_str, qual_pass = "", None

            rows.append({
                "Experiment": inj.experiment_name,
                "Folder": inj.folder_name,
                "Sample": inj.sample_name,
                "Peak_Number": 1,
                "Compound": compound.name,
                "Channel": cfg.quant_channel,
                "Quantifier_mz": compound.quantifier_mz,
                "Retention_Time_min": round(res.rt_min, 3) if res.rt_min is not None else None,
                "Area_Integral": round(res.area, 2),
                "Height": round(res.height, 2),
                "Width_min": round(res.width_min, 4) if res.width_min is not None else None,
                "Qualifier_Ratios": qual_str,
                "Ion_Ratio_Pass": qual_pass,
                "Detected": res.detected,
                "Integration_Method": f"GC_{cfg.quant_channel.upper()}_Trapezoid",
            })
        print(f"  [{i}/{len(injections)}] {inj.folder_name} '{inj.sample_name}' "
              f"-> {n_detected}/{len(compounds)} compounds detected")

    df = pd.DataFrame(rows)
    print(f"\n  Total rows: {len(df)}  ({df['Detected'].sum() if not df.empty else 0} detected)")
    return df


def process_gc_experiment(data_root: Path, compounds_csv: Path,
                          cfg: GCProcessingConfig) -> pd.DataFrame:
    """Convenience: discover injections + compounds and build the peak table."""
    injections = find_gc_injections(data_root)
    if not injections:
        raise ValueError(f"No GC-MS .D folders (with data.ms) found in {data_root}")
    print(f"Found {len(injections)} GC-MS injections in {data_root}")
    compounds = load_gc_compounds(compounds_csv)
    return build_gc_peak_table(injections, compounds, cfg)
