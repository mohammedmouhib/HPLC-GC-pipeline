"""
Consolidation: sum sub-peaks within a compound window and zero-fill (Stage 2).

Produces one row per (injection x compound). Every compound appears for every
injection even when undetected (area 0), so downstream stats and plots operate
on a complete grid. This is the canonical plot input (consolidated_peaks.csv).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def consolidate_compound_areas(assigned_df: pd.DataFrame, compounds_df: pd.DataFrame, all_samples_df: pd.DataFrame) -> pd.DataFrame:
    """Sum sub-peaks per (injection, compound), then zero-fill the full grid."""
    print("\n" + "=" * 70)
    print("CONSOLIDATING COMPOUND AREAS (SUMMATION + ZERO-FILL)")
    print("=" * 70)

    meta_cols = [
        c for c in [
            "Experiment", "Folder", "Sample", "Strain",
            "Time_Point_h", "Replicate", "Replicate_Group", "Wavelength_nm",
        ]
        if c in assigned_df.columns
    ]

    group_sizes = assigned_df.groupby(["Folder", "Compound"]).size()
    multi_peak = group_sizes[group_sizes > 1]
    if len(multi_peak) > 0:
        print(f"\n  Multi-peak summations ({len(multi_peak)} cases):")
        for (folder, compound), count in multi_peak.items():
            rows = assigned_df[(assigned_df["Folder"] == folder) & (assigned_df["Compound"] == compound)]
            rts = rows["Retention_Time_min"].dropna().round(3).tolist()
            total = rows["Area_Integral"].sum()
            print(f"    {folder} | {compound}: {count} peaks at RT {rts} -> summed area = {total:.2f}")
    else:
        print("\n  No multi-peak summations detected.")

    consolidated = (
        assigned_df
        .groupby(meta_cols + ["Compound"], as_index=False)
        .agg(
            Area_Integral=("Area_Integral", "sum"),
            Height_mAU=("Height_mAU", "max"),
            n_peaks_summed=("Area_Integral", "count"),
        )
    )

    compounds_list = sorted(compounds_df["Compound"].tolist())
    sample_meta = all_samples_df[meta_cols].drop_duplicates().reset_index(drop=True)
    sample_meta["_key"] = 1
    compounds_frame = pd.DataFrame({"Compound": compounds_list, "_key": 1})
    full_grid = sample_meta.merge(compounds_frame, on="_key").drop(columns="_key")

    merge_keys = meta_cols + ["Compound"]
    result = full_grid.merge(
        consolidated[merge_keys + ["Area_Integral", "Height_mAU", "n_peaks_summed"]],
        on=merge_keys,
        how="left",
    )

    result["Area_Integral"] = result["Area_Integral"].fillna(0.0)
    result["Height_mAU"] = result["Height_mAU"].fillna(0.0)
    result["n_peaks_summed"] = result["n_peaks_summed"].fillna(0).astype(int)

    sample_totals = result.groupby("Folder")["Area_Integral"].transform("sum")
    result["Area_Percent"] = np.where(
        sample_totals > 0,
        (result["Area_Integral"] / sample_totals * 100).round(2),
        0.0,
    )

    n_zeros = (result["Area_Integral"] == 0).sum()
    n_summed = (result["n_peaks_summed"] > 1).sum()
    print(
        f"\n  Output: {len(result)} entries "
        f"({result['Folder'].nunique()} samples x {len(compounds_list)} compounds)"
    )
    if n_summed > 0:
        print(f"  Entries with summed sub-peaks: {n_summed}")
    print(f"  Zero-filled (not detected)   : {n_zeros}")

    return (
        result.sort_values(["Strain", "Time_Point_h", "Replicate", "Compound"])
        .reset_index(drop=True)
    )
