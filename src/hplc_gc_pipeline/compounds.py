"""
Compound definitions and RT-window assignment (Stage 2).

Reads ``compounds.csv`` (Compound, RT_low, RT_high [, Notes]) and tags each
peak with the compound whose retention-time window contains it.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_compounds(compounds_csv: Path) -> pd.DataFrame:
    """Load and validate the compound RT-window table."""
    compounds_csv = Path(compounds_csv)
    if not compounds_csv.exists():
        raise FileNotFoundError(
            f"compounds.csv not found: {compounds_csv}\n"
            f"Columns required: Compound, RT_low, RT_high (RT in minutes)"
        )

    df = pd.read_csv(compounds_csv)
    missing = {"Compound", "RT_low", "RT_high"} - set(df.columns)
    if missing:
        raise ValueError(f"compounds.csv missing required columns: {missing}")

    print("\n" + "=" * 70)
    print("COMPOUND DEFINITIONS")
    print("=" * 70)
    for _, r in df.iterrows():
        note = f" [{r['Notes']}]" if "Notes" in df.columns and pd.notna(r.get("Notes")) else ""
        print(f"  {r['Compound']:<20} RT {r['RT_low']:.2f}--{r['RT_high']:.2f} min{note}")
    return df


def assign_compounds(df: pd.DataFrame, compounds_df: pd.DataFrame, exclude_unknown: bool = True) -> pd.DataFrame:
    """Assign a compound name to each peak from its retention time."""
    print("\n" + "=" * 70)
    print("ASSIGNING COMPOUNDS BY RETENTION TIME")
    print("=" * 70)

    def find_compound(rt):
        if pd.isna(rt):
            return "Unknown"
        for _, c in compounds_df.iterrows():
            if c["RT_low"] <= rt <= c["RT_high"]:
                return c["Compound"]
        return "Unknown"

    df = df.copy()
    df["Compound"] = df["Retention_Time_min"].apply(find_compound)
    n_unknown = (df["Compound"] == "Unknown").sum()

    print(f"\n  Assigned: {len(df) - n_unknown} / {len(df)} peaks")
    if n_unknown:
        print(f"  Warning: {n_unknown} peaks outside all RT windows (Unknown)")
    for c in compounds_df["Compound"]:
        print(f"    - {c:<20} {(df['Compound'] == c).sum()} peaks")

    if exclude_unknown:
        df = df[df["Compound"] != "Unknown"].copy()
        if n_unknown:
            print(f"\n  Unknown peaks excluded ({n_unknown} dropped)")
    return df
