"""
Sample-name parsing and replicate pooling (Stage 2).

Sample names are encoded as ``"<time_h>, <strain>, <replicate>"`` (e.g.
``"24, s1, A"``). The regex is configurable so a different naming convention can
be plugged in without touching code. Injections that don't match (blanks,
standards) are dropped with a warning.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


def parse_sample_name(sample_name: str, pattern: str) -> dict:
    """Parse one sample name into time/strain/replicate fields.

    ``Time_Point_h`` is None when the name doesn't match; the caller drops
    those rows.
    """
    match = re.match(pattern, str(sample_name).strip(), re.IGNORECASE)
    if match:
        strain = match.group("strain").strip()
        time = int(match.group("time"))
        rep = match.group("replicate").upper()
        return {
            "Time_Point_h": time,
            "Strain": strain,
            "Replicate": rep,
            "Replicate_Group": f"{strain}_t{time}h",
        }
    print(f"  Warning: could not parse sample name '{sample_name}' -- row skipped")
    return {"Time_Point_h": None, "Strain": sample_name, "Replicate": None, "Replicate_Group": None}


def load_peak_data(peak_results_csv: Path, pattern: str) -> pd.DataFrame:
    """Load peak_results.csv, parse sample names, drop unparseable rows."""
    print("=" * 70)
    print("LOADING PEAK DATA")
    print("=" * 70)

    df = pd.read_csv(peak_results_csv)
    print(f"  Loaded {len(df)} peaks from {Path(peak_results_csv).name}")

    parsed = df["Sample"].apply(lambda s: parse_sample_name(s, pattern))
    for key in ["Strain", "Time_Point_h", "Replicate", "Replicate_Group"]:
        df[key] = parsed.apply(lambda x, k=key: x[k])

    n_before = len(df)
    df_valid = df[df["Time_Point_h"].notna()].copy()
    df_valid["Time_Point_h"] = df_valid["Time_Point_h"].astype(int)
    n_dropped = n_before - len(df_valid)

    if n_dropped > 0:
        dropped = df[df["Time_Point_h"].isna()]["Sample"].unique()
        print(f"\n  Skipped {n_dropped} rows from {len(dropped)} unparseable sample(s):")
        for s in dropped:
            print(f"    - '{s}' ({(df['Sample'] == s).sum()} peaks)")
        print("    (likely blank/standard injections -- safe to ignore)\n")

    print(
        f"  Retained: {len(df_valid)} peaks | "
        f"{df_valid['Sample'].nunique()} samples | "
        f"{df_valid['Strain'].nunique()} strains | "
        f"{df_valid['Time_Point_h'].nunique()} time points"
    )
    return df_valid


def pool_replicates_across_experiments(df: pd.DataFrame, group_cols=None) -> pd.DataFrame:
    """Give each injection a unique replicate letter within its group.

    All peaks from one Folder (injection) share a letter. When the same letter
    appears more than once within a (Strain, Time_Point_h) group (e.g. the same
    replicate label reused across two runs), duplicates are renamed to the next
    unused letter.
    """
    df = df.copy()
    letters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    if group_cols is None:
        group_cols = ["Strain", "Time_Point_h"]

    folder_to_rep = {}
    for _, g in df.groupby(group_cols, sort=False):
        sample_info = (
            g[["Folder", "Replicate"]]
            .drop_duplicates(subset=["Folder"])
            .sort_values("Replicate")
        )
        used = set()
        counts = {}
        for _, row in sample_info.iterrows():
            folder = row["Folder"]
            rep = str(row["Replicate"])
            counts[rep] = counts.get(rep, 0) + 1
            if counts[rep] == 1 and rep not in used:
                new_rep = rep
            else:
                new_rep = next((c for c in letters if c not in used), None)
                if new_rep is None:
                    raise ValueError("Ran out of replicate letters A-Z")
            used.add(new_rep)
            folder_to_rep[folder] = new_rep

    df["Replicate"] = df["Folder"].map(folder_to_rep).fillna(df["Replicate"])
    df["Replicate_Group"] = df.apply(
        lambda r: f"{r['Strain']}_t{int(r['Time_Point_h'])}h", axis=1
    )
    return df
