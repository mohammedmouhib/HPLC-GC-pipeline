"""
Replicate statistics (Stage 2).

Mean / SD / CV% per (strain x time point x compound) for both peak area and
concentration. Groups with CV above a configurable threshold are flagged.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_replicate_stats(df: pd.DataFrame, cv_warning_threshold: float = 15.0) -> pd.DataFrame:
    """Aggregate replicates into mean/SD/CV per (strain, time, compound)."""
    print("\n" + "=" * 70)
    print("COMPUTING REPLICATE STATISTICS")
    print("=" * 70)

    stats = []
    for (strain, time, compound), group in df.groupby(["Strain", "Time_Point_h", "Compound"]):
        areas = group["Area_Integral"].dropna()
        heights = group["Height_mAU"].dropna()
        area_pcts = group["Area_Percent"].dropna()
        concs = group["Concentration_uM"].dropna() if "Concentration_uM" in group.columns else pd.Series([], dtype=float)
        reps = sorted(group["Replicate"].dropna().astype(str).tolist())

        n = len(areas)
        mean_a = areas.mean() if n > 0 else np.nan
        std_a = areas.std() if n > 1 else 0.0
        cv = (std_a / mean_a * 100) if (n > 1 and mean_a > 0) else np.nan
        mean_c = concs.mean() if len(concs) > 0 else np.nan
        std_c = concs.std() if len(concs) > 1 else 0.0

        stats.append({
            "Strain": strain,
            "Time_Point_h": int(time),
            "Replicate_Group": f"{strain}_t{int(time)}h",
            "Compound": compound,
            "n_replicates": n,
            "Replicates": ",".join(reps),
            "Mean_Area": round(mean_a, 2) if pd.notna(mean_a) else np.nan,
            "Std_Area": round(std_a, 2),
            "CV_Percent": round(cv, 1) if pd.notna(cv) else np.nan,
            "Mean_Height_mAU": round(heights.mean(), 2) if len(heights) > 0 else np.nan,
            "Std_Height_mAU": round(heights.std(), 2) if len(heights) > 1 else 0.0,
            "Mean_Area_Percent": round(area_pcts.mean(), 2) if len(area_pcts) > 0 else np.nan,
            "Std_Area_Percent": round(area_pcts.std(), 2) if len(area_pcts) > 1 else 0.0,
            "Mean_Conc_uM": round(mean_c, 3) if pd.notna(mean_c) else np.nan,
            "Std_Conc_uM": round(std_c, 3),
        })

    stats_df = pd.DataFrame(stats).sort_values(["Compound", "Strain", "Time_Point_h"]).reset_index(drop=True)

    print(
        f"\n  {len(stats_df)} groups "
        f"({stats_df['Compound'].nunique()} compounds x "
        f"{stats_df['Strain'].nunique()} strains x "
        f"{stats_df['Time_Point_h'].nunique()} time points)"
    )

    high_cv = stats_df[stats_df["CV_Percent"] > cv_warning_threshold]
    if len(high_cv):
        print(f"\n  Warning: CV > {cv_warning_threshold}% in {len(high_cv)} group(s):")
        for _, r in high_cv.iterrows():
            print(f"    - {r['Replicate_Group']} | {r['Compound']}: "
                  f"CV = {r['CV_Percent']:.1f}% (n={r['n_replicates']})")
    return stats_df
