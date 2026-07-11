"""
Stage 2 Bokeh dashboard + CSV export.

Eight plot types (summary, raw replicates, combined, bar charts) in two units
(peak area and concentration uM), written to a single analysis_plots.html.
Aesthetic constants live here as documented module defaults; only operational
settings (palette, strain order/filter, bar time points, grid columns) come
from config.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import bokeh.plotting
import bokeh.io
from bokeh.layouts import column, gridplot
from bokeh.models import ColumnDataSource

# --- Aesthetic defaults (see config decision: styling stays in code) --------
DEFAULT_PALETTE = [
    "#377eb8", "#e41a1c", "#984ea3", "#4daf4a",
    "#ff7f00", "#000000", "#a65628", "#f781bf",
]
DASH_STYLES = {"A": "solid", "B": "dashed", "C": "dotted", "D": "dashdot"}
PUB_STYLE = dict(background_fill_color="white", outline_line_color="black", outline_line_width=1.5)
AXIS_FONT = "Times New Roman"
AXIS_LABEL_SIZE = "24pt"
TICK_LABEL_SIZE = "14pt"
LINE_WIDTH = 2
MARKER_SIZE = 6
BAND_ALPHA = 0.07
PLOT_WIDTH = 900
PLOT_HEIGHT = 750
BAR_ALPHA = 0.8
ERR_BAR_WIDTH = 0.03
ERR_BAR_COLOR = "black"
ERR_BAR_LW = 1.5


# ---------------------------------------------------------------------------
# Styling helpers
# ---------------------------------------------------------------------------

def _style_axes(p, x_label, y_label, grid=False):
    for axis in [p.xaxis, p.yaxis]:
        axis.axis_label_text_font = AXIS_FONT
        axis.axis_label_text_font_size = AXIS_LABEL_SIZE
        axis.major_label_text_font = AXIS_FONT
        axis.major_label_text_font_size = TICK_LABEL_SIZE
    p.xaxis.axis_label = x_label
    p.yaxis.axis_label = y_label
    p.xgrid.visible = False
    p.ygrid.visible = grid
    if grid:
        p.ygrid.grid_line_alpha = 0.3


def _style_legend(p):
    p.legend.location = "top_left"
    p.legend.label_text_font = AXIS_FONT
    p.legend.label_text_font_size = TICK_LABEL_SIZE
    p.legend.background_fill_alpha = 0.0
    p.legend.border_line_color = None
    p.legend.click_policy = "hide"


def _style_legend_replicate(p, n_entries):
    if not p.legend:
        return
    legend = p.legend[0]
    legend.label_text_font = AXIS_FONT
    legend.label_text_font_size = TICK_LABEL_SIZE
    legend.background_fill_alpha = 0.9
    legend.border_line_color = "black"
    legend.border_line_width = 1
    legend.click_policy = "hide"
    legend.ncols = max(1, (n_entries - 1) // 12 + 1)
    if legend in p.center:
        p.center.remove(legend)
    p.add_layout(legend, "right")


def _add_error_bars(p, x, y_mean, y_std):
    for xc, ym, ye in zip(x, y_mean, y_std):
        if ye <= 0 or np.isnan(ye):
            continue
        p.segment(x0=xc, y0=max(0, ym - ye), x1=xc, y1=ym + ye, color=ERR_BAR_COLOR, line_width=ERR_BAR_LW)
        p.segment(x0=xc - ERR_BAR_WIDTH, y0=ym + ye, x1=xc + ERR_BAR_WIDTH, y1=ym + ye, color=ERR_BAR_COLOR, line_width=ERR_BAR_LW)
        bottom = max(0, ym - ye)
        if bottom > 0:
            p.segment(x0=xc - ERR_BAR_WIDTH, y0=bottom, x1=xc + ERR_BAR_WIDTH, y1=bottom, color=ERR_BAR_COLOR, line_width=ERR_BAR_LW)


def _make_figure(title, x_label, y_label, grid=False, x_range=None):
    kwargs = dict(title=title, width=PLOT_WIDTH, height=PLOT_HEIGHT, sizing_mode="fixed", **PUB_STYLE)
    if x_range is not None:
        kwargs["x_range"] = x_range
    p = bokeh.plotting.figure(**kwargs)
    _style_axes(p, x_label, y_label, grid=grid)
    return p


def _grid_layout(plots, n_cols=2):
    rows = []
    for i in range(0, len(plots), n_cols):
        chunk = plots[i:i + n_cols]
        while len(chunk) < n_cols:
            chunk.append(None)
        rows.append(chunk)
    return gridplot(rows, sizing_mode="fixed")


def _combos(consolidated_df):
    return (
        consolidated_df[["Strain", "Replicate"]]
        .drop_duplicates()
        .sort_values(["Strain", "Replicate"])
        .values.tolist()
    )


# ---------------------------------------------------------------------------
# Plot builders. Each takes the palette so colouring follows config.
# `value` selects area vs concentration ("Area" or "Conc").
# ---------------------------------------------------------------------------

def _mean_cols(value):
    return ("Mean_Area", "Std_Area") if value == "Area" else ("Mean_Conc_uM", "Std_Conc_uM")


def _raw_col(value):
    return "Area_Integral" if value == "Area" else "Concentration_uM"


def _make_summary_plots(stats_df, compounds, strains, palette, value, y_label, title_suffix):
    mean_c, std_c = _mean_cols(value)
    plots = []
    for compound in compounds:
        p = _make_figure(f"{compound} -- {title_suffix}", "Time (h)", y_label)
        comp_df = stats_df[stats_df["Compound"] == compound]
        for s_idx, strain in enumerate(strains):
            sub = comp_df[comp_df["Strain"] == strain].sort_values("Time_Point_h")
            if sub.empty:
                continue
            clr = palette[s_idx % len(palette)]
            src = ColumnDataSource({
                "x": sub["Time_Point_h"].values,
                "y": sub[mean_c].values,
                "y1": (sub[mean_c] - sub[std_c].fillna(0)).values,
                "y2": (sub[mean_c] + sub[std_c].fillna(0)).values,
            })
            p.varea("x", "y1", "y2", source=src, fill_color=clr, fill_alpha=BAND_ALPHA, legend_label=strain)
            p.line("x", "y", source=src, color=clr, line_width=LINE_WIDTH, legend_label=strain)
            p.scatter("x", "y", source=src, color=clr, size=MARKER_SIZE, marker="circle", legend_label=strain)
        _style_legend(p)
        plots.append(p)
    return plots


def _make_raw_replicate_plots(consolidated_df, compounds, strains, palette, value, y_label, title_suffix):
    raw_c = _raw_col(value)
    combos = _combos(consolidated_df)
    plots = []
    for compound in compounds:
        p = _make_figure(f"{compound} -- {title_suffix}", "Time (h)", y_label)
        comp_df = consolidated_df[consolidated_df["Compound"] == compound]
        for strain, rep in combos:
            sub = comp_df[(comp_df["Strain"] == strain) & (comp_df["Replicate"] == rep)].sort_values("Time_Point_h")
            if sub.empty:
                continue
            s_idx = strains.index(strain) if strain in strains else 0
            clr = palette[s_idx % len(palette)]
            dash = DASH_STYLES.get(rep, "solid")
            src = ColumnDataSource({"x": sub["Time_Point_h"].values, "y": sub[raw_c].values})
            p.line("x", "y", source=src, color=clr, line_width=LINE_WIDTH, line_dash=dash, legend_label=f"{strain}-{rep}")
        _style_legend_replicate(p, len(combos))
        plots.append(p)
    return plots


def _make_combined_plots(consolidated_df, stats_df, compounds, strains, palette, value, y_label, title_suffix, mean_suffix):
    raw_c = _raw_col(value)
    mean_c, std_c = _mean_cols(value)
    combos = _combos(consolidated_df)
    plots = []
    for compound in compounds:
        p = _make_figure(f"{compound} -- {title_suffix}", "Time (h)", y_label)
        comp_raw = consolidated_df[consolidated_df["Compound"] == compound]
        comp_stat = stats_df[stats_df["Compound"] == compound]

        for strain, rep in combos:
            sub = comp_raw[(comp_raw["Strain"] == strain) & (comp_raw["Replicate"] == rep)].sort_values("Time_Point_h")
            if sub.empty:
                continue
            s_idx = strains.index(strain) if strain in strains else 0
            clr = palette[s_idx % len(palette)]
            dash = DASH_STYLES.get(rep, "solid")
            src = ColumnDataSource({"x": sub["Time_Point_h"].values, "y": sub[raw_c].values})
            p.line("x", "y", source=src, color=clr, line_width=LINE_WIDTH - 0.6, line_dash=dash, alpha=0.7, legend_label=f"{strain}-{rep}")

        for s_idx, strain in enumerate(strains):
            sub = comp_stat[comp_stat["Strain"] == strain].sort_values("Time_Point_h")
            if sub.empty:
                continue
            clr = palette[s_idx % len(palette)]
            src = ColumnDataSource({
                "x": sub["Time_Point_h"].values,
                "y": sub[mean_c].values,
                "y1": (sub[mean_c] - sub[std_c].fillna(0)).values,
                "y2": (sub[mean_c] + sub[std_c].fillna(0)).values,
            })
            p.varea("x", "y1", "y2", source=src, fill_color=clr, fill_alpha=BAND_ALPHA, legend_label=f"{strain} {mean_suffix}")
            p.line("x", "y", source=src, color=clr, line_width=LINE_WIDTH + 1, legend_label=f"{strain} {mean_suffix}")
            p.scatter("x", "y", source=src, color=clr, size=MARKER_SIZE + 2, marker="circle", legend_label=f"{strain} {mean_suffix}")

        _style_legend_replicate(p, len(combos) + len(strains))
        plots.append(p)
    return plots


def _make_bar_charts(stats_df, compounds, strains, time_points, palette, value, y_label, title_prefix):
    mean_c, std_c = _mean_cols(value)
    plots = []
    n_strains = len(strains)
    bar_width = 0.75 / max(n_strains, 1)
    for tp in time_points:
        tp_df = stats_df[stats_df["Time_Point_h"] == tp]
        if tp_df.empty:
            continue
        p = _make_figure(f"{title_prefix} -- t = {tp} h", "Compound", y_label, grid=True, x_range=compounds)
        for s_idx, strain in enumerate(strains):
            strain_df = tp_df[tp_df["Strain"] == strain]
            clr = palette[s_idx % len(palette)]
            x_coords, y_vals, y_errs = [], [], []
            for c_idx, compound in enumerate(compounds):
                row = strain_df[strain_df["Compound"] == compound]
                mean_v = float(row[mean_c].values[0]) if len(row) > 0 else 0.0
                std_v = float(row[std_c].fillna(0).values[0]) if len(row) > 0 else 0.0
                x_coords.append(c_idx + (s_idx - (n_strains - 1) / 2) * bar_width)
                y_vals.append(mean_v)
                y_errs.append(std_v)
            p.vbar(x=x_coords, top=y_vals, width=bar_width * 0.9, color=clr, alpha=BAR_ALPHA, legend_label=strain)
            _add_error_bars(p, x_coords, y_vals, y_errs)
        _style_legend(p)
        plots.append(p)
    return plots


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def export_analysis_results(assigned_df, consolidated_df, stats_df, analysis_dir: Path) -> Path:
    """Write assigned/consolidated/summary CSVs and per-strain time series."""
    print("\n" + "=" * 70)
    print("EXPORTING ANALYSIS RESULTS")
    print("=" * 70)

    analysis_dir = Path(analysis_dir)
    analysis_dir.mkdir(exist_ok=True, parents=True)

    p1 = analysis_dir / "assigned_peaks.csv"
    assigned_df.to_csv(p1, index=False)
    print(f"\n  Assigned peaks     : {p1.name} ({len(assigned_df)} rows)")

    p2 = analysis_dir / "consolidated_peaks.csv"
    consolidated_df.to_csv(p2, index=False)
    print(f"  Consolidated peaks : {p2.name} ({len(consolidated_df)} rows) <- plot input")

    p3 = analysis_dir / "replicate_summary.csv"
    stats_df.to_csv(p3, index=False)
    print(f"  Replicate summary  : {p3.name} ({len(stats_df)} rows)")

    for strain in sorted(stats_df["Strain"].unique()):
        pivot = (
            stats_df[stats_df["Strain"] == strain]
            .pivot_table(index="Time_Point_h", columns="Compound", values="Mean_Area", aggfunc="first")
        )
        pivot.index.name = "Time_Point_h"
        pivot.columns.name = None
        ts_path = analysis_dir / f"time_series_{strain.replace(' ', '_')}.csv"
        pivot.to_csv(ts_path)
        print(f"  Time series ({strain:<12}): {ts_path.name}")

    return analysis_dir


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def generate_plots(consolidated_df, stats_df, analysis_dir: Path,
                   palette=None, strain_order=None, plot_strains=None, n_cols=2) -> Path:
    """Build the area + concentration plot sets and save a single HTML."""
    print("\n" + "=" * 70)
    print("GENERATING BOKEH PLOTS")
    print("=" * 70)

    palette = palette or DEFAULT_PALETTE
    analysis_dir = Path(analysis_dir)
    analysis_dir.mkdir(exist_ok=True, parents=True)

    compounds = sorted(stats_df["Compound"].unique().tolist())
    time_points = sorted(consolidated_df["Time_Point_h"].unique().tolist())

    known = set(stats_df["Strain"].unique())
    if strain_order is None:
        all_strains = sorted(known)
    else:
        all_strains = [s for s in strain_order if s in known]
        missing = [s for s in strain_order if s not in known]
        if missing:
            print(f"  Warning: strains in strain_order not found in data: {missing}")

    if plot_strains is not None:
        not_found = [s for s in plot_strains if s not in all_strains]
        if not_found:
            print(f"  Warning: plot_strains not found in data: {not_found}")
        strains = [s for s in all_strains if s in plot_strains]
        plot_consolidated = consolidated_df[consolidated_df["Strain"].isin(strains)].copy()
        plot_stats = stats_df[stats_df["Strain"].isin(strains)].copy()
        print(f"  Plot filter: {strains} ({len(all_strains) - len(strains)} strains hidden)")
    else:
        strains = all_strains
        plot_consolidated = consolidated_df
        plot_stats = stats_df

    print(f"  Compounds  : {compounds}")
    print(f"  Strains    : {strains}")
    print(f"  Time points: {time_points} h")

    html_path = analysis_dir / "analysis_plots.html"
    bokeh.io.output_file(str(html_path))

    # Area-based
    summary_area = _make_summary_plots(plot_stats, compounds, strains, palette, "Area", "Peak Area", "Mean +/- SD")
    raw_area = _make_raw_replicate_plots(plot_consolidated, compounds, strains, palette, "Area", "Peak Area", "Raw Replicates")
    combined_area = _make_combined_plots(plot_consolidated, plot_stats, compounds, strains, palette, "Area", "Peak Area", "Combined", "(mean)")
    bars_area = _make_bar_charts(plot_stats, compounds, strains, time_points, palette, "Area", "Mean Peak Area", "Peak Areas")

    # Concentration-based
    summary_conc = _make_summary_plots(plot_stats, compounds, strains, palette, "Conc", "Concentration (uM)", "Mean +/- SD (uM)")
    raw_conc = _make_raw_replicate_plots(plot_consolidated, compounds, strains, palette, "Conc", "Concentration (uM)", "Raw Replicates (uM)")
    combined_conc = _make_combined_plots(plot_consolidated, plot_stats, compounds, strains, palette, "Conc", "Concentration (uM)", "Combined (uM)", "(mean uM)")
    bars_conc = _make_bar_charts(plot_stats, compounds, strains, time_points, palette, "Conc", "Mean Concentration (uM)", "Concentration")

    grid_summary = _grid_layout(summary_area + summary_conc, n_cols)
    grid_raw = _grid_layout(raw_area + raw_conc, n_cols)
    grid_combined = _grid_layout(combined_area + combined_conc, n_cols)
    grid_bars = _grid_layout(bars_area + bars_conc, n_cols)

    layout = column(grid_summary, grid_raw, grid_combined, grid_bars, sizing_mode="fixed")
    bokeh.io.save(layout)

    print(f"\n  Saved: {html_path}")
    print(f"   area: {len(summary_area)} summary, {len(raw_area)} raw, {len(combined_area)} combined, {len(bars_area)} bars")
    print(f"   uM  : {len(summary_conc)} summary, {len(raw_conc)} raw, {len(combined_conc)} combined, {len(bars_conc)} bars")
    return html_path
