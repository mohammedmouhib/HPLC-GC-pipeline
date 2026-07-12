"""
GC-MS Stage 1 visualisation.

Per-injection PNGs (the quantifier-ion EIC for each compound with its integrated
area shaded), a reused HTML gallery, a static overlay, and an interactive Bokeh
TIC overlay. Integration is recomputed from the loaded traces via the same
functions the peak table uses, so the plots always match the reported areas.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .gc_agilent import GCInjection, load_gc_traces
from .gc_processing import GCCompound, GCProcessingConfig, als_baseline, integrate_peak, _channel_signal
from .plotting_chromatograms import export_html_gallery


# ---------------------------------------------------------------------------
# Per-injection plot
# ---------------------------------------------------------------------------

def plot_gc_chromatogram(inj: GCInjection, compounds: list[GCCompound],
                         cfg: GCProcessingConfig, save_path: Path) -> None:
    """Render one GC injection: the quantification trace per compound with the
    integrated peak shaded and annotated."""
    traces = load_gc_traces(inj)
    fig, ax = plt.subplots(figsize=(14, 6))
    colors = plt.cm.tab10.colors

    # Faint TIC for context. The solvent peak can dwarf the compound ions, so
    # the y-axis is scaled to the compound EICs below rather than the TIC.
    ax.plot(traces.time_min, traces.tic, color="0.85", lw=0.7, label="TIC", zorder=1)

    compound_ymax = 0.0
    for i, compound in enumerate(compounds):
        color = colors[i % len(colors)]
        try:
            t, raw = _channel_signal(traces, compound, cfg.quant_channel)
        except Exception:
            continue
        corrected = raw - als_baseline(raw, cfg.baseline.lam, cfg.baseline.p, cfg.baseline.niter)
        chan = (f"m/z {int(compound.quantifier_mz)}" if cfg.quant_channel == "eic"
                else cfg.quant_channel.upper())
        ax.plot(t, corrected, color=color, lw=1.0, alpha=0.9,
                label=f"{compound.name} ({chan})", zorder=3)

        res = integrate_peak(t, corrected, compound.rt_low, compound.rt_high, cfg.peak_detection)
        if res.detected and res.left_idx is not None:
            lo, hi = res.left_idx, res.right_idx + 1
            ax.fill_between(t[lo:hi], 0, corrected[lo:hi], color=color, alpha=0.35, zorder=4)
            ax.axvline(res.rt_min, color=color, lw=0.6, ls="--", alpha=0.6, zorder=2)
            ax.annotate(f"{compound.name}\n{res.rt_min:.2f} min\narea {res.area:.0f}",
                        xy=(res.rt_min, res.height), xytext=(0, 8), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8, color=color,
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor=color, alpha=0.85),
                        zorder=6)
        # Track the compound-signal scale within each RT window (excludes the
        # solvent front, which lives outside every window).
        win = (t >= compound.rt_low) & (t <= compound.rt_high)
        if win.any():
            compound_ymax = max(compound_ymax, float(np.max(corrected[win])))

    if compound_ymax > 0:
        ax.set_ylim(bottom=-0.08 * compound_ymax, top=compound_ymax * 1.25)

    ax.set_xlabel("Retention Time (min)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Intensity (ion counts)", fontsize=11, fontweight="bold")
    ax.set_title(f"{inj.sample_name}  ({inj.folder_name})  |  GC-MS {cfg.quant_channel.upper()}",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3, ls="--")
    plt.tight_layout()
    plt.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Static overlay
# ---------------------------------------------------------------------------

def plot_gc_overlay(injections: list[GCInjection], output_dir: Path) -> Path:
    """Static matplotlib overlay of every injection's TIC."""
    fig, ax = plt.subplots(figsize=(16, 8))
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(injections), 1)))
    for i, inj in enumerate(injections):
        traces = load_gc_traces(inj)
        label = f"{inj.folder_name}: {inj.sample_name}"
        ax.plot(traces.time_min, traces.tic, lw=1, alpha=0.7, color=colors[i],
                label=label[:40] + ("..." if len(label) > 40 else ""))
    ax.set_xlabel("Retention Time (min)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Total Ion Current", fontsize=12, fontweight="bold")
    ax.set_title("GC-MS TIC Overlay", fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    path = Path(output_dir) / "gc_tic_overlay.png"
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  GC TIC overlay: {path.name}")
    return path


# ---------------------------------------------------------------------------
# Interactive Bokeh TIC overlay
# ---------------------------------------------------------------------------

def plot_gc_overlay_interactive(injections: list[GCInjection], output_dir: Path) -> Path:
    """Interactive Bokeh overlay of TICs with a per-injection checkbox panel."""
    from bokeh.plotting import figure
    from bokeh.models import ColumnDataSource, CheckboxGroup, CustomJS, Button, HoverTool
    from bokeh.layouts import row, column
    from bokeh.palettes import Category20_20
    import bokeh.io

    n = len(injections)
    if n <= 20:
        colors = list(Category20_20[:max(n, 1)])
    else:
        import matplotlib.cm as cm
        colors = ["#%02x%02x%02x" % tuple(int(c * 255) for c in cm.tab20(i / n)[:3]) for i in range(n)]

    p = figure(title="GC-MS TIC Overlay", x_axis_label="Retention Time (min)",
               y_axis_label="Total Ion Current", width=1050, height=620, sizing_mode="fixed",
               background_fill_color="white", outline_line_color="black", outline_line_width=1.5)
    p.xgrid.visible = False
    p.ygrid.grid_line_alpha = 0.3
    p.add_tools(HoverTool(tooltips=[("Sample", "@sample"), ("RT", "@x{0.00} min"), ("TIC", "@y{0.0}")]))

    renderers, labels = [], []
    for i, inj in enumerate(injections):
        traces = load_gc_traces(inj)
        src = ColumnDataSource({"x": traces.time_min, "y": traces.tic,
                                "sample": [inj.sample_name] * len(traces.time_min)})
        renderers.append(p.line("x", "y", source=src, color=colors[i], line_width=2, alpha=0.85, visible=False))
        labels.append(f"{inj.folder_name}  |  {inj.sample_name}")

    checkbox = CheckboxGroup(labels=labels, active=[], width=320)
    args = dict(renderers=renderers, cb=checkbox)
    checkbox.js_on_change("active", CustomJS(args=args, code="""
        for (let i = 0; i < renderers.length; i++) { renderers[i].visible = cb.active.includes(i); }
    """))
    show_all = Button(label="Show All", button_type="success", width=155)
    hide_all = Button(label="Hide All", button_type="danger", width=155)
    show_all.js_on_click(CustomJS(args=dict(**args, n=n), code="""
        cb.active = Array.from({length: n}, (_, i) => i);
        for (const r of renderers) r.visible = true;"""))
    hide_all.js_on_click(CustomJS(args=args, code="""
        cb.active = []; for (const r of renderers) r.visible = false;"""))

    side = column(row(show_all, hide_all), checkbox, width=340)
    layout = row(p, side)
    path = Path(output_dir) / "gc_tic_overlay_interactive.html"
    bokeh.io.output_file(str(path))
    bokeh.io.save(layout)
    print(f"  GC interactive overlay: {path.name} ({n} traces)")
    return path


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def render_gc_plots(injections: list[GCInjection], compounds: list[GCCompound],
                    cfg: GCProcessingConfig, results_dir: Path) -> None:
    """Per-injection PNGs + HTML gallery + static & interactive TIC overlays."""
    results_dir = Path(results_dir)
    plots_dir = results_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n  Rendering {len(injections)} GC per-injection plots...")
    for inj in injections:
        safe = inj.sample_name.replace(" ", "_").replace("/", "-")
        out = plots_dir / f"{inj.experiment_name}_{inj.folder_name.replace('.D', '')}_{safe}.png"
        try:
            plot_gc_chromatogram(inj, compounds, cfg, out)
        except Exception as exc:
            print(f"    plot failed for {inj.sample_name}: {exc}")

    export_html_gallery(plots_dir, results_dir, title="GC Chromatogram Gallery")
    plot_gc_overlay(injections, results_dir)
    plot_gc_overlay_interactive(injections, results_dir)
