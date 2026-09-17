"""
Stage 1 chromatogram visualisation.

Per-injection PNGs, a self-contained HTML gallery, a static overlay, and an
interactive Bokeh overlay. The original code probed a list of candidate
attribute names to find the time/signal arrays; MOCCA2's interface is known, so
this uses ``chrom.time`` / ``chrom.data`` / ``chrom.baseline`` directly via the
``_time`` and ``_corrected_signal`` helpers.
"""

from __future__ import annotations

import base64
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless: never require a display to save PNGs
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# MOCCA2 array accessors (single source of truth for both matplotlib + Bokeh)
# ---------------------------------------------------------------------------

def _time(chrom) -> np.ndarray:
    return np.asarray(chrom.time).flatten()


def _raw_signal(chrom) -> np.ndarray:
    raw = np.asarray(chrom.data)
    if raw.ndim == 2:
        raw = raw[0]
    elif raw.ndim > 2:
        raw = raw.reshape(raw.shape[0], -1)[0]
    return raw


def _corrected_signal(chrom) -> np.ndarray:
    raw = _raw_signal(chrom)
    if getattr(chrom, "baseline", None) is not None:
        baseline = np.asarray(chrom.baseline)
        if baseline.ndim == 2:
            baseline = baseline[0]
        return raw - baseline
    return raw.copy()


# ---------------------------------------------------------------------------
# Per-injection plot
# ---------------------------------------------------------------------------

def plot_chromatogram(chrom, wavelength: int, save_path: Path) -> None:
    """Render one chromatogram: raw + baseline + corrected + deconvolved comps."""
    fig, ax = plt.subplots(figsize=(14, 6))
    time_data = _time(chrom)
    raw = _raw_signal(chrom)

    baseline = None
    if getattr(chrom, "baseline", None) is not None:
        baseline = np.asarray(chrom.baseline)
        if baseline.ndim == 2:
            baseline = baseline[0]
    corrected = (raw - baseline) if baseline is not None else raw.copy()

    ax.plot(time_data, raw, color="0.75", lw=0.8, label="Raw signal", zorder=1)
    if baseline is not None:
        ax.plot(time_data, baseline, "b--", lw=1.0, label="Baseline", zorder=2)
    ax.plot(time_data, corrected, "k-", lw=1.5, label="Corrected signal", zorder=3)

    tab_colors = plt.cm.tab10.colors
    comp_global = 0
    if getattr(chrom, "peaks", None):
        for peak in chrom.peaks:
            left = max(0, int(peak.left))
            right = min(len(time_data) - 1, int(peak.right))
            x_win = time_data[left:right]
            peak_components = getattr(peak, "components", None) or []

            if not peak_components:
                ax.fill_between(
                    x_win, 0, corrected[left:right], alpha=0.20, color="grey",
                    label="Peak (no deconv)" if comp_global == 0 else "", zorder=4,
                )
                continue

            for comp in peak_components:
                color = tab_colors[comp_global % len(tab_colors)]
                conc = np.asarray(comp.concentration, dtype=float).flatten()
                spec = np.asarray(comp.spectrum, dtype=float).flatten()
                scale = float(spec[0]) if spec.size > 0 else 1.0
                y_comp = conc * scale

                n = min(len(x_win), len(y_comp))
                x_p, y_p = x_win[:n], y_comp[:n]

                area = getattr(comp, "integral", None)
                label = f"Comp {comp_global + 1}"
                if area is not None:
                    label += f"  (integral = {area:.1f})"

                ax.fill_between(x_p, 0, y_p, alpha=0.40, color=color, label=label, zorder=4)
                ax.plot(x_p, y_p, color=color, lw=1.4, zorder=5)

                apex_i = int(np.argmax(y_p))
                rt = float(x_p[apex_i])
                ax.text(
                    rt, float(y_p[apex_i]) * 1.03, f"C{comp_global + 1}\n{rt:.2f} min",
                    ha="center", va="bottom", fontsize=8, color=color,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor=color, alpha=0.85),
                    zorder=6,
                )
                comp_global += 1

    ax.set_xlabel("Retention Time (min)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Intensity (mAU)", fontsize=11, fontweight="bold")
    ax.set_title(
        f"{chrom.sample_name}  ({chrom.folder_name})  @ {wavelength} nm  "
        f"|  {comp_global} deconvolved components",
        fontsize=12, fontweight="bold",
    )
    handles, labels_ = ax.get_legend_handles_labels()
    by_label = dict(zip(labels_, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3, ls="--")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_manual_bounds_preview(
    time_data: np.ndarray,
    signal_raw: np.ndarray,
    segments: list,
    sample_name: str,
    folder_name: str,
    wavelength: int,
    save_path: Path,
) -> None:
    """Chromatogram + shaded manual-integration windows + area/RT annotations.

    signal_raw is expected to be already baseline-corrected (MOCCA corrects in-place).
    """
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(time_data, signal_raw, color="0.35", lw=1.2, label="Corrected signal", zorder=2)

    tab_colors = plt.cm.tab10.colors
    for i, (name, left_rt, right_rt, seg_t, seg_corr, area, apex_rt, height) in enumerate(segments):
        color = tab_colors[i % len(tab_colors)]
        ax.fill_between(seg_t, 0, seg_corr, alpha=0.40, color=color, zorder=3, label=name)
        ax.axvline(left_rt,  color=color, lw=1.2, ls="--", zorder=4)
        ax.axvline(right_rt, color=color, lw=1.2, ls="--", zorder=4)
        ax.text(
            apex_rt, height * 1.04,
            f"{name}\nArea: {area:.0f}\n@ {apex_rt:.3f} min",
            ha="center", va="bottom", fontsize=9, color=color,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor=color, alpha=0.90),
            zorder=5,
        )

    ax.set_xlabel("Retention Time (min)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Intensity (mAU)", fontsize=11, fontweight="bold")
    ax.set_title(
        f"{sample_name}  ({folder_name})  @ {wavelength} nm  —  Manual bounds integration",
        fontsize=12, fontweight="bold",
    )
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3, ls="--")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Static overlay
# ---------------------------------------------------------------------------

def plot_overlay(dataset, output_dir: Path) -> Path:
    """Static matplotlib overlay of all baseline-corrected traces."""
    fig, ax = plt.subplots(figsize=(16, 8))
    n = len(dataset.chromatograms)
    colors = plt.cm.tab20(np.linspace(0, 1, max(n, 1)))

    for i, chrom in enumerate(dataset.chromatograms):
        time_data = _time(chrom)
        signal = _corrected_signal(chrom)
        label = f"{chrom.folder_name}: {chrom.sample_name}"
        if len(label) > 40:
            label = label[:37] + "..."
        ax.plot(time_data, signal, linewidth=1, label=label, color=colors[i], alpha=0.7)

    ax.set_xlabel("Retention Time (min)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Intensity (mAU)", fontsize=12, fontweight="bold")
    ax.set_title(f"Overlay of All Chromatograms ({dataset.wavelength} nm)", fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    overlay_path = Path(output_dir) / "overlay_plot.png"
    plt.savefig(overlay_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Overlay plot: {overlay_path.name}")
    return overlay_path


# ---------------------------------------------------------------------------
# Long-format trace CSV
# ---------------------------------------------------------------------------

def export_chromatogram_traces_csv(dataset, output_dir: Path) -> Path:
    """Export baseline-corrected traces in long format for downstream use."""
    import pandas as pd

    output_dir = Path(output_dir)
    rows, n_ok, n_failed = [], 0, 0
    for chrom in dataset.chromatograms:
        try:
            for t, s in zip(_time(chrom), _corrected_signal(chrom)):
                rows.append({
                    "Folder": chrom.folder_name,
                    "Sample": chrom.sample_name,
                    "Time_min": round(float(t), 5),
                    "Signal_mAU": round(float(s), 4),
                })
            n_ok += 1
        except Exception as exc:
            print(f"  Warning: could not export trace for '{chrom.sample_name}': {exc}")
            n_failed += 1

    df = pd.DataFrame(rows)
    csv_path = output_dir / "chromatogram_traces.csv"
    df.to_csv(csv_path, index=False)
    print(f"  Chromatogram traces: {csv_path.name} ({n_ok} traces, {len(df)} points"
          + (f", {n_failed} failed" if n_failed else "") + ")")
    return csv_path


# ---------------------------------------------------------------------------
# HTML gallery (base64-embedded PNGs)
# ---------------------------------------------------------------------------

def export_html_gallery(plots_dir: Path, output_dir: Path, title: str = "Chromatogram Gallery") -> Path | None:
    """Self-contained HTML viewer with sidebar nav and arrow-key navigation."""
    plots_dir, output_dir = Path(plots_dir), Path(output_dir)
    png_files = sorted(plots_dir.glob("*.png"))
    if not png_files:
        print(f"  Warning: no PNG files in {plots_dir}")
        return None

    images = []
    for png_path in png_files:
        with open(png_path, "rb") as f:
            images.append({"name": png_path.stem, "data": base64.b64encode(f.read()).decode("utf-8")})
    n = len(images)

    nav_links = "\n        ".join(
        f'<a href="#" onclick="goTo({i}); return false;" id="nav_{i}">{i+1}. {img["name"][:40]}</a>'
        for i, img in enumerate(images)
    )
    image_divs = "\n        ".join(
        f'<div class="slide" id="slide_{i}" style="display:none;">'
        f'<img src="data:image/png;base64,{img["data"]}" alt="{img["name"]}"></div>'
        for i, img in enumerate(images)
    )

    html = _GALLERY_TEMPLATE.format(title=title, n=n, nav_links=nav_links, image_divs=image_divs)
    gallery_path = output_dir / "chromatogram_gallery.html"
    with open(gallery_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"  HTML gallery: {gallery_path.name} ({n} plots)")
    return gallery_path


# ---------------------------------------------------------------------------
# Interactive Bokeh overlay
# ---------------------------------------------------------------------------

def plot_overlay_interactive_bokeh(dataset, output_dir: Path, palette: list[str] | None = None) -> Path:
    """Interactive Bokeh overlay: scrollable checkbox panel, dynamic legend."""
    from bokeh.plotting import figure
    from bokeh.models import (ColumnDataSource, CheckboxGroup, CustomJS,
                              Button, Div, HoverTool, InlineStyleSheet)
    from bokeh.layouts import row, column
    from bokeh.palettes import Category20_20
    import bokeh.io

    n_chroms = len(dataset.chromatograms)
    if n_chroms <= 20:
        colors = list(Category20_20[:max(n_chroms, 1)])
    else:
        import matplotlib.cm as cm
        colors = ["#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))
                  for r, g, b, _ in [cm.tab20(i / n_chroms) for i in range(n_chroms)]]

    pub_style = dict(background_fill_color="white", outline_line_color="black", outline_line_width=1.5)
    p = figure(
        title=f"Chromatogram Overlay @ {dataset.wavelength} nm",
        x_axis_label="Retention Time (min)", y_axis_label="Intensity (mAU)",
        width=1050, height=620, sizing_mode="fixed", **pub_style,
    )
    p.xgrid.visible = False
    p.ygrid.grid_line_alpha = 0.3
    p.xaxis.axis_label_text_font_size = "13pt"
    p.yaxis.axis_label_text_font_size = "13pt"
    p.title.text_font_size = "14pt"
    p.add_tools(HoverTool(tooltips=[
        ("Sample", "@sample"), ("Folder", "@folder"),
        ("RT", "@x{0.000} min"), ("Intensity", "@y{0.0} mAU"),
    ]))

    renderers, cb_labels = [], []
    for i, chrom in enumerate(dataset.chromatograms):
        time_data = _time(chrom)
        signal = _corrected_signal(chrom)
        source = ColumnDataSource({
            "x": time_data, "y": signal,
            "sample": [chrom.sample_name] * len(time_data),
            "folder": [chrom.folder_name] * len(time_data),
        })
        renderers.append(p.line("x", "y", source=source, color=colors[i], line_width=2, alpha=0.85, visible=False))
        cb_labels.append(f"{chrom.folder_name}  |  {chrom.sample_name}")

    scroll_css = InlineStyleSheet(css="""
        :host { overflow-y:auto; max-height:420px; display:block;
                border:1px solid #d0d0d0; border-radius:4px; padding:6px 4px; background:#fafafa; }
    """)
    checkbox = CheckboxGroup(labels=cb_labels, active=[], width=300, stylesheets=[scroll_css])
    legend_div = Div(
        text="<div style='padding:8px; color:#999; font-size:11px; font-style:italic;'>No samples selected</div>",
        width=300,
        stylesheets=[InlineStyleSheet(css="""
            :host { overflow-y:auto; max-height:160px; display:block;
                    border:1px solid #b0c4de; border-radius:4px; background:#f0f4f8; }
        """)],
    )

    legend_update_js = """
    (function updateLegend() {
        if (cb.active.length === 0) {
            legend_div.text = "<div style='padding:8px; color:#999; font-size:11px; font-style:italic;'>No samples selected</div>";
            return;
        }
        let html = "<div style='padding:4px 0;'>";
        for (const i of cb.active) {
            html += "<div style='display:flex; align-items:center; margin:3px 4px; gap:6px;'>"
                  + "<span style='display:inline-block; width:14px; height:14px; flex-shrink:0; border-radius:2px; background:" + colors[i] + ";'></span>"
                  + "<span style='font-size:11px; font-family:Arial; word-break:break-all;'>" + labels[i] + "</span></div>";
        }
        legend_div.text = html + "</div>";
    })();
    """
    cb_args = dict(renderers=renderers, cb=checkbox, legend_div=legend_div, colors=colors, labels=cb_labels)

    checkbox.js_on_change("active", CustomJS(args=cb_args, code="""
        for (let i = 0; i < renderers.length; i++) { renderers[i].visible = cb.active.includes(i); }
    """ + legend_update_js))

    btn_show = Button(label="Show All", button_type="success", width=145)
    btn_hide = Button(label="Hide All", button_type="danger", width=145)
    btn_show.js_on_click(CustomJS(args=dict(**cb_args, n=len(renderers)), code="""
        cb.active = Array.from({length: n}, (_, i) => i);
        for (const r of renderers) r.visible = true;
    """ + legend_update_js))
    btn_hide.js_on_click(CustomJS(args=cb_args, code="""
        cb.active = [];
        for (const r of renderers) r.visible = false;
    """ + legend_update_js))

    header = Div(text=f"<b style='font-size:13px; font-family:Arial;'>All samples ({len(renderers)}) — scroll &amp; click:</b>", width=300)
    legend_header = Div(text="<b style='font-size:13px; font-family:Arial;'>Selected curves:</b>", width=300)
    hint = Div(text="<span style='font-size:11px; color:#888; font-family:Arial;'>Hover plot for RT &amp; intensity values.</span>", width=300)

    # No sizing_mode="fixed" on the row/column: their children (the figure and
    # the fixed-width panel) already carry sizes, and forcing "fixed" without
    # explicit width/height triggers a harmless Bokeh W-1005 warning.
    side_panel = column(header, row(btn_show, btn_hide), checkbox, legend_header, legend_div, hint, width=320)
    layout = row(p, side_panel)

    html_path = Path(output_dir) / "chromatogram_overlay_interactive.html"
    bokeh.io.output_file(str(html_path))
    bokeh.io.save(layout)
    print(f"  Interactive overlay: {html_path.name} ({len(renderers)} traces)")
    return html_path


_GALLERY_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>{title}</title>
  <style>
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ font-family: Arial, sans-serif; background: #e0e0e0; display: flex; height: 100vh; overflow: hidden; }}
    #sidebar {{ width: 240px; min-width: 240px; background: #2c3e50; color: #fff; display: flex; flex-direction: column; height: 100vh; }}
    #sidebar h2 {{ padding: 14px 12px; font-size: 13px; background: #1a252f; border-bottom: 1px solid #3d5166; flex-shrink: 0; }}
    #sidebar .subtitle {{ padding: 7px 12px; font-size: 11px; color: #aaa; background: #1a252f; border-bottom: 1px solid #3d5166; flex-shrink: 0; }}
    #nav-list {{ flex: 1; overflow-y: auto; padding: 4px 0; }}
    #nav-list a {{ display: block; padding: 5px 12px; font-size: 11px; color: #ccc; text-decoration: none; border-left: 3px solid transparent; word-break: break-all; }}
    #nav-list a:hover {{ background: #3d5166; color: #fff; }}
    #nav-list a.active {{ background: #2980b9; color: #fff; border-left-color: #5dade2; }}
    #main {{ flex: 1; display: flex; flex-direction: column; overflow: hidden; }}
    #toolbar {{ background: #fff; padding: 10px 18px; border-bottom: 1px solid #ddd; display: flex; align-items: center; gap: 12px; flex-shrink: 0; }}
    #toolbar h1 {{ font-size: 14px; color: #333; flex: 1; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }}
    #counter {{ font-size: 13px; color: #666; white-space: nowrap; }}
    button {{ padding: 6px 14px; border: 1px solid #aaa; border-radius: 4px; background: #f0f0f0; cursor: pointer; font-size: 12px; }}
    button:hover {{ background: #ddd; }}
    button:disabled {{ opacity: 0.35; cursor: default; }}
    #viewer {{ flex: 1; overflow-y: auto; padding: 20px; display: flex; justify-content: center; align-items: flex-start; }}
    .slide img {{ width: 100%; max-width: 1500px; height: auto; border: 1px solid #ccc; border-radius: 4px; box-shadow: 2px 2px 10px rgba(0,0,0,0.18); background: #fff; display: block; }}
  </style>
</head>
<body>
<div id="sidebar">
  <h2>{title}</h2>
  <div class="subtitle">{n} chromatograms — click to jump</div>
  <div id="nav-list">
        {nav_links}
  </div>
</div>
<div id="main">
  <div id="toolbar">
    <h1 id="sample-name"></h1>
    <span id="counter">— / {n}</span>
    <button id="btn-prev" onclick="prev()">&#9664; Prev</button>
    <button id="btn-next" onclick="next()">Next &#9654;</button>
  </div>
  <div id="viewer">
        {image_divs}
  </div>
</div>
<script>
  const N = {n};
  let current = 0;
  function goTo(i) {{
    if (i < 0 || i >= N) return;
    document.getElementById('slide_' + current).style.display = 'none';
    document.getElementById('nav_' + current).classList.remove('active');
    current = i;
    document.getElementById('slide_' + current).style.display = 'block';
    const navEl = document.getElementById('nav_' + current);
    navEl.classList.add('active');
    navEl.scrollIntoView({{ block: 'nearest' }});
    document.getElementById('counter').textContent = (current + 1) + ' / ' + N;
    document.getElementById('sample-name').textContent = navEl.textContent.trim();
    document.getElementById('btn-prev').disabled = (current === 0);
    document.getElementById('btn-next').disabled = (current === N - 1);
    document.getElementById('viewer').scrollTop = 0;
  }}
  function prev() {{ goTo(current - 1); }}
  function next() {{ goTo(current + 1); }}
  document.addEventListener('keydown', function(e) {{
    if (e.key === 'ArrowRight' || e.key === 'ArrowDown') next();
    if (e.key === 'ArrowLeft' || e.key === 'ArrowUp') prev();
  }});
  goTo(0);
</script>
</body>
</html>"""
