"""HPLC Results page — chromatogram gallery · overlay · analysis dashboard."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import _bokeh_content_height, embed_html, result_header

exp          = Path(st.session_state.loaded_exp)
analysis_dir = st.session_state.get("w_analysis_dir", "analysis")
results_dir  = st.session_state.get("w_results_dir",  "results")
n_cols       = st.session_state.get("w_ncols", 3)

overlay   = exp / results_dir  / "chromatogram_overlay_interactive.html"
dashboard = exp / analysis_dir / "analysis_plots.html"
gallery   = exp / results_dir  / "chromatogram_gallery.html"

# ── 1. Per-injection chromatogram gallery ─────────────────────────────────
if gallery.exists():
    st.subheader("Per-injection chromatograms")
    result_header(gallery)
    embed_html(gallery, height=900, scrolling=True)

# ── 2. Chromatogram overlay ───────────────────────────────────────────────
if overlay.exists():
    if gallery.exists():
        st.divider()
    st.subheader("Chromatogram overlay")
    result_header(overlay)
    embed_html(overlay, height=700)

# ── 3. Analysis dashboard (full content height, no inner scroll) ──────────
if dashboard.exists():
    if overlay.exists() or gallery.exists():
        st.divider()
    st.subheader("Analysis")
    result_header(dashboard)
    embed_html(dashboard,
               height=_bokeh_content_height(dashboard, n_cols=n_cols),
               scrolling=False)
