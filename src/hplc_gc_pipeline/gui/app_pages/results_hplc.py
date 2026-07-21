"""HPLC Results page."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import embed_html, result_header

exp          = Path(st.session_state.loaded_exp)
analysis_dir = st.session_state.get("w_analysis_dir", "analysis")
results_dir  = st.session_state.get("w_results_dir",  "results")

dashboard = exp / analysis_dir / "analysis_plots.html"
overlay   = exp / results_dir  / "chromatogram_overlay_interactive.html"
gallery   = exp / results_dir  / "chromatogram_gallery.html"

if dashboard.exists():
    result_header(dashboard)
    embed_html(dashboard, height=900)

if overlay.exists():
    st.divider()
    st.subheader("Chromatogram overlay")
    result_header(overlay)
    embed_html(overlay, height=700)

if gallery.exists():
    st.divider()
    st.subheader("Per-injection chromatograms")
    result_header(gallery)
    embed_html(gallery, height=900, scrolling=True)
