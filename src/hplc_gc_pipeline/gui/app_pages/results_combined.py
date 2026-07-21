"""Combined HPLC + GC-MS Results page."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import embed_html, result_header

exp = Path(st.session_state.loaded_exp)
dashboard = exp / "combined_analysis" / "analysis_plots.html"

if dashboard.exists():
    result_header(dashboard)
    embed_html(dashboard, height=900)
