"""Results page — view analysis dashboards and download outputs."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import embed_html, output_row

exp = Path(st.session_state.loaded_exp)
analysis_dir = st.session_state.get("w_analysis_dir", "analysis")
results_dir = st.session_state.get("w_results_dir", "results")

dashboard = exp / analysis_dir / "analysis_plots.html"
overlay = exp / results_dir / "chromatogram_overlay_interactive.html"
gallery = exp / results_dir / "chromatogram_gallery.html"
gc_dashboard = exp / "gc_analysis" / "analysis_plots.html"
gc_overlay = exp / "gc_results" / "gc_tic_overlay_interactive.html"
gc_gallery = exp / "gc_results" / "chromatogram_gallery.html"
combined = exp / "combined_analysis" / "analysis_plots.html"

hplc_files = (dashboard, overlay, gallery)
gc_files = (gc_dashboard, gc_overlay, gc_gallery)

if not any(p.exists() for p in (*hplc_files, *gc_files, combined)):
    st.space("large")
    with st.container(horizontal_alignment="center"):
        st.markdown(":material/analytics:")
        st.subheader("No results yet", text_alignment="center")
        st.caption(
            "Go to **Run** and click **Run both** to process your injections "
            "and generate calibrated concentration plots.",
            text_alignment="center",
        )
    st.space("large")
else:
    if combined.exists():
        st.subheader("Combined HPLC + GC dashboard")
        output_row("Combined dashboard", combined)
        embed_html(combined, height=900)

    if any(p.exists() for p in hplc_files):
        if combined.exists():
            st.divider()
        st.subheader("HPLC — analysis dashboard")
        if dashboard.exists():
            output_row("HPLC dashboard", dashboard)
            embed_html(dashboard, height=900)
        if overlay.exists():
            st.markdown("**Chromatogram overlay (interactive)**")
            output_row("HPLC interactive overlay", overlay)
            embed_html(overlay, height=650)
        if gallery.exists():
            st.markdown("**Per-injection chromatogram gallery**")
            output_row("HPLC chromatogram gallery", gallery)
            with st.expander("Preview gallery inline"):
                embed_html(gallery, height=820)

    if any(p.exists() for p in gc_files):
        st.divider()
        st.subheader("GC-MS — analysis dashboard")
        if gc_dashboard.exists():
            output_row("GC dashboard", gc_dashboard)
            embed_html(gc_dashboard, height=900)
        if gc_overlay.exists():
            st.markdown("**TIC overlay (interactive)**")
            output_row("GC TIC overlay", gc_overlay)
            embed_html(gc_overlay, height=650)
        if gc_gallery.exists():
            st.markdown("**Per-injection EIC gallery**")
            output_row("GC chromatogram gallery", gc_gallery)
            with st.expander("Preview GC gallery inline"):
                embed_html(gc_gallery, height=820)
