"""HPLC Results page — chromatogram gallery · overlay · analysis dashboard."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from hplc_gc_pipeline.gui.common import (
    _bokeh_content_height, embed_html, result_header, run_stage, load_into_state,
)

exp          = Path(st.session_state.loaded_exp)
analysis_dir = st.session_state.get("w_analysis_dir", "analysis")
results_dir  = st.session_state.get("w_results_dir",  "results")
n_cols       = st.session_state.get("w_ncols", 3)

overlay   = exp / results_dir  / "chromatogram_overlay_interactive.html"
dashboard = exp / analysis_dir / "analysis_plots.html"
gallery   = exp / results_dir  / "chromatogram_gallery.html"

# Cache key incremented after each reintegration to bust the browser iframe cache.
_cache_v = st.session_state.get("ri_cache_v", 0)

# ── 1. Per-injection chromatogram gallery ─────────────────────────────────
if gallery.exists():
    st.subheader("Per-injection chromatograms")
    result_header(gallery)
    embed_html(gallery, height=900, scrolling=True, cache_key=_cache_v)

# ── 2. Chromatogram overlay ───────────────────────────────────────────────
if overlay.exists():
    if gallery.exists():
        st.divider()
    st.subheader("Chromatogram overlay")
    result_header(overlay)
    embed_html(overlay, height=700, cache_key=_cache_v)

# ── 3. Re-integrate a single injection ───────────────────────────────────
_ro = not st.session_state.get("edit_mode", True)
st.divider()

# ── Preview from last reintegration (outside the expander so it survives rerun) ──
_ri_prev = st.session_state.get("ri_preview_png")
if _ri_prev and Path(_ri_prev).exists():
    _ri_prev_folder = st.session_state.get("ri_preview_folder", "")
    st.markdown(f"**Last re-integration result — `{_ri_prev_folder}`**")
    st.image(str(_ri_prev), use_container_width=True)
    if st.button("Clear preview", icon=":material/close:", key="ri_clear"):
        st.session_state.pop("ri_preview_png", None)
        st.session_state.pop("ri_preview_folder", None)
        st.rerun()
    st.divider()

with st.expander("Re-integrate one injection", icon=":material/build:"):
    if _ro:
        st.info(
            "Enable **Edit mode** in the sidebar to use this tool.",
            icon=":material/visibility:",
        )
    else:
        st.caption(
            "Re-process one injection in-place without re-running the full pipeline. "
            "Patches `peak_results.csv`, refreshes the plot in the gallery above, "
            "and re-runs the analysis automatically."
        )

        # Build folder list from peak_results.csv (already processed injections)
        peak_csv = exp / results_dir / "peak_results.csv"
        _folders: list[str] = []
        if peak_csv.exists():
            try:
                _folders = sorted(
                    pd.read_csv(peak_csv, usecols=["Folder"])["Folder"].dropna().unique().tolist()
                )
            except Exception:
                pass
        if not _folders:
            _data_root_str = st.session_state.get("w_datadir", "")
            _data_root = (exp / _data_root_str) if _data_root_str else exp
            _folders = sorted(p.name for p in _data_root.glob("*.D"))

        if not _folders:
            st.warning("No injections found. Run Stage 1 (Process) first.")
        else:
            ri_folder = st.selectbox(
                "Injection to re-integrate",
                _folders,
                help="Pick the folder that had a bad deconvolution in the gallery above.",
                key="ri_folder",
            )

            with st.expander("Override MOCCA2 parameters for this injection"):
                st.caption(
                    "Leave at defaults to re-run with the same settings as the full pipeline. "
                    "Increase **min R²** to force MOCCA2 to split merged peaks. "
                    "Enable **No deconvolution** to fall back to manual trapezoid integration."
                )
                ri_col1, ri_col2 = st.columns(2)
                with ri_col1:
                    ri_min_r2 = st.number_input(
                        "min R²",
                        min_value=0.50, max_value=1.00,
                        value=float(st.session_state.get("w_minr2", 0.95)),
                        step=0.01, format="%.2f",
                        help="Higher → more components tried before accepting the fit.",
                        key="ri_min_r2",
                    )
                    ri_max_comps = st.number_input(
                        "max components",
                        min_value=1, max_value=20,
                        value=int(st.session_state.get("w_maxcomps", 5)),
                        step=1,
                        key="ri_max_comps",
                    )
                with ri_col2:
                    ri_min_height = st.text_input(
                        "min height (blank = use config)",
                        value="",
                        placeholder="e.g. 500",
                        key="ri_min_height",
                    )
                    ri_min_prom = st.text_input(
                        "min prominence (blank = use config)",
                        value="",
                        placeholder="e.g. 200",
                        key="ri_min_prom",
                    )
                ri_no_deconv = st.checkbox(
                    "No deconvolution (manual trapezoid fallback)",
                    value=False,
                    key="ri_no_deconv",
                )

            st.caption(
                "**Manual RT bounds** — use when MOCCA2 and auto-trapezoid both fail "
                "(e.g. two overlapping peaks that aren't split correctly). "
                "For each compound, specify the left and right RT limits to integrate. "
                "The pipeline integrates the baseline-corrected signal in that exact window "
                "and replaces any existing peaks there. Leave at 0 to skip a compound."
            )
            _compounds_path = exp / st.session_state.get("w_compounds_file", "compounds.csv")
            _bound_compounds: list[str] = []
            if _compounds_path.exists():
                try:
                    _bound_compounds = pd.read_csv(_compounds_path)["Compound"].dropna().tolist()
                except Exception:
                    pass

            ri_bounds_active = st.checkbox(
                "Specify manual RT bounds for one or more compounds",
                value=False,
                key="ri_bounds_active",
            )
            ri_bounds: dict[str, tuple[float, float]] = {}
            if ri_bounds_active and _bound_compounds:
                for _cname in _bound_compounds:
                    _ckey = _cname.replace(" ", "_").replace("/", "_")
                    _bcol1, _bcol2 = st.columns(2)
                    with _bcol1:
                        _bleft = st.number_input(
                            f"{_cname}  left RT (min)",
                            min_value=0.0, max_value=60.0, value=0.0,
                            step=0.01, format="%.3f",
                            key=f"ri_bound_left_{_ckey}",
                        )
                    with _bcol2:
                        _bright = st.number_input(
                            f"{_cname}  right RT (min)",
                            min_value=0.0, max_value=60.0, value=0.0,
                            step=0.01, format="%.3f",
                            key=f"ri_bound_right_{_ckey}",
                        )
                    if _bleft > 0 and _bright > _bleft:
                        ri_bounds[_cname] = (_bleft, _bright)

            ri_skip_analysis = st.checkbox(
                "Skip analysis re-run (patch peak_results.csv only)",
                value=False,
                key="ri_skip_analysis",
                help="Use when re-integrating multiple injections; run Analyze once at the end.",
            )

            if st.button(
                f"Re-integrate  {ri_folder}",
                type="primary",
                icon=":material/refresh:",
                key="ri_run_btn",
            ):
                extra: list[str] = ["--injection", ri_folder]

                cfg_min_r2 = float(st.session_state.get("w_minr2", 0.95))
                cfg_max_comps = int(st.session_state.get("w_maxcomps", 5))
                if abs(ri_min_r2 - cfg_min_r2) > 1e-6:
                    extra += ["--min-r2", str(ri_min_r2)]
                if ri_max_comps != cfg_max_comps:
                    extra += ["--max-comps", str(ri_max_comps)]
                _h = ri_min_height.strip()
                if _h:
                    extra += ["--min-height", _h]
                _p = ri_min_prom.strip()
                if _p:
                    extra += ["--min-prominence", _p]
                if ri_no_deconv:
                    extra.append("--no-deconvolution")
                for _bc, (_bl, _br) in ri_bounds.items():
                    extra += ["--bounds", f"{_bc}={_bl}:{_br}"]
                if ri_skip_analysis:
                    extra.append("--skip-analysis")

                with st.spinner(f"Re-integrating {ri_folder} …"):
                    code = run_stage(exp, "reintegrate", extra)

                if code == 0:
                    _stem = ri_folder.replace(".D", "")
                    _plots_dir = exp / results_dir / "plots"
                    _preview_path = None
                    # Prefer the manual-bounds preview when bounds were specified.
                    if ri_bounds:
                        _mb = sorted(_plots_dir.glob(f"*{_stem}*manual_preview*.png"))
                        if _mb:
                            _preview_path = str(_mb[-1])
                    if _preview_path is None:
                        _all = [
                            p for p in sorted(_plots_dir.glob(f"*{_stem}*.png"))
                            if "manual_preview" not in p.name
                        ]
                        if _all:
                            _preview_path = str(_all[0])
                    if _preview_path:
                        st.session_state.ri_preview_png = _preview_path
                        st.session_state.ri_preview_folder = ri_folder
                    # Bump cache key so all embedded iframes reload on next render.
                    st.session_state.ri_cache_v = _cache_v + 1
                    load_into_state(exp)
                    st.rerun()
                else:
                    st.error(
                        f"Re-integration failed (exit {code}). See the terminal log above.",
                        icon=":material/error:",
                    )

# ── 4. Analysis dashboard (full content height, no inner scroll) ──────────
if dashboard.exists():
    st.divider()
    st.subheader("Analysis")
    result_header(dashboard)
    embed_html(dashboard,
               height=_bokeh_content_height(dashboard, n_cols=n_cols),
               scrolling=False,
               cache_key=_cache_v)
