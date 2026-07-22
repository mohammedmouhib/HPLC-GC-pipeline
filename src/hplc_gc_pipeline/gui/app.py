"""Entry point for the HPLC / GC-MS pipeline GUI.

Handles page config, the sidebar (folder picker + status), and navigation.
All shared helpers live in common.py; page content lives in app_pages/.
"""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import (
    _count_d_folders,
    _has_any_results,
    _pick_folder,
    initial_exp_dir,
    load_into_state,
)

# --------------------------------------------------------------------------
# Page config  (must be the very first Streamlit call)
# --------------------------------------------------------------------------
st.set_page_config(
    page_title="HPLC / GC-MS pipeline",
    page_icon=":material/science:",
    layout="wide",
)

# --------------------------------------------------------------------------
# Sidebar — folder picker + status
# --------------------------------------------------------------------------
with st.sidebar:
    # Apply any pending edit-mode default BEFORE the toggle widget renders.
    if "_reset_edit_mode" in st.session_state:
        st.session_state.edit_mode = st.session_state.pop("_reset_edit_mode")

    # Edit mode toggle — first sidebar item once a folder is loaded,
    # positioned just below the navigation links.
    if "doc" in st.session_state:
        st.toggle(
            "Edit mode",
            key="edit_mode",
            help="Off: browse parameters and results without changes. "
                 "On: modify parameters and re-run the pipeline.",
        )
        st.divider()

    st.header("Experiment folder")

    if st.button("Browse for folder", icon=":material/folder_open:", width="stretch"):
        picked = _pick_folder()
        if picked and Path(picked).is_dir():
            st.session_state.exp_path_input = picked
            load_into_state(Path(picked))
            st.rerun()

    exp_str = st.text_input(
        "Experiment folder path",
        key="exp_path_input",
        value=st.session_state.get("loaded_exp", initial_exp_dir()),
        label_visibility="collapsed",
        placeholder="or paste a path here…",
        help="Folder containing your Agilent .D injection files (or subfolders of them).",
    )
    if st.button("Load", type="primary", width="stretch", icon=":material/refresh:"):
        if exp_str and Path(exp_str).is_dir():
            load_into_state(Path(exp_str))
            st.rerun()
        else:
            st.error("Not a valid directory.")

    # Status panel — only shown once a folder is loaded
    if "doc" in st.session_state:
        st.divider()
        _exp = Path(st.session_state.loaded_exp)
        _s = st.session_state
        _config_ok = (_exp / "hplc_config.yaml").exists()
        _data_root = _exp / _s.get("w_datadir") if _s.get("w_datadir") else _exp
        _n = _count_d_folders(_data_root)
        _results_ok = _has_any_results(_exp)

        st.caption("**Status**")
        st.badge(
            "Config ready" if _config_ok else "No config yet",
            icon=":material/check:" if _config_ok else ":material/warning:",
            color="green" if _config_ok else "orange",
        )
        st.badge(
            f"{_n} injections found" if _n else "No .D files found",
            icon=":material/science:",
            color="green" if _n else "red",
        )
        if _results_ok:
            st.badge("Results ready", icon=":material/analytics:", color="blue")

        if _s.get("dirty"):
            st.warning("Unsaved edits — apply in Parameters.", icon=":material/edit:")

        st.space("small")
        if not _config_ok:
            st.info("Go to **Initialize** to write a starter config.",
                    icon=":material/arrow_forward:")
        elif not _results_ok:
            st.info("Go to **Run** → click **Run both** to process your data.",
                    icon=":material/arrow_forward:")
        else:
            st.caption(":material/check_circle: Ready — see Results in the sidebar.")

# --------------------------------------------------------------------------
# Auto-load when a path was pre-supplied on the command line, or re-init
# if doc exists but w_* keys were dropped (e.g. after a navigation rerun).
# --------------------------------------------------------------------------
if exp_str and Path(exp_str).is_dir() and (
    "doc" not in st.session_state or "w_gc_enabled" not in st.session_state
):
    load_into_state(Path(exp_str))

# --------------------------------------------------------------------------
# Navigation
# --------------------------------------------------------------------------

if "doc" not in st.session_state:
    pages = [st.Page("app_pages/home.py", title="Get started", icon=":material/home:", default=True)]
else:
    exp = Path(st.session_state.loaded_exp)
    analysis_dir = st.session_state.get("w_analysis_dir", "analysis")
    results_dir  = st.session_state.get("w_results_dir",  "results")

    result_pages = []
    if any((exp / "gc_analysis" / f).exists() or (exp / "gc_results" / f).exists()
           for f in ("analysis_plots.html",
                     "gc_tic_overlay_interactive.html",
                     "chromatogram_gallery.html")):
        result_pages.append(
            st.Page("app_pages/results_gc.py", title="GC-MS",
                    icon=":material/biotech:")
        )
    if any(p.exists() for p in (
        exp / analysis_dir / "analysis_plots.html",
        exp / results_dir  / "chromatogram_overlay_interactive.html",
        exp / results_dir  / "chromatogram_gallery.html",
    )):
        result_pages.append(
            st.Page("app_pages/results_hplc.py", title="HPLC",
                    icon=":material/water_drop:")
        )
    if (exp / "combined_analysis" / "analysis_plots.html").exists():
        result_pages.append(
            st.Page("app_pages/results_combined.py", title="Combined",
                    icon=":material/merge:")
        )

    core_pages = [
        st.Page("app_pages/home.py", title="Guide", icon=":material/book:"),
        st.Page("app_pages/initialize.py", title="Initialize", icon=":material/auto_awesome:"),
        st.Page("app_pages/parameters.py", title="Parameters",
                icon=":material/tune:", default=True),
        st.Page("app_pages/run.py", title="Run", icon=":material/play_arrow:"),
    ]
    pages = {"": core_pages, "Results": result_pages} if result_pages else core_pages

pg = st.navigation(pages)
pg.run()
