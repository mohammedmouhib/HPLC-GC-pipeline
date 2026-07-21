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
# Sidebar — folder picker + status (appears above the navigation links)
# --------------------------------------------------------------------------
with st.sidebar:
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
            st.success("Loaded.")
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
            st.info("Go to **Run → Initialize** to generate a starter config.",
                    icon=":material/arrow_forward:")
        elif not _results_ok:
            st.info("Go to **Run → Run both** to process your data.",
                    icon=":material/arrow_forward:")
        else:
            st.caption(":material/check_circle: Ready — open **Results**.")

# --------------------------------------------------------------------------
# Auto-load when a path was pre-supplied on the command line
# --------------------------------------------------------------------------
if "doc" not in st.session_state and exp_str and Path(exp_str).is_dir():
    load_into_state(Path(exp_str))

# --------------------------------------------------------------------------
# Navigation
# --------------------------------------------------------------------------
def _welcome() -> None:
    st.space("large")
    with st.container(horizontal_alignment="center"):
        st.markdown(":material/science:")
        st.title("HPLC / GC-MS pipeline", text_alignment="center")
        st.caption(
            "Process raw Agilent .D injections → calibrated concentrations → interactive plots",
            text_alignment="center",
        )

    st.space("large")
    col1, col2, col3 = st.columns(3)
    with col1:
        with st.container(border=True):
            st.markdown(":material/folder_open: **Step 1 — Load**")
            st.caption(
                "Click **Browse for folder** in the sidebar to pick your experiment folder. "
                "The folder should contain Agilent `.D` injection subfolders."
            )
    with col2:
        with st.container(border=True):
            st.markdown(":material/auto_awesome: **Step 2 — Initialize** *(new experiments)*")
            st.caption(
                "Open **Run** and click **Initialize**. "
                "The pipeline scans your `.D` files, detects peaks and m/z values, "
                "and writes a starter config automatically."
            )
    with col3:
        with st.container(border=True):
            st.markdown(":material/play_arrow: **Step 3 — Run & explore**")
            st.caption(
                "Click **Run both** to process all injections and build calibration curves. "
                "Results appear in **Results** as interactive Bokeh dashboards."
            )

    st.space("large")
    st.caption(
        "Already have a configured experiment? Just enter the folder path in the sidebar — "
        "the pipeline will pick up any existing config automatically.",
        text_alignment="center",
    )


if "doc" not in st.session_state:
    pages = [st.Page(_welcome, title="Get started", icon=":material/home:", default=True)]
else:
    pages = [
        st.Page("app_pages/parameters.py", title="Parameters",
                icon=":material/tune:", default=True),
        st.Page("app_pages/run.py", title="Run",
                icon=":material/play_arrow:"),
        st.Page("app_pages/results.py", title="Results",
                icon=":material/analytics:"),
    ]

pg = st.navigation(pages)
pg.run()
