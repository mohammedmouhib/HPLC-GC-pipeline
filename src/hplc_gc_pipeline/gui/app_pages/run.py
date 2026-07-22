"""Run page — run the pipeline stages (Process / Analyze / Run both)."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import run_stage

exp = Path(st.session_state.loaded_exp)

_ro = not st.session_state.get("edit_mode", True)

st.caption(f"Experiment: `{exp}`")

if _ro:
    st.info(
        "View only — toggle **Edit mode** in the sidebar to run the pipeline.",
        icon=":material/visibility:",
    )

with st.container(border=True):
    st.markdown(":material/play_arrow: **Run the pipeline**")
    st.caption(
        "**Process** (Stage 1) — integrate peaks and generate chromatogram plots.  "
        "**Analyze** (Stage 2) — calibrate, quantify, and build the Bokeh dashboard.  "
        "**Run both** — full pipeline end-to-end."
    )
    if st.session_state.get("dirty"):
        st.warning(
            "Unsaved parameter edits — go to **Parameters** and click **Apply changes** first.",
            icon=":material/edit:",
        )
    with st.container(horizontal=True):
        go_process = st.button("Process", icon=":material/table_chart:", width="stretch", disabled=_ro)
        go_analyze = st.button("Analyze", icon=":material/analytics:", width="stretch", disabled=_ro)
        go_run = st.button("Run both", type="primary", icon=":material/play_arrow:", width="stretch", disabled=_ro)

if go_process or go_analyze or go_run:
    sub = "process" if go_process else "analyze" if go_analyze else "run"
    with st.spinner(f"Running `hplc {sub}` …"):
        code = run_stage(exp, sub)
    (st.success if code == 0 else st.error)(
        f"`hplc {sub}` finished (exit {code}).",
        icon=":material/check:" if code == 0 else ":material/error:",
    )
    st.caption(f"Log saved to `{st.session_state.last_log_path}`")
elif st.session_state.get("last_log"):
    with st.expander("Last run log", icon=":material/terminal:"):
        st.code(st.session_state.last_log, language="text")

if st.session_state.get("last_log_path") and Path(st.session_state.last_log_path).exists():
    lp = Path(st.session_state.last_log_path)
    st.download_button(
        "Download last log", data=lp.read_bytes(),
        file_name=lp.name, mime="text/plain", icon=":material/download:",
    )
