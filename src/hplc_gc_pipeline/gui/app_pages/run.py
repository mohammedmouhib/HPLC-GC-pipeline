"""Run page — initialize config and run the pipeline stages."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import load_into_state, run_stage

exp = Path(st.session_state.loaded_exp)

st.caption(f"Experiment: `{exp}`")

# ── Step 1: Initialize ────────────────────────────────────────────────────
config_path = exp / "hplc_config.yaml"
gc_compounds_path = exp / st.session_state.get("w_gc_compounds_file", "gc_compounds.csv")
files_exist = config_path.exists() or gc_compounds_path.exists()

with st.container(border=True):
    hdr_col, badge_col = st.columns([5, 1], vertical_alignment="center")
    with hdr_col:
        st.markdown(":material/auto_awesome: **Step 1 · Initialize** — auto-detect from GC-MS data")
    with badge_col:
        if config_path.exists():
            st.badge("Config ready", icon=":material/check:", color="green")
        else:
            st.badge("No config", icon=":material/warning:", color="orange")
    st.caption(
        "Scans your `.D` files, identifies standard injections, detects peak RT windows "
        "and quantifier ions, and writes a starter `hplc_config.yaml` + `gc_compounds.csv`. "
        "Run once per new experiment, then review compound names in **Parameters**."
    )
    with st.expander("Init options", icon=":material/settings:"):
        ic1, ic2, ic3 = st.columns(3)
        init_rt_min = ic1.number_input(
            "RT scan start (min)", value=1.5, min_value=0.0, step=0.5,
            help="Peaks before this time are ignored (skips solvent front).",
        )
        init_rt_max = ic2.number_input("RT scan end (min)", value=30.0, min_value=1.0, step=5.0)
        init_rt_margin = ic3.number_input(
            "RT window margin (min)", value=0.08, min_value=0.0, step=0.01, format="%.2f",
            help="Added to each side of a detected peak edge when writing the RT window.",
        )
    if files_exist:
        st.warning(
            "Config or compounds file already exists and will be overwritten.",
            icon=":material/warning:",
        )
    if st.button("Initialize", icon=":material/auto_awesome:", type="secondary"):
        extra = [
            "--force",
            "--rt-min", str(init_rt_min),
            "--rt-max", str(init_rt_max),
            "--rt-margin", str(init_rt_margin),
        ]
        with st.spinner("Running `hplc init` …"):
            code = run_stage(exp, "init", extra_args=extra)
        if code == 0:
            st.success("Config files written. Reloading parameters…", icon=":material/check:")
            load_into_state(exp)
            st.rerun()
        else:
            st.error(f"`hplc init` failed (exit {code}). See log above.")

# ── Step 2: Run pipeline ──────────────────────────────────────────────────
with st.container(border=True):
    st.markdown(":material/play_arrow: **Step 2 · Run the pipeline**")
    st.caption(
        "**Process** (Stage 1) — integrate peaks and generate chromatogram plots.  "
        "**Analyze** (Stage 2) — calibrate, quantify, and build the Bokeh dashboard.  "
        "**Run both** — full pipeline end-to-end."
    )
    if st.session_state.get("dirty"):
        st.warning("Unsaved parameter edits — go to **Parameters** and click **Apply changes** first.",
                   icon=":material/edit:")
    with st.container(horizontal=True):
        go_process = st.button("Process", icon=":material/table_chart:", width="stretch")
        go_analyze = st.button("Analyze", icon=":material/analytics:", width="stretch")
        go_run = st.button("Run both", type="primary", icon=":material/play_arrow:", width="stretch")

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
    st.download_button("Download last log", data=lp.read_bytes(),
                       file_name=lp.name, mime="text/plain", icon=":material/download:")
