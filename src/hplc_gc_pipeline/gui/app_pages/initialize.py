"""Initialize page — create config files for a new experiment."""

from __future__ import annotations

import shutil
from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import (
    PKG_TEMPLATE,
    _pick_folder,
    load_into_state,
    run_stage,
)

exp = Path(st.session_state.loaded_exp)
_ro = not st.session_state.get("edit_mode", True)

config_path = exp / "hplc_config.yaml"
gc_compounds_path = exp / st.session_state.get("w_gc_compounds_file", "gc_compounds.csv")

st.caption(f"Experiment: `{exp}`")

if _ro:
    st.info(
        "View only — toggle **Edit mode** in the sidebar to run Initialize.",
        icon=":material/visibility:",
    )

col_hdr, col_badge = st.columns([5, 1], vertical_alignment="center")
with col_hdr:
    st.markdown(":material/auto_awesome: **Initialize** — create config files for this experiment")
with col_badge:
    if config_path.exists():
        st.badge("Config ready", icon=":material/check:", color="green")
    else:
        st.badge("No config", icon=":material/warning:", color="orange")

st.caption(
    "Run once per new experiment to generate `hplc_config.yaml` and (for GC-MS) "
    "`gc_compounds.csv`. Then review and refine compound names and RT windows in **Parameters**."
)

# ── Step A: Modality ───────────────────────────────────────────────────────
with st.container(border=True):
    st.markdown("**Step A · Experiment modality**")
    mode = st.segmented_control(
        "Modality",
        options=["HPLC only", "GC only", "HPLC + GC"],
        default="HPLC + GC",
        key="init_mode",
        label_visibility="collapsed",
        help="HPLC only: writes a template config. GC only / HPLC + GC: scans .D files to "
             "auto-detect peaks and quantifier ions.",
        disabled=_ro,
    )

# ── Step B: Config source ──────────────────────────────────────────────────
with st.container(border=True):
    st.markdown("**Step B · Config source**")
    source = st.segmented_control(
        "Config source",
        options=["Auto-detect from .D files", "Copy from a previous project"],
        default="Auto-detect from .D files",
        key="init_source",
        label_visibility="collapsed",
        help="Auto-detect scans your .D files and writes starter config. "
             "Copy reuses config files from a previous experiment with the same compounds.",
        disabled=_ro,
    )

# ── Option A: Auto-detect ─────────────────────────────────────────────────
if source == "Auto-detect from .D files":
    with st.container(border=True):
        if mode == "HPLC only":
            st.info(
                "For HPLC-only experiments, Initialize writes a **template `hplc_config.yaml`** "
                "with placeholder values. Edit compound names and RT windows in **Parameters** "
                "after initializing.",
                icon=":material/info:",
            )
        else:
            st.markdown("**GC peak detection settings**")
            ic1, ic2, ic3 = st.columns(3)
            init_rt_min = ic1.number_input(
                "RT scan start (min)", value=1.5, min_value=0.0, step=0.5,
                help="Peaks before this time are ignored (skips solvent front).",
                disabled=_ro,
            )
            init_rt_max = ic2.number_input(
                "RT scan end (min)", value=30.0, min_value=1.0, step=5.0, disabled=_ro,
            )
            init_rt_margin = ic3.number_input(
                "RT window margin (min)", value=0.08, min_value=0.0, step=0.01, format="%.2f",
                help="Added to each side of a detected peak edge when writing the RT window.",
                disabled=_ro,
            )

        files_exist = config_path.exists() or gc_compounds_path.exists()
        if files_exist:
            st.warning(
                "Config or compounds file already exists and will be overwritten.",
                icon=":material/warning:",
            )

        if st.button("Initialize", icon=":material/auto_awesome:", type="primary", disabled=_ro):
            if mode == "HPLC only":
                shutil.copy2(PKG_TEMPLATE, config_path)
                st.success(
                    "Template config written. Edit compound names and RT windows in **Parameters**.",
                    icon=":material/check:",
                )
                load_into_state(exp)
                st.rerun()
            else:
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

# ── Option B: Copy from previous project ─────────────────────────────────
elif source == "Copy from a previous project":
    with st.container(border=True):
        st.caption(
            "Copy `hplc_config.yaml` and CSV files from a previous experiment that used "
            "the same compounds and processing parameters."
        )
        c_browse, c_input, c_scan = st.columns([1, 4, 1])
        if c_browse.button("Browse", icon=":material/folder_open:", key="cp_browse"):
            picked = _pick_folder()
            if picked:
                st.session_state.cp_src = picked
        src_str = c_input.text_input(
            "Source folder",
            key="cp_src",
            label_visibility="collapsed",
            placeholder="path to previous experiment folder…",
        )
        c_scan.button("Load", icon=":material/search:", key="cp_scan_btn")

        if src_str:
            src = Path(src_str)
            if not src.is_dir():
                st.error("Not a valid directory.")
            else:
                known = ["hplc_config.yaml"]
                known += sorted(p.name for p in src.glob("*.csv"))
                candidates = [f for f in known if (src / f).exists()]

                if not candidates:
                    st.info("No recognised config files found in that folder.")
                else:
                    will_overwrite = [f for f in candidates if (exp / f).exists()]
                    if will_overwrite:
                        st.warning(
                            "These files already exist and will be overwritten: "
                            + ", ".join(f"`{f}`" for f in will_overwrite),
                            icon=":material/warning:",
                        )

                    st.caption(f"Files found in `{src.name}/` — select which to copy:")
                    cols = st.columns(3)
                    selected = [
                        f for i, f in enumerate(candidates)
                        if cols[i % 3].checkbox(f, value=True, key=f"cp_cb_{f}")
                    ]

                    if selected and st.button(
                        "Copy files", icon=":material/content_copy:", type="primary",
                        disabled=_ro,
                    ):
                        for fname in selected:
                            shutil.copy2(src / fname, exp / fname)
                        st.toast(
                            f"Copied {len(selected)} file(s) from `{src.name}`.",
                            icon=":material/check:",
                        )
                        load_into_state(exp)
                        st.rerun()
