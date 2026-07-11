"""Streamlit GUI for the HPLC pipeline.

Launched by ``hplc gui`` (see launch.py). It is a thin front-end over the same
three input files and the same CLI the pipeline already uses:

    1. Edit hplc_config.yaml + compounds.csv + standard.csv in the browser.
    2. "Apply" writes them back to the experiment folder (YAML comments kept).
    3. Run buttons call `hplc process/analyze/run` and stream the log.
    4. Result links open the output HTML files in the browser.

Nothing about the processing backend changes: the GUI only touches the files
and invokes the CLI.
"""

from __future__ import annotations

import subprocess
import sys
import webbrowser
from io import StringIO
from pathlib import Path

import pandas as pd
import streamlit as st
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

yaml = YAML()  # round-trip loader: preserves comments and formatting
yaml.preserve_quotes = True

PKG_TEMPLATE = Path(__file__).with_name("config_template.yaml")
COMPOUND_COLS = ["Compound", "RT_low", "RT_high", "Notes"]
STANDARD_COLS = ["Compound", "Area_Integral", "concentration"]


# --------------------------------------------------------------------------
# File helpers
# --------------------------------------------------------------------------
def initial_exp_dir() -> str:
    """Experiment dir passed after `--` on the streamlit command line, if any."""
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    return args[0] if args else ""


def load_config_doc(exp: Path) -> CommentedMap:
    """Round-trip-load the experiment's YAML, or the packaged template if none."""
    src = exp / "hplc_config.yaml"
    if not src.exists():
        src = PKG_TEMPLATE
    with open(src, "r", encoding="utf-8") as fh:
        return yaml.load(fh)


def cfg_get(doc, path, default=None):
    cur = doc
    for key in path:
        if not isinstance(cur, dict) or key not in cur or cur[key] is None:
            return default
        cur = cur[key]
    return cur


def cfg_set(doc, path, value):
    cur = doc
    for key in path[:-1]:
        if key not in cur or not isinstance(cur[key], dict):
            cur[key] = CommentedMap()
        cur = cur[key]
    cur[path[-1]] = value


def read_csv_or_empty(path: Path, cols: list[str]) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path)
    return pd.DataFrame(columns=cols)


def _csv_list(text: str):
    """'a, b, c' -> ['a','b','c']; '' -> None (meaning 'all/omit')."""
    items = [t.strip() for t in str(text).split(",") if t.strip()]
    return items or None


def _int_list(text: str):
    items = [t.strip() for t in str(text).split(",") if t.strip()]
    return [int(t) for t in items] or None


# --------------------------------------------------------------------------
# State: load the three files into widget-backing session state
# --------------------------------------------------------------------------
def load_into_state(exp: Path) -> None:
    doc = load_config_doc(exp)
    st.session_state.doc = doc
    st.session_state.w_wavelength = int(cfg_get(doc, ["processing", "wavelength_nm"], 262))
    st.session_state.w_blank = cfg_get(doc, ["processing", "blank_folder_name"], "") or ""
    st.session_state.w_datadir = cfg_get(doc, ["processing", "data_dir"], "") or ""
    st.session_state.w_deconv = bool(cfg_get(doc, ["processing", "deconvolution", "enabled"], True))
    st.session_state.w_minr2 = float(cfg_get(doc, ["processing", "deconvolution", "min_r2"], 0.95))
    st.session_state.w_maxcomps = int(cfg_get(doc, ["processing", "deconvolution", "max_comps"], 5))
    st.session_state.w_excl = bool(cfg_get(doc, ["analysis", "exclude_unknown"], True))
    st.session_state.w_cv = float(cfg_get(doc, ["analysis", "cv_warning_threshold"], 15.0))
    st.session_state.w_clamp = bool(cfg_get(doc, ["analysis", "calibration", "clamp_negative_to_zero"], True))
    st.session_state.w_zero = bool(cfg_get(doc, ["analysis", "calibration", "zero_area_zero_conc"], True))
    st.session_state.w_pattern = cfg_get(doc, ["analysis", "sample_name", "pattern"], "")
    st.session_state.w_ncols = int(cfg_get(doc, ["plots", "n_cols"], 2))
    st.session_state.w_strains = ", ".join(cfg_get(doc, ["plots", "plot_strains"], []) or [])
    st.session_state.w_order = ", ".join(cfg_get(doc, ["plots", "strain_order"], []) or [])
    st.session_state.w_bartp = ", ".join(str(x) for x in (cfg_get(doc, ["plots", "bar_time_points"], []) or []))
    st.session_state.compounds_df = read_csv_or_empty(exp / "compounds.csv", COMPOUND_COLS)
    st.session_state.standard_df = read_csv_or_empty(exp / "standard.csv", STANDARD_COLS)
    # Reset the data-editor widgets so they re-seed from the reloaded frames.
    st.session_state.pop("ed_compounds", None)
    st.session_state.pop("ed_standard", None)
    st.session_state.loaded_exp = str(exp)


def apply_changes(exp: Path, compounds_df, standard_df) -> tuple[bool, str]:
    """Write all three files, then validate the config. Returns (ok, message)."""
    doc = st.session_state.doc
    cfg_set(doc, ["processing", "wavelength_nm"], int(st.session_state.w_wavelength))
    cfg_set(doc, ["processing", "blank_folder_name"], st.session_state.w_blank or None)
    cfg_set(doc, ["processing", "data_dir"], st.session_state.w_datadir or None)
    cfg_set(doc, ["processing", "deconvolution", "enabled"], bool(st.session_state.w_deconv))
    cfg_set(doc, ["processing", "deconvolution", "min_r2"], float(st.session_state.w_minr2))
    cfg_set(doc, ["processing", "deconvolution", "max_comps"], int(st.session_state.w_maxcomps))
    cfg_set(doc, ["analysis", "exclude_unknown"], bool(st.session_state.w_excl))
    cfg_set(doc, ["analysis", "cv_warning_threshold"], float(st.session_state.w_cv))
    cfg_set(doc, ["analysis", "calibration", "clamp_negative_to_zero"], bool(st.session_state.w_clamp))
    cfg_set(doc, ["analysis", "calibration", "zero_area_zero_conc"], bool(st.session_state.w_zero))
    cfg_set(doc, ["analysis", "sample_name", "pattern"], st.session_state.w_pattern)
    cfg_set(doc, ["plots", "n_cols"], int(st.session_state.w_ncols))
    cfg_set(doc, ["plots", "plot_strains"], _csv_list(st.session_state.w_strains))
    cfg_set(doc, ["plots", "strain_order"], _csv_list(st.session_state.w_order))
    cfg_set(doc, ["plots", "bar_time_points"], _int_list(st.session_state.w_bartp))

    exp.mkdir(parents=True, exist_ok=True)
    buf = StringIO()
    yaml.dump(doc, buf)
    (exp / "hplc_config.yaml").write_text(buf.getvalue(), encoding="utf-8")
    compounds_df.dropna(how="all").to_csv(exp / "compounds.csv", index=False)
    standard_df.dropna(how="all").to_csv(exp / "standard.csv", index=False)

    # Validate by loading through the real config machinery.
    try:
        from hplc_gc_pipeline.config import load_config
        load_config(exp)
    except Exception as exc:  # noqa: BLE001 - surface any validation error verbatim
        return False, f"Files written, but config is invalid: {exc}"
    return True, "All three files written and the config validates."


# --------------------------------------------------------------------------
# Running the pipeline (stream the CLI log into the page)
# --------------------------------------------------------------------------
NOISE = ("Sparse", "spsolve", "flatfit", "UserWarning", "RuntimeWarning")


def run_stage(exp: Path, subcommand: str) -> int:
    cmd = [sys.executable, "-m", "hplc_gc_pipeline.cli", subcommand, str(exp)]
    box = st.empty()
    lines: list[str] = []
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    )
    assert proc.stdout is not None
    for raw in proc.stdout:
        if any(n in raw for n in NOISE):
            continue
        lines.append(raw.rstrip())
        box.code("\n".join(lines[-400:]), language="text")
    proc.wait()
    st.session_state.last_log = "\n".join(lines)
    return proc.returncode


# --------------------------------------------------------------------------
# UI
# --------------------------------------------------------------------------
st.set_page_config(page_title="HPLC pipeline", layout="wide")
st.title("HPLC pipeline")

with st.sidebar:
    st.header("Experiment folder")
    exp_str = st.text_input(
        "Path", value=st.session_state.get("loaded_exp", initial_exp_dir()),
        help="Folder holding your .D data (or a Data/ subfolder) plus the config and CSVs.",
    )
    if st.button("Load / Reload", type="primary", use_container_width=True):
        if exp_str and Path(exp_str).is_dir():
            load_into_state(Path(exp_str))
            st.success("Loaded.")
        else:
            st.error("Not a directory.")

# Auto-load once if a valid folder was passed on the command line.
if "doc" not in st.session_state and exp_str and Path(exp_str).is_dir():
    load_into_state(Path(exp_str))

if "doc" not in st.session_state:
    st.info("Enter an experiment folder in the sidebar and click **Load / Reload**.")
    st.stop()

exp = Path(st.session_state.loaded_exp)

tab_params, tab_run = st.tabs(["1 · Parameters", "2 · Run & results"])

with tab_params:
    st.subheader("Processing (Stage 1)")
    c1, c2, c3 = st.columns(3)
    c1.number_input("Wavelength (nm)", key="w_wavelength", step=1)
    c2.text_input("Blank .D folder", key="w_blank", placeholder="e.g. 091-0202.D  (blank = none)")
    c3.text_input("Data subfolder", key="w_datadir", placeholder="e.g. Data  (blank = experiment root)")
    c4, c5, c6 = st.columns(3)
    c4.checkbox("Deconvolution enabled", key="w_deconv")
    c5.number_input("Min R² to accept fit", key="w_minr2", step=0.01, format="%.2f")
    c6.number_input("Max components / peak", key="w_maxcomps", step=1)

    st.subheader("Analysis (Stage 2)")
    a1, a2 = st.columns(2)
    a1.checkbox("Exclude peaks outside all RT windows", key="w_excl")
    a2.number_input("CV%% warning threshold", key="w_cv", step=1.0)
    a3, a4 = st.columns(2)
    a3.checkbox("Clamp negative concentrations to 0", key="w_clamp")
    a4.checkbox("Zero area → 0 µM (recommended fix)", key="w_zero")
    st.text_input("Sample-name regex", key="w_pattern")

    st.subheader("Plots")
    p1, p2, p3, p4 = st.columns(4)
    p1.number_input("Grid columns", key="w_ncols", step=1, min_value=1)
    p2.text_input("Plot strains", key="w_strains", placeholder="all (comma-separated)")
    p3.text_input("Strain order", key="w_order", placeholder="alphabetical")
    p4.text_input("Bar time points", key="w_bartp", placeholder="all (e.g. 0, 24, 50)")

    st.subheader("compounds.csv — retention-time windows")
    compounds_df = st.data_editor(
        st.session_state.compounds_df, num_rows="dynamic",
        use_container_width=True, key="ed_compounds",
    )
    st.subheader("standard.csv — calibration points")
    standard_df = st.data_editor(
        st.session_state.standard_df, num_rows="dynamic",
        use_container_width=True, key="ed_standard",
    )

    if st.button("💾  Apply — write all three files", type="primary"):
        ok, msg = apply_changes(exp, compounds_df, standard_df)
        (st.success if ok else st.error)(msg)

with tab_run:
    st.caption(f"Experiment folder: `{exp}`")
    st.warning("Click **Apply** on the Parameters tab first if you changed anything.")
    r1, r2, r3 = st.columns(3)
    go_process = r1.button("▶  Process (Stage 1)", use_container_width=True)
    go_analyze = r2.button("▶  Analyze (Stage 2)", use_container_width=True)
    go_run = r3.button("▶  Run both", type="primary", use_container_width=True)

    if go_process or go_analyze or go_run:
        sub = "process" if go_process else "analyze" if go_analyze else "run"
        with st.spinner(f"Running `hplc {sub}` …"):
            code = run_stage(exp, sub)
        (st.success if code == 0 else st.error)(
            f"`hplc {sub}` finished (exit {code})."
        )

    st.subheader("Results")
    results_dir = cfg_get(st.session_state.doc, ["output", "results_dirname"], "results")
    analysis_dir = cfg_get(st.session_state.doc, ["output", "analysis_dirname"], "analysis")
    outputs = {
        "Analysis dashboard": exp / analysis_dir / "analysis_plots.html",
        "Chromatogram gallery": exp / results_dir / "chromatogram_gallery.html",
        "Interactive overlay": exp / results_dir / "chromatogram_overlay_interactive.html",
    }
    any_found = False
    for label, path in outputs.items():
        if path.exists():
            any_found = True
            cols = st.columns([1, 4])
            if cols[0].button(f"Open ↗", key=f"open_{label}"):
                webbrowser.open(path.resolve().as_uri())
            cols[1].caption(f"{label} — `{path}`")
    if not any_found:
        st.info("No result files yet — run a stage above.")
