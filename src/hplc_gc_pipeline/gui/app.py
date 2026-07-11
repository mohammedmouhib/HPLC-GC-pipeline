"""Streamlit GUI for the HPLC pipeline.

Launched by ``hplc gui`` (see launch.py). A thin front-end over the same three
input files and the same CLI the pipeline already uses:

    * Parameters tab - edit hplc_config.yaml + compounds.csv + standard.csv and
      write them back ("Apply"); YAML comments are preserved.
    * Run tab        - run process / analyze / run, stream the log, save it to
      a .txt file in the experiment folder.
    * Results tab    - view the analysis dashboard and chromatogram overlay
      inline, plus open / download every output.

The processing backend is untouched: the GUI only edits the files and invokes
the CLI.
"""

from __future__ import annotations

import subprocess
import sys
import webbrowser
from datetime import datetime
from io import StringIO
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

yaml = YAML()  # round-trip loader: preserves comments and formatting
yaml.preserve_quotes = True

HERE = Path(__file__).parent
PKG_TEMPLATE = HERE / "config_template.yaml"
PKG_COMPOUNDS = HERE / "compounds_default.csv"
PKG_STANDARD = HERE / "standard_default.csv"
COMPOUND_COLS = ["Compound", "RT_low", "RT_high", "Notes"]
STANDARD_COLS = ["Compound", "Area_Integral", "concentration"]


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def initial_exp_dir() -> str:
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    return args[0] if args else ""


def load_config_doc(exp: Path) -> CommentedMap:
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


def read_csv_or_default(path: Path, pkg_default: Path, cols: list[str]) -> tuple[pd.DataFrame, bool]:
    """Return (frame, seeded). `seeded` = True when it came from the packaged
    example rather than the experiment folder."""
    if path.exists():
        return pd.read_csv(path), False
    if pkg_default.exists():
        return pd.read_csv(pkg_default), True
    return pd.DataFrame(columns=cols), True


def _csv_list(text):
    items = [t.strip() for t in str(text).split(",") if t.strip()]
    return items or None


def _lines_list(text):
    items = [t.strip() for t in str(text).splitlines() if t.strip()]
    return items or None


def _num_or_none(text, cast):
    s = str(text).strip()
    if s == "" or s.lower() == "null":
        return None
    try:
        return cast(s)
    except ValueError:
        return None


def _human_size(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def _file_meta(path: Path) -> str:
    stat = path.stat()
    when = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M")
    return f"{_human_size(stat.st_size)} · {when}"


def _count_d_folders(root: Path) -> int:
    if not root.is_dir():
        return 0
    return len(list(root.glob("*.D"))) + len(list(root.glob("*/*.D")))


def _dimension_options(exp: Path, analysis_dir: str, col: str):
    """Distinct values of a column in consolidated_peaks.csv, or None."""
    p = exp / analysis_dir / "consolidated_peaks.csv"
    if not p.exists():
        return None
    try:
        df = pd.read_csv(p)
    except Exception:
        return None
    if col not in df.columns:
        return None
    vals = df[col].dropna().unique().tolist()
    return sorted(vals) if vals else None


# --------------------------------------------------------------------------
# State
# --------------------------------------------------------------------------
def _mark_dirty() -> None:
    st.session_state.dirty = True


def load_into_state(exp: Path) -> None:
    doc = load_config_doc(exp)
    st.session_state.doc = doc
    g = lambda p, d=None: cfg_get(doc, p, d)  # noqa: E731

    # Processing
    st.session_state.w_wavelength = int(g(["processing", "wavelength_nm"], 262))
    st.session_state.w_blank = g(["processing", "blank_folder_name"], "") or ""
    st.session_state.w_datadir = g(["processing", "data_dir"], "") or ""
    st.session_state.w_deconv = bool(g(["processing", "deconvolution", "enabled"], True))
    st.session_state.w_minr2 = float(g(["processing", "deconvolution", "min_r2"], 0.95))
    st.session_state.w_maxcomps = int(g(["processing", "deconvolution", "max_comps"], 5))
    st.session_state.w_models = ", ".join(g(["processing", "deconvolution", "models"], []) or [])
    st.session_state.w_min_height = "" if g(["processing", "peak_detection", "min_height"]) is None else str(g(["processing", "peak_detection", "min_height"]))
    st.session_state.w_min_prom = "" if g(["processing", "peak_detection", "min_prominence"]) is None else str(g(["processing", "peak_detection", "min_prominence"]))
    st.session_state.w_min_width = "" if g(["processing", "peak_detection", "min_width"]) is None else str(g(["processing", "peak_detection", "min_width"]))
    st.session_state.w_distance = "" if g(["processing", "peak_detection", "distance"]) is None else str(g(["processing", "peak_detection", "distance"]))

    # Analysis
    st.session_state.w_excl = bool(g(["analysis", "exclude_unknown"], True))
    st.session_state.w_cv = float(g(["analysis", "cv_warning_threshold"], 15.0))
    st.session_state.w_clamp = bool(g(["analysis", "calibration", "clamp_negative_to_zero"], True))
    st.session_state.w_zero = bool(g(["analysis", "calibration", "zero_area_zero_conc"], True))
    st.session_state.w_pattern = g(["analysis", "sample_name", "pattern"], "")
    st.session_state.w_compounds_file = g(["analysis", "compounds_file"], "compounds.csv")
    st.session_state.w_standard_file = g(["analysis", "standard_file"], "standard.csv")

    # Plots
    st.session_state.w_ncols = int(g(["plots", "n_cols"], 2))
    st.session_state.w_strains_txt = ", ".join(g(["plots", "plot_strains"], []) or [])
    st.session_state.w_order = ", ".join(g(["plots", "strain_order"], []) or [])
    st.session_state.w_bartp_txt = ", ".join(str(x) for x in (g(["plots", "bar_time_points"], []) or []))
    st.session_state.w_palette = "\n".join(g(["plots", "palette"], []) or [])

    # Output
    st.session_state.w_results_dir = g(["output", "results_dirname"], "results")
    st.session_state.w_analysis_dir = g(["output", "analysis_dirname"], "analysis")

    # CSVs (seed from packaged examples if the folder has none)
    cfile = st.session_state.w_compounds_file
    sfile = st.session_state.w_standard_file
    st.session_state.compounds_df, st.session_state.compounds_seeded = read_csv_or_default(exp / cfile, PKG_COMPOUNDS, COMPOUND_COLS)
    st.session_state.standard_df, st.session_state.standard_seeded = read_csv_or_default(exp / sfile, PKG_STANDARD, STANDARD_COLS)

    for k in ("ed_compounds", "ed_standard", "w_strains_ms", "w_bartp_ms"):
        st.session_state.pop(k, None)
    st.session_state.loaded_exp = str(exp)
    st.session_state.dirty = False


def apply_changes(exp: Path, compounds_df, standard_df, plot_strains, bar_time_points) -> tuple[bool, str]:
    doc = st.session_state.doc
    s = st.session_state
    cfg_set(doc, ["processing", "wavelength_nm"], int(s.w_wavelength))
    cfg_set(doc, ["processing", "blank_folder_name"], s.w_blank or None)
    cfg_set(doc, ["processing", "data_dir"], s.w_datadir or None)
    cfg_set(doc, ["processing", "deconvolution", "enabled"], bool(s.w_deconv))
    cfg_set(doc, ["processing", "deconvolution", "min_r2"], float(s.w_minr2))
    cfg_set(doc, ["processing", "deconvolution", "max_comps"], int(s.w_maxcomps))
    cfg_set(doc, ["processing", "deconvolution", "models"], _csv_list(s.w_models))
    cfg_set(doc, ["processing", "peak_detection", "min_height"], _num_or_none(s.w_min_height, float))
    cfg_set(doc, ["processing", "peak_detection", "min_prominence"], _num_or_none(s.w_min_prom, float))
    cfg_set(doc, ["processing", "peak_detection", "min_width"], _num_or_none(s.w_min_width, int))
    cfg_set(doc, ["processing", "peak_detection", "distance"], _num_or_none(s.w_distance, int))
    cfg_set(doc, ["analysis", "exclude_unknown"], bool(s.w_excl))
    cfg_set(doc, ["analysis", "cv_warning_threshold"], float(s.w_cv))
    cfg_set(doc, ["analysis", "calibration", "clamp_negative_to_zero"], bool(s.w_clamp))
    cfg_set(doc, ["analysis", "calibration", "zero_area_zero_conc"], bool(s.w_zero))
    cfg_set(doc, ["analysis", "sample_name", "pattern"], s.w_pattern)
    cfg_set(doc, ["analysis", "compounds_file"], s.w_compounds_file)
    cfg_set(doc, ["analysis", "standard_file"], s.w_standard_file)
    cfg_set(doc, ["plots", "n_cols"], int(s.w_ncols))
    cfg_set(doc, ["plots", "plot_strains"], plot_strains or None)
    cfg_set(doc, ["plots", "strain_order"], _csv_list(s.w_order))
    cfg_set(doc, ["plots", "bar_time_points"], bar_time_points or None)
    cfg_set(doc, ["plots", "palette"], _lines_list(s.w_palette))
    cfg_set(doc, ["output", "results_dirname"], s.w_results_dir or "results")
    cfg_set(doc, ["output", "analysis_dirname"], s.w_analysis_dir or "analysis")

    exp.mkdir(parents=True, exist_ok=True)
    buf = StringIO()
    yaml.dump(doc, buf)
    (exp / "hplc_config.yaml").write_text(buf.getvalue(), encoding="utf-8")
    compounds_df.dropna(how="all").to_csv(exp / s.w_compounds_file, index=False)
    standard_df.dropna(how="all").to_csv(exp / s.w_standard_file, index=False)

    try:
        from hplc_gc_pipeline.config import load_config
        load_config(exp)
    except Exception as exc:  # noqa: BLE001
        return False, f"Files written, but config is invalid: {exc}"
    s.dirty = False
    s.compounds_seeded = s.standard_seeded = False
    return True, "All three files written and the config validates."


# --------------------------------------------------------------------------
# Running the pipeline
# --------------------------------------------------------------------------
NOISE = ("Sparse", "spsolve", "flatfit", "UserWarning", "RuntimeWarning")


def run_stage(exp: Path, subcommand: str) -> int:
    cmd = [sys.executable, "-m", "hplc_gc_pipeline.cli", subcommand, str(exp)]
    box = st.empty()
    lines: list[str] = []
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
    assert proc.stdout is not None
    for raw in proc.stdout:
        if any(n in raw for n in NOISE):
            continue
        lines.append(raw.rstrip())
        box.code("\n".join(lines[-500:]), language="text")
    proc.wait()
    log = "\n".join(lines)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logs_dir = exp / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / f"run_{subcommand}_{ts}.txt"
    header = f"# hplc {subcommand}  ·  {datetime.now():%Y-%m-%d %H:%M:%S}  ·  exit {proc.returncode}\n\n"
    log_path.write_text(header + log, encoding="utf-8")
    st.session_state.last_log = log
    st.session_state.last_log_path = str(log_path)
    return proc.returncode


def embed_html(path: Path, height: int = 820) -> None:
    components.html(path.read_text(encoding="utf-8"), height=height, scrolling=True)


def output_row(label: str, path: Path) -> None:
    st.markdown(f"**{label}**")
    st.caption(f"{path} · {_file_meta(path)}")
    with st.container(horizontal=True):
        if st.button("Open", icon=":material/open_in_new:", key=f"open_{label}"):
            webbrowser.open(path.resolve().as_uri())
        st.download_button("Download", data=path.read_bytes(), file_name=path.name,
                           mime="text/html", icon=":material/download:", key=f"dl_{label}")


# --------------------------------------------------------------------------
# UI
# --------------------------------------------------------------------------
st.set_page_config(page_title="HPLC pipeline", page_icon=":material/science:", layout="wide")
st.title("HPLC pipeline")

with st.sidebar:
    st.header("Experiment folder")
    exp_str = st.text_input(
        "Experiment folder path", value=st.session_state.get("loaded_exp", initial_exp_dir()),
        label_visibility="collapsed",
        help="Folder holding your .D data (or a Data/ subfolder) plus the config and CSVs.",
    )
    if st.button("Load / reload", type="primary", width="stretch", icon=":material/refresh:"):
        if exp_str and Path(exp_str).is_dir():
            load_into_state(Path(exp_str))
            st.success("Loaded.")
        else:
            st.error("Not a directory.")

if "doc" not in st.session_state and exp_str and Path(exp_str).is_dir():
    load_into_state(Path(exp_str))

if "doc" not in st.session_state:
    st.info("Enter an experiment folder in the sidebar and click **Load / Reload**.")
    st.stop()

exp = Path(st.session_state.loaded_exp)
analysis_dir = st.session_state.get("w_analysis_dir", "analysis")
results_dir = st.session_state.get("w_results_dir", "results")

with st.sidebar:
    st.divider()
    st.subheader("Status")
    data_root = exp / st.session_state.w_datadir if st.session_state.w_datadir else exp
    for fname in (st.session_state.w_compounds_file, st.session_state.w_standard_file, "hplc_config.yaml"):
        mark = ":material/check_circle:" if (exp / fname).exists() else ":material/radio_button_unchecked:"
        st.write(f"{mark} {fname}")
    st.write(f":material/science: `.D` injections found: **{_count_d_folders(data_root)}**")
    if st.session_state.get("dirty"):
        st.warning("Unsaved edits — click **Apply changes**.")

tab_params, tab_run, tab_results = st.tabs(
    ["Parameters", "Run", "Results"], on_change="rerun"
)

# ==========================================================================
# PARAMETERS
# ==========================================================================
with tab_params:
    with st.container(border=True):
        st.subheader("Processing · Stage 1")
        c1, c2, c3 = st.columns(3)
        c1.number_input("Wavelength (nm)", key="w_wavelength", step=1, on_change=_mark_dirty)
        c2.text_input("Blank .D folder", key="w_blank",
                      placeholder="e.g. 091-0202.D (blank = none)", on_change=_mark_dirty)
        c3.text_input("Data subfolder", key="w_datadir",
                      placeholder="e.g. Data (blank = root)", on_change=_mark_dirty)
        c4, c5, c6 = st.columns(3)
        c4.checkbox("Deconvolution enabled", key="w_deconv", on_change=_mark_dirty)
        c5.number_input("Min R² to accept fit", key="w_minr2", step=0.01, format="%.2f", on_change=_mark_dirty)
        c6.number_input("Max components / peak", key="w_maxcomps", step=1, on_change=_mark_dirty)

    with st.container(border=True):
        st.subheader("Analysis · Stage 2")
        a1, a2, a3, a4 = st.columns(4)
        a1.checkbox("Exclude peaks outside RT windows", key="w_excl", on_change=_mark_dirty)
        a2.number_input("CV%% warning threshold", key="w_cv", step=1.0, on_change=_mark_dirty)
        a3.checkbox("Clamp negative conc. to 0", key="w_clamp", on_change=_mark_dirty)
        a4.checkbox("Zero area → 0 µM (fix)", key="w_zero", on_change=_mark_dirty)
        st.text_input("Sample-name regex", key="w_pattern", on_change=_mark_dirty,
                      help='Named groups: time, strain, replicate. Default: "24, s1, A".')

    with st.container(border=True):
        st.subheader("Plots · samples to include")
        strain_opts = _dimension_options(exp, analysis_dir, "Strain")
        time_opts = _dimension_options(exp, analysis_dir, "Time_Point_h")
        pc1, pc2 = st.columns(2)
        with pc1:
            if strain_opts:
                default = [x for x in _csv_list(st.session_state.w_strains_txt) or [] if x in strain_opts]
                plot_strains = st.multiselect("Strains to plot (empty = all)", strain_opts,
                                              default=default, key="w_strains_ms", on_change=_mark_dirty)
            else:
                plot_strains = _csv_list(st.text_input("Strains to plot", key="w_strains_txt",
                                          placeholder="all (comma-separated)", on_change=_mark_dirty,
                                          help="Run Stage 2 once to get a checklist of detected strains."))
        with pc2:
            if time_opts:
                t_default = [int(x) for x in _csv_list(st.session_state.w_bartp_txt) or [] if str(x).isdigit()]
                t_default = [x for x in t_default if x in [int(o) for o in time_opts]]
                bar_time_points = st.multiselect("Bar-chart time points (empty = all)",
                                                 [int(o) for o in time_opts], default=t_default,
                                                 key="w_bartp_ms", on_change=_mark_dirty)
            else:
                txt = st.text_input("Bar-chart time points", key="w_bartp_txt",
                                    placeholder="all (e.g. 0, 24, 50)", on_change=_mark_dirty)
                bar_time_points = [int(x) for x in (_csv_list(txt) or []) if str(x).strip().isdigit()] or None
        pc3, pc4 = st.columns(2)
        pc3.number_input("Grid columns", key="w_ncols", step=1, min_value=1, on_change=_mark_dirty)
        pc4.text_input("Strain order (colour order)", key="w_order",
                       placeholder="alphabetical", on_change=_mark_dirty)

    with st.expander("Advanced — peak detection, models, palette, file names & output folders"):
        st.caption("Leave a peak-detection field blank to use MOCCA2's own default.")
        d1, d2, d3, d4 = st.columns(4)
        d1.text_input("min_height", key="w_min_height", on_change=_mark_dirty)
        d2.text_input("min_prominence", key="w_min_prom", on_change=_mark_dirty)
        d3.text_input("min_width", key="w_min_width", on_change=_mark_dirty)
        d4.text_input("distance", key="w_distance", on_change=_mark_dirty)
        st.text_input("Deconvolution models (in order)", key="w_models",
                      placeholder="FraserSuzuki, Gaussian", on_change=_mark_dirty)
        f1, f2 = st.columns(2)
        f1.text_input("compounds file name", key="w_compounds_file", on_change=_mark_dirty)
        f2.text_input("standard file name", key="w_standard_file", on_change=_mark_dirty)
        o1, o2 = st.columns(2)
        o1.text_input("results folder name", key="w_results_dir", on_change=_mark_dirty)
        o2.text_input("analysis folder name", key="w_analysis_dir", on_change=_mark_dirty)
        st.text_area("Colour palette (one hex per line)", key="w_palette", height=120, on_change=_mark_dirty)

    with st.container(border=True):
        st.subheader("compounds.csv — retention-time windows")
        if st.session_state.get("compounds_seeded"):
            st.caption("No compounds file in the folder yet — showing example values; apply to save.")
        compounds_df = st.data_editor(st.session_state.compounds_df, num_rows="dynamic",
                                      key="ed_compounds", on_change=_mark_dirty)
        st.subheader("standard.csv — calibration points")
        if st.session_state.get("standard_seeded"):
            st.caption("No standard file in the folder yet — showing example values; apply to save.")
        standard_df = st.data_editor(st.session_state.standard_df, num_rows="dynamic",
                                     key="ed_standard", on_change=_mark_dirty)

    if st.button("Apply changes", type="primary", icon=":material/save:"):
        ok, msg = apply_changes(exp, compounds_df, standard_df, plot_strains, bar_time_points)
        (st.success if ok else st.error)(msg)

# ==========================================================================
# RUN
# ==========================================================================
with tab_run:
    st.caption(f"Experiment folder: `{exp}`")
    if st.session_state.get("dirty"):
        st.warning("Unsaved parameter edits — go to **Parameters** and click **Apply changes** first.")
    with st.container(horizontal=True):
        go_process = st.button("Process (Stage 1)", width="stretch", icon=":material/play_arrow:")
        go_analyze = st.button("Analyze (Stage 2)", width="stretch", icon=":material/play_arrow:")
        go_run = st.button("Run both", type="primary", width="stretch", icon=":material/play_arrow:")

    if go_process or go_analyze or go_run:
        sub = "process" if go_process else "analyze" if go_analyze else "run"
        with st.spinner(f"Running `hplc {sub}` …"):
            code = run_stage(exp, sub)
        (st.success if code == 0 else st.error)(f"`hplc {sub}` finished (exit {code}).")
        st.caption(f"Log saved to `{st.session_state.last_log_path}`")
    elif st.session_state.get("last_log"):
        with st.expander("Last run log"):
            st.code(st.session_state.last_log, language="text")

    if st.session_state.get("last_log_path") and Path(st.session_state.last_log_path).exists():
        lp = Path(st.session_state.last_log_path)
        st.download_button("Download last log (.txt)", data=lp.read_bytes(),
                           file_name=lp.name, mime="text/plain", icon=":material/download:")

# ==========================================================================
# RESULTS
# ==========================================================================
# Only render (and read/embed the large HTML files) when the Results tab is
# actually open — avoids loading multi-MB files on every rerun of other tabs.
if tab_results.open:
  with tab_results:
    dashboard = exp / analysis_dir / "analysis_plots.html"
    overlay = exp / results_dir / "chromatogram_overlay_interactive.html"
    gallery = exp / results_dir / "chromatogram_gallery.html"

    if not any(p.exists() for p in (dashboard, overlay, gallery)):
        st.info("No result files yet — run a stage on the **Run** tab.")
    else:
        st.subheader("Analysis dashboard")
        if dashboard.exists():
            output_row("Analysis dashboard", dashboard)
            embed_html(dashboard, height=900)
        else:
            st.caption("Run **Analyze** to generate the dashboard.")

        st.divider()
        st.subheader("Chromatogram overlay (interactive)")
        if overlay.exists():
            output_row("Interactive overlay", overlay)
            embed_html(overlay, height=650)
        else:
            st.caption("Run **Process** to generate the overlay.")

        if gallery.exists():
            st.divider()
            st.subheader("Per-injection chromatogram gallery")
            output_row("Chromatogram gallery", gallery)
            with st.expander("Preview gallery inline"):
                embed_html(gallery, height=820)
