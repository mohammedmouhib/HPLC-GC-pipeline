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
PKG_GC_COMPOUNDS = HERE / "gc_compounds_default.csv"
COMPOUND_COLS = ["Compound", "RT_low", "RT_high", "Notes"]
STANDARD_COLS = ["Compound", "Area_Integral", "concentration"]
GC_COMPOUND_COLS = ["Compound", "RT_low", "RT_high", "quantifier_mz",
                    "qualifier_mz", "ion_ratio_tol", "calibrate_as", "Notes"]
GC_SAMPLE_PATTERN_DEFAULT = r"^\s*(?P<time>\d+)\s*_\s*(?P<strain>[^_]+?)\s*_\s*(?P<replicate>[A-Za-z])\s*$"
GC_STANDARD_PATTERN_DEFAULT = r"^(?P<conc>\d+(?:\.\d+)?)\s*(?P<unit>[a-zA-Zµ]*M)_(?P<compound>[^_]+)"


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


def _has_any_results(exp: Path) -> bool:
    candidates = [
        exp / "gc_analysis" / "analysis_plots.html",
        exp / "analysis" / "analysis_plots.html",
        exp / "combined_analysis" / "analysis_plots.html",
    ]
    return any(p.exists() for p in candidates)


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

    # GC modality (only present when the config has a `gc:` block)
    st.session_state.w_gc_enabled = "gc" in doc
    st.session_state.w_gc_datadir = g(["gc", "data_dir"], "") or ""
    st.session_state.w_gc_compounds_file = g(["gc", "compounds_file"], "gc_compounds.csv")
    st.session_state.w_gc_quant_channel = g(["gc", "processing", "quant_channel"], "eic")
    st.session_state.w_gc_cal_source = g(["gc", "calibration", "source"], "injections")
    st.session_state.w_gc_std_pattern = g(["gc", "calibration", "standard_pattern"], GC_STANDARD_PATTERN_DEFAULT)
    st.session_state.w_gc_std_file = g(["gc", "calibration", "standard_file"], "gc_standard.csv")
    st.session_state.w_gc_dilution = float(g(["gc", "calibration", "dilution_factor"], 1.0))
    st.session_state.w_gc_force_origin = bool(g(["gc", "calibration", "force_through_origin"], False))
    st.session_state.w_gc_clamp = bool(g(["gc", "calibration", "clamp_negative_to_zero"], True))
    st.session_state.w_gc_require_ratio = bool(g(["gc", "calibration", "require_ion_ratio_pass"], False))
    st.session_state.w_gc_pattern = g(["gc", "sample_name", "pattern"], GC_SAMPLE_PATTERN_DEFAULT)
    st.session_state.w_gc_cv = float(g(["gc", "cv_warning_threshold"], 15.0))

    # CSVs (seed from packaged examples if the folder has none)
    cfile = st.session_state.w_compounds_file
    sfile = st.session_state.w_standard_file
    st.session_state.compounds_df, st.session_state.compounds_seeded = read_csv_or_default(exp / cfile, PKG_COMPOUNDS, COMPOUND_COLS)
    st.session_state.standard_df, st.session_state.standard_seeded = read_csv_or_default(exp / sfile, PKG_STANDARD, STANDARD_COLS)
    gcfile = st.session_state.w_gc_compounds_file
    st.session_state.gc_compounds_df, st.session_state.gc_compounds_seeded = read_csv_or_default(exp / gcfile, PKG_GC_COMPOUNDS, GC_COMPOUND_COLS)

    for k in ("ed_compounds", "ed_standard", "ed_gc_compounds", "w_strains_ms", "w_bartp_ms"):
        st.session_state.pop(k, None)
    st.session_state.loaded_exp = str(exp)
    st.session_state.dirty = False


def apply_changes(exp: Path, compounds_df, standard_df, gc_compounds_df,
                  plot_strains, bar_time_points) -> tuple[bool, str]:
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

    # GC modality: write the `gc:` block when enabled, drop it when not.
    if s.w_gc_enabled:
        cfg_set(doc, ["gc", "data_dir"], s.w_gc_datadir or None)
        cfg_set(doc, ["gc", "compounds_file"], s.w_gc_compounds_file)
        cfg_set(doc, ["gc", "processing", "quant_channel"], s.w_gc_quant_channel)
        cfg_set(doc, ["gc", "calibration", "source"], s.w_gc_cal_source)
        cfg_set(doc, ["gc", "calibration", "standard_pattern"], s.w_gc_std_pattern)
        cfg_set(doc, ["gc", "calibration", "standard_file"], s.w_gc_std_file)
        cfg_set(doc, ["gc", "calibration", "dilution_factor"], float(s.w_gc_dilution))
        cfg_set(doc, ["gc", "calibration", "force_through_origin"], bool(s.w_gc_force_origin))
        cfg_set(doc, ["gc", "calibration", "clamp_negative_to_zero"], bool(s.w_gc_clamp))
        cfg_set(doc, ["gc", "calibration", "require_ion_ratio_pass"], bool(s.w_gc_require_ratio))
        cfg_set(doc, ["gc", "sample_name", "pattern"], s.w_gc_pattern)
        cfg_set(doc, ["gc", "cv_warning_threshold"], float(s.w_gc_cv))
    elif "gc" in doc:
        del doc["gc"]

    exp.mkdir(parents=True, exist_ok=True)
    buf = StringIO()
    yaml.dump(doc, buf)
    (exp / "hplc_config.yaml").write_text(buf.getvalue(), encoding="utf-8")
    compounds_df.dropna(how="all").to_csv(exp / s.w_compounds_file, index=False)
    standard_df.dropna(how="all").to_csv(exp / s.w_standard_file, index=False)
    if s.w_gc_enabled:
        gc_compounds_df.dropna(how="all").to_csv(exp / s.w_gc_compounds_file, index=False)

    try:
        from hplc_gc_pipeline.config import load_config
        load_config(exp)
    except Exception as exc:  # noqa: BLE001
        return False, f"Files written, but config is invalid: {exc}"
    s.dirty = False
    s.compounds_seeded = s.standard_seeded = s.gc_compounds_seeded = False
    return True, "Files written and the config validates."


# --------------------------------------------------------------------------
# Running the pipeline
# --------------------------------------------------------------------------
NOISE = ("Sparse", "spsolve", "flatfit", "UserWarning", "RuntimeWarning")


def run_stage(exp: Path, subcommand: str, extra_args: list[str] | None = None) -> int:
    cmd = [sys.executable, "-m", "hplc_gc_pipeline.cli", subcommand, str(exp)] + (extra_args or [])
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
# Page config
# --------------------------------------------------------------------------
st.set_page_config(page_title="HPLC / GC-MS pipeline", page_icon=":material/science:", layout="wide")

# --------------------------------------------------------------------------
# Sidebar
# --------------------------------------------------------------------------
with st.sidebar:
    st.header("Experiment folder")
    exp_str = st.text_input(
        "Experiment folder path", value=st.session_state.get("loaded_exp", initial_exp_dir()),
        label_visibility="collapsed",
        placeholder="/path/to/your/experiment",
        help="Folder containing your Agilent .D injection files (or subfolders of them).",
    )
    if st.button("Load / reload", type="primary", width="stretch", icon=":material/refresh:"):
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

        # Contextual next-step hint
        st.space("small")
        if not _config_ok:
            st.info("Go to **Run** → **Initialize** to generate a starter config.", icon=":material/arrow_forward:")
        elif not _results_ok:
            st.info("Go to **Run** → **Run both** to process your data.", icon=":material/arrow_forward:")
        else:
            st.caption(":material/check_circle: Ready — open the **Results** tab.")

# --------------------------------------------------------------------------
# Welcome screen (no folder loaded yet)
# --------------------------------------------------------------------------
if "doc" not in st.session_state and exp_str and Path(exp_str).is_dir():
    load_into_state(Path(exp_str))

if "doc" not in st.session_state:
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
                "Enter your experiment folder path in the sidebar and click **Load / Reload**. "
                "The folder should contain Agilent `.D` injection subfolders."
            )
    with col2:
        with st.container(border=True):
            st.markdown(":material/auto_awesome: **Step 2 — Initialize** *(new experiments)*")
            st.caption(
                "Open the **Run** tab and click **Initialize**. "
                "The pipeline scans your `.D` files, detects peaks and m/z values, "
                "and writes a starter config automatically."
            )
    with col3:
        with st.container(border=True):
            st.markdown(":material/play_arrow: **Step 3 — Run & explore**")
            st.caption(
                "Click **Run both** to process all injections and build calibration curves. "
                "Results appear in the **Results** tab as interactive Bokeh dashboards."
            )

    st.space("large")
    st.caption(
        "Already have a configured experiment? Just enter the folder path above — "
        "the pipeline will pick up any existing config automatically.",
        text_alignment="center",
    )
    st.stop()

exp = Path(st.session_state.loaded_exp)
analysis_dir = st.session_state.get("w_analysis_dir", "analysis")
results_dir = st.session_state.get("w_results_dir", "results")

# --------------------------------------------------------------------------
# Tabs
# --------------------------------------------------------------------------
tab_params, tab_run, tab_results = st.tabs(
    [":material/tune: Parameters", ":material/play_arrow: Run", ":material/analytics: Results"],
)

# ==========================================================================
# PARAMETERS
# ==========================================================================
with tab_params:

    # ── GC-MS modality ────────────────────────────────────────────────────
    with st.container(border=True):
        gc_hdr, gc_badge = st.columns([5, 1], vertical_alignment="center")
        with gc_hdr:
            st.subheader("GC-MS modality")
        with gc_badge:
            if st.session_state.w_gc_enabled:
                st.badge("Enabled", color="green")
            else:
                st.badge("Disabled", color="gray")

        st.toggle(
            "Enable GC-MS processing",
            key="w_gc_enabled",
            on_change=_mark_dirty,
            help="Adds a `gc:` block to hplc_config.yaml. Enable when your experiment includes GC-MS .D injections.",
        )
        gc_compounds_df = st.session_state.gc_compounds_df

        if st.session_state.w_gc_enabled:
            g1, g2, g3 = st.columns(3)
            g1.text_input(
                "GC data subfolder", key="w_gc_datadir",
                placeholder="blank = experiment root",
                on_change=_mark_dirty,
                help="Subfolder holding GC .D files. Leave blank to scan the experiment folder directly.",
            )
            g2.segmented_control(
                "Quantification channel", ["eic", "tic", "fid"],
                key="w_gc_quant_channel", on_change=_mark_dirty,
                help="eic = quantifier-ion extracted chromatogram (most selective, recommended). "
                     "tic = total ion current. fid = flame-ionisation detector.",
            )
            g3.number_input(
                "Sample dilution factor", key="w_gc_dilution",
                min_value=0.0, step=1.0, on_change=_mark_dirty,
                help="Multiply sample concentrations by this to recover the undiluted value. "
                     "1:10 dilution → 10. Standards are never scaled.",
            )
            c1, c2 = st.columns(2)
            c1.segmented_control(
                "Calibration source", ["injections", "csv"],
                key="w_gc_cal_source", on_change=_mark_dirty,
                help="injections = build the curve from standard .D folders; "
                     "csv = read a gc_standard.csv you supply.",
            )
            c2.toggle(
                "Force fit through origin",
                key="w_gc_force_origin", on_change=_mark_dirty,
                help="Fit area = slope · conc (no intercept). Recommended for trace GC-MS; "
                     "avoids a non-zero concentration at zero area.",
            )

            with st.expander("Advanced GC settings", icon=":material/settings:"):
                st.text_input(
                    "Sample-name pattern", key="w_gc_pattern", on_change=_mark_dirty,
                    help='Regex with named groups: time, strain, replicate. '
                         'Default matches "50_s11_A" (underscore) or '
                         '"24, s11, A" (comma, as auto-detected by init).',
                )
                st.text_input(
                    "Standard-name pattern", key="w_gc_std_pattern", on_change=_mark_dirty,
                    help="Groups: conc, unit, compound. Matches names like '250uM_34DMS_hexane'.",
                )
                e1, e2 = st.columns(2)
                e1.text_input("GC compounds file", key="w_gc_compounds_file", on_change=_mark_dirty)
                e2.text_input("GC standard CSV (source = csv)", key="w_gc_std_file", on_change=_mark_dirty)
                x1, x2, x3 = st.columns(3)
                x1.toggle("Clamp negative conc. to 0", key="w_gc_clamp", on_change=_mark_dirty)
                x2.toggle("Require ion-ratio confirmation", key="w_gc_require_ratio", on_change=_mark_dirty,
                           help="Report qualifier-ratio failures as 0 µM instead of the fitted value.")
                x3.number_input("CV% warning threshold", key="w_gc_cv", step=1.0, on_change=_mark_dirty)

            if st.session_state.get("gc_compounds_seeded"):
                st.caption("No gc_compounds.csv in this folder yet — showing example values. Run **Initialize** or edit below and click **Apply**.")
            st.markdown("**gc_compounds.csv** — RT windows + quantifier / qualifier m/z")
            st.caption(
                'qualifier_mz format: `"149=0.40;91=0.36"` (m/z=expected_ratio). '
                "ion_ratio_tol = relative tolerance on those ratios. "
                "calibrate_as = another compound whose standards are used to calibrate this row "
                "(valid when the two share the quantifier ion)."
            )
            gc_compounds_df = st.data_editor(
                st.session_state.gc_compounds_df, num_rows="dynamic",
                key="ed_gc_compounds", on_change=_mark_dirty,
            )
        else:
            st.caption(
                "Enable to process GC-MS `.D` data alongside (or instead of) HPLC data. "
                "When both modalities are active, a combined dashboard is written to `combined_analysis/`."
            )

    # ── HPLC modality ─────────────────────────────────────────────────────
    with st.expander(
        "HPLC settings" + (" — not active" if st.session_state.w_gc_enabled and not any(
            k in st.session_state.get("doc", {}) for k in ("processing", "analysis", "hplc")
        ) else ""),
        icon=":material/water_drop:",
        expanded=not st.session_state.w_gc_enabled,
    ):
        st.caption("Configure these when your experiment includes HPLC / DAD data.")
        with st.container(border=True):
            st.subheader("Processing · Stage 1")
            c1, c2, c3 = st.columns(3)
            c1.number_input("Detection wavelength (nm)", key="w_wavelength", step=1, on_change=_mark_dirty)
            c2.text_input(
                "Blank injection (.D folder name)", key="w_blank",
                placeholder="e.g. 091-0202.D (blank = none)", on_change=_mark_dirty,
            )
            c3.text_input(
                "Data subfolder", key="w_datadir",
                placeholder="e.g. Data (blank = experiment root)", on_change=_mark_dirty,
            )
            c4, c5, c6 = st.columns(3)
            c4.toggle("Peak deconvolution", key="w_deconv", on_change=_mark_dirty,
                      help="If disabled, peaks are integrated by trapezoid fallback only.")
            c5.number_input("Min R² to accept fit", key="w_minr2", step=0.01, format="%.2f", on_change=_mark_dirty)
            c6.number_input("Max components per peak", key="w_maxcomps", step=1, on_change=_mark_dirty)

        with st.container(border=True):
            st.subheader("Analysis · Stage 2")
            a1, a2, a3, a4 = st.columns(4)
            a1.toggle("Exclude peaks outside RT windows", key="w_excl", on_change=_mark_dirty)
            a2.number_input("CV% warning threshold", key="w_cv", step=1.0, on_change=_mark_dirty)
            a3.toggle("Clamp negative conc. to 0", key="w_clamp", on_change=_mark_dirty)
            a4.toggle("Zero area → 0 µM", key="w_zero", on_change=_mark_dirty,
                      help="Report zero-area compounds as 0 µM rather than the calibration intercept.")
            st.text_input(
                "Sample-name pattern", key="w_pattern", on_change=_mark_dirty,
                help='Named groups: time, strain, replicate. Default matches "24, s1, A".',
            )

        with st.container(border=True):
            st.subheader("compounds.csv — retention-time windows")
            if st.session_state.get("compounds_seeded"):
                st.caption("No compounds file yet — showing example values; apply to save.")
            compounds_df = st.data_editor(st.session_state.compounds_df, num_rows="dynamic",
                                          key="ed_compounds", on_change=_mark_dirty)
            st.subheader("standard.csv — calibration points")
            if st.session_state.get("standard_seeded"):
                st.caption("No standard file yet — showing example values; apply to save.")
            standard_df = st.data_editor(st.session_state.standard_df, num_rows="dynamic",
                                         key="ed_standard", on_change=_mark_dirty)

    # ── Plots ──────────────────────────────────────────────────────────────
    with st.expander("Plot options", icon=":material/bar_chart:"):
        strain_opts = _dimension_options(exp, analysis_dir, "Strain")
        time_opts = _dimension_options(exp, analysis_dir, "Time_Point_h")
        pc1, pc2 = st.columns(2)
        with pc1:
            if strain_opts:
                default = [x for x in _csv_list(st.session_state.w_strains_txt) or [] if x in strain_opts]
                plot_strains = st.multiselect("Strains to plot (empty = all)", strain_opts,
                                              default=default, key="w_strains_ms", on_change=_mark_dirty)
            else:
                plot_strains = _csv_list(st.text_input(
                    "Strains to plot", key="w_strains_txt",
                    placeholder="all (comma-separated)", on_change=_mark_dirty,
                    help="Run Stage 2 once to get a checklist of detected strains.",
                ))
        with pc2:
            if time_opts:
                t_default = [int(x) for x in _csv_list(st.session_state.w_bartp_txt) or [] if str(x).isdigit()]
                t_default = [x for x in t_default if x in [int(o) for o in time_opts]]
                bar_time_points = st.multiselect(
                    "Bar-chart time points (empty = all)",
                    [int(o) for o in time_opts], default=t_default,
                    key="w_bartp_ms", on_change=_mark_dirty,
                )
            else:
                txt = st.text_input("Bar-chart time points", key="w_bartp_txt",
                                    placeholder="all (e.g. 0, 24, 50)", on_change=_mark_dirty)
                bar_time_points = [int(x) for x in (_csv_list(txt) or []) if str(x).strip().isdigit()] or None
        pc3, pc4 = st.columns(2)
        pc3.number_input("Grid columns", key="w_ncols", step=1, min_value=1, on_change=_mark_dirty)
        pc4.text_input("Strain colour order", key="w_order",
                       placeholder="alphabetical", on_change=_mark_dirty,
                       help="Comma-separated list controls which strain gets which colour.")

    with st.expander("Advanced — peak detection, deconvolution models, palette, file names", icon=":material/build:"):
        st.caption("Leave peak-detection fields blank to use MOCCA2's own defaults.")
        d1, d2, d3, d4 = st.columns(4)
        d1.text_input("min_height", key="w_min_height", on_change=_mark_dirty)
        d2.text_input("min_prominence", key="w_min_prom", on_change=_mark_dirty)
        d3.text_input("min_width", key="w_min_width", on_change=_mark_dirty)
        d4.text_input("distance", key="w_distance", on_change=_mark_dirty)
        st.text_input("Deconvolution models (in order)", key="w_models",
                      placeholder="FraserSuzuki, Gaussian", on_change=_mark_dirty)
        f1, f2 = st.columns(2)
        f1.text_input("HPLC compounds file name", key="w_compounds_file", on_change=_mark_dirty)
        f2.text_input("HPLC standard file name", key="w_standard_file", on_change=_mark_dirty)
        o1, o2 = st.columns(2)
        o1.text_input("Results folder name", key="w_results_dir", on_change=_mark_dirty)
        o2.text_input("Analysis folder name", key="w_analysis_dir", on_change=_mark_dirty)
        st.text_area("Colour palette (one hex per line)", key="w_palette", height=120, on_change=_mark_dirty)

    st.space("small")
    if st.button("Apply changes", type="primary", icon=":material/save:"):
        ok, msg = apply_changes(exp, compounds_df, standard_df, gc_compounds_df,
                                plot_strains, bar_time_points)
        (st.success if ok else st.error)(msg)

# ==========================================================================
# RUN
# ==========================================================================
with tab_run:
    st.caption(f"Experiment: `{exp}`")

    # ── Step 1: Initialize ────────────────────────────────────────────────
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

    # ── Step 2: Run pipeline ──────────────────────────────────────────────
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

# ==========================================================================
# RESULTS
# ==========================================================================
with tab_results:
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
                "Go to the **Run** tab and click **Run both** to process your injections "
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
