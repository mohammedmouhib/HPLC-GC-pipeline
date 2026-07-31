"""Shared helpers, state functions, and pipeline runners for the GUI pages."""

from __future__ import annotations

import http.server
import os
import posixpath
import subprocess
import sys
import threading
import urllib.parse
import webbrowser
from datetime import datetime
from io import StringIO
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

yaml = YAML()
yaml.preserve_quotes = True

HERE = Path(__file__).parent
PKG_TEMPLATE = HERE / "config_template.yaml"
PKG_COMPOUNDS = HERE / "compounds_default.csv"
PKG_STANDARD = HERE / "standard_default.csv"
PKG_GC_COMPOUNDS = HERE / "gc_compounds_default.csv"

COMPOUND_COLS = ["Compound", "RT_low", "RT_high", "calibrate_as", "Notes"]
STANDARD_COLS = ["Compound", "Area_Integral", "concentration"]
GC_COMPOUND_COLS = ["Compound", "RT_low", "RT_high", "quantifier_mz",
                    "qualifier_mz", "ion_ratio_tol", "calibrate_as", "Notes"]

GC_SAMPLE_PATTERN_DEFAULT = r"^\s*(?P<time>\d+)\s*_\s*(?P<strain>[^_]+?)\s*_\s*(?P<replicate>[A-Za-z])\s*$"
GC_STANDARD_PATTERN_DEFAULT = r"^(?P<conc>\d+(?:\.\d+)?)\s*(?P<unit>[a-zA-Zµ]*M)_(?P<compound>[^_]+)"

NOISE = ("Sparse", "spsolve", "flatfit", "UserWarning", "RuntimeWarning")


# --------------------------------------------------------------------------
# Local static-file server (serves experiment directory for HTML previews)
# --------------------------------------------------------------------------
class _RootedHandler(http.server.SimpleHTTPRequestHandler):
    """Serves files from a dynamically-updated root (the experiment folder)."""
    root: str = ""

    def translate_path(self, path: str) -> str:
        path = path.split("?", 1)[0].split("#", 1)[0]
        path = urllib.parse.unquote(path, errors="surrogatepass")
        path = posixpath.normpath(path)
        parts = [p for p in path.split("/") if p and p not in (".", "..")]
        return os.path.join(self.root or os.getcwd(), *parts)

    def log_message(self, *args) -> None:
        pass  # suppress per-request logs


_fs: dict = {"httpd": None, "port": 0}


def _ensure_file_server(root: str) -> int:
    """Start the daemon file server once; update its root on each call.

    Returns the bound port, or 0 if startup failed.
    """
    _RootedHandler.root = root
    if _fs["httpd"] is not None:
        return _fs["port"]

    from http.server import ThreadingHTTPServer
    for candidate in (8502, 0):   # prefer fixed port, fall back to OS-assigned
        try:
            httpd = ThreadingHTTPServer(("127.0.0.1", candidate), _RootedHandler)
            break
        except OSError:
            continue
    else:
        return 0

    port = httpd.server_address[1]
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    _fs["httpd"] = httpd
    _fs["port"] = port
    return port


# --------------------------------------------------------------------------
# Misc helpers
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


def read_csv_or_default(
    path: Path, pkg_default: Path, cols: list[str]
) -> tuple[pd.DataFrame, bool]:
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


def _pick_folder() -> str | None:
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes("-topmost", 1)
        folder = filedialog.askdirectory(title="Select experiment folder")
        root.destroy()
        return folder or None
    except Exception:
        return None


def _has_any_results(exp: Path) -> bool:
    return any(
        (exp / d / "analysis_plots.html").exists()
        for d in ("gc_analysis", "analysis", "combined_analysis")
    )


# --------------------------------------------------------------------------
# Session-state management
# --------------------------------------------------------------------------
def _mark_dirty() -> None:
    st.session_state.dirty = True


def load_into_state(exp: Path) -> None:
    doc = load_config_doc(exp)
    st.session_state.doc = doc
    g = lambda p, d=None: cfg_get(doc, p, d)  # noqa: E731

    st.session_state.w_wavelength = int(g(["processing", "wavelength_nm"], 262))
    st.session_state.w_blank = g(["processing", "blank_folder_name"], "") or ""
    _dd = g(["processing", "data_dir"], "") or ""
    st.session_state.w_datadir = ", ".join(_dd) if isinstance(_dd, list) else _dd
    st.session_state.w_deconv = bool(g(["processing", "deconvolution", "enabled"], True))
    st.session_state.w_minr2 = float(g(["processing", "deconvolution", "min_r2"], 0.95))
    st.session_state.w_maxcomps = int(g(["processing", "deconvolution", "max_comps"], 5))
    st.session_state.w_models = ", ".join(g(["processing", "deconvolution", "models"], []) or [])
    st.session_state.w_min_height = "" if g(["processing", "peak_detection", "min_height"]) is None else str(g(["processing", "peak_detection", "min_height"]))
    st.session_state.w_min_prom = "" if g(["processing", "peak_detection", "min_prominence"]) is None else str(g(["processing", "peak_detection", "min_prominence"]))
    st.session_state.w_min_width = "" if g(["processing", "peak_detection", "min_width"]) is None else str(g(["processing", "peak_detection", "min_width"]))
    st.session_state.w_distance = "" if g(["processing", "peak_detection", "distance"]) is None else str(g(["processing", "peak_detection", "distance"]))

    st.session_state.w_hplc_enabled = any(k in doc for k in ("processing", "analysis", "hplc"))
    st.session_state.w_excl = bool(g(["analysis", "exclude_unknown"], True))
    st.session_state.w_cv = float(g(["analysis", "cv_warning_threshold"], 15.0))
    st.session_state.w_clamp = bool(g(["analysis", "calibration", "clamp_negative_to_zero"], True))
    st.session_state.w_zero = bool(g(["analysis", "calibration", "zero_area_zero_conc"], True))
    st.session_state.w_cal_source = g(["analysis", "calibration", "source"], "csv")
    st.session_state.w_std_pattern = g(["analysis", "calibration", "standard_pattern"], GC_STANDARD_PATTERN_DEFAULT)
    st.session_state.w_dilution = float(g(["analysis", "calibration", "dilution_factor"], 1.0))
    st.session_state.w_force_origin = bool(g(["analysis", "calibration", "force_through_origin"], False))
    st.session_state.w_pattern = g(["analysis", "sample_name", "pattern"], "")
    st.session_state.w_compounds_file = g(["analysis", "compounds_file"], "compounds.csv") or "compounds.csv"
    st.session_state.w_standard_file = g(["analysis", "standard_file"], "standard.csv") or "standard.csv"

    st.session_state.w_ncols = int(g(["plots", "n_cols"], 2))
    st.session_state.w_strains_txt = ", ".join(g(["plots", "plot_strains"], []) or [])
    st.session_state.w_order = ", ".join(g(["plots", "strain_order"], []) or [])
    st.session_state.w_bartp_txt = ", ".join(str(x) for x in (g(["plots", "bar_time_points"], []) or []))
    st.session_state.w_palette = "\n".join(g(["plots", "palette"], []) or [])

    st.session_state.w_results_dir = g(["output", "results_dirname"], "results")
    st.session_state.w_analysis_dir = g(["output", "analysis_dirname"], "analysis")

    st.session_state.w_gc_enabled = "gc" in doc
    st.session_state.w_gc_datadir = g(["gc", "data_dir"], "") or ""
    st.session_state.w_gc_compounds_file = g(["gc", "compounds_file"], "gc_compounds.csv") or "gc_compounds.csv"
    st.session_state.w_gc_quant_channel = g(["gc", "processing", "quant_channel"], "eic")
    st.session_state.w_gc_cal_source = g(["gc", "calibration", "source"], "injections")
    st.session_state.w_gc_std_pattern = g(["gc", "calibration", "standard_pattern"], GC_STANDARD_PATTERN_DEFAULT)
    st.session_state.w_gc_std_file = g(["gc", "calibration", "standard_file"], "gc_standard.csv") or "gc_standard.csv"
    st.session_state.w_gc_dilution = float(g(["gc", "calibration", "dilution_factor"], 1.0))
    st.session_state.w_gc_force_origin = bool(g(["gc", "calibration", "force_through_origin"], False))
    st.session_state.w_gc_clamp = bool(g(["gc", "calibration", "clamp_negative_to_zero"], True))
    st.session_state.w_gc_require_ratio = bool(g(["gc", "calibration", "require_ion_ratio_pass"], False))
    st.session_state.w_gc_pattern = g(["gc", "sample_name", "pattern"], GC_SAMPLE_PATTERN_DEFAULT)
    st.session_state.w_gc_cv = float(g(["gc", "cv_warning_threshold"], 15.0))

    cfile = st.session_state.w_compounds_file
    sfile = st.session_state.w_standard_file
    st.session_state.compounds_df, st.session_state.compounds_seeded = read_csv_or_default(
        exp / cfile, PKG_COMPOUNDS, COMPOUND_COLS)
    st.session_state.standard_df, st.session_state.standard_seeded = read_csv_or_default(
        exp / sfile, PKG_STANDARD, STANDARD_COLS)
    gcfile = st.session_state.w_gc_compounds_file
    st.session_state.gc_compounds_df, st.session_state.gc_compounds_seeded = read_csv_or_default(
        exp / gcfile, PKG_GC_COMPOUNDS, GC_COMPOUND_COLS)

    for k in ("ed_compounds", "ed_standard", "ed_gc_compounds", "w_strains_ms", "w_bartp_ms"):
        st.session_state.pop(k, None)
    # Stage the mode default via a non-widget key so it can be applied before the
    # toggle widget renders in app.py.  Only reset when the folder actually changes
    # so that Copy / Initialize on the same folder preserve the user's choice.
    if str(exp) != st.session_state.get("loaded_exp"):
        st.session_state._reset_edit_mode = not _has_any_results(exp)
    st.session_state.loaded_exp = str(exp)
    st.session_state.dirty = False
    _ensure_file_server(str(exp))


def apply_changes(
    exp: Path,
    compounds_df: pd.DataFrame,
    standard_df: pd.DataFrame,
    gc_compounds_df: pd.DataFrame,
    plot_strains,
    bar_time_points,
) -> tuple[bool, str]:
    doc = st.session_state.doc
    s = st.session_state

    if s.get("w_hplc_enabled", True):
        cfile = s.w_compounds_file or "compounds.csv"
        sfile = s.w_standard_file or "standard.csv"
        cfg_set(doc, ["processing", "wavelength_nm"], int(s.w_wavelength))
        cfg_set(doc, ["processing", "blank_folder_name"], s.w_blank or None)
        _dd_parts = [p.strip() for p in (s.w_datadir or "").split(",") if p.strip()]
        cfg_set(doc, ["processing", "data_dir"],
                None if not _dd_parts else
                (_dd_parts[0] if len(_dd_parts) == 1 else _dd_parts))
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
        cfg_set(doc, ["analysis", "calibration", "source"], s.w_cal_source or "csv")
        cfg_set(doc, ["analysis", "calibration", "standard_pattern"], s.w_std_pattern or "")
        cfg_set(doc, ["analysis", "calibration", "dilution_factor"], float(s.w_dilution))
        cfg_set(doc, ["analysis", "calibration", "force_through_origin"], bool(s.w_force_origin))
        cfg_set(doc, ["analysis", "sample_name", "pattern"], s.w_pattern)
        cfg_set(doc, ["analysis", "compounds_file"], cfile)
        cfg_set(doc, ["analysis", "standard_file"], sfile)
    else:
        cfile = s.get("w_compounds_file") or "compounds.csv"
        sfile = s.get("w_standard_file") or "standard.csv"
        for k in ("processing", "analysis"):
            if k in doc:
                del doc[k]

    cfg_set(doc, ["plots", "n_cols"], int(s.w_ncols))
    cfg_set(doc, ["plots", "plot_strains"], plot_strains or None)
    cfg_set(doc, ["plots", "strain_order"], _csv_list(s.w_order))
    cfg_set(doc, ["plots", "bar_time_points"], bar_time_points or None)
    cfg_set(doc, ["plots", "palette"], _lines_list(s.w_palette))
    cfg_set(doc, ["output", "results_dirname"], s.w_results_dir or "results")
    cfg_set(doc, ["output", "analysis_dirname"], s.w_analysis_dir or "analysis")

    if s.w_gc_enabled:
        gcfile = s.w_gc_compounds_file or "gc_compounds.csv"
        gcstdfile = s.w_gc_std_file or "gc_standard.csv"
        cfg_set(doc, ["gc", "data_dir"], s.w_gc_datadir or None)
        cfg_set(doc, ["gc", "compounds_file"], gcfile)
        cfg_set(doc, ["gc", "processing", "quant_channel"], s.w_gc_quant_channel)
        cfg_set(doc, ["gc", "calibration", "source"], s.w_gc_cal_source)
        cfg_set(doc, ["gc", "calibration", "standard_pattern"], s.w_gc_std_pattern)
        cfg_set(doc, ["gc", "calibration", "standard_file"], gcstdfile)
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
    if s.get("w_hplc_enabled", True):
        compounds_df.dropna(how="all").to_csv(exp / cfile, index=False)
        standard_df.dropna(how="all").to_csv(exp / sfile, index=False)
    if s.w_gc_enabled:
        gc_compounds_df.dropna(how="all").to_csv(exp / gcfile, index=False)

    try:
        from hplc_gc_pipeline.config import load_config
        load_config(exp)
    except Exception as exc:  # noqa: BLE001
        return False, f"Files written, but config is invalid: {exc}"
    s.dirty = False
    s.compounds_seeded = s.standard_seeded = s.gc_compounds_seeded = False
    return True, "Files written and the config validates."


# --------------------------------------------------------------------------
# Pipeline runner
# --------------------------------------------------------------------------
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


# --------------------------------------------------------------------------
# Result helpers
# --------------------------------------------------------------------------
def _bokeh_content_height(path: Path, n_cols: int = 3) -> int:
    """Estimate the pixel height of a Bokeh dashboard HTML file.

    Counts ``"height": N`` entries in the JSON payload (one per figure),
    divides by n_cols to get the number of rows, and multiplies by the
    plot height. Falls back to 2000 if parsing fails.
    """
    import math, re
    try:
        html = path.read_text(encoding="utf-8", errors="replace")
        heights = [int(h) for h in re.findall(r'"height":(\d+)', html)]
        if not heights:
            return 2000
        plot_h = max(set(heights), key=heights.count)   # modal value
        n_plots = len(heights)
        n_rows = math.ceil(n_plots / max(n_cols, 1))
        return n_rows * plot_h + 300    # +300 for toolbars, titles, spacing
    except Exception:
        return 2000


def embed_html(path: Path, height: int = 900, scrolling: bool = False) -> None:
    """Embed an HTML result file. Uses the local file server when available
    (proper Bokeh sizing); falls back to inline injection."""
    port = _fs.get("port", 0)
    if port:
        try:
            exp = Path(st.session_state.loaded_exp)
            rel = path.relative_to(exp).as_posix()
            components.iframe(
                src=f"http://127.0.0.1:{port}/{rel}",
                height=height,
                scrolling=scrolling,
            )
            return
        except (ValueError, KeyError):
            pass
    components.html(path.read_text(encoding="utf-8"), height=height, scrolling=True)


def result_header(path: Path) -> None:
    """Open-in-browser + download buttons with file metadata."""
    st.caption(f"`{path.name}` · {_file_meta(path)}")
    with st.container(horizontal=True):
        if st.button("Open in browser", icon=":material/open_in_new:",
                     key=f"open_{path.name}"):
            webbrowser.open(path.resolve().as_uri())
        st.download_button("Download", data=path.read_bytes(), file_name=path.name,
                           mime="text/html", icon=":material/download:",
                           key=f"dl_{path.name}")
