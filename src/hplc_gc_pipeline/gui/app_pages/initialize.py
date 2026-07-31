"""Initialize page — create config files for a new experiment."""

from __future__ import annotations

import shutil
from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import (
    GC_STANDARD_PATTERN_DEFAULT,
    PKG_COMPOUNDS,
    PKG_STANDARD,
    PKG_GC_COMPOUNDS,
    PKG_TEMPLATE,
    _file_meta,
    _pick_folder,
    load_into_state,
    run_stage,
)

exp = Path(st.session_state.loaded_exp)
_ro = not st.session_state.get("edit_mode", True)

config_path = exp / "hplc_config.yaml"

st.caption(f"Experiment: `{exp}`")

if _ro:
    st.info(
        "View only — loading existing config files still works in read-only mode. "
        "Toggle **Edit mode** in the sidebar to auto-detect or copy config files.",
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
    "Run once per new experiment to generate `hplc_config.yaml` and all supporting CSV files. "
    "Then review and refine compound names and RT windows in **Parameters**."
)

# ---------------------------------------------------------------------------
# Pre-scan: inspect the experiment folder without touching any files
# ---------------------------------------------------------------------------

def _scan_experiment(folder: Path) -> dict:
    """Fast pathlib scan — no GC/HPLC loading, just folder structure."""
    hplc_dirs: dict[str, int] = {}  # subdir name -> count of HPLC .D folders
    gc_dirs: dict[str, int] = {}

    def _classify(d: Path) -> str:
        """'gc' | 'hplc' | 'other' for a .D folder."""
        if not (d.is_dir() and d.name.endswith(".D")):
            return "other"
        return "gc" if (d / "data.ms").exists() else "hplc"

    def _scan_dir(directory: Path, label: str) -> None:
        hplc_n = sum(1 for p in directory.iterdir() if _classify(p) == "hplc")
        gc_n = sum(1 for p in directory.iterdir() if _classify(p) == "gc")
        if hplc_n:
            hplc_dirs[label] = hplc_n
        if gc_n:
            gc_dirs[label] = gc_n

    try:
        _scan_dir(folder, "(root)")
    except PermissionError:
        pass

    for item in sorted(folder.iterdir()):
        if not item.is_dir() or item.name.startswith(".") or item.name.endswith(".D"):
            continue
        try:
            _scan_dir(item, item.name)
        except PermissionError:
            pass

    return {"hplc_dirs": hplc_dirs, "gc_dirs": gc_dirs}


scan = _scan_experiment(exp)
hplc_dirs = scan["hplc_dirs"]
gc_dirs = scan["gc_dirs"]
has_hplc = bool(hplc_dirs)
has_gc = bool(gc_dirs)

# Display detected structure
with st.container(border=True):
    st.markdown("**Detected experiment structure**")

    if not has_hplc and not has_gc:
        st.warning(
            "No `.D` injection folders found in this experiment folder. "
            "Make sure you loaded the correct folder, or use **Copy from previous project** below.",
            icon=":material/warning:",
        )
    else:
        rows = []
        for name, count in sorted(hplc_dirs.items()):
            display = "(experiment root)" if name == "(root)" else f"`{name}/`"
            rows.append(f"| HPLC (DAD) | {display} | {count} injections |")
        for name, count in sorted(gc_dirs.items()):
            display = "(experiment root)" if name == "(root)" else f"`{name}/`"
            rows.append(f"| GC-MS | {display} | {count} injections |")

        table = "| Type | Location | Count |\n|------|----------|-------|\n" + "\n".join(rows)
        st.markdown(table)


# ---------------------------------------------------------------------------
# Step A: Modality
# ---------------------------------------------------------------------------

auto_mode = (
    "HPLC + GC" if (has_hplc and has_gc) else
    "GC only" if has_gc else
    "HPLC only"
)

with st.container(border=True):
    st.markdown("**Step A · Experiment modality**")
    mode = st.segmented_control(
        "Modality",
        options=["HPLC only", "GC only", "HPLC + GC"],
        default=auto_mode,
        key="init_mode",
        label_visibility="collapsed",
        help="Auto-selected from detected data. Override if needed.",
        disabled=_ro,
    )

include_hplc = mode in ("HPLC only", "HPLC + GC")
include_gc = mode in ("GC only", "HPLC + GC")

# ---------------------------------------------------------------------------
# Step B: Source
# ---------------------------------------------------------------------------

with st.container(border=True):
    st.markdown("**Step B · Config source**")
    source = st.segmented_control(
        "Config source",
        options=["Auto-detect from .D files", "Copy from a previous project",
                 "Read from this project"],
        default="Auto-detect from .D files",
        key="init_source",
        label_visibility="collapsed",
        help="Auto-detect scans .D files and writes a ready-to-run config. "
             "Copy reuses files from another experiment. "
             "Read from this project loads config files already in this folder into the form.",
    )

# ---------------------------------------------------------------------------
# Auto-detect flow
# ---------------------------------------------------------------------------
if source == "Auto-detect from .D files":

    with st.container(border=True):
        if include_gc:
            st.markdown("**GC-MS peak detection settings**")
            ic1, ic2, ic3 = st.columns(3)
            init_rt_min = ic1.number_input(
                "RT scan start (min)", value=1.5, min_value=0.0, step=0.5,
                help="Peaks before this time are ignored (skips solvent front).",
                disabled=_ro,
            )
            init_rt_max = ic2.number_input(
                "RT scan end (min)", value=30.0, min_value=1.0, step=5.0,
                disabled=_ro,
            )
            init_rt_margin = ic3.number_input(
                "RT window margin (min)", value=0.08, min_value=0.0, step=0.01, format="%.2f",
                help="Added to each side of a detected peak edge when writing the RT window.",
                disabled=_ro,
            )
        else:
            st.info(
                "HPLC-only: a template `hplc_config.yaml` will be written with the detected "
                "data subfolder. Edit compound RT windows in **Parameters** afterwards.",
                icon=":material/info:",
            )

        # Show files that will be written
        st.markdown("**Files that will be written:**")
        file_rows = [
            ("`hplc_config.yaml`",
             "Main config with HPLC + GC blocks" if (include_hplc and include_gc) else
             "Main config (HPLC)" if include_hplc else "Main config (GC)"),
        ]
        if include_gc:
            file_rows.append(("`gc_compounds.csv`",
                               "GC compound table — auto-detected peaks; rename Peak_N to real names"))
        if include_hplc:
            exists_c = (exp / "compounds.csv").exists()
            exists_s = (exp / "standard.csv").exists()
            file_rows.append((
                "`compounds.csv`",
                "Seeded with placeholder HPLC compounds — edit RT windows for your compounds"
                + (" *(already exists — will not overwrite)*" if exists_c else ""),
            ))
            file_rows.append((
                "`standard.csv`",
                "Seeded with example HPLC calibration data — replace with your standards"
                + (" *(already exists — will not overwrite)*" if exists_s else ""),
            ))

        for fname, desc in file_rows:
            st.markdown(f"- {fname} — {desc}")

        files_exist = config_path.exists() or (include_gc and (exp / "gc_compounds.csv").exists())
        if files_exist:
            st.warning(
                "`hplc_config.yaml`" +
                (" and `gc_compounds.csv`" if include_gc and (exp / "gc_compounds.csv").exists() else "") +
                " already exist and will be overwritten.",
                icon=":material/warning:",
            )

        if st.button("Initialize", icon=":material/auto_awesome:", type="primary", disabled=_ro):
            mode_arg = "both" if (include_hplc and include_gc) else ("hplc" if include_hplc else "gc")
            extra = [
                "--force",
                "--mode", mode_arg,
                "--rt-min", str(init_rt_min if include_gc else 1.5),
                "--rt-max", str(init_rt_max if include_gc else 30.0),
                "--rt-margin", str(init_rt_margin if include_gc else 0.08),
            ]
            with st.spinner("Running `hplc init` …"):
                code = run_stage(exp, "init", extra_args=extra)

            if code == 0:
                # Seed HPLC CSVs if not written by the CLI (extra safety)
                for src, fname in ((PKG_COMPOUNDS, "compounds.csv"), (PKG_STANDARD, "standard.csv")):
                    dst = exp / fname
                    if not dst.exists() and src.exists():
                        shutil.copy2(src, dst)

                written = ["`hplc_config.yaml`"]
                if include_gc:
                    written.append("`gc_compounds.csv`")
                if include_hplc:
                    written += ["`compounds.csv`", "`standard.csv`"]

                st.toast(
                    "Config files written: " + ", ".join(written),
                    icon=":material/check:",
                )
                load_into_state(exp)
                st.rerun()
            else:
                st.error(f"`hplc init` failed (exit {code}). See log above.")

# ---------------------------------------------------------------------------
# Copy from previous project
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Read from this project
# ---------------------------------------------------------------------------
elif source == "Read from this project":
    with st.container(border=True):
        st.caption(
            "Load config files already present in this experiment folder into the form. "
            "No files are written or overwritten — use this after editing `hplc_config.yaml` "
            "or CSV files externally (text editor, Excel, etc.)."
        )

        known = ["hplc_config.yaml"] + sorted(p.name for p in exp.glob("*.csv"))
        present = [f for f in known if (exp / f).exists()]

        if not present:
            st.info(
                "No config files found in this experiment folder. "
                "Use **Auto-detect** to generate them, or **Copy from a previous project** "
                "to reuse an existing configuration.",
                icon=":material/info:",
            )
        else:
            st.markdown("**Config files found in this folder:**")
            for fname in present:
                st.markdown(f"- `{fname}` — {_file_meta(exp / fname)}")

            if st.button(
                "Load into form", icon=":material/upload_file:", type="primary",
            ):
                load_into_state(exp)
                st.toast(
                    f"Loaded {len(present)} config file(s) into the form.",
                    icon=":material/check:",
                )
                st.rerun()

# ---------------------------------------------------------------------------
# Calibration quick-setup (shown once a config exists)
# ---------------------------------------------------------------------------
if config_path.exists() and not _ro:
    with st.expander(
        "Optional: calibration quick-setup",
        icon=":material/science:",
        expanded=False,
    ):
        from io import StringIO
        from hplc_gc_pipeline.gui.common import cfg_set
        from ruamel.yaml import YAML as _YAML

        _yaml = _YAML()
        _yaml.preserve_quotes = True

        st.caption(
            "Answer a few questions to pre-configure calibration in `hplc_config.yaml`. "
            "You can also change these settings at any time in **Parameters**."
        )

        has_stds = st.toggle(
            "This experiment includes standard/calibration injections in the .D data",
            key="init_has_stds",
            value=False,
        )

        if has_stds:
            st.markdown(
                "Standard injections must follow a naming convention so the pipeline can "
                "detect them automatically. The default pattern matches names like "
                "`250uM_pCA_hexane` — concentration, unit, compound name, separated by `_`."
            )
            std_pattern = st.text_input(
                "Standard-name pattern (regex)",
                value=GC_STANDARD_PATTERN_DEFAULT,
                key="init_std_pattern",
                help="Named groups required: conc, unit, compound.",
            )
            cal_source = st.segmented_control(
                "Calibration source",
                options=["injections", "both"],
                default="injections",
                key="init_cal_source",
                help="injections = use only these .D standard injections.  "
                     "both = pool with any pre-existing standard.csv as well.",
            )
        else:
            cal_source = "csv"
            std_pattern = GC_STANDARD_PATTERN_DEFAULT

        dilution = st.number_input(
            "Sample dilution factor",
            value=1.0, min_value=0.0, step=1.0,
            key="init_dilution",
            help="If samples were diluted before injection, enter the inverse dilution "
                 "(e.g. 1:10 dilution → 10). Standards are never scaled.",
        )

        if st.button(
            "Apply calibration settings to config",
            icon=":material/save:",
            type="primary",
        ):
            try:
                with open(config_path, "r", encoding="utf-8") as fh:
                    doc = _yaml.load(fh)
                cfg_set(doc, ["analysis", "calibration", "source"], cal_source)
                cfg_set(doc, ["analysis", "calibration", "dilution_factor"], float(dilution))
                if has_stds:
                    cfg_set(doc, ["analysis", "calibration", "standard_pattern"], std_pattern)
                buf = StringIO()
                _yaml.dump(doc, buf)
                config_path.write_text(buf.getvalue(), encoding="utf-8")
                load_into_state(exp)
                st.toast("Calibration settings saved to hplc_config.yaml.", icon=":material/check:")
            except Exception as exc:
                st.error(f"Failed to write config: {exc}")
