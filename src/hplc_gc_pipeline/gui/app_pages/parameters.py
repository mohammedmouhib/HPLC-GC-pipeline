"""Parameters page — edit hplc_config.yaml + compounds + standards."""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from hplc_gc_pipeline.gui.common import (
    GC_SAMPLE_PATTERN_DEFAULT,
    GC_STANDARD_PATTERN_DEFAULT,
    _csv_list,
    _dimension_options,
    _mark_dirty,
    apply_changes,
    load_into_state,
)

exp = Path(st.session_state.loaded_exp)
analysis_dir = st.session_state.get("w_analysis_dir", "analysis")

# ── GC-MS modality ────────────────────────────────────────────────────────
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
                     'Default matches "50_s11_A" (underscore) or "24, s11, A" (comma).',
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
            st.caption("No gc_compounds.csv yet — showing example values. Run **Initialize** or edit below.")
        st.markdown("**gc_compounds.csv** — RT windows + quantifier / qualifier m/z")
        st.caption(
            'qualifier_mz format: `"149=0.40;91=0.36"` (m/z=expected_ratio). '
            "ion_ratio_tol = relative tolerance. "
            "calibrate_as = compound whose standards calibrate this row."
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

# ── HPLC modality ─────────────────────────────────────────────────────────
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

# ── Plots ──────────────────────────────────────────────────────────────────
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
