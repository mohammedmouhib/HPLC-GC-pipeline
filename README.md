# HPLC_GC_pipeline

A reusable, config-driven pipeline for HPLC data from metabolic engineering
experiments. It processes Agilent Chemstation `.D` injection folders, detects
and integrates chromatographic peaks with [MOCCA2](https://github.com/bayer-group/MOCCA),
assigns peaks to named compounds by retention-time window, applies linear
calibration curves to convert areas to µM concentrations, and generates
publication-quality plots.

The program is fully separated from experiment data: you point it at any
experiment folder, and every tunable value lives in that folder's
`hplc_config.yaml`. Nothing is hardcoded to a particular dataset or machine.

> The name reserves room for future GC support; today the implemented pipeline
> is the HPLC/MOCCA2 path described here.

## Install

Requires Python ≥ 3.10.

```bash
git clone <your-repo-url> HPLC_GC_pipeline
cd HPLC_GC_pipeline
python -m venv .venv && source .venv/bin/activate
pip install -e .          # installs deps and the `hplc` command
```

## Quick start

1. Create an **experiment folder** containing your Agilent `.D` folders (either
   directly, or grouped one level deep in run subfolders). If your raw data
   lives in a nested subfolder (e.g. `Data/<run>/*.D`) that you'd rather not
   move, keep the config/CSVs/outputs at the experiment root and set
   `processing.data_dir: Data` in the config (see below).
2. Copy the config template and reference tables into it:

   ```bash
   cp config/config.example.yaml   /path/to/experiment/hplc_config.yaml
   cp examples/compounds.csv       /path/to/experiment/
   cp examples/standard.csv        /path/to/experiment/
   ```

3. Edit `hplc_config.yaml` (at minimum, set `wavelength_nm` and, if you have
   one, `blank_folder_name`) and edit `compounds.csv` / `standard.csv` for your
   compounds and calibration.
4. Run the pipeline:

   ```bash
   hplc run /path/to/experiment      # both stages
   # or run them separately:
   hplc process /path/to/experiment  # Stage 1 only
   hplc analyze /path/to/experiment  # Stage 2 only (needs Stage 1 output)
   ```

Use `--config /some/other.yaml` to point at a config file outside the
experiment folder. On a mistake (missing/invalid config, bad path) the command
prints a single-line `Error: …`; add `--debug` to any subcommand to see the
full traceback.

## Web GUI (optional)

A browser GUI wraps the whole workflow — edit the config and both CSVs, write
them with one button, run any stage, and open the result HTMLs — without
touching the backend (it edits the same three files and calls the same CLI).

```bash
pip install -e ".[gui]"          # one-time: installs streamlit + ruamel.yaml
hplc gui /path/to/experiment     # opens http://localhost:8501 in your browser
```

The experiment folder is optional (`hplc gui`) — you can set it in the sidebar.
`--port` changes the port. It binds to `localhost` only. YAML comments are
preserved when the GUI saves the config.

## Input interface

Everything the pipeline needs lives in the experiment folder:

| File | Purpose |
|---|---|
| `hplc_config.yaml` | All parameters (see `config/config.example.yaml`; every field documented) |
| `compounds.csv` | `Compound, RT_low, RT_high [, Notes]` — retention-time windows (min) |
| `standard.csv` | `Compound, Area_Integral, concentration` — calibration points (µM) |
| `*.D/` folders | Agilent injections (with `sample.xml` and `DAD*.ch`) |

**Sample-name convention.** Sample names (read from each `sample.xml`) are
expected as `"<time_h>, <strain>, <replicate>"`, e.g. `"24, s1, A"`. Injections
that don't match (blanks, standards) are skipped with a warning. The pattern is
a config field, so a different convention can be plugged in without code changes.

## Output interface

Outputs are written **back into the experiment folder**, never into the program
directory:

```
<experiment>/
  results/                                  # Stage 1
    peak_results.csv                        # one row per detected peak
    sample_summary.csv                      # peak count + blank flag per injection
    chromatogram_traces.csv                 # baseline-corrected signal, long format
    plots/*.png                             # one plot per injection
    chromatogram_gallery.html               # self-contained PNG viewer
    overlay_plot.png                        # static overlay of all traces
    chromatogram_overlay_interactive.html   # interactive Bokeh overlay
  analysis/                                 # Stage 2
    assigned_peaks.csv                      # peaks tagged with a compound
    consolidated_peaks.csv                  # (injection × compound) grid, zero-filled, + Concentration_uM
    replicate_summary.csv                   # mean ± SD, CV% per (strain × time × compound)
    time_series_<strain>.csv                # pivoted time course per strain
    analysis_plots.html                     # full Bokeh dashboard (area + µM)
```

## How it works

**Stage 1 — `hplc process`**
1. Discover `.D` folders (flat or one level of run subfolders) and read each
   `sample.xml` for the sample name.
2. Load each chromatogram via MOCCA2 at the configured wavelength.
3. Baseline-correct → detect peaks → deconvolve (MOCCA2 auto, then the
   configured models as fallbacks).
4. Integrate: MOCCA2 component integrals when deconvolution succeeded, otherwise
   a manual trapezoid over the baseline-corrected signal. The
   `Integration_Method` column records which path was used.
5. Write the peak table, traces, and plots.

**Stage 2 — `hplc analyze`**
1. Parse sample names → strain / time / replicate; pool duplicate replicate
   letters across runs.
2. Assign peaks to compounds by RT window.
3. Consolidate: sum sub-peaks within a window, then zero-fill so every compound
   appears for every injection.
4. Build per-compound linear calibration curves and convert area → µM.
5. Compute replicate mean ± SD and CV%, flagging high-CV groups.
6. Export CSVs and the Bokeh dashboard.

## Configuration

See `config/config.example.yaml` — every field is commented with what it
controls and its default. The most common knobs:

| Field | Effect |
|---|---|
| `processing.wavelength_nm` | Detector wavelength to integrate |
| `processing.blank_folder_name` | `.D` folder of the blank (`null` if none) |
| `processing.data_dir` | Subfolder holding raw `.D` data (`null` = experiment root) |
| `analysis.exclude_unknown` | Drop peaks outside all RT windows |
| `analysis.cv_warning_threshold` | CV% above which a warning prints |
| `plots.plot_strains` | Subset of strains to plot (`null` = all) |
| `plots.strain_order` | Controls colour/order assignment |
| `plots.bar_time_points` | Restrict bar charts to certain time points |

## Known fixes vs. the original scripts

This project is a refactor of two standalone scripts. Behaviour is preserved
except for these deliberate corrections (documented so results are traceable):

- **Zero-area concentration.** Undetected (zero-filled) compounds now report
  `0 µM` instead of the calibration intercept `b` (which produced artefacts like
  pCA ≈ 15 µM at zero area). Toggle with `analysis.calibration.zero_area_zero_conc`;
  set `false` to reproduce the original numbers exactly.
- **Duplicate `process_single_chrom`.** The original defined this function twice
  and Python silently used the second (deconvolution) version. Only that
  behaviour is kept.
- **Hardcoded paths removed.** All absolute paths in the old `__main__` blocks
  are replaced by the experiment-folder argument and the YAML config.
- **Robust array access.** Overlay plots use MOCCA2's known API
  (`chrom.time` / `chrom.data` / `chrom.baseline`) instead of probing candidate
  attribute names.

### Data caveats (not code — check your inputs)

- `34DMS` and `3M4HS` share identical calibration points in the example
  `standard.csv` (likely a placeholder). Verify or replace with real measurements.
- Compounds with an RT window but no calibration data (e.g. `34DHS`, `4HS`,
  `4MS`, and the slash-named `3M4HS / 3H4MS`, which doesn't match the `3M4HS`
  calibration key) produce `NaN` concentrations. Stage 2 warns about these at
  startup.

## Optional tool

`tools/agilent_to_csv.py` bulk-converts binary `.ch` files to CSV (UTF-16, for
Excel). It is **not** needed by the pipeline — MOCCA2 reads `.D` folders
natively — and takes the dataset path as an argument:

```bash
python tools/agilent_to_csv.py /path/to/folder/with/.D/subfolders
```

## Project layout

```
HPLC_GC_pipeline/
├── src/hplc_gc_pipeline/     # the package
│   ├── cli.py                # `hplc` entry point
│   ├── config.py             # YAML load + validation
│   ├── pipeline.py           # stage orchestration
│   ├── agilent.py            # .D discovery + loading
│   ├── processing.py         # MOCCA2 processing + peak table (Stage 1)
│   ├── sample_names.py       # name parsing + replicate pooling (Stage 2)
│   ├── compounds.py          # RT-window assignment
│   ├── consolidate.py        # summation + zero-fill
│   ├── calibration.py        # area → concentration
│   ├── stats.py              # replicate statistics
│   ├── plotting_chromatograms.py  # Stage 1 plots
│   └── plotting_dashboard.py      # Stage 2 dashboard
├── config/config.example.yaml
├── examples/                 # sample compounds.csv + standard.csv
├── tools/agilent_to_csv.py   # optional standalone converter
├── tests/
├── requirements.txt
└── pyproject.toml
```
