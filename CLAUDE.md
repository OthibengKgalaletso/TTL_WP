# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

TTL_pipeline is an in-house untargeted LC-MS metabolomics pipeline. Each stage is a standalone, independently testable Python module under `modules/<module_name>/`. The end goal is to wire the modules together as a Nextflow/nf-core pipeline with a GUI front end. Design decisions should keep that in mind: every module has a CLI with an exit code, plus a programmatic API (an inputs dataclass and a class with `run()`) that a GUI or workflow engine can call directly.

Current modules:
- `modules/Upload_and_validator_module/`: `metabolomics_validator.py` checks inputs before anything runs downstream.
- `modules/visualization/`: `visualization_module.py` handles QC, normalization/scaling and multivariate/univariate stats for **unannotated** data.

## Environment and commands

All modules share one environment, Python 3.11. The root `requirements.txt` covers every module. When a module needs a new dependency, append it there instead of creating a per-module requirements file.

```bash
conda create -n TTL_pipeline python=3.11 pandas pytest -y && conda activate TTL_pipeline
pip install -r requirements.txt
```

Commands are run from **inside the module folder**. The `run_check.sh` scripts `cd` to their own directory, and tests insert the module's parent folder into `sys.path`.

```bash
cd modules/Upload_and_validator_module
./run_check.sh [data_folder]           # default: sample_data/valid_files
pytest -v                              # whole suite
pytest tests/test_validator.py::test_name -v   # single test

cd modules/visualization
./run_check.sh path/to/folder [ATTRIBUTE_col] [none|sum|median] [none|pareto|auto]
# --venn-min-fraction is only on the Python CLI (visualization_module.py), not run_check.sh
pytest -v
```

For example, `./run_check.sh ../Upload_and_validator_module/sample_data/valid_files ATTRIBUTE_Tissue` runs the visualization module on the real sample dataset. Output goes to `<folder>/../viz_output/`.

No linter or formatter is configured.

## Architecture

### Shared conventions across modules
- **Module independence is intentional.** Small helpers such as `_norm`, `_find_col`, `_find_col_containing`, `_filename_variants` and `_read_table` are copied into each module rather than imported between them. Keep each module runnable on its own, and don't add cross-module imports.
- **Flexible column matching.** Columns are found by normalized name (lowercase, alphanumeric only) against a list of accepted aliases, because real exports differ (MZmine, GNPS, XCMS). Examples: `feature_id`/`id`/`peak_id`/`row ID`, and `filename` for the metadata column.
- **Loose filename matching.** `_filename_variants` compares names with and without the extension, so metadata `filename` values, mzML basenames and data-matrix sample columns can match whether or not they include `.mzML`.
- **Grouping columns follow the GNPS `ATTRIBUTE_*` convention.** For example, `ATTRIBUTE_Tissue` in the sample data has groups MCL and MCR.
- **Tests build fixtures on the fly.** They write small mzML/CSV/MGF files in `tmp_path` and call the Python API, not the CLI. `sample_data/` is for running the module by hand, not for the tests.

### Upload & Validation module
`MetabolomicsValidator(InputSet).run()` runs four steps in order. Each step adds `Issue(step, severity, file, message)` records to a list:
1. Upload
2. File types
3. Integrity & format
4. Compatibility (do the files agree with each other)

Severity works as follows:
- `ERROR` fails validation.
- `WARNING` does not fail it. For example, a blank or NA matrix cell is a warning, but a non-numeric string is an error.

`ValidationReport` renders as text or JSON, and the CLI exits 0 or 1 so a workflow can stop on bad input. Inputs:
- **Required:** mzML files, sample metadata, data matrix, and the **full** MZmine 3/4 feature table ("Export to CSV file (modular)", `id`/`mz`/`rt`/`datafile:<file>:area|height` columns).
- **Optional:** MGF file (MZmine export only), annotation list. These are checked only when supplied.
- The feature table and MGF are MZmine-only by design (SPEC VAL-EX-01). The GNPS/IIMN quant table is rejected on purpose, so don't add aliases to make it pass.

**`modules/Upload_and_validator_module/SPEC.md` is the source of truth for this module's behaviour.** It lists every rule with a stable ID (`VAL-xx-nn`) and its severity, the API/JSON/exit-code contract, the mapping from tests to requirements, and the known limitations. Before changing validator behaviour, follow its §9 change-control steps: spec, then tests, then code, then docs, then verify, with a semver bump and a changelog entry. Also keep the docstring at the top of `metabolomics_validator.py` in sync with SPEC §5.

### Visualization module
**`modules/visualization/SPEC.md` is the source of truth for this module's behaviour** (rule IDs `VIZ-xx-nn`). Follow its §9 change control, the same as the validator, and keep `__version__` in `visualization_module.py` in step with the spec version.

`VisualizationModule(VisualizationInputs).run()` calls its stages in a fixed order:
1. feature summary
2. chromatograms
3. normalization summary
4. distribution shape
5. PCA (2D, plus 3D when there are ≥4 samples) and correlation
6. PLS-DA
7. OPLS-DA
8. volcano plot
9. feature overlap (Venn for 2–3 groups, UpSet-style chart for 4+)

Each stage saves PNGs and adds a text section to the report. `build_report()` then writes everything to `analysis_report.txt` in `output_dir`. Matplotlib uses the headless `Agg` backend.

Key data-handling rules. They're deliberate, so don't "fix" them:
- **Raw data stays raw for QC.** The feature summary and chromatograms use unmodified values.
- **Normalization is per sample** (`sum`/`median`/`none`, default `sum`). It applies to every cross-sample analysis.
- **Scaling is per feature** (`pareto`/`auto`/`none`, default `pareto`). It applies **only** to PCA, PLS-DA and OPLS-DA. It is **not** applied to the volcano plot, so fold changes stay real intensity ratios, or to the correlation matrix, which is Pearson on log1p of normalized, unscaled values.
- **Missing values.** NA is treated as 0 for PCA, correlation, PLS-DA and OPLS-DA. The volcano plot instead skips NAs (`nanmean`, `nan_policy="omit"`).
- **Group requirements.** OPLS-DA and the volcano plot need exactly 2 groups. Every supervised analysis needs at least 2 matched samples per group. Feature overlap needs ≥2 groups but only 1 sample per group. If these conditions aren't met, the analysis is skipped with a note in the report rather than raising an error.
- **Grouping column resolution:**
  - `--group-column` or `VisualizationInputs.group_column`, if set, wins (any metadata column is accepted).
  - Otherwise, a single `ATTRIBUTE_*` column is used automatically.
  - With several candidates, an interactive run (TTY) asks you to pick one. A non-interactive run skips the supervised analyses. Nothing should ever block waiting for input when there's no TTY.
  - The prompt is a stand-in for a future GUI picker.
  - PCA colouring and the feature-summary group counts use the same resolved column, falling back to a looser `Class`/`Group`/`Condition`/`*attribute*` match only when nothing resolves.
- **Group colours** come from one tab10 map in sorted group order (`_group_color_map`), so a group has the same colour in every plot. The volcano plot colours significant features by the group they're higher in.
- **Feature overlap.** A feature is "captured" by a group when it's detected (> 0, not NA, raw values) in ≥ `venn_min_fraction` (default 0.5) of that group's samples. The Venn layout is fixed, not area-proportional. Membership goes to `feature_overlap.csv`.
- **Data matrix shape.** It accepts either a bare `feature_id, sample...` matrix or the MZmine shape `row ID, row m/z, row retention time, samples...`. The m/z and RT columns are excluded from the per-sample statistics.
- **Chromatograms.** TIC/BPC use MS1 scans only. They prefer the precomputed cvParams `MS:1000285` and `MS:1000505`, and fall back to decoding the binary arrays (32- or 64-bit, zlib-compressed or not) when either is missing. RT is reported in seconds; minute-valued scan start times are converted.
- **Statistics:**
  - PLS-DA: `PLSRegression` on one-hot group labels.
  - OPLS-DA: one orthogonal component (Trygg & Wold 2002).
  - Volcano plot: Welch's t-test, Benjamini-Hochberg q-values, and a default significance overlay of p < 0.05 and |log2FC| > 1.
- **Deferred:** an interactive MS1/MS2 spectra viewer, and interactive versions of the Venn/UpSet and 3D PCA plots (static PNGs for now), all planned for the GUI phase.

## Sample data

`modules/Upload_and_validator_module/sample_data/valid_files/` holds a real dataset:
- 12 MeOH positive-mode mzML runs: tissues MCL and MCR × biological reps R1–R3 × technical reps 1a/2a. Each file is about 55 MB.
- Matching metadata, data matrix, GNPS/IIMN quant table, MGF and annotation files.

`invalid_files_examples/` holds files that are meant to fail validation. The visualization module also reads this folder.

## Known doc drift

The READMEs still refer to `modules/upload_validation` and `sample_data/valid/` / `sample_data/invalid/`. The actual paths are `modules/Upload_and_validator_module` and `sample_data/valid_files/` / `sample_data/invalid_files_examples/`. The root README's module table doesn't list the visualization module yet.
