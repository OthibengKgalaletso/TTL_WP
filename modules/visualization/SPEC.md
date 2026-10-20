# Visualization (unannotated) module: specification

| | |
|---|---|
| Spec version | 1.2.0 (matches `__version__` in the module) |
| Module file | `visualization_module.py` |
| Status | Describes the code as it stands on 2026-10-20 (module not yet committed) |
| Owner | TTL_pipeline |

This document is the source of truth for the behaviour of the visualization
module. Each rule has a stable ID (`VIZ-xx-nn`) so tests, issues and future
GUI/Nextflow wiring can point at it. If the code and this spec disagree, treat
it as a bug in one of them and resolve it via §9.

---

## 1. Purpose and scope

The module takes a validated, **unannotated** LC-MS metabolomics dataset and
produces QC plots, normalized/scaled multivariate views, supervised group
comparisons and a plain-text report.

In scope:
- Feature summary (counts, missingness, intensity distribution).
- TIC/BPC chromatograms from mzML.
- Sample-wise normalization and feature-wise scaling, with diagnostic plots.
- PCA and sample-sample correlation.
- PLS-DA, OPLS-DA and a volcano plot, grouped by a metadata column.
- 3D PCA (static).
- Feature overlap between groups: a Venn diagram for 2–3 groups, an
  UpSet-style chart for 4+ groups.

Out of scope:
- Input validation. The module assumes its inputs already passed the
  Upload & Validation module. It does not re-check them, and malformed
  input may raise or degrade silently (see §8).
- Annotation-aware analyses. The feature table, MGF and annotation list are
  not read.
- Interactive spectra viewing. The MS1/MS2 viewer is deferred to the GUI
  phase (VIZ-DF-01).
- Interactive plots. The Venn/UpSet and 3D PCA plots are static PNGs for
  now (VIZ-DF-03).

## 2. Inputs

`VisualizationInputs` dataclass:

| Field | Type | Default | Required | Notes |
|---|---|---|---|---|
| `data_matrix_file` | `str` | `None` | **yes** | Feature × sample intensity table. |
| `metadata_file` | `str` | `None` | no | Sample metadata. Without it, supervised analyses are skipped and PCA is uncoloured. |
| `mzml_files` | `list[str]` | `[]` | no | One chromatogram per existing file. |
| `output_dir` | `str` | `"viz_output"` | no | Created if missing. Existing files are overwritten. |
| `group_column` | `str` | `None` | no | Grouping column for PLS-DA/OPLS-DA/volcano. See VIZ-GR-*. |
| `normalization` | `str` | `"sum"` | no | `none` \| `sum` \| `median`. |
| `scaling` | `str` | `"pareto"` | no | `none` \| `pareto` \| `auto`. |
| `venn_min_fraction` | `float` | `0.5` | no | Feature-overlap capture threshold, in (0, 1]. See VIZ-OV-03. *(1.2.0)* |

### 2.1 Data matrix

| ID | Rule |
|---|---|
| VIZ-IN-01 | Read with delimiter sniffing (comma, semicolon or tab), BOM-safe, all cells as strings first. If both read attempts fail, `ValueError` is raised. |
| VIZ-IN-02 | The identifier column is the first column whose normalized name matches `feature_id`, `id`, `peak_id` or `row ID`. If none matches, `ValueError` is raised. |
| VIZ-IN-03 | Descriptive columns matching `row m/z`/`m/z`/`mz` and `row retention time`/`retention time`/`rt` are excluded from the sample columns. |
| VIZ-IN-04 | Every remaining column is a sample column, kept in file order. Values are coerced to numeric. Non-numeric cells become NaN without a warning. |
| VIZ-IN-05 | Both the bare shape (`feature_id, s1, s2, …`) and the MZmine shape (`row ID, row m/z, row retention time, s1, …`) are accepted. |

### 2.2 Metadata

| ID | Rule |
|---|---|
| VIZ-IN-10 | Read with the same reader as VIZ-IN-01. If parsing fails, `ValueError` is raised. |
| VIZ-IN-11 | The filename column matches `filename` or `file_name`. If it is missing, the metadata loads, but every analysis that needs sample-to-group matching treats metadata as unavailable. |
| VIZ-IN-12 | **Display group column** (used for PCA colouring and the feature-summary group counts): the supervised grouping column from VIZ-GR-01…05 when one resolves, so all plots are grouped the same way. Only if nothing resolves does it fall back to the first of `Class`/`Group`/`Condition`, otherwise the first column whose name *contains* `attribute`. *(Changed in 1.1.0; previously always used the fallback rule.)* |
| VIZ-IN-13 | Sample-to-metadata matching uses `_filename_variants`: case-insensitive, with or without extension. A sample column with no match gets the label `None`. |

### 2.3 mzML

| ID | Rule |
|---|---|
| VIZ-IN-20 | Paths that are not existing files are skipped without a note. |
| VIZ-IN-21 | Files are stream-parsed (`iterparse`) and spectra are cleared after use, so memory stays bounded on large files (~55 MB samples). |

## 3. Outputs

All outputs are written to `output_dir`.

| File | Produced by | Condition |
|---|---|---|
| `analysis_report.txt` | `run()` | Always (on a successful run). |
| `missing_values_per_sample.png` | feature summary | Always. |
| `intensity_distribution.png` | feature summary | Always (empty axes if there are no positive values). |
| `chromatogram_<stem>.png` | chromatograms | Per mzML file with ≥1 usable MS1 scan. |
| `normalization_effect.png` | normalization summary | Always. |
| `sample_distribution_gaussian.png` | distribution shape | Always. |
| `feature_distribution_gaussian.png` | distribution shape | Always. |
| `pca.png` | PCA | ≥2 samples. |
| `pca_3d.png` | PCA | ≥4 samples and ≥3 features (3 components). *(1.2.0)* |
| `correlation_heatmap.png`, `sample_correlation_matrix.csv` | correlation | ≥2 samples. |
| `pls_da.png` | PLS-DA | Supervised preconditions met (VIZ-SU-01). |
| `opls_da.png` | OPLS-DA | Supervised preconditions met, plus exactly 2 groups. |
| `volcano_plot.png`, `volcano_plot_results.csv` | volcano | Supervised preconditions met, plus exactly 2 groups. |
| `feature_overlap.csv` | feature overlap | Overlap preconditions met (VIZ-OV-01). *(1.2.0)* |
| `venn_diagram.png` | feature overlap | Same, with 2 or 3 groups. *(1.2.0)* |
| `upset_plot.png` | feature overlap | Same, with 4 or more groups. *(1.2.0)* |

| ID | Rule |
|---|---|
| VIZ-OU-01 | PNGs are saved at 150 dpi with `bbox_inches="tight"` through the headless `Agg` backend. No display is required. |
| VIZ-OU-02 | The report has a fixed header (title, ISO timestamp), then one section per stage in run order, then a list of every generated plot path. |
| VIZ-OU-03 | The feature-summary section lists per-sample intensity stats for at most 10 samples, followed by "… and N more sample(s)". |
| VIZ-OU-04 | `volcano_plot_results.csv` has one row per feature (all features, not only significant ones) with columns `<id column>, log2fc, p_value, q_value`. |
| VIZ-OU-05 | `sample_correlation_matrix.csv` is a square sample × sample matrix, indexed by sample column name. |
| VIZ-OU-06 | `feature_overlap.csv` has one row per feature: `<id column>`, one boolean `captured_<group>` column per group (sorted group order), and `overlap`. *(1.2.0)* |

## 4. Processing pipeline

`VisualizationModule.run()` executes the stages in this fixed order
(VIZ-PL-01):

1. `load()`: check `venn_min_fraction`, then matrix, then metadata, then create `output_dir`, then normalize.
2. `feature_summary()`
3. `chromatograms()`
4. `normalization_summary()`
5. `distribution_shape_summary()`
6. `pca_and_correlation()`
7. `pls_da(group_column)`
8. `opls_da(group_column)`
9. `volcano_plot(group_column)`
10. `feature_overlap(group_column)` *(1.2.0)*
11. `build_report()`, then write `analysis_report.txt`.

Each stage method is public and can be called on its own after `load()`.
Each one returns a `dict` and appends one section to the report.

### 4.1 Which data each stage sees

| Stage | Raw | Normalized | Scaled | NA handling |
|---|---|---|---|---|
| Feature summary | ✔ | | | Counted as missing; excluded from stats. |
| Chromatograms | mzML | | | n/a |
| Normalization summary | ✔ (before) | ✔ (after) | | `skipna` sums. |
| Distribution: samples | ✔ log1p | ✔ log1p | | NaN dropped; negatives clipped to 0. |
| Distribution: features | | ✔ (before) | ✔ (after) | NaN → 0. |
| PCA | | ✔ | ✔ | NaN → 0. |
| Correlation | | ✔ log1p | **no** | NaN → 0. |
| PLS-DA / OPLS-DA | | ✔ | ✔ | NaN → 0. |
| Volcano | | ✔ | **no** | NaN omitted (`nanmean`, `nan_policy="omit"`). |
| Feature overlap | ✔ | | | NaN and 0 both count as not detected. |

## 5. Requirements

### 5.1 Normalization (sample-wise): `VIZ-NO`

| ID | Rule |
|---|---|
| VIZ-NO-01 | Normalization is computed once in `load()` and reused by every cross-sample stage (PCA, correlation, PLS-DA, OPLS-DA, volcano, distribution diagnostics). |
| VIZ-NO-02 | `sum`: each sample is divided by its NaN-skipping column sum, then multiplied by the mean of those sums. After this, all sample totals are equal and within-sample proportions are preserved. |
| VIZ-NO-03 | `median`: same as VIZ-NO-02, using the column median instead of the sum. |
| VIZ-NO-04 | `none`: an unchanged copy is returned. |
| VIZ-NO-05 | A sample whose sum/median is 0 has its divisor treated as NaN, so the whole column becomes NaN. If the mean target is non-finite or 0, the data is returned unchanged. |
| VIZ-NO-06 | Any other method raises `ValueError`. |
| VIZ-NO-07 | Raw (unnormalized) values are kept for the feature summary and chromatograms (QC must show the real data). |

### 5.2 Scaling (feature-wise): `VIZ-SC`

| ID | Rule |
|---|---|
| VIZ-SC-01 | Scaling takes a samples × features array and is applied after NA→0 filling. |
| VIZ-SC-02 | `pareto`: `(x − mean) / sqrt(std)`, with population std (`ddof=0`). A zero std is replaced with 1. |
| VIZ-SC-03 | `auto`: `StandardScaler` (mean-centre, unit variance). |
| VIZ-SC-04 | `none`: mean-centre only. |
| VIZ-SC-05 | Any other method raises `ValueError`. |
| VIZ-SC-06 | Scaling applies to PCA, PLS-DA, OPLS-DA and the feature-distribution diagnostic **only**. It is never applied to the volcano plot (fold changes must stay real intensity ratios) or to the correlation matrix. |
| VIZ-SC-07 | In PLS-DA and OPLS-DA, scaling is fitted on the matched samples only, not on all samples. |

### 5.3 Feature summary: `VIZ-FS`

| ID | Rule |
|---|---|
| VIZ-FS-01 | Reports the feature count (matrix rows), the total missing cells, missing % of all sample cells, and the missing count per sample. |
| VIZ-FS-02 | Per-sample `min`, `max`, `mean`, `median`, `std` come from `DataFrame.describe()` on raw values. |
| VIZ-FS-03 | The intensity histogram is `log10` of all non-NaN values > 0 across all samples, with 60 bins. |
| VIZ-FS-04 | When metadata and a display group column (VIZ-IN-12) exist, the group counts are reported. They count **metadata rows**, not matched samples. |

### 5.4 Chromatograms: `VIZ-CH`

| ID | Rule |
|---|---|
| VIZ-CH-01 | Only spectra with `MS:1000511` (ms level) = 1 and a `MS:1000016` (scan start time) are included. |
| VIZ-CH-02 | TIC comes from `MS:1000285` and BPC from `MS:1000505` when present on the spectrum. |
| VIZ-CH-03 | Fallback: if TIC **or** BPC is absent *(1.1.0; previously only when TIC was absent)*, the intensity array (`MS:1000515`) is decoded from base64 as 32-bit float, or 64-bit if `MS:1000523` is present, zlib-decompressed if `MS:1000574` is present. Missing TIC is then the array sum and missing BPC the array max. |
| VIZ-CH-04 | A scan is kept only if both TIC and BPC end up defined. Points are sorted by RT. |
| VIZ-CH-05 | RT is reported in seconds. A scan start time with `unitAccession="UO:0000031"` or `unitName="minute"` is multiplied by 60. Any other unit (normally `UO:0000010`, second) is taken as seconds. *(1.1.0; previously the unit was ignored.)* |
| VIZ-CH-06 | A file with 0 usable MS1 scans gets no plot. Its report line shows `0 MS1 scans -> None`. |

### 5.5 PCA and correlation: `VIZ-PC`

| ID | Rule |
|---|---|
| VIZ-PC-01 | Requires ≥2 samples. Otherwise the section reads "Fewer than 2 samples — PCA/correlation skipped." |
| VIZ-PC-02 | PCA is fitted on all samples (normalized, NA→0, scaled) with up to 3 components (VIZ-PC-05). The explained-variance ratio is reported. *(Was `min(2, …)` before 1.2.0.)* |
| VIZ-PC-03 | Points are coloured by the display group column (VIZ-IN-12) when it exists. Samples with no metadata match are drawn as grey `x` markers labelled "unmatched". The result dict carries `group_column` and `n_unmatched`, and the report names the column and, when non-zero, the unmatched count. *(1.1.0; previously unmatched samples were silently not drawn.)* |
| VIZ-PC-04 | The correlation matrix is Pearson (`np.corrcoef`) between samples on `log1p` of normalized, NA→0, **unscaled** values. The heatmap uses a fixed colour range of [−1, 1] (`coolwarm`). |
| VIZ-PC-05 | The PCA is fitted once with `n_components = min(3, n_samples − 1, n_features)`. The 2D plot uses PC1 and PC2, which are the same components a 2-component fit would give. `explained_variance_ratio` has up to 3 entries. *(1.2.0)* |
| VIZ-PC-06 | With 3 components, `pca_3d.png` is a static 3D scatter of PC1–PC3, with axes labelled by % variance. It uses the colouring and unmatched-sample marking from VIZ-PC-03, and the result dict carries `plot_3d`. With fewer components there is no 3D plot. The result carries `note_3d`, and the report says "3D PCA skipped: needs ≥4 samples and ≥3 features." The 2D plot is unaffected. *(1.2.0)* |
| VIZ-PC-07 | Group colours in every plot (PCA, 3D PCA, PLS-DA, OPLS-DA, volcano, Venn/UpSet) come from one map: tab10, assigned in `sorted(…, key=str)` group order. *(1.2.0)* |

### 5.6 Grouping column resolution: `VIZ-GR`

Applies to PLS-DA, OPLS-DA and the volcano plot.

| ID | Rule |
|---|---|
| VIZ-GR-01 | If there is no metadata or no filename column, the analysis is skipped: "No sample metadata available." |
| VIZ-GR-02 | An explicit `group_column` (API field or `--group-column`) wins. It matches any metadata column exactly or by normalized name. It does **not** have to be `ATTRIBUTE_*`-prefixed. If nothing matches, the analysis is skipped: "Column '<x>' not found in metadata." |
| VIZ-GR-03 | Otherwise, the candidates are columns whose normalized name starts with `attribute` (GNPS convention). Exactly one candidate is used automatically. |
| VIZ-GR-04 | With more than one candidate and no explicit choice, the module skips with a note that lists them and says to use `--group-column`. |
| VIZ-GR-05 | With zero candidates, it skips: "No ATTRIBUTE_ grouping column found in metadata." |
| VIZ-GR-06 | CLI only: with more than one candidate, no `--group-column`, and stdin a TTY, the user is prompted to pick by number. Enter or an invalid choice leaves the column unset, so VIZ-GR-04 applies. **Without a TTY there is never a prompt.** The module API itself never prompts. |

### 5.7 Supervised preconditions: `VIZ-SU`

| ID | Rule |
|---|---|
| VIZ-SU-01 | After resolution, samples are matched to labels (VIZ-IN-13). Unmatched samples are excluded. The analysis is skipped (returning `{"skipped": True, "reason": …}`, with a "Skipped: …" report section) when: fewer than 2 samples matched; fewer than 2 distinct groups; or any group has fewer than 2 samples. |
| VIZ-SU-02 | Groups are ordered by `sorted(…, key=str)`. That order defines "group A" and "group B" in two-group analyses. |
| VIZ-SU-03 | A skipped analysis never raises. The remaining stages still run. |

### 5.8 PLS-DA: `VIZ-PL`

| ID | Rule |
|---|---|
| VIZ-PL-10 | `PLSRegression` on scaled X against one-hot group labels (`pd.get_dummies`), with `n_components = max(1, min(2, n_samples − 1, n_features))`. Supports ≥2 groups. |
| VIZ-PL-11 | The plot shows LV1 vs LV2 (LV2 = 0 if there is only one component), coloured by group. No cross-validation, R²/Q², VIP or permutation test is computed. |

### 5.9 OPLS-DA: `VIZ-OP`

| ID | Rule |
|---|---|
| VIZ-OP-01 | Requires exactly 2 groups. Otherwise it skips with a note recommending PLS-DA. |
| VIZ-OP-02 | Single orthogonal component (Trygg & Wold 2002). y is coded 0 for group A and 1 for group B, then mean-centred. Steps: `w = Xᵀy/‖·‖`; `t = Xw`; `p = Xᵀt/(tᵀt)`; `w_orth = p − (wᵀp/wᵀw)w`, normalized; `t_orth = Xw_orth`; remove `t_orth p_orthᵀ`; then the predictive score `t_pred` comes from the filtered X. |
| VIZ-OP-03 | If `‖w_orth‖ < 1e-12`, it skips: "Could not separate an orthogonal component…" |
| VIZ-OP-04 | The plot shows the predictive component (x) vs the orthogonal component (y), with zero lines. No model statistics, S-plot or VIP are produced. |

### 5.10 Volcano plot: `VIZ-VO`

| ID | Rule |
|---|---|
| VIZ-VO-01 | Requires exactly 2 groups. |
| VIZ-VO-02 | Uses normalized, **unscaled**, non-imputed values. |
| VIZ-VO-03 | `log2FC = log2((mean_B + ε) / (mean_A + ε))`, with group means from `nanmean`. ε is one tenth of the smallest positive value in either group (`1e-9` if there is none). The direction is B over A, per VIZ-SU-02. |
| VIZ-VO-04 | p-values come from Welch's t-test (`ttest_ind`, `equal_var=False`, `nan_policy="omit"`) per feature. q-values are Benjamini-Hochberg, with NaN p passed through as NaN q and q forced monotone and clipped to [0, 1]. |
| VIZ-VO-05 | A feature is significant if **raw** `p < 0.05` and `|log2FC| > 1`. The thresholds are fixed. q-values are reported but do not drive the overlay. |
| VIZ-VO-06 | `-log10(p)` is computed with p clipped at 1e-300. |
| VIZ-VO-07 | Significant points are coloured by direction. log2FC > 0 means higher in group B and uses B's colour; log2FC < 0 means higher in group A and uses A's colour (VIZ-PC-07). Non-significant points are grey. Legend entries read "higher in <group> (n)". The result dict adds `n_higher_in = {A: n, B: n}`, which sums to `n_significant`, and the report gives both counts. *(1.2.0; previously all significant points were one colour.)* |

### 5.11 Distribution diagnostics: `VIZ-DS`

| ID | Rule |
|---|---|
| VIZ-DS-01 | Sample panel: all non-NaN raw values vs all non-NaN normalized values, pooled, clipped at 0, `log1p`-transformed. |
| VIZ-DS-02 | Feature panel: the exact normalized NA→0 matrix fed to `scale_features()`, before and after scaling, pooled and not log-transformed. |
| VIZ-DS-03 | Each panel is a 60-bin density histogram with a fitted normal curve (when std > 0). The stats are `mean`, `std`, `skew` (only when n > 2) and `n`. An empty input hides the axes and returns `None`. |
| VIZ-DS-04 | These diagnostics are visual only. No normality test is performed. |

### 5.12 Feature overlap: `VIZ-OV` *(1.2.0)*

| ID | Rule |
|---|---|
| VIZ-OV-01 | The grouping column is resolved with VIZ-GR-01…05, then samples are matched with VIZ-IN-13, and unmatched samples are excluded. The stage needs ≥2 groups with ≥1 matched sample each. Unlike VIZ-SU-01, it does **not** need ≥2 samples per group. Otherwise it is skipped with `{"skipped": True, "reason": …}` and a "Skipped: …" report section. It never raises. |
| VIZ-OV-02 | Uses raw values. A cell is *detected* if it is > 0 and not NA. |
| VIZ-OV-03 | A feature is *captured* by a group if (detected samples in the group / group size) ≥ `venn_min_fraction` (default 0.5, with a 1e-12 float tolerance). A value outside (0, 1] raises `ValueError` in `load()` (API) or is an argparse error (CLI, exit 2). |
| VIZ-OV-04 | Each feature belongs to exactly one *exclusive intersection*: the set of groups that capture it. It is labelled `<g> only`, `<g1> & <g2> …` (sorted group order), or `none`. |
| VIZ-OV-05 | 2–3 groups: `venn_diagram.png`, with fixed-layout circles (radius 1) in the group colours. Each of the 3 or 7 regions shows its exclusive count, and each label shows `<group> (<n> captured)`. The circles are **not** area-proportional, and the subtitle says so. |
| VIZ-OV-06 | 4+ groups: `upset_plot.png`. Top: bars for each non-empty exclusive intersection, largest first (ties broken by fewer groups, then name), capped at the 30 largest, with a note giving the number omitted. Bottom: a dot matrix (filled dots in group colours, joined by a line). Left: a bar per group of its total captured. |
| VIZ-OV-07 | The report section "Feature overlap (Venn/UpSet)" lists: the column, the threshold, each group's captured and unique counts, the number shared by all groups, the number captured by no group, a truncation note if any, and the plot and CSV paths. |
| VIZ-OV-08 | The result dict has: `group_column`, `groups`, `min_fraction`, `captured_per_group`, `unique_per_group`, `shared_all`, `none`, `intersections` (label → count, in plot order, excluding `none`), `plot` and `csv`. `sum(intersections) + none` equals the feature count. |

### 5.13 Deferred: `VIZ-DF`

| ID | Item |
|---|---|
| VIZ-DF-01 | Interactive MS1/MS2 spectra viewer and export. Planned for the GUI phase. |
| VIZ-DF-02 | GUI picker for the grouping column. It will set `VisualizationInputs.group_column` directly and replace the VIZ-GR-06 prompt. |
| VIZ-DF-03 | Interactive Venn/UpSet (e.g. clicking a region to list its features) and an interactive (rotatable) 3D PCA. Both are static PNGs for now; the underlying data is already exposed (`feature_overlap.csv`, the result dicts). Planned for the GUI phase. |

## 6. Interface contract

### 6.1 Python API

```python
from visualization_module import VisualizationInputs, VisualizationModule

module = VisualizationModule(VisualizationInputs(data_matrix_file=..., metadata_file=..., ...))
report_text: str = module.run()
```

- `run()` returns the report text and writes `analysis_report.txt`.
- Fatal problems raise `ValueError`: an unparseable matrix or metadata, no
  identifier column, an invalid normalization/scaling method, or a
  `venn_min_fraction` outside (0, 1]. An invalid
  scaling method is only detected at stage 5, after the earlier plots have
  already been written.
- Each stage returns a `dict`. Skipped supervised stages return
  `{"skipped": True, "reason": str}`.
- Public helpers that tests use directly: `normalize_matrix`,
  `pareto_scale`, `scale_features`, `extract_chromatogram`.

### 6.2 CLI

```
python visualization_module.py --data-matrix FILE [--metadata FILE] [--mzml FILE ...]
       [--output-dir DIR] [--group-column COL]
       [--normalization {none,sum,median}] [--scaling {none,pareto,auto}]
       [--venn-min-fraction FRACTION]
```

| Exit code | Meaning |
|---|---|
| 0 | Run completed. Supervised analyses may still have been skipped (see the report). |
| 1 | Uncaught exception (e.g. VIZ-IN-01/02/10), with a Python traceback. |
| 2 | argparse usage error (missing `--data-matrix`, invalid choice, `--venn-min-fraction` outside (0, 1]). |

The report is also printed to stdout.

`run_check.sh <folder> [group_col] [normalization] [scaling]` picks the first
`*.csv`/`*.tsv` file whose name contains `metadata` and the first containing
`matrix`, plus every `*.mzML` in the folder, and writes to
`<folder>/../viz_output`. It exits 1 if the folder or the matrix is missing.

## 7. Test traceability

Tests are in `tests/test_visualization.py`.

| Test | Requirements |
|---|---|
| `test_extract_chromatogram_uses_precomputed_tic_and_bpc` | VIZ-CH-01, CH-02, CH-04 |
| `test_extract_chromatogram_falls_back_to_decoding_binary_array_uncompressed` | VIZ-CH-03 (32-bit) |
| `test_extract_chromatogram_falls_back_to_decoding_binary_array_zlib_compressed` | VIZ-CH-03 (zlib) |
| `test_feature_summary_counts_and_excludes_descriptive_columns` | VIZ-IN-03, IN-05, FS-01, FS-04 |
| `test_pca_and_correlation_outputs` | VIZ-PC-01, PC-02, PC-04, OU-05 |
| `test_benjamini_hochberg_monotonic_and_bounded` | VIZ-VO-04 |
| `test_pls_da_two_groups` | VIZ-PL-10, GR-03 |
| `test_opls_da_two_groups_separates_predictive_component` | VIZ-OP-01, OP-02 |
| `test_volcano_plot_flags_the_shifted_features` | VIZ-VO-03, VO-04, VO-05, OU-04 |
| `test_supervised_analyses_skip_gracefully_with_only_one_group` | VIZ-SU-01, SU-03 |
| `test_multiple_attribute_columns_without_explicit_choice_skips_with_note` | VIZ-GR-02, GR-04 |
| `test_opls_da_and_volcano_require_exactly_two_groups` | VIZ-OP-01, VO-01, PL-10 |
| `test_normalize_matrix_sum_equalizes_sample_totals` | VIZ-NO-02 |
| `test_normalize_matrix_median_and_none` | VIZ-NO-03, NO-04, NO-06 |
| `test_pareto_scale_matches_hand_computed_formula` | VIZ-SC-02 |
| `test_scale_features_none_only_centers` | VIZ-SC-04, SC-05 |
| `test_sum_normalization_corrects_a_pure_dilution_difference_between_groups` | VIZ-NO-01, NO-02, VO-02 |
| `test_normalization_summary_reports_method_and_saves_plot` | VIZ-NO-03, SC-03 (config only) |
| `test_hist_with_gaussian_fit_returns_stats_and_handles_empty` | VIZ-DS-03 |
| `test_distribution_shape_summary_produces_plots_and_reports_skew_reduction` | VIZ-DS-01, DS-02 |
| `test_full_run_produces_report_and_plots` | VIZ-PL-01, OU-02, OV-07 |
| `test_extract_chromatogram_converts_minutes_to_seconds` | VIZ-CH-05 |
| `test_extract_chromatogram_decodes_array_when_only_bpc_missing` | VIZ-CH-03 (partial cvParams) |
| `test_pca_and_feature_summary_use_the_supervised_group_column` | VIZ-IN-12, FS-04, PC-03, GR-03 |
| `test_pca_reports_unmatched_samples` | VIZ-PC-03 |
| `test_feature_overlap_two_groups_counts_and_csv` | VIZ-OV-02, OV-03, OV-04, OV-05, OV-08, OU-06 |
| `test_feature_overlap_three_groups_region_counts` | VIZ-OV-04, OV-05, OV-08 |
| `test_feature_overlap_min_fraction_threshold` | VIZ-OV-03 |
| `test_feature_overlap_four_groups_uses_upset` | VIZ-OV-06, OV-08 |
| `test_feature_overlap_skips_with_one_group_or_no_metadata` | VIZ-OV-01 |
| `test_invalid_venn_min_fraction_raises` | VIZ-OV-03 |
| `test_volcano_counts_significant_features_per_group` | VIZ-VO-07 |
| `test_pca_3d_written_with_enough_samples` | VIZ-PC-05, PC-06 |
| `test_pca_3d_skipped_with_three_samples` | VIZ-PC-06 |

**Requirements with no test:** VIZ-IN-01/02 (error paths), IN-10/11, IN-13
(extension-less matching), IN-20; CH-03 (64-bit), CH-06; GR-01, GR-02
(not-found branch), GR-05, GR-06 (CLI prompt / no-TTY); SU-01 (group < 2
samples branch); OP-03; NO-05; SC-03 numerical check; SC-06
(correlation/volcano unscaled); PC-07 and the actual colours in VO-07 (both
checked by eye only); OV-06 truncation at 30 bars; the CLI and its exit
codes, including `--venn-min-fraction`; and `run_check.sh`.

## 8. Known limitations and discrepancies

### 8.1 Limitations (current behaviour, by design or not yet addressed)

1. **Binary array encodings (VIZ-CH-03).** Only 32- and 64-bit float
   arrays are supported. Integer arrays and MS-Numpress are not.
2. **Chromatogram filename collisions.** Two mzML files with the same stem
   in different folders overwrite each other's PNG.
3. **Display-column fallback (VIZ-IN-12).** When no supervised column
   resolves (e.g. several `ATTRIBUTE_` columns, none chosen), PCA colouring
   still uses the looser fallback rule, while the supervised analyses are
   skipped.
4. **A zero-total sample becomes all-NaN after normalization (VIZ-NO-05)**,
   and then all-zero in multivariate stages, with no warning.
5. **Fixed significance thresholds on raw p (VIZ-VO-05).** They cannot be
   configured, and q-values are not used for the overlay.
6. **No model validation** for PLS-DA or OPLS-DA (no CV, R²/Q², permutation
   test or VIP). The scores plots alone can overfit when there are few
   samples.
7. **Silent coercion (VIZ-IN-04).** Non-numeric matrix cells become NaN with
   no warning. The module relies on the validator having run first.
8. **Partial outputs on a late failure.** An invalid `scaling` passed
    through the API fails at stage 5, after the earlier plots have been
    written, and no report is produced.
9. **Group-column resolution runs several times.** `_resolve_group_column`
   runs for the feature summary, PCA, each supervised analysis and the
   feature overlap. This is harmless but redundant.
10. **The Venn diagram is not area-proportional (VIZ-OV-05).** Only the
    numbers carry quantitative meaning.
11. **The UpSet chart is capped at 30 bars (VIZ-OV-06).** Smaller
    intersections appear only in the CSV.
12. **Static 3D PCA (VIZ-PC-06).** The view angle is fixed (matplotlib
    default), so some separation may be hidden until the interactive
    version (VIZ-DF-03).
13. **Capture ignores intensity (VIZ-OV-02).** Any value > 0 counts as
    detected, so gap-filled noise-level values count. Raise
    `venn_min_fraction` to compensate.

### 8.2 Documentation discrepancies

None known as of 1.1.0. The discrepancies found in 1.0.0 (scaling of the
correlation matrix, `--group-column` help, the volcano "original values"
wording, the per-module `requirements.txt`, and the missing spec
reference/version in the docstring) were resolved in 1.1.0.

## 9. Change control

Behaviour changes follow the same order as the validator module:

1. **Spec.** Edit this file. Add, change or retire rule IDs. Never reuse a
   retired ID; mark it *Retired in x.y.z* instead.
2. **Tests.** Add or update tests in `tests/test_visualization.py` and the
   §7 mapping.
3. **Code.** Change `visualization_module.py`.
4. **Docs.** Update the module docstring, the module `README.md`, and the
   root `CLAUDE.md` if a "deliberate rule" there is affected.
5. **Verify.** `pytest -v` passes, and `./run_check.sh
   ../Upload_and_validator_module/sample_data/valid_files ATTRIBUTE_Tissue`
   completes on the real dataset.
6. **Version.** Bump the semver and add a changelog entry (§10):
   - **Major**: a changed output file name or format, a changed default
     (normalization, scaling, thresholds), changed API fields or return keys,
     or a changed exit-code meaning.
   - **Minor**: a new stage, output, option or rule that keeps existing
     behaviour.
   - **Patch**: bug fixes that bring code in line with the spec, or doc-only
     changes.

The data-handling rules in VIZ-NO-07, VIZ-SC-06, the VIZ-PC/VIZ-VO NA
handling and VIZ-GR-06 (never block without a TTY) are deliberate. Changing
any of them is a major change and needs explicit sign-off.

## 10. Changelog

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-10-20 | First spec, written from the existing implementation and test suite. No code changes. |
| 1.1.0 | 2026-10-20 | VIZ-IN-12/FS-04: PCA colouring and group counts use the supervised grouping column. VIZ-PC-03: unmatched samples drawn and counted. VIZ-CH-03: decode the array when either TIC or BPC is missing. VIZ-CH-05: minute RTs converted to seconds. Added `__version__`. Docs fixed (correlation is unscaled, `--group-column` accepts any column, volcano uses normalized values). Removed the duplicate `modules/visualization/requirements.txt` (identical to root). 4 new tests. |
| 1.2.0 | 2026-10-20 | New feature-overlap stage (VIZ-OV-01…08): Venn for 2–3 groups, UpSet for 4+, `feature_overlap.csv`, and a `venn_min_fraction` input / `--venn-min-fraction` flag. 3D PCA (VIZ-PC-05/06). Volcano significant points coloured by the group they are higher in, with `n_higher_in` added (VIZ-VO-07). Shared group-colour map (VIZ-PC-07). VIZ-DF-03 added. 10 new tests. |
