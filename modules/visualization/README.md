# Visualization (unannotated) module

Phase 2 of TTL_pipeline. Takes the same core inputs as the upload/validation
module — mzML files, sample metadata, and a data matrix (no feature table,
MGF, or annotation list needed here, hence "unannotated") — and produces a
feature summary, chromatograms, sample-wise normalization + feature-wise
scaling, PCA, a sample-sample correlation heatmap, and supervised group
comparisons (PLS-DA, OPLS-DA, a volcano plot), all saved as PNG plots plus
a text analysis report. It also draws a 3D PCA and a feature-overlap
view (Venn diagram / UpSet chart) of which features each group captures.
This covers the full design diagram except the interactive MS1/MS2
spectra viewer and interactive plots — see "Deferred" below.

## What it produces

Given `--data-matrix` (required) and optionally `--metadata` / `--mzml` /
`--group-column` / `--normalization` / `--scaling` / `--venn-min-fraction`, it writes to
`--output-dir` (default `viz_output/`):

- `analysis_report.txt` — a text report covering everything below.
- `missing_values_per_sample.png` — bar chart of missing (NA) feature
  counts per sample.
- `intensity_distribution.png` — histogram of log10(intensity) across all
  samples.
- `chromatogram_<filename>.png` — one per mzML file: TIC and BPC vs
  retention time, extracted from MS1 scans.
- `normalization_effect.png` — per-sample total intensity before vs after
  normalization, so a systematic per-sample scale difference (and whether
  it got corrected) is visible at a glance.
- `sample_distribution_gaussian.png` — pooled, log1p-transformed sample
  intensities before vs after normalization, each panel a histogram with
  a fitted Gaussian curve overlaid, plus mean/std/skewness in the report.
- `feature_distribution_gaussian.png` — pooled feature values (the exact
  input scale_features() receives) before vs after scaling, same
  histogram-plus-Gaussian-fit treatment — shows visually and numerically
  (skewness) how scaling pulls a handful of dominant high-abundance
  features back into a tighter, more symmetric distribution.
- `pca.png` — PCA scatter of samples (PC1 vs PC2, on normalized + scaled
  data), colored by the same grouping column the supervised analyses use
  (see below). Samples with no metadata match are drawn as grey
  "unmatched" points and counted in the report.
- `pca_3d.png` — static 3D PCA (PC1–PC3), same colouring. Needs at least
  4 samples and 3 features; otherwise skipped with a note.
- `correlation_heatmap.png` + `sample_correlation_matrix.csv` —
  sample-sample Pearson correlation (on normalized, log1p-transformed
  intensities).
- `pls_da.png` — PLS-DA scores plot (LV1 vs LV2, on normalized + scaled
  data), colored by group. Handles two or more groups.
- `opls_da.png` — OPLS-DA scores plot (predictive vs orthogonal
  component, on normalized + scaled data). Two groups only (see
  "Two-group-only analyses" below).
- `volcano_plot.png` + `volcano_plot_results.csv` — log2 fold change vs
  -log10(p-value) per feature (on normalized but **unscaled** data — see
  "Normalization and scaling" below for why), with Benjamini-Hochberg
  q-values in the CSV. Two groups only. Significant features are coloured
  by the group they are higher in (same colour as that group in the PCA /
  PLS-DA / OPLS-DA plots); non-significant ones are grey.
- `venn_diagram.png` (2–3 groups) or `upset_plot.png` (4+ groups) +
  `feature_overlap.csv` — which features each group captures, and which
  are shared or unique. See "Feature overlap" below.

## Choosing the grouping column (PLS-DA / OPLS-DA / volcano plot)

Only metadata columns following GNPS's `ATTRIBUTE_*` naming convention are
considered automatically — e.g. `ATTRIBUTE_Tissue`, `ATTRIBUTE_Group`.
`--group-column` accepts any metadata column. The PCA colouring and the
feature-summary group counts use the same resolved column; only when no
column resolves do they fall back to a looser `Class`/`Group`/`Condition`/
`*attribute*` match.

- **Exactly one `ATTRIBUTE_*` column** in the metadata: used automatically.
- **More than one**: running the script interactively prompts you to pick
  one from a numbered list. Running non-interactively (no TTY — e.g. from
  another script or CI) instead skips these three analyses with a note in
  the report listing the candidates, so nothing hangs waiting for input.
  Pass `--group-column ATTRIBUTE_Whatever` to skip the prompt/ambiguity
  entirely.
- **None**: PLS-DA/OPLS-DA/volcano plot are skipped with a note; the rest
  of the report (feature summary, chromatograms, PCA, correlation) still
  runs normally.

This interactive prompt is a script-era stand-in for the eventual GUI: a
future GUI will present the same list of `ATTRIBUTE_*` columns as a picker
and set `VisualizationInputs.group_column` directly, bypassing the prompt
— the underlying `VisualizationModule`/`VisualizationInputs` API is
already GUI-ready in that sense.

## Feature overlap (Venn / UpSet)

A feature counts as **captured** by a group when it is detected (value
> 0, not NA, in the raw matrix) in at least `--venn-min-fraction` of that
group's samples (default `0.5`, i.e. at least half). Each feature is then
assigned to exactly one exclusive intersection, such as `MCL only`,
`MCL & MCR` or `none`.

- **2–3 groups**: a Venn diagram showing the count in every region, with
  each group's total in its label. The circles are a fixed layout and
  **not** proportional to the counts.
- **4+ groups**: an UpSet-style chart (Venn diagrams stop being readable).
  The top bars are intersection sizes, largest first, showing the 30
  largest. The dot matrix below shows which groups form each intersection,
  and the left bars show each group's total.
- `feature_overlap.csv` lists every feature's per-group `captured_<group>`
  flags and its `overlap` label, so any region can be pulled out for
  follow-up.

The grouping column is resolved exactly as for PLS-DA (see below). It needs
at least 2 groups but, unlike PLS-DA, only 1 sample per group.

## Two-group-only analyses

OPLS-DA and the volcano plot are inherently pairwise comparisons and
require the resolved grouping column to have **exactly 2** groups (after
excluding any samples that couldn't be matched to a metadata row). If it
has more, both are skipped with a note explaining the group count — use
PLS-DA instead, which handles any number of groups. Each group also needs
at least 2 matched samples; smaller groups are skipped with a clear reason
rather than silently producing meaningless (NaN) statistics.

## Normalization and scaling

Two distinct preprocessing steps are applied before cross-sample
comparisons (normalization for all of PCA/correlation/PLS-DA/OPLS-DA/volcano;
scaling for PCA/PLS-DA/OPLS-DA only) — the raw
feature-summary stats/plots and chromatograms are untouched, since those
are meant to show the real, unmodified data for QC purposes.

**Normalization** (`--normalization`, default `sum`) is sample-wise: it
corrects for differences in each sample's overall amount/dilution/
instrument response, so samples become comparable before anything else.
This is what fixes the exact issue flagged from the real dataset — two
groups running at systematically different overall intensity (e.g. from
loading amount or instrument drift, not real biology) can otherwise
inflate the apparent number of "significant" features across the board.
Options: `sum` (total ion current normalization — each sample divided by
its own total, rescaled to keep the overall magnitude comparable),
`median` (same idea, more robust to a few very large peaks dominating the
sum), or `none`.

**Scaling** (`--scaling`, default `pareto`) is feature-wise, applied on
top of normalization, but **only** for PCA/PLS-DA/OPLS-DA — deliberately
**not** for the volcano plot or the correlation matrix (which is Pearson on
log1p of normalized, unscaled intensities). Pareto scaling (van den Berg et
al., 2006 — mean-center each feature, divide by the square root of its
standard deviation) is the metabolomics-standard choice for multivariate
methods: it prevents a handful of very high-abundance features from
dominating PCA/PLS-DA purely by magnitude, while being less aggressive
than full autoscaling. Options: `pareto`, `auto` (standard autoscaling /
unit-variance, equivalent to `StandardScaler`), or `none` (mean-center
only).

The volcano plot uses normalized-but-**unscaled** values for its fold
change and t-test on purpose: scaling mean-centers each feature
independently (so values can go negative) and rescales it relative to its
own variance — after that, "fold change" no longer means anything as a
ratio of real, comparable quantities. Normalization alone (which just
corrects each *sample's* overall scale, not each feature's) keeps the
per-feature values meaningful as real intensities, so the fold change
stays interpretable.

### Seeing it as a distribution, not just numbers

`sample_distribution_gaussian.png` and `feature_distribution_gaussian.png`
plot these effects as actual histograms with a fitted Gaussian (bell)
curve overlaid, rather than just before/after summary numbers:

- **Samples**: pooled, log1p-transformed intensities, before vs after
  normalization. Untransformed intensities are strongly right-skewed and
  won't look bell-shaped at all regardless of normalization — the log
  transform is what makes the comparison meaningful here. Normalization's
  effect (correcting each *sample's* overall scale) is really the
  per-sample view in `normalization_effect.png`; this pooled view mainly
  shows whether the overall intensity distribution shape is reasonable.
- **Features**: pooled values in their normalized, **not** log-transformed
  units — the exact input `scale_features()` receives — before vs after
  scaling. Before scaling, a handful of very high-abundance features
  dominate the pooled range and the distribution is typically extremely
  right-skewed (nowhere near Gaussian, often looking like a single spike
  near zero with a long tail). After Pareto/auto scaling, it collapses
  into a much tighter, more symmetric distribution centered on 0.

Both panels report the fitted mean, standard deviation, and skewness in
`analysis_report.txt` — skewness near 0 is what "Gaussian-shaped"
corresponds to numerically. These are visual/diagnostic aids, not a
formal normality test (e.g. Shapiro-Wilk).

## Statistical methodology

- **PLS-DA**: `sklearn.cross_decomposition.PLSRegression` against a
  one-hot-encoded group matrix (standard multi-class PLS-DA formulation).
- **OPLS-DA**: single-orthogonal-component OPLS (Trygg & Wold, 2002) —
  removes the variation in the (scaled) feature matrix that's orthogonal
  to the group label, then fits one predictive PLS component on what's
  left. This is the standard formulation behind the classic OPLS-DA
  scores plot (predictive vs orthogonal axis).
- **Volcano plot**: log2 fold change of group means (with a small
  pseudocount so it's defined at zero), Welch's t-test (unequal variance)
  per feature, Benjamini-Hochberg FDR correction for the q-value column.
  The default significance overlay is `p < 0.05` and `|log2FC| > 1`.

## Running it

```bash
cd modules/visualization
python visualization_module.py \
  --mzml path/to/*.mzML \
  --metadata path/to/metadata.csv \
  --data-matrix path/to/data_matrix.csv \
  --group-column ATTRIBUTE_Tissue \
  --normalization sum \
  --scaling pareto \
  --output-dir viz_output
```

`--group-column`, `--normalization`, and `--scaling` are all optional —
omit `--group-column` and the script auto-resolves a single `ATTRIBUTE_*`
column, or prompts/skips as described above; omit `--normalization`/
`--scaling` and it defaults to `sum`/`pareto`.

Or use the convenience wrapper, which detects the metadata/matrix files by
keyword in the filename (same approach as the upload_validation module's
`run_check.sh`) and optionally forwards a group column, normalization, and
scaling method:

```bash
./run_check.sh path/to/folder [ATTRIBUTE_column_name] [normalization] [scaling]
```

## Data-matrix shape

Accepts either a bare `feature_id, sample1, sample2, ...` matrix, or
MZmine's full feature-table shape `row ID, row m/z, row retention time,
<sample columns>` — the `row m/z` / `row retention time` columns, when
present, are excluded from the per-sample statistics automatically.

## mzML chromatogram extraction

Prefers each spectrum's precomputed `total ion current` (MS:1000285) and
`base peak intensity` (MS:1000505) cvParams, which most converters
(msconvert, MZmine, etc.) already write out — this is fast even on large
files. Falls back to decoding the intensity binary array itself (handling
both 32-bit and 64-bit floats, zlib-compressed or not) only when those
cvParams are missing — including when only one of the two is present.
Only MS1 scans are included, matching standard TIC/BPC chromatogram
convention. Retention times are plotted in seconds; scan start times stored
in minutes (`UO:0000031`) are converted.

## Assumptions worth knowing about

- Missing/NA values in the data matrix are imputed as 0 for PCA, the
  correlation matrix, and PLS-DA/OPLS-DA specifically (not for the
  missing-value count/report, which reports them as-is). This is a
  common but not the only valid choice for exploratory-stage untargeted
  metabolomics data — revisit if a different imputation strategy is
  needed downstream.
- The volcano plot's fold change and t-test use the normalized (but
  unscaled, non-imputed) values per feature, via `nanmean`/`nan_policy="omit"` —
  missing values are excluded per-comparison rather than treated as zero,
  since zero-imputation would badly distort a ratio-based statistic like
  fold change.

## Deferred to a later increment

- **Interactive MS1/MS2 spectra viewer + file export** — this implies an
  actual UI, which is out of scope for a plain script; revisit once the
  pipeline's GUI phase is underway (confirmed as the plan for this and
  the grouping-column picker above).
- **Interactive Venn/UpSet and 3D PCA** — static PNGs for now; the GUI
  phase will make them interactive (e.g. rotating the PCA, clicking a Venn
  region to list its features from `feature_overlap.csv`).

## Specification

`SPEC.md` is the source of truth for this module's behaviour (rule IDs
`VIZ-xx-nn`, API/CLI contract, test mapping, known limitations, change
control). Change behaviour spec-first, as described in its §9.

## Running tests

```bash
cd modules/visualization
pytest -v
```
