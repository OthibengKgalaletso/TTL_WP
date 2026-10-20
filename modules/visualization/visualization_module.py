"""
Visualization (unannotated) module — Phase 2 of TTL_pipeline.

Consumes the same core inputs as the upload/validation module (mzML files,
sample metadata, data matrix — no feature table / MGF / annotation list,
hence "unannotated") and produces:

  - A feature summary: feature count, missing values, intensity
    distribution, sample summary (each also saved as a plot where it
    makes sense).
  - Chromatograms (TIC/BPC) extracted from the mzML files.
  - Sample-wise normalization (sum/median/none) and a before/after plot,
    applied ahead of every cross-sample comparison below.
  - Distribution-shape plots (pooled histogram + fitted Gaussian curve)
    for samples (raw vs normalized) and features (normalized-but-
    unscaled vs scaled), so the effect of each step is visible as an
    actual distribution shape rather than just summary numbers.
  - PCA (on normalized + feature-wise-scaled data — pareto scaling by
    default, the metabolomics standard), plus a sample-sample correlation
    matrix and heatmap (Pearson on log1p of normalized, unscaled data).
    A 3D PCA (PC1-PC3) is also drawn when there are enough samples.
  - PLS-DA and OPLS-DA scores plots (also on normalized + scaled data),
    and a volcano plot (on normalized-but-UNscaled data, so its fold
    change stays interpretable; significant features coloured by the
    group they are higher in) — all grouped by a metadata column
    following GNPS's ATTRIBUTE_ naming convention.
  - A feature-overlap view: which features each group captures (detected
    in >= venn_min_fraction of its samples), as a Venn diagram for 2-3
    groups or an UpSet-style chart for 4+, plus a membership CSV.
  - A text analysis report tying all of the above together.

All plots are saved as PNG files into an output directory (static images,
not a live/interactive viewer — see the module README for why, and for
what's deferred).

Deferred to a later increment:
  - An interactive MS1/MS2 spectra viewer — a UI/GUI decision, out of
    scope for a plain script; revisit once the pipeline's GUI phase is
    underway.
  - Interactive versions of the Venn/UpSet and 3D PCA plots (static PNGs
    for now), for the same reason.

SPEC.md in this folder is the source of truth for behaviour; keep
__version__ in step with its spec version.
"""

from __future__ import annotations

import argparse
import base64
import os
import re
import sys
import zlib
from dataclasses import dataclass, field
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")  # headless/script-safe backend — no display needed
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from scipy import stats as scipy_stats
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

__version__ = "1.2.0"


# --------------------------------------------------------------------- #
# Small shared helpers (kept local rather than importing the validator,
# so this module stays independently runnable — same philosophy as the
# upload_validation module).
# --------------------------------------------------------------------- #

def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _local_tag(tag: str) -> str:
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _find_col(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    norm_candidates = {_norm(c) for c in candidates}
    for col in df.columns:
        if _norm(col) in norm_candidates:
            return col
    return None


def _find_col_containing(df: pd.DataFrame, substrings: list[str]) -> Optional[str]:
    subs = [_norm(s) for s in substrings]
    for col in df.columns:
        nc = _norm(col)
        if any(s in nc for s in subs):
            return col
    return None


def _read_table(path: str) -> Optional[pd.DataFrame]:
    """Delimiter-sniffing table reader (comma/semicolon/tab), BOM-safe."""
    try:
        return pd.read_csv(path, sep=None, engine="python", dtype=str, encoding="utf-8-sig")
    except Exception:
        try:
            return pd.read_csv(path, dtype=str, encoding="utf-8-sig")
        except Exception:
            return None


def _filename_variants(name: str) -> set[str]:
    name = str(name).strip()
    variants = {name.lower()}
    variants.add(os.path.splitext(name)[0].lower())
    return variants


def _benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values (NaNs pass through as NaN)."""
    p = np.asarray(pvalues, dtype=float)
    n = np.sum(~np.isnan(p))
    q = np.full_like(p, np.nan)
    if n == 0:
        return q
    valid_idx = np.where(~np.isnan(p))[0]
    order = valid_idx[np.argsort(p[valid_idx])]
    ranked = p[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p-value downward
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    q[order] = adjusted
    return q


def _group_color_map(unique_groups) -> dict:
    """tab10 colour per group, assigned in sorted order, so a given group
    gets the same colour in every plot (PCA, PLS-DA, OPLS-DA, volcano,
    Venn)."""
    cmap = plt.get_cmap("tab10")
    ordered = sorted(unique_groups, key=str)
    return {g: cmap(i % 10) for i, g in enumerate(ordered)}


# --------------------------------------------------------------------- #
# Normalization (sample-wise) and scaling (feature-wise)
# --------------------------------------------------------------------- #

def normalize_matrix(values: pd.DataFrame, method: str) -> pd.DataFrame:
    """Sample-wise normalization: corrects each sample column for overall
    differences in total amount/dilution/instrument response, so samples
    become comparable before any cross-sample comparison. Operates on a
    features x samples DataFrame (rows=features, columns=samples).

    - "sum": total sum (a.k.a. total ion current / TIC) normalization —
      each sample is divided by its own total, then rescaled by the mean
      total across samples so the overall magnitude stays comparable to
      the input (a convenience for reading plots — the ratio between
      samples is what matters, not the absolute scale).
    - "median": same idea, using each sample's median instead of its sum
      (more robust to a handful of very large peaks dominating the sum).
    - "none": returned unchanged.
    """
    if method == "none":
        return values.copy()
    if method == "sum":
        totals = values.sum(axis=0, skipna=True)
    elif method == "median":
        totals = values.median(axis=0, skipna=True)
    else:
        raise ValueError(f"Unknown normalization method: {method!r} (expected none/sum/median)")
    totals_safe = totals.replace(0, np.nan)
    target = totals_safe.mean(skipna=True)
    if not np.isfinite(target) or target == 0:
        return values.copy()
    return values.div(totals_safe, axis=1) * target


def pareto_scale(sample_matrix: np.ndarray) -> np.ndarray:
    """Pareto scaling (van den Berg et al., 2006): mean-center each feature
    (column) and divide by the square root of its standard deviation,
    rather than the standard deviation itself (as autoscaling/StandardScaler
    would). This is the standard scaling method in metabolomics for
    multivariate analyses (PCA/PLS-DA/OPLS-DA) — it moderates the influence
    of high-variance features less aggressively than autoscaling, keeping
    the data closer to its original measurement scale.

    Expects sample_matrix as samples x features (rows=samples).
    """
    mean = sample_matrix.mean(axis=0)
    std = sample_matrix.std(axis=0, ddof=0)
    std_safe = np.where(std == 0, 1.0, std)
    return (sample_matrix - mean) / np.sqrt(std_safe)


def scale_features(sample_matrix: np.ndarray, method: str) -> np.ndarray:
    """Feature-wise scaling ahead of a multivariate method. "pareto" uses
    pareto_scale(); "auto" is standard autoscaling (mean-center + unit
    variance, i.e. StandardScaler); "none" only mean-centers (no variance
    scaling), which is rarely what you want for PCA/PLS-DA/OPLS-DA but is
    available for comparison."""
    if method == "pareto":
        return pareto_scale(sample_matrix)
    if method == "auto":
        return StandardScaler().fit_transform(sample_matrix)
    if method == "none":
        return sample_matrix - sample_matrix.mean(axis=0)
    raise ValueError(f"Unknown scaling method: {method!r} (expected none/pareto/auto)")


def _hist_with_gaussian_fit(ax, values: np.ndarray, color: str, bins: int = 60) -> Optional[dict]:
    """Plots a density histogram of `values` on `ax` with a fitted Gaussian
    (normal) curve overlaid, so the distribution's shape can be read off
    visually. Returns the fitted mean/std/skewness (or None if there's
    nothing finite to plot) — skewness near 0 is what "looks Gaussian"
    actually corresponds to, since a histogram alone can be misleading at
    a glance.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        ax.set_visible(False)
        return None

    ax.hist(values, bins=bins, density=True, color=color, alpha=0.6)
    mu = float(values.mean())
    sigma = float(values.std())
    if sigma > 0:
        x = np.linspace(values.min(), values.max(), 200)
        pdf = scipy_stats.norm.pdf(x, mu, sigma)
        ax.plot(x, pdf, color="black", linewidth=1.5, label=f"Gaussian fit (μ={mu:.2f}, σ={sigma:.2f})")
        ax.legend(fontsize=7)
    skew = float(scipy_stats.skew(values)) if values.size > 2 else None
    return {"mean": mu, "std": sigma, "skew": skew, "n": int(values.size)}


# --------------------------------------------------------------------- #
# mzML chromatogram extraction
# --------------------------------------------------------------------- #

def _decode_binary_array(binary_text: str, is_64bit: bool, is_zlib: bool) -> np.ndarray:
    raw = base64.b64decode(binary_text)
    if is_zlib:
        raw = zlib.decompress(raw)
    dtype = np.float64 if is_64bit else np.float32
    return np.frombuffer(raw, dtype=dtype)


def extract_chromatogram(path: str) -> dict:
    """Streams an mzML file and returns MS1 retention times, TIC, and BPC.

    Prefers precomputed per-spectrum cvParams (total ion current / base
    peak intensity), which almost every mzML converter (e.g. msconvert,
    MZmine) writes out — this avoids decoding every binary intensity
    array. Falls back to decoding the intensity array itself only when
    either cvParam is missing. Retention times are returned in seconds
    (minute-valued scan start times are converted).
    """
    rt: list[float] = []
    tic: list[float] = []
    bpc: list[float] = []

    current: Optional[dict] = None
    in_binary_array = False
    array_is_intensity = False
    array_is_64bit = False
    array_is_zlib = False

    for event, elem in ET.iterparse(path, events=("start", "end")):
        tag = _local_tag(elem.tag)

        if event == "start":
            if tag == "spectrum":
                current = {
                    "ms_level": None,
                    "rt": None,
                    "tic": None,
                    "bpc": None,
                    "intensity_array": None,
                }
            elif tag == "binaryDataArray" and current is not None:
                in_binary_array = True
                array_is_intensity = False
                array_is_64bit = False
                array_is_zlib = False
            elif tag == "cvParam" and current is not None:
                accession = elem.get("accession")
                value = elem.get("value")
                if not in_binary_array:
                    if accession == "MS:1000511":  # ms level
                        try:
                            current["ms_level"] = int(value)
                        except (TypeError, ValueError):
                            pass
                    elif accession == "MS:1000016":  # scan start time
                        try:
                            rt_val = float(value)
                        except (TypeError, ValueError):
                            rt_val = None
                        if rt_val is not None:
                            # Report RT in seconds regardless of how the
                            # file stores it (MZmine/msconvert often write
                            # minutes, UO:0000031).
                            unit_acc = elem.get("unitAccession")
                            unit_name = (elem.get("unitName") or "").lower()
                            if unit_acc == "UO:0000031" or unit_name == "minute":
                                rt_val *= 60.0
                            current["rt"] = rt_val
                    elif accession == "MS:1000285":  # total ion current
                        try:
                            current["tic"] = float(value)
                        except (TypeError, ValueError):
                            pass
                    elif accession == "MS:1000505":  # base peak intensity
                        try:
                            current["bpc"] = float(value)
                        except (TypeError, ValueError):
                            pass
                else:
                    if accession == "MS:1000515":  # intensity array
                        array_is_intensity = True
                    elif accession == "MS:1000523":  # 64-bit float
                        array_is_64bit = True
                    elif accession == "MS:1000574":  # zlib compression
                        array_is_zlib = True

        elif event == "end":
            if tag == "binary" and in_binary_array and current is not None and array_is_intensity:
                if elem.text and (current["tic"] is None or current["bpc"] is None):
                    try:
                        arr = _decode_binary_array(elem.text.strip(), array_is_64bit, array_is_zlib)
                        current["intensity_array"] = arr
                    except Exception:
                        pass
            elif tag == "binaryDataArray":
                in_binary_array = False
            elif tag == "spectrum" and current is not None:
                if current["ms_level"] == 1 and current["rt"] is not None:
                    tic_val = current["tic"]
                    bpc_val = current["bpc"]
                    if (tic_val is None or bpc_val is None) and current["intensity_array"] is not None:
                        arr = current["intensity_array"]
                        if arr.size:
                            if tic_val is None:
                                tic_val = float(arr.sum())
                            if bpc_val is None:
                                bpc_val = float(arr.max())
                    if tic_val is not None and bpc_val is not None:
                        rt.append(current["rt"])
                        tic.append(tic_val)
                        bpc.append(bpc_val)
                current = None
                elem.clear()

    order = np.argsort(rt)
    return {
        "rt": np.array(rt)[order],
        "tic": np.array(tic)[order],
        "bpc": np.array(bpc)[order],
    }


# --------------------------------------------------------------------- #
# Inputs / module
# --------------------------------------------------------------------- #

@dataclass
class VisualizationInputs:
    mzml_files: list[str] = field(default_factory=list)
    metadata_file: Optional[str] = None
    data_matrix_file: Optional[str] = None
    output_dir: str = "viz_output"
    # Which metadata column defines the groups for PLS-DA/OPLS-DA/volcano
    # plot. Only ATTRIBUTE_-prefixed columns (GNPS convention) are
    # considered automatically. Leave unset to auto-resolve: if metadata
    # has exactly one ATTRIBUTE_ column it's used automatically; if it has
    # more than one, the CLI prompts interactively (a future GUI would
    # instead present the same list of candidates as a picker and set
    # this field directly, skipping the prompt).
    group_column: Optional[str] = None
    # Sample-wise normalization applied before any cross-sample comparison
    # (PCA/correlation/PLS-DA/OPLS-DA/volcano plot) — "sum" (total ion
    # current normalization), "median", or "none". Does NOT affect the raw
    # feature-summary stats/chromatograms, which intentionally show the
    # unmodified data for QC purposes.
    normalization: str = "sum"
    # Feature-wise scaling applied on top of the normalized data, but only
    # for the multivariate methods (PCA/correlation/PLS-DA/OPLS-DA) — NOT
    # for the volcano plot, whose fold-change/t-test need normalized-but-
    # unscaled values to stay interpretable as real fold changes. "pareto"
    # (the metabolomics standard), "auto" (autoscaling/unit variance), or
    # "none" (mean-center only).
    scaling: str = "pareto"
    # Feature-overlap (Venn/UpSet) capture rule: a feature counts as
    # captured by a group when it is detected (> 0, not NA) in at least
    # this fraction of the group's samples. Must be in (0, 1].
    venn_min_fraction: float = 0.5


class VisualizationModule:
    def __init__(self, inputs: VisualizationInputs) -> None:
        self.inputs = inputs
        self._matrix_df: Optional[pd.DataFrame] = None
        self._normalized_df: Optional[pd.DataFrame] = None
        self._feature_id_col: Optional[str] = None
        self._sample_cols: list[str] = []
        self._metadata_df: Optional[pd.DataFrame] = None
        self._metadata_filename_col: Optional[str] = None
        self._metadata_group_col: Optional[str] = None
        self._plots: list[str] = []
        self._sections: list[str] = []
        self._group_column_used: Optional[str] = None

    # -- loading ---------------------------------------------------- #

    def _load_matrix(self) -> None:
        df = _read_table(self.inputs.data_matrix_file)
        if df is None:
            raise ValueError(f"Could not parse data matrix: {self.inputs.data_matrix_file}")

        feature_id_col = _find_col(df, ["feature_id", "id", "peak_id", "row ID"])
        if feature_id_col is None:
            raise ValueError("Data matrix: no recognizable identifier column (feature_id/id/peak_id/row ID).")

        mz_col = _find_col(df, ["row m/z", "m/z", "mz"])
        rt_col = _find_col(df, ["row retention time", "retention time", "rt"])
        descriptive = {c for c in (mz_col, rt_col) if c is not None}

        sample_cols = [c for c in df.columns if c != feature_id_col and c not in descriptive]
        df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce")

        self._matrix_df = df
        self._feature_id_col = feature_id_col
        self._sample_cols = sample_cols

    def _load_metadata(self) -> None:
        if not self.inputs.metadata_file:
            return
        df = _read_table(self.inputs.metadata_file)
        if df is None:
            raise ValueError(f"Could not parse metadata: {self.inputs.metadata_file}")
        filename_col = _find_col(df, ["filename", "file_name"])
        group_col = _find_col(df, ["Class", "Group", "Condition"]) or _find_col_containing(df, ["attribute"])
        self._metadata_df = df
        self._metadata_filename_col = filename_col
        self._metadata_group_col = group_col

    def load(self) -> None:
        if not (0 < self.inputs.venn_min_fraction <= 1):
            raise ValueError(
                f"venn_min_fraction must be in (0, 1], got {self.inputs.venn_min_fraction!r}"
            )
        self._load_matrix()
        self._load_metadata()
        os.makedirs(self.inputs.output_dir, exist_ok=True)
        self._normalized_df = normalize_matrix(self._matrix_df[self._sample_cols], self.inputs.normalization)

    def normalization_summary(self) -> dict:
        """Reports the normalization method applied and plots per-sample
        total intensity before vs after, so a systematic per-sample scale
        difference (e.g. one group loaded/ionizing more overall) — and
        whether normalization corrected it — is visible at a glance."""
        raw_totals = self._matrix_df[self._sample_cols].sum(axis=0, skipna=True)
        norm_totals = self._normalized_df.sum(axis=0, skipna=True)

        fig, axes = plt.subplots(1, 2, figsize=(max(8, len(self._sample_cols) * 0.6), 4), sharey=False)
        for ax, totals, title in ((axes[0], raw_totals, "Before"), (axes[1], norm_totals, "After")):
            ax.bar(range(len(self._sample_cols)), totals.to_numpy(), color="#4C72B0")
            ax.set_xticks(range(len(self._sample_cols)))
            ax.set_xticklabels(self._sample_cols, rotation=90, fontsize=6)
            ax.set_title(title)
        axes[0].set_ylabel("Total intensity per sample")
        fig.suptitle(f"Normalization effect (method: {self.inputs.normalization})")
        plot_path = self._save_plot(fig, "normalization_effect.png")

        self._sections.append(
            f"--- Normalization / scaling ---\n"
            f"Normalization: {self.inputs.normalization} "
            f"(applied to PCA/correlation/PLS-DA/OPLS-DA/volcano; raw values are kept "
            f"for the feature-summary stats/plots above and chromatograms below)\n"
            f"Scaling: {self.inputs.scaling} "
            f"(applied on top of normalization for PCA/PLS-DA/OPLS-DA only — "
            f"the correlation matrix uses log1p of normalized-but-unscaled values, and "
            f"the volcano plot's fold-change/t-test use normalized-but-unscaled values so "
            f"they stay interpretable as real fold changes)\n"
            f"Plot: {plot_path}\n"
        )
        return {"normalization": self.inputs.normalization, "scaling": self.inputs.scaling, "plot": plot_path}

    def distribution_shape_summary(self) -> dict:
        """Plots pooled intensity distributions as histograms with a fitted
        Gaussian curve overlaid — one figure for samples (raw vs
        normalized), one for features (normalized-but-unscaled vs scaled).

        Sample panel: raw and normalized values are log1p-transformed
        first. Untransformed intensities are strongly right-skewed and
        won't resemble a bell curve at all regardless of normalization —
        the log transform is what makes a Gaussian-shaped comparison
        meaningful here; normalization's own job is correcting each
        sample's overall scale, which shows up as the "after" curve
        sitting more consistently across samples (see
        normalization_effect.png), not necessarily as a visibly different
        bell shape in this pooled view.

        Feature panel: reuses the exact normalized matrix that
        PCA/PLS-DA/OPLS-DA feed into scale_features() as "before scaling",
        and the exact scaled matrix (pareto/auto/none) they actually
        consume afterward — so this shows precisely what scaling does to
        the data those methods see. Before scaling, a handful of very
        high-abundance features dominate the pooled range and the
        distribution is typically extremely right-skewed (nowhere near
        Gaussian); scaling mean-centers and rescales each feature
        individually, so pooled across all features it collapses into a
        much tighter, more symmetric distribution around 0.

        Neither panel is a formal normality test (e.g. Shapiro-Wilk) —
        it's a visual/diagnostic aid. Skewness near 0 is what "Gaussian-
        shaped" corresponds to numerically; the fitted curve's mean/std/
        skew are recorded in the report for exactly that reason.
        """
        sample_cols = self._sample_cols
        raw_vals = self._matrix_df[sample_cols].to_numpy(dtype=float)
        norm_vals = self._normalized_df[sample_cols].to_numpy(dtype=float)

        # -- Samples: pooled raw vs normalized, log1p-transformed -- #
        raw_flat = raw_vals[~np.isnan(raw_vals)]
        raw_log = np.log1p(np.clip(raw_flat, 0, None))
        norm_flat = norm_vals[~np.isnan(norm_vals)]
        norm_log = np.log1p(np.clip(norm_flat, 0, None))

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        stats_sample_before = _hist_with_gaussian_fit(axes[0], raw_log, "#4C72B0")
        axes[0].set_title("Before normalization")
        stats_sample_after = _hist_with_gaussian_fit(axes[1], norm_log, "#55A868")
        axes[1].set_title("After normalization")
        for ax in axes:
            ax.set_xlabel("log1p(intensity)")
        axes[0].set_ylabel("Density")
        fig.suptitle("Sample intensity distribution (pooled across all samples)")
        sample_plot_path = self._save_plot(fig, "sample_distribution_gaussian.png")

        # -- Features: normalized-but-unscaled vs scaled, pooled -- #
        # Uses the exact same matrix scale_features() receives in
        # pca_and_correlation()/_prepare_supervised() (raw normalized
        # units, NOT log-transformed) so this is a true "before/after"
        # of what PCA/PLS-DA/OPLS-DA actually consume.
        filled = np.nan_to_num(norm_vals, nan=0.0)
        sample_matrix = filled.T  # samples x features
        before_flat = sample_matrix.ravel()
        scaled_matrix = scale_features(sample_matrix, self.inputs.scaling)
        after_flat = scaled_matrix.ravel()

        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        stats_feature_before = _hist_with_gaussian_fit(axes[0], before_flat, "#4C72B0")
        axes[0].set_title("Normalized, before scaling")
        stats_feature_after = _hist_with_gaussian_fit(axes[1], after_flat, "#C44E52")
        axes[1].set_title(f"After {self.inputs.scaling} scaling")
        axes[0].set_xlabel("Normalized intensity")
        axes[1].set_xlabel("Scaled value")
        axes[0].set_ylabel("Density")
        fig.suptitle("Feature-value distribution (pooled across all features)")
        feature_plot_path = self._save_plot(fig, "feature_distribution_gaussian.png")

        self._sections.append(_format_distribution_shape_text(
            stats_sample_before, stats_sample_after,
            stats_feature_before, stats_feature_after,
            sample_plot_path, feature_plot_path, self.inputs.scaling,
        ))

        return {
            "sample_plot": sample_plot_path,
            "feature_plot": feature_plot_path,
            "sample_before": stats_sample_before,
            "sample_after": stats_sample_after,
            "feature_before": stats_feature_before,
            "feature_after": stats_feature_after,
        }

    def _save_plot(self, fig, filename: str) -> str:
        path = os.path.join(self.inputs.output_dir, filename)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        self._plots.append(path)
        return path

    def _lookup_group_labels(self, group_column: str) -> list[Optional[str]]:
        """Maps each sample column to its metadata group label (or None if
        unmatched), via the same loose filename matching used elsewhere."""
        lookup: dict[str, str] = {}
        for _, row in self._metadata_df.iterrows():
            for variant in _filename_variants(str(row[self._metadata_filename_col])):
                lookup[variant] = str(row[group_column])
        labels: list[Optional[str]] = []
        for c in self._sample_cols:
            match = None
            for variant in _filename_variants(c):
                if variant in lookup:
                    match = lookup[variant]
                    break
            labels.append(match)
        return labels

    def attribute_columns(self) -> list[str]:
        """Metadata columns following GNPS's ATTRIBUTE_ naming convention —
        the only columns considered for PLS-DA/OPLS-DA/volcano grouping."""
        if self._metadata_df is None:
            return []
        return [c for c in self._metadata_df.columns if _norm(c).startswith("attribute")]

    def _resolve_group_column(self, explicit: Optional[str]) -> tuple[Optional[str], Optional[str]]:
        """Resolves the grouping column for the supervised analyses.

        Returns (column_name, note). column_name is None when no usable
        column could be resolved, in which case note explains why (and,
        for the multiple-ATTRIBUTE_-columns case, lists the candidates so
        a caller — the CLI's interactive prompt today, a GUI picker later
        — can ask the user to choose one via --group-column / the
        `group_column` field).
        """
        if self._metadata_df is None or self._metadata_filename_col is None:
            return None, "No sample metadata available."

        if explicit:
            for c in self._metadata_df.columns:
                if c == explicit or _norm(c) == _norm(explicit):
                    return c, None
            return None, f"Column '{explicit}' not found in metadata."

        attribute_cols = self.attribute_columns()
        if len(attribute_cols) == 1:
            return attribute_cols[0], None
        if len(attribute_cols) > 1:
            return None, (
                "Multiple ATTRIBUTE_ grouping columns found in metadata "
                f"({', '.join(attribute_cols)}); specify one with --group-column."
            )
        return None, "No ATTRIBUTE_ grouping column found in metadata."

    def _display_group_column(self) -> Optional[str]:
        """Column used to colour the PCA and count groups in the feature
        summary. Uses the same resolution as the supervised analyses so
        every plot is grouped by the same column; falls back to the looser
        Class/Group/Condition/ATTRIBUTE match only when that resolution
        fails (e.g. several ATTRIBUTE_ columns and none chosen)."""
        if self._metadata_df is None:
            return None
        resolved, _ = self._resolve_group_column(self.inputs.group_column)
        return resolved or self._metadata_group_col

    # -- feature summary --------------------------------------------- #

    def feature_summary(self) -> dict:
        df = self._matrix_df
        sample_cols = self._sample_cols
        values = df[sample_cols]

        feature_count = int(df.shape[0])

        total_cells = values.size
        missing_count = int(values.isna().to_numpy().sum())
        missing_pct = (missing_count / total_cells * 100) if total_cells else 0.0
        missing_per_sample = values.isna().sum().to_dict()

        desc = values.describe().T  # per-sample count/mean/std/min/.../max
        intensity_stats = desc[["min", "max", "mean", "50%", "std"]].rename(
            columns={"50%": "median"}
        ).to_dict(orient="index")

        sample_summary: dict = {"samples": list(sample_cols), "n_samples": len(sample_cols)}
        display_col = self._display_group_column()
        if display_col is not None:
            counts = self._metadata_df[display_col].value_counts().to_dict()
            sample_summary["group_column"] = display_col
            sample_summary["group_counts"] = counts

        # Missing values per sample — bar chart
        fig, ax = plt.subplots(figsize=(max(6, len(sample_cols) * 0.4), 4))
        ax.bar(range(len(sample_cols)), [missing_per_sample[c] for c in sample_cols], color="#4C72B0")
        ax.set_xticks(range(len(sample_cols)))
        ax.set_xticklabels(sample_cols, rotation=90, fontsize=6)
        ax.set_ylabel("Missing (NA) feature count")
        ax.set_title("Missing values per sample")
        self._save_plot(fig, "missing_values_per_sample.png")

        # Intensity distribution — log10 histogram of all non-missing, >0 values
        fig, ax = plt.subplots(figsize=(6, 4))
        flat = values.to_numpy().astype(float).ravel()
        flat = flat[~np.isnan(flat) & (flat > 0)]
        if flat.size:
            ax.hist(np.log10(flat), bins=60, color="#55A868")
        ax.set_xlabel("log10(intensity)")
        ax.set_ylabel("Cell count")
        ax.set_title("Intensity distribution (all samples)")
        self._save_plot(fig, "intensity_distribution.png")

        self._sections.append(_format_feature_summary_text(
            feature_count, missing_count, missing_pct, intensity_stats, sample_summary
        ))

        return {
            "feature_count": feature_count,
            "missing_count": missing_count,
            "missing_pct": missing_pct,
            "missing_per_sample": missing_per_sample,
            "intensity_stats": intensity_stats,
            "sample_summary": sample_summary,
        }

    # -- chromatograms ------------------------------------------------ #

    def chromatograms(self) -> list[dict]:
        results = []
        for path in self.inputs.mzml_files:
            if not os.path.isfile(path):
                continue
            data = extract_chromatogram(path)
            name = Path(path).stem
            if data["rt"].size:
                fig, ax = plt.subplots(figsize=(7, 4))
                ax.plot(data["rt"], data["tic"], label="TIC", color="#4C72B0", linewidth=0.8)
                ax.plot(data["rt"], data["bpc"], label="BPC", color="#C44E52", linewidth=0.8)
                ax.set_xlabel("Retention time (s)")
                ax.set_ylabel("Intensity")
                ax.set_title(f"Chromatogram (MS1): {name}")
                ax.legend()
                plot_path = self._save_plot(fig, f"chromatogram_{name}.png")
            else:
                plot_path = None
            results.append({
                "file": path,
                "n_ms1_scans": int(data["rt"].size),
                "plot": plot_path,
            })
        self._sections.append(_format_chromatogram_text(results))
        return results

    # -- PCA + correlation --------------------------------------------- #

    def pca_and_correlation(self) -> dict:
        sample_cols = self._sample_cols
        values = self._normalized_df[sample_cols].to_numpy(dtype=float)

        # Missing values imputed with 0 for PCA/correlation purposes —
        # an explicit, documented assumption (non-detected features are
        # commonly treated as absent/zero at this exploratory stage).
        filled = np.nan_to_num(values, nan=0.0)
        # samples-as-rows for both PCA and correlation
        sample_matrix = filled.T

        n_samples = sample_matrix.shape[0]
        result: dict = {"n_samples": n_samples}

        groups = None
        display_col = self._display_group_column()
        if display_col is not None and self._metadata_filename_col is not None:
            groups = self._lookup_group_labels(display_col)
            result["group_column"] = display_col
            result["n_unmatched"] = sum(g is None for g in groups)

        if n_samples >= 2:
            scaled = scale_features(sample_matrix, self.inputs.scaling)
            # Fit up to 3 components once; the 2D plot uses PC1/PC2 (the
            # same components a 2-component fit would give) and the 3D
            # plot adds PC3 when there are enough samples/features.
            n_components = min(3, n_samples - 1, scaled.shape[1]) or 1
            pca = PCA(n_components=n_components)
            coords = pca.fit_transform(scaled)
            explained = pca.explained_variance_ratio_.tolist()
            result["explained_variance_ratio"] = explained

            # Pad to 3 columns so 1-/2-component fits plot on a flat axis.
            coords3 = np.zeros((n_samples, 3))
            coords3[:, :coords.shape[1]] = coords

            def scatter_by_group(ax, dims):
                def pts(idx):
                    return [coords3[idx, d] for d in dims]

                if groups and any(g is not None for g in groups):
                    color_map = _group_color_map({g for g in groups if g is not None})
                    for g, color in color_map.items():
                        idx = [i for i, gg in enumerate(groups) if gg == g]
                        ax.scatter(*pts(idx), label=str(g), color=color)
                    # Samples with no metadata match still shape the PCA
                    # fit, so draw them too rather than silently dropping them.
                    unmatched_idx = [i for i, gg in enumerate(groups) if gg is None]
                    if unmatched_idx:
                        ax.scatter(*pts(unmatched_idx), label="unmatched", color="#999999", marker="x")
                    ax.legend(title=display_col, fontsize=7)
                else:
                    ax.scatter(*pts(list(range(n_samples))))

            def pc_label(i):
                return f"PC{i+1} ({explained[i]*100:.1f}%)" if len(explained) > i else f"PC{i+1}"

            fig, ax = plt.subplots(figsize=(6, 5))
            scatter_by_group(ax, (0, 1))
            ax.set_xlabel(pc_label(0))
            ax.set_ylabel(pc_label(1))
            ax.set_title("PCA (samples)")
            self._save_plot(fig, "pca.png")

            if coords.shape[1] >= 3:
                fig = plt.figure(figsize=(7, 6))
                ax = fig.add_subplot(projection="3d")
                scatter_by_group(ax, (0, 1, 2))
                ax.set_xlabel(pc_label(0))
                ax.set_ylabel(pc_label(1))
                ax.set_zlabel(pc_label(2))
                # shrink the cube slightly so the z-axis label isn't clipped
                ax.set_box_aspect(None, zoom=0.85)
                ax.set_title("PCA 3D (samples)")
                result["plot_3d"] = self._save_plot(fig, "pca_3d.png")
            else:
                result["note_3d"] = "3D PCA skipped: needs ≥4 samples and ≥3 features."

            # Sample-sample correlation matrix (Pearson, on log1p intensities)
            log_matrix = np.log1p(sample_matrix)
            corr = np.corrcoef(log_matrix)
            corr_df = pd.DataFrame(corr, index=sample_cols, columns=sample_cols)
            corr_df.to_csv(os.path.join(self.inputs.output_dir, "sample_correlation_matrix.csv"))

            fig, ax = plt.subplots(figsize=(max(6, len(sample_cols) * 0.35), max(5, len(sample_cols) * 0.35)))
            im = ax.imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
            ax.set_xticks(range(len(sample_cols)))
            ax.set_xticklabels(sample_cols, rotation=90, fontsize=6)
            ax.set_yticks(range(len(sample_cols)))
            ax.set_yticklabels(sample_cols, fontsize=6)
            ax.set_title("Sample-sample correlation (heatmap)")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            self._save_plot(fig, "correlation_heatmap.png")

            result["correlation_matrix_csv"] = os.path.join(self.inputs.output_dir, "sample_correlation_matrix.csv")
        else:
            result["note"] = "Fewer than 2 samples — PCA/correlation skipped."

        self._sections.append(_format_pca_text(result))
        return result

    # -- shared setup for the supervised analyses ------------------------ #

    def _prepare_supervised(self, group_column: Optional[str]):
        """Resolves the group column and returns (X, labels, unique_groups,
        resolved_column) for samples with a matched group label, or None
        (after recording a skip note) when preconditions aren't met."""
        resolved, note = self._resolve_group_column(group_column)
        if resolved is None:
            return None, note
        self._group_column_used = resolved

        labels_all = self._lookup_group_labels(resolved)
        sample_cols = self._sample_cols
        valid = [i for i, l in enumerate(labels_all) if l is not None]
        if len(valid) < 2:
            return None, f"Fewer than 2 samples matched to a group via '{resolved}'."

        values = self._normalized_df[sample_cols].to_numpy(dtype=float)
        filled = np.nan_to_num(values, nan=0.0)
        sample_matrix = filled.T  # samples x features

        X = sample_matrix[valid]
        labels = [labels_all[i] for i in valid]
        matched_cols = [sample_cols[i] for i in valid]
        unique_groups = sorted(set(labels), key=str)
        if len(unique_groups) < 2:
            return None, f"Only one group ('{unique_groups[0] if unique_groups else '?'}') present via '{resolved}' — need at least 2."

        counts = {g: labels.count(g) for g in unique_groups}
        too_small = [g for g, n in counts.items() if n < 2]
        if too_small:
            return None, (
                f"Group(s) {too_small} via '{resolved}' have fewer than 2 samples — "
                "not enough to compute meaningful statistics."
            )

        return (X, labels, matched_cols, unique_groups, resolved), None

    # -- PLS-DA ----------------------------------------------------------- #

    def pls_da(self, group_column: Optional[str] = None) -> dict:
        prepared, note = self._prepare_supervised(group_column)
        if prepared is None:
            self._sections.append(f"--- PLS-DA ---\nSkipped: {note}\n")
            return {"skipped": True, "reason": note}

        X, labels, matched_cols, unique_groups, resolved = prepared
        X_scaled = scale_features(X, self.inputs.scaling)
        Y = pd.get_dummies(labels).to_numpy(dtype=float)

        n_components = max(1, min(2, X_scaled.shape[0] - 1, X_scaled.shape[1]))
        pls = PLSRegression(n_components=n_components)
        scores, _ = pls.fit_transform(X_scaled, Y)

        fig, ax = plt.subplots(figsize=(6, 5))
        color_map = _group_color_map(unique_groups)
        for g in unique_groups:
            idx = [i for i, l in enumerate(labels) if l == g]
            y_vals = scores[idx, 1] if scores.shape[1] > 1 else np.zeros(len(idx))
            ax.scatter(scores[idx, 0], y_vals, label=str(g), color=color_map[g])
        ax.legend(title=resolved, fontsize=7)
        ax.set_xlabel("LV1")
        ax.set_ylabel("LV2" if scores.shape[1] > 1 else "")
        ax.set_title(f"PLS-DA (grouped by {resolved})")
        plot_path = self._save_plot(fig, "pls_da.png")

        result = {
            "group_column": resolved,
            "groups": unique_groups,
            "n_samples_used": len(labels),
            "n_components": n_components,
            "plot": plot_path,
        }
        self._sections.append(
            f"--- PLS-DA ---\nGrouping column: {resolved} ({', '.join(unique_groups)})\n"
            f"Samples used: {len(labels)} / {len(self._sample_cols)}\nPlot: {plot_path}\n"
        )
        return result

    # -- OPLS-DA (single orthogonal component, two-group only) ----------- #

    def opls_da(self, group_column: Optional[str] = None) -> dict:
        prepared, note = self._prepare_supervised(group_column)
        if prepared is None:
            self._sections.append(f"--- OPLS-DA ---\nSkipped: {note}\n")
            return {"skipped": True, "reason": note}

        X, labels, matched_cols, unique_groups, resolved = prepared
        if len(unique_groups) != 2:
            note = (
                f"OPLS-DA needs exactly 2 groups; '{resolved}' has {len(unique_groups)} "
                f"({', '.join(unique_groups)}). Use PLS-DA for more than 2 groups."
            )
            self._sections.append(f"--- OPLS-DA ---\nSkipped: {note}\n")
            return {"skipped": True, "reason": note}

        group_a, group_b = unique_groups
        y = np.array([0.0 if l == group_a else 1.0 for l in labels])
        y = y - y.mean()

        X_scaled = scale_features(X, self.inputs.scaling)

        # Single-orthogonal-component OPLS (Trygg & Wold, 2002): remove the
        # variation in X that's orthogonal to y, then fit one predictive
        # PLS component on what's left. Gives the classic OPLS-DA scores
        # plot: predictive component (x) vs orthogonal component (y).
        w = X_scaled.T @ y
        w = w / np.linalg.norm(w)
        t = X_scaled @ w
        p = X_scaled.T @ t / (t @ t)
        w_orth = p - ((w @ p) / (w @ w)) * w
        w_orth_norm = np.linalg.norm(w_orth)
        if w_orth_norm < 1e-12:
            note = "Could not separate an orthogonal component (data too collinear with the group labels)."
            self._sections.append(f"--- OPLS-DA ---\nSkipped: {note}\n")
            return {"skipped": True, "reason": note}
        w_orth = w_orth / w_orth_norm
        t_orth = X_scaled @ w_orth
        p_orth = X_scaled.T @ t_orth / (t_orth @ t_orth)
        X_filtered = X_scaled - np.outer(t_orth, p_orth)

        w_pred = X_filtered.T @ y
        w_pred = w_pred / np.linalg.norm(w_pred)
        t_pred = X_filtered @ w_pred

        fig, ax = plt.subplots(figsize=(6, 5))
        color_map = _group_color_map(unique_groups)
        for g in unique_groups:
            idx = [i for i, l in enumerate(labels) if l == g]
            ax.scatter(t_pred[idx], t_orth[idx], label=str(g), color=color_map[g])
        ax.axhline(0, color="grey", linewidth=0.5)
        ax.axvline(0, color="grey", linewidth=0.5)
        ax.legend(title=resolved, fontsize=7)
        ax.set_xlabel("Predictive component")
        ax.set_ylabel("Orthogonal component")
        ax.set_title(f"OPLS-DA: {group_a} vs {group_b}")
        plot_path = self._save_plot(fig, "opls_da.png")

        result = {
            "group_column": resolved,
            "groups": unique_groups,
            "n_samples_used": len(labels),
            "plot": plot_path,
        }
        self._sections.append(
            f"--- OPLS-DA ---\nGrouping column: {resolved} ({group_a} vs {group_b})\n"
            f"Samples used: {len(labels)} / {len(self._sample_cols)}\nPlot: {plot_path}\n"
        )
        return result

    # -- Volcano plot (two-group only) ------------------------------------ #

    def volcano_plot(self, group_column: Optional[str] = None) -> dict:
        prepared, note = self._prepare_supervised(group_column)
        if prepared is None:
            self._sections.append(f"--- Volcano plot ---\nSkipped: {note}\n")
            return {"skipped": True, "reason": note}

        _, labels, matched_cols, unique_groups, resolved = prepared
        if len(unique_groups) != 2:
            note = (
                f"Volcano plot needs exactly 2 groups; '{resolved}' has {len(unique_groups)} "
                f"({', '.join(unique_groups)})."
            )
            self._sections.append(f"--- Volcano plot ---\nSkipped: {note}\n")
            return {"skipped": True, "reason": note}

        group_a, group_b = unique_groups
        a_cols = [c for c, l in zip(matched_cols, labels) if l == group_a]
        b_cols = [c for c, l in zip(matched_cols, labels) if l == group_b]

        # Normalized-but-unscaled values: normalization corrects for overall
        # per-sample loading differences (the whole point, given fold-change
        # should reflect real biology, not a batch-scale artifact), while
        # scaling (pareto/auto) is deliberately NOT applied here — it would
        # make "fold change" meaningless (scaling can shift values negative
        # and rescales each feature independently of its original units).
        a_vals = self._normalized_df[a_cols].to_numpy(dtype=float)
        b_vals = self._normalized_df[b_cols].to_numpy(dtype=float)

        mean_a = np.nanmean(a_vals, axis=1)
        mean_b = np.nanmean(b_vals, axis=1)
        # small pseudocount (relative to the smallest nonzero value seen)
        # so log2 fold change is defined even for zero/near-zero means.
        nonzero = np.concatenate([a_vals[a_vals > 0], b_vals[b_vals > 0]])
        eps = float(nonzero.min()) / 10 if nonzero.size else 1e-9
        log2fc = np.log2((mean_b + eps) / (mean_a + eps))

        with np.errstate(invalid="ignore"):
            t_stat, p_val = scipy_stats.ttest_ind(a_vals, b_vals, axis=1, equal_var=False, nan_policy="omit")
        p_val = np.asarray(p_val, dtype=float)
        q_val = _benjamini_hochberg(p_val)
        neg_log10_p = -np.log10(np.clip(p_val, 1e-300, None))

        feature_ids = self._matrix_df[self._feature_id_col].to_numpy()
        results_df = pd.DataFrame({
            self._feature_id_col: feature_ids,
            "log2fc": log2fc,
            "p_value": p_val,
            "q_value": q_val,
        })
        csv_path = os.path.join(self.inputs.output_dir, "volcano_plot_results.csv")
        results_df.to_csv(csv_path, index=False)

        sig = (p_val < 0.05) & (np.abs(log2fc) > 1)
        # Colour significant features by the group they are higher in,
        # using the same per-group colours as the PCA/PLS-DA/OPLS-DA plots.
        higher_b = sig & (log2fc > 0)
        higher_a = sig & (log2fc < 0)
        n_higher_in = {group_a: int(higher_a.sum()), group_b: int(higher_b.sum())}
        color_map = _group_color_map(unique_groups)
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(log2fc[~sig], neg_log10_p[~sig], color="#999999", s=10, label="not significant")
        ax.scatter(log2fc[higher_b], neg_log10_p[higher_b], color=color_map[group_b], s=10,
                   label=f"higher in {group_b} ({n_higher_in[group_b]})")
        ax.scatter(log2fc[higher_a], neg_log10_p[higher_a], color=color_map[group_a], s=10,
                   label=f"higher in {group_a} ({n_higher_in[group_a]})")
        ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.7)
        ax.axvline(1, color="grey", linestyle="--", linewidth=0.7)
        ax.axvline(-1, color="grey", linestyle="--", linewidth=0.7)
        ax.set_xlabel(f"log2 fold change ({group_b} / {group_a})")
        ax.set_ylabel("-log10(p-value)")
        ax.set_title(f"Volcano plot: {group_b} vs {group_a}")
        ax.legend(fontsize=7)
        plot_path = self._save_plot(fig, "volcano_plot.png")

        result = {
            "group_column": resolved,
            "groups": unique_groups,
            "n_significant": int(sig.sum()),
            "n_higher_in": n_higher_in,
            "csv": csv_path,
            "plot": plot_path,
        }
        self._sections.append(
            f"--- Volcano plot ---\nGrouping column: {resolved} ({group_b} vs {group_a})\n"
            f"Significant features (p<0.05 & |log2FC|>1): {int(sig.sum())} / {len(feature_ids)} "
            f"(higher in {group_b}: {n_higher_in[group_b]}, higher in {group_a}: {n_higher_in[group_a]})\n"
            f"Test: Welch's t-test (unequal variance), Benjamini-Hochberg q-values in {csv_path}\n"
            f"Plot: {plot_path}\n"
        )
        return result

    # -- Feature overlap (Venn for 2-3 groups, UpSet for 4+) ------------- #

    def feature_overlap(self, group_column: Optional[str] = None) -> dict:
        """Which features each group captures, and which are shared or
        unique. A feature is captured by a group when it is detected
        (> 0, not NA) in at least `venn_min_fraction` of that group's
        samples. Uses raw values: detection is a QC property, and
        normalization doesn't change zeros/NAs anyway."""
        def skip(note):
            self._sections.append(f"--- Feature overlap (Venn/UpSet) ---\nSkipped: {note}\n")
            return {"skipped": True, "reason": note}

        resolved, note = self._resolve_group_column(group_column)
        if resolved is None:
            return skip(note)
        labels = self._lookup_group_labels(resolved)
        groups = sorted({l for l in labels if l is not None}, key=str)
        if len(groups) < 2:
            return skip(
                f"Need at least 2 groups with matched samples via '{resolved}'; found {len(groups)}."
            )

        min_fraction = self.inputs.venn_min_fraction
        raw = self._matrix_df[self._sample_cols].to_numpy(dtype=float)
        detected = np.nan_to_num(raw, nan=0.0) > 0  # features x samples
        captured: dict = {}
        for g in groups:
            idx = [i for i, l in enumerate(labels) if l == g]
            # small tolerance so e.g. 7/10 >= 0.7 isn't lost to float error
            captured[g] = detected[:, idx].mean(axis=1) >= min_fraction - 1e-12

        n_features = raw.shape[0]
        memberships = [frozenset(g for g in groups if captured[g][i]) for i in range(n_features)]
        set_counts: dict = {}
        for m in memberships:
            if m:
                set_counts[m] = set_counts.get(m, 0) + 1
        ordered_sets = sorted(set_counts, key=lambda s: (-set_counts[s], len(s), sorted(s, key=str)))

        def label(s):
            if not s:
                return "none"
            members = [g for g in groups if g in s]
            return f"{members[0]} only" if len(members) == 1 else " & ".join(members)

        overlap_df = pd.DataFrame({self._feature_id_col: self._matrix_df[self._feature_id_col].to_numpy()})
        for g in groups:
            overlap_df[f"captured_{g}"] = captured[g]
        overlap_df["overlap"] = [label(m) for m in memberships]
        csv_path = os.path.join(self.inputs.output_dir, "feature_overlap.csv")
        overlap_df.to_csv(csv_path, index=False)

        captured_per_group = {g: int(captured[g].sum()) for g in groups}
        unique_per_group = {g: set_counts.get(frozenset([g]), 0) for g in groups}
        shared_all = set_counts.get(frozenset(groups), 0)
        n_none = sum(1 for m in memberships if not m)
        color_map = _group_color_map(groups)

        if len(groups) <= 3:
            fig = _draw_venn(groups, set_counts, captured_per_group, color_map)
            fig.suptitle(
                f"Feature overlap by {resolved} (detected in ≥{min_fraction:.0%} of group samples)\n"
                "circle sizes are not proportional to counts",
                fontsize=9,
            )
            plot_path = self._save_plot(fig, "venn_diagram.png")
            truncated = 0
        else:
            fig, truncated = _draw_upset(groups, ordered_sets, set_counts, captured_per_group, color_map)
            fig.suptitle(
                f"Feature overlap by {resolved} (detected in ≥{min_fraction:.0%} of group samples)",
                fontsize=9,
            )
            plot_path = self._save_plot(fig, "upset_plot.png")

        result = {
            "group_column": resolved,
            "groups": groups,
            "min_fraction": min_fraction,
            "captured_per_group": captured_per_group,
            "unique_per_group": unique_per_group,
            "shared_all": shared_all,
            "none": n_none,
            "intersections": {label(s): set_counts[s] for s in ordered_sets},
            "plot": plot_path,
            "csv": csv_path,
        }

        lines = [
            "--- Feature overlap (Venn/UpSet) ---",
            f"Grouping column: {resolved} ({', '.join(groups)})",
            f"Captured = detected (>0, not NA) in ≥{min_fraction:.0%} of a group's samples (raw values)",
        ]
        for g in groups:
            lines.append(f"  {g}: {captured_per_group[g]} captured, {unique_per_group[g]} unique")
        lines.append(f"Shared by all groups: {shared_all}")
        lines.append(f"Captured by no group: {n_none}")
        if truncated:
            lines.append(f"UpSet plot shows the 30 largest intersections ({truncated} smaller ones omitted; all are in the CSV)")
        lines.append(f"Plot: {plot_path}")
        lines.append(f"Membership table: {csv_path}")
        lines.append("")
        self._sections.append("\n".join(lines))
        return result

    # -- report --------------------------------------------------------- #

    def build_report(self) -> str:
        lines = [
            "=" * 70,
            "METABOLOMICS DATA — VISUALIZATION (UNANNOTATED) MODULE REPORT",
            "=" * 70,
            f"Run at: {datetime.now().isoformat(timespec='seconds')}",
            "",
        ]
        lines.extend(self._sections)
        lines.append("-" * 70)
        lines.append("Generated plots:")
        for p in self._plots:
            lines.append(f"  - {p}")
        lines.append("=" * 70)
        return "\n".join(lines)

    def run(self) -> str:
        self.load()
        self.feature_summary()
        self.chromatograms()
        self.normalization_summary()
        self.distribution_shape_summary()
        self.pca_and_correlation()
        self.pls_da(self.inputs.group_column)
        self.opls_da(self.inputs.group_column)
        self.volcano_plot(self.inputs.group_column)
        self.feature_overlap(self.inputs.group_column)
        report = self.build_report()
        report_path = os.path.join(self.inputs.output_dir, "analysis_report.txt")
        Path(report_path).write_text(report + "\n")
        return report


# --------------------------------------------------------------------- #
# Feature-overlap drawing (static; interactive versions deferred to GUI)
# --------------------------------------------------------------------- #

# Fixed (not area-proportional) layouts: circle centres, radius 1, and the
# text position for each exclusive region, keyed by member-index tuple.
_VENN_LAYOUTS = {
    2: {
        "centers": [(-0.55, 0.0), (0.55, 0.0)],
        "labels": [(-0.9, 1.2), (0.9, 1.2)],
        "regions": {(0,): (-0.95, 0.0), (1,): (0.95, 0.0), (0, 1): (0.0, 0.0)},
        "ylim": (-1.2, 1.5),
    },
    3: {
        "centers": [(-0.55, 0.35), (0.55, 0.35), (0.0, -0.6)],
        "labels": [(-1.25, 1.55), (1.25, 1.55), (0.0, -1.85)],
        "regions": {
            (0,): (-1.0, 0.65), (1,): (1.0, 0.65), (2,): (0.0, -1.15),
            (0, 1): (0.0, 0.85), (0, 2): (-0.62, -0.35), (1, 2): (0.62, -0.35),
            (0, 1, 2): (0.0, 0.05),
        },
        "ylim": (-2.1, 1.9),
    },
}


def _draw_venn(groups, set_counts, captured_per_group, color_map):
    layout = _VENN_LAYOUTS[len(groups)]
    ymin, ymax = layout["ylim"]
    fig, ax = plt.subplots(figsize=(6, 6 * (ymax - ymin) / 4.0))
    for g, (cx, cy) in zip(groups, layout["centers"]):
        ax.add_patch(Circle((cx, cy), 1.0, facecolor=color_map[g], edgecolor=color_map[g], alpha=0.35, linewidth=1.5))
    for g, (lx, ly) in zip(groups, layout["labels"]):
        ax.text(lx, ly, f"{g}\n({captured_per_group[g]} captured)", ha="center", va="center",
                fontsize=9, color=color_map[g], fontweight="bold")
    for r in range(1, len(groups) + 1):
        for members in combinations(range(len(groups)), r):
            count = set_counts.get(frozenset(groups[i] for i in members), 0)
            x, y = layout["regions"][members]
            ax.text(x, y, str(count), ha="center", va="center", fontsize=11)
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig


def _draw_upset(groups, ordered_sets, set_counts, captured_per_group, color_map, max_bars: int = 30):
    """UpSet-style chart: bars of exclusive-intersection sizes (largest
    first, capped at max_bars), a dot matrix of which groups form each
    intersection, and per-group totals. Returns (fig, n_omitted)."""
    shown = ordered_sets[:max_bars]
    n_omitted = len(ordered_sets) - len(shown)
    n_rows = len(groups)

    fig = plt.figure(figsize=(max(7, 0.4 * len(shown) + 3), 4 + 0.35 * n_rows))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[3, max(1.2, 0.35 * n_rows)],
                          wspace=0.05, hspace=0.05)
    ax_bars = fig.add_subplot(gs[0, 1])
    ax_dots = fig.add_subplot(gs[1, 1], sharex=ax_bars)
    ax_sets = fig.add_subplot(gs[1, 0], sharey=ax_dots)

    xs = np.arange(len(shown))
    heights = [set_counts[s] for s in shown]
    ax_bars.bar(xs, heights, color="#4C4C4C")
    for x, h in zip(xs, heights):
        ax_bars.text(x, h, str(h), ha="center", va="bottom", fontsize=7)
    ax_bars.set_ylabel("Features in intersection")
    ax_bars.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in ("top", "right"):
        ax_bars.spines[spine].set_visible(False)

    ys = np.arange(n_rows)
    for x, s in zip(xs, shown):
        filled = [i for i, g in enumerate(groups) if g in s]
        ax_dots.scatter([x] * n_rows, ys, color="#DDDDDD", s=40, zorder=1)
        ax_dots.scatter([x] * len(filled), filled, color=[color_map[groups[i]] for i in filled], s=40, zorder=3)
        if len(filled) > 1:
            ax_dots.plot([x, x], [min(filled), max(filled)], color="#4C4C4C", linewidth=1.5, zorder=2)
    ax_dots.set_yticks(ys)
    ax_dots.set_yticklabels(groups, fontsize=8)
    ax_dots.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_dots.tick_params(axis="y", left=False, labelleft=False, labelright=True)
    ax_dots.set_ylim(-0.5, n_rows - 0.5)
    ax_dots.invert_yaxis()
    for spine in ax_dots.spines.values():
        spine.set_visible(False)

    ax_sets.barh(ys, [captured_per_group[g] for g in groups], color=[color_map[g] for g in groups])
    ax_sets.invert_xaxis()
    ax_sets.set_xlabel("Captured", fontsize=8)
    ax_sets.tick_params(axis="y", left=False, labelleft=False)
    for spine in ("top", "left"):
        ax_sets.spines[spine].set_visible(False)

    if n_omitted:
        ax_bars.set_title(f"{n_omitted} smaller intersection(s) not shown", fontsize=7)
    return fig, n_omitted


# --------------------------------------------------------------------- #
# Text formatting helpers
# --------------------------------------------------------------------- #

def _format_feature_summary_text(feature_count, missing_count, missing_pct, intensity_stats, sample_summary) -> str:
    lines = ["--- Feature summary ---", f"Feature count: {feature_count}"]
    lines.append(f"Missing values: {missing_count} ({missing_pct:.2f}% of all peak cells)")
    lines.append(f"Samples: {sample_summary['n_samples']}")
    if "group_column" in sample_summary:
        lines.append(f"Grouping column: {sample_summary['group_column']}")
        for g, n in sample_summary["group_counts"].items():
            lines.append(f"  {g}: {n}")
    lines.append("Intensity distribution (per sample min/max/mean/median/std):")
    for sample, stats in list(intensity_stats.items())[:10]:
        lines.append(
            f"  {sample}: min={stats['min']:.1f} max={stats['max']:.1f} "
            f"mean={stats['mean']:.1f} median={stats['median']:.1f} std={stats['std']:.1f}"
        )
    if len(intensity_stats) > 10:
        lines.append(f"  ... and {len(intensity_stats) - 10} more sample(s)")
    lines.append("")
    return "\n".join(lines)


def _format_chromatogram_text(results: list[dict]) -> str:
    lines = ["--- Chromatograms (TIC/BPC) ---"]
    for r in results:
        lines.append(f"{r['file']}: {r['n_ms1_scans']} MS1 scans -> {r['plot']}")
    lines.append("")
    return "\n".join(lines)


def _format_distribution_shape_text(
    stats_sample_before, stats_sample_after,
    stats_feature_before, stats_feature_after,
    sample_plot_path, feature_plot_path, scaling,
) -> str:
    def fmt(stats):
        if not stats:
            return "n/a"
        skew = f"{stats['skew']:.2f}" if stats["skew"] is not None else "n/a"
        return f"mean={stats['mean']:.2f}, std={stats['std']:.2f}, skew={skew} (n={stats['n']})"

    lines = ["--- Distribution shape (Gaussian fit) ---"]
    lines.append("Sample intensities, pooled and log1p-transformed (see normalization_effect.png for the per-sample view):")
    lines.append(f"  Before normalization: {fmt(stats_sample_before)}")
    lines.append(f"  After normalization:  {fmt(stats_sample_after)}")
    lines.append(f"  Plot: {sample_plot_path}")
    lines.append("Feature values, pooled (input to PCA/PLS-DA/OPLS-DA, not log-transformed):")
    lines.append(f"  Normalized, before {scaling} scaling: {fmt(stats_feature_before)}")
    lines.append(f"  After {scaling} scaling:              {fmt(stats_feature_after)}")
    lines.append(f"  Plot: {feature_plot_path}")
    lines.append(
        "Note: skewness near 0 is what a Gaussian/bell shape corresponds to numerically; "
        "these are visual/diagnostic fits, not a formal normality test (e.g. Shapiro-Wilk)."
    )
    lines.append("")
    return "\n".join(lines)


def _format_pca_text(result: dict) -> str:
    lines = ["--- PCA / correlation ---"]
    if "note" in result:
        lines.append(result["note"])
    else:
        ev = result.get("explained_variance_ratio", [])
        ev_str = ", ".join(f"PC{i+1}={v*100:.1f}%" for i, v in enumerate(ev))
        lines.append(f"Explained variance: {ev_str}")
        if "plot_3d" in result:
            lines.append(f"3D PCA: {result['plot_3d']}")
        elif "note_3d" in result:
            lines.append(result["note_3d"])
        if "group_column" in result:
            lines.append(f"PCA coloured by: {result['group_column']}")
            if result.get("n_unmatched"):
                lines.append(
                    f"  {result['n_unmatched']} sample(s) not matched to metadata "
                    "(shown as grey 'unmatched' points)"
                )
        lines.append(
            "Correlation: Pearson on log1p of normalized, unscaled intensities "
            f"-> {result.get('correlation_matrix_csv')}"
        )
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------- #

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Feature summary, chromatograms, PCA and correlation heatmap "
                    "for a validated metabolomics dataset (mzML + metadata + data matrix)."
    )
    parser.add_argument("--mzml", nargs="+", metavar="FILE", default=[], help="One or more .mzML files.")
    parser.add_argument("--metadata", metavar="FILE", help="Sample metadata table.")
    parser.add_argument("--data-matrix", metavar="FILE", required=True, help="Feature x sample intensity matrix.")
    parser.add_argument("--output-dir", metavar="DIR", default="viz_output", help="Where plots/report are written.")
    parser.add_argument(
        "--group-column", metavar="COLUMN",
        help="Metadata column defining groups for PCA colouring and "
             "PLS-DA/OPLS-DA/volcano plot (any column name; normally an "
             "ATTRIBUTE_ column). If omitted, a single ATTRIBUTE_ column is used "
             "automatically; with more than one you'll be prompted to choose.",
    )
    parser.add_argument(
        "--normalization", choices=["none", "sum", "median"], default="sum",
        help="Sample-wise normalization before any cross-sample comparison "
             "(default: sum, a.k.a. total ion current normalization).",
    )
    parser.add_argument(
        "--scaling", choices=["none", "pareto", "auto"], default="pareto",
        help="Feature-wise scaling for PCA/PLS-DA/OPLS-DA only, "
             "applied on top of normalization (default: pareto, the "
             "metabolomics-standard choice). Not applied to the correlation "
             "matrix or the volcano plot, which use normalized-but-unscaled values.",
    )
    parser.add_argument(
        "--venn-min-fraction", type=_fraction, default=0.5, metavar="FRACTION",
        help="Feature-overlap (Venn/UpSet) capture rule: a feature counts as "
             "captured by a group when detected (>0, not NA) in at least this "
             "fraction of the group's samples, in (0, 1] (default: 0.5).",
    )
    return parser


def _fraction(value: str) -> float:
    try:
        f = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {value!r}")
    if not (0 < f <= 1):
        raise argparse.ArgumentTypeError(f"must be in (0, 1], got {f}")
    return f


def _prompt_for_group_column(metadata_file: Optional[str], explicit: Optional[str]) -> Optional[str]:
    """Interactively asks which ATTRIBUTE_ column to use when metadata has
    more than one and none was given via --group-column. This is the
    script-era stand-in for what a future GUI would show as a picker —
    the GUI would just set VisualizationInputs.group_column directly and
    never reach this prompt."""
    if explicit or not metadata_file:
        return explicit
    df = _read_table(metadata_file)
    if df is None:
        return explicit
    attribute_cols = [c for c in df.columns if _norm(c).startswith("attribute")]
    if len(attribute_cols) <= 1:
        return explicit
    if not sys.stdin.isatty():
        return explicit  # non-interactive context — let the module report the ambiguity
    print("Multiple ATTRIBUTE_ grouping columns found in metadata:")
    for i, c in enumerate(attribute_cols, 1):
        print(f"  {i}. {c}")
    choice = input(
        f"Select a grouping column for PLS-DA/OPLS-DA/volcano plot [1-{len(attribute_cols)}], "
        "or press Enter to skip those analyses: "
    ).strip()
    if not choice:
        return explicit
    try:
        return attribute_cols[int(choice) - 1]
    except (ValueError, IndexError):
        print("Invalid selection — skipping PLS-DA/OPLS-DA/volcano plot.")
        return explicit


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    group_column = _prompt_for_group_column(args.metadata, args.group_column)

    inputs = VisualizationInputs(
        mzml_files=args.mzml,
        metadata_file=args.metadata,
        data_matrix_file=args.data_matrix,
        output_dir=args.output_dir,
        group_column=group_column,
        normalization=args.normalization,
        scaling=args.scaling,
        venn_min_fraction=args.venn_min_fraction,
    )
    module = VisualizationModule(inputs)
    report = module.run()
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
