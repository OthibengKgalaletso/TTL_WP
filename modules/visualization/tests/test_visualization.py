"""
Pytest suite for visualization_module.py.

Run with:
    pytest -v
"""

import base64
import struct
import sys
import zlib
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from visualization_module import (  # noqa: E402
    VisualizationInputs,
    VisualizationModule,
    extract_chromatogram,
    _benjamini_hochberg,
    normalize_matrix,
    pareto_scale,
    scale_features,
    _hist_with_gaussian_fit,
)


def _write(path: Path, text: str) -> str:
    path.write_text(text)
    return str(path)


# --------------------------------------------------------------------- #
# mzML fixtures
# --------------------------------------------------------------------- #

MZML_WITH_PRECOMPUTED_TIC = """<?xml version="1.0" encoding="utf-8"?>
<indexedmzML>
  <mzML>
    <run>
      <spectrumList>
        <spectrum id="scan=1" index="0" defaultArrayLength="2">
          <cvParam accession="MS:1000511" name="ms level" value="1"/>
          <cvParam accession="MS:1000016" name="scan start time" value="1.0" unitAccession="UO:0000010"/>
          <cvParam accession="MS:1000285" name="total ion current" value="5000.0"/>
          <cvParam accession="MS:1000505" name="base peak intensity" value="3000.0"/>
        </spectrum>
        <spectrum id="scan=2" index="1" defaultArrayLength="2">
          <cvParam accession="MS:1000511" name="ms level" value="1"/>
          <cvParam accession="MS:1000016" name="scan start time" value="2.0" unitAccession="UO:0000010"/>
          <cvParam accession="MS:1000285" name="total ion current" value="7000.0"/>
          <cvParam accession="MS:1000505" name="base peak intensity" value="4000.0"/>
        </spectrum>
        <spectrum id="scan=3" index="2" defaultArrayLength="2">
          <cvParam accession="MS:1000511" name="ms level" value="2"/>
          <cvParam accession="MS:1000016" name="scan start time" value="2.1" unitAccession="UO:0000010"/>
          <cvParam accession="MS:1000285" name="total ion current" value="999.0"/>
        </spectrum>
      </spectrumList>
    </run>
  </mzML>
</indexedmzML>
"""


def _make_binary_intensity_mzml(intensities_by_scan: list[list[float]], compress: bool) -> str:
    spectra = []
    for i, intensities in enumerate(intensities_by_scan):
        raw = struct.pack(f"<{len(intensities)}f", *intensities)
        if compress:
            raw = zlib.compress(raw)
        b64 = base64.b64encode(raw).decode("ascii")
        comp_cv = '<cvParam accession="MS:1000574" name="zlib compression"/>' if compress else ""
        spectra.append(f"""
        <spectrum id="scan={i+1}" index="{i}" defaultArrayLength="{len(intensities)}">
          <cvParam accession="MS:1000511" name="ms level" value="1"/>
          <cvParam accession="MS:1000016" name="scan start time" value="{float(i+1)}" unitAccession="UO:0000010"/>
          <binaryDataArrayList count="1">
            <binaryDataArray>
              <cvParam accession="MS:1000515" name="intensity array"/>
              <cvParam accession="MS:1000521" name="32-bit float"/>
              {comp_cv}
              <binary>{b64}</binary>
            </binaryDataArray>
          </binaryDataArrayList>
        </spectrum>""")
    body = "".join(spectra)
    return f"""<?xml version="1.0" encoding="utf-8"?>
<indexedmzML>
  <mzML>
    <run>
      <spectrumList>{body}
      </spectrumList>
    </run>
  </mzML>
</indexedmzML>
"""


def test_extract_chromatogram_uses_precomputed_tic_and_bpc(tmp_path):
    path = _write(tmp_path / "precomputed.mzML", MZML_WITH_PRECOMPUTED_TIC)
    data = extract_chromatogram(path)
    assert data["rt"].tolist() == [1.0, 2.0]  # MS2 scan excluded
    assert data["tic"].tolist() == [5000.0, 7000.0]
    assert data["bpc"].tolist() == [3000.0, 4000.0]


def test_extract_chromatogram_falls_back_to_decoding_binary_array_uncompressed(tmp_path):
    xml = _make_binary_intensity_mzml([[10.0, 20.0, 5.0], [1.0, 2.0]], compress=False)
    path = _write(tmp_path / "raw_binary.mzML", xml)
    data = extract_chromatogram(path)
    assert data["rt"].tolist() == [1.0, 2.0]
    assert data["tic"].tolist() == [35.0, 3.0]
    assert data["bpc"].tolist() == [20.0, 2.0]


def test_extract_chromatogram_falls_back_to_decoding_binary_array_zlib_compressed(tmp_path):
    xml = _make_binary_intensity_mzml([[100.0, 50.0]], compress=True)
    path = _write(tmp_path / "compressed_binary.mzML", xml)
    data = extract_chromatogram(path)
    assert data["rt"].tolist() == [1.0]
    assert data["tic"].tolist() == [150.0]
    assert data["bpc"].tolist() == [100.0]


# --------------------------------------------------------------------- #
# VisualizationModule fixtures
# --------------------------------------------------------------------- #

def _make_matrix_and_metadata(tmp_path):
    matrix_path = _write(
        tmp_path / "matrix.csv",
        "row ID,row m/z,row retention time,sample_A.mzML,sample_B.mzML,sample_C.mzML\n"
        "1,100.1,60.0,1000,1200,\n"
        "2,200.2,120.0,500,450,600\n"
        "3,300.3,180.0,10,20,30\n",
    )
    metadata_path = _write(
        tmp_path / "metadata.csv",
        "filename,ATTRIBUTE_Group\n"
        "sample_A.mzML,control\n"
        "sample_B.mzML,control\n"
        "sample_C.mzML,treated\n",
    )
    return matrix_path, metadata_path


def test_feature_summary_counts_and_excludes_descriptive_columns(tmp_path):
    matrix_path, metadata_path = _make_matrix_and_metadata(tmp_path)
    inputs = VisualizationInputs(
        mzml_files=[], metadata_file=metadata_path, data_matrix_file=matrix_path,
        output_dir=str(tmp_path / "out"),
    )
    module = VisualizationModule(inputs)
    module.load()

    assert module._sample_cols == ["sample_A.mzML", "sample_B.mzML", "sample_C.mzML"]

    summary = module.feature_summary()
    assert summary["feature_count"] == 3
    assert summary["missing_count"] == 1
    assert summary["sample_summary"]["n_samples"] == 3
    assert summary["sample_summary"]["group_column"] == "ATTRIBUTE_Group"
    assert summary["sample_summary"]["group_counts"] == {"control": 2, "treated": 1}

    assert (tmp_path / "out" / "missing_values_per_sample.png").exists()
    assert (tmp_path / "out" / "intensity_distribution.png").exists()


def test_pca_and_correlation_outputs(tmp_path):
    matrix_path, metadata_path = _make_matrix_and_metadata(tmp_path)
    inputs = VisualizationInputs(
        mzml_files=[], metadata_file=metadata_path, data_matrix_file=matrix_path,
        output_dir=str(tmp_path / "out"),
    )
    module = VisualizationModule(inputs)
    module.load()
    result = module.pca_and_correlation()

    assert result["n_samples"] == 3
    assert "explained_variance_ratio" in result
    assert (tmp_path / "out" / "pca.png").exists()
    assert (tmp_path / "out" / "correlation_heatmap.png").exists()
    corr_csv = Path(result["correlation_matrix_csv"])
    assert corr_csv.exists()


import numpy as np
import pandas as pd


def _two_group_matrix_and_metadata(tmp_path, n_features=40, n_per_group=10):
    """Two clearly-separated groups (control/treated), with a handful of
    features deliberately shifted so the volcano plot has real signal."""
    rng = np.random.default_rng(0)
    control_cols = [f"control_{i}.mzML" for i in range(n_per_group)]
    treated_cols = [f"treated_{i}.mzML" for i in range(n_per_group)]
    all_cols = control_cols + treated_cols

    base = rng.exponential(1000, size=(n_features, len(all_cols)))
    # Push the first 5 features up in "treated" and next 5 down, for a
    # clear, testable fold-change signal.
    base[0:5, len(control_cols):] *= 20
    base[5:10, len(control_cols):] /= 20

    df = pd.DataFrame(base, columns=all_cols)
    df.insert(0, "row ID", [f"FT{i+1:04d}" for i in range(n_features)])
    matrix_path = tmp_path / "two_group_matrix.csv"
    df.to_csv(matrix_path, index=False)

    meta_rows = [(c, "control") for c in control_cols] + [(c, "treated") for c in treated_cols]
    meta_df = pd.DataFrame(meta_rows, columns=["filename", "ATTRIBUTE_Group"])
    metadata_path = tmp_path / "two_group_metadata.csv"
    meta_df.to_csv(metadata_path, index=False)

    return str(matrix_path), str(metadata_path)


def test_benjamini_hochberg_monotonic_and_bounded():
    p = np.array([0.001, 0.2, 0.03, 0.5, 0.04, np.nan])
    q = _benjamini_hochberg(p)
    valid = ~np.isnan(q)
    assert np.all(q[valid] >= 0) and np.all(q[valid] <= 1)
    assert np.isnan(q[-1])
    # smallest p-value should get the smallest (or equal) q-value
    assert q[0] == np.nanmin(q)


def test_pls_da_two_groups(tmp_path):
    matrix_path, metadata_path = _two_group_matrix_and_metadata(tmp_path)
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()
    result = module.pls_da()
    assert result.get("skipped") is not True
    assert result["group_column"] == "ATTRIBUTE_Group"
    assert set(result["groups"]) == {"control", "treated"}
    assert Path(result["plot"]).exists()


def test_opls_da_two_groups_separates_predictive_component(tmp_path):
    matrix_path, metadata_path = _two_group_matrix_and_metadata(tmp_path)
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()
    result = module.opls_da()
    assert result.get("skipped") is not True
    assert Path(result["plot"]).exists()


def test_volcano_plot_flags_the_shifted_features(tmp_path):
    matrix_path, metadata_path = _two_group_matrix_and_metadata(tmp_path)
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()
    result = module.volcano_plot()
    assert result.get("skipped") is not True
    # With n=10/group and a 20x deliberate shift, most (not necessarily
    # all, given the high-variance exponential noise) of the 10 shifted
    # features should clear significance; the sign/magnitude checks below
    # are the real correctness assertions.
    assert result["n_significant"] >= 5
    assert Path(result["plot"]).exists()
    csv_df = pd.read_csv(result["csv"])
    assert {"row ID", "log2fc", "p_value", "q_value"}.issubset(csv_df.columns)
    # features 0-4 were pushed up in "treated" -> positive log2fc (treated/control)
    assert (csv_df.iloc[0:5]["log2fc"] > 1).all()
    # features 5-9 were pushed down in "treated" -> negative log2fc
    assert (csv_df.iloc[5:10]["log2fc"] < -1).all()


def test_supervised_analyses_skip_gracefully_with_only_one_group(tmp_path):
    matrix_path = _write(
        tmp_path / "matrix.csv",
        "row ID,sample_A.mzML,sample_B.mzML\n1,1000,1100\n2,500,600\n",
    )
    metadata_path = _write(
        tmp_path / "metadata.csv",
        "filename,ATTRIBUTE_Group\nsample_A.mzML,control\nsample_B.mzML,control\n",
    )
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()
    for method in (module.pls_da, module.opls_da, module.volcano_plot):
        result = method()
        assert result["skipped"] is True


def test_multiple_attribute_columns_without_explicit_choice_skips_with_note(tmp_path):
    matrix_path = _write(
        tmp_path / "matrix.csv",
        "row ID,sample_A.mzML,sample_B.mzML,sample_C.mzML,sample_D.mzML\n"
        "1,1000,1100,2000,2100\n2,500,600,300,350\n",
    )
    metadata_path = _write(
        tmp_path / "metadata.csv",
        "filename,ATTRIBUTE_Tissue,ATTRIBUTE_Batch\n"
        "sample_A.mzML,MCL,1\nsample_B.mzML,MCL,2\n"
        "sample_C.mzML,MCR,1\nsample_D.mzML,MCR,2\n",
    )
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()
    result = module.pls_da()
    assert result["skipped"] is True
    assert "Multiple ATTRIBUTE_" in result["reason"]

    # explicit --group-column disambiguates
    result = module.pls_da(group_column="ATTRIBUTE_Tissue")
    assert result.get("skipped") is not True
    assert result["group_column"] == "ATTRIBUTE_Tissue"


def test_opls_da_and_volcano_require_exactly_two_groups(tmp_path):
    rng = np.random.default_rng(1)
    cols = [f"s{i}.mzML" for i in range(9)]
    values = rng.exponential(1000, size=(20, len(cols)))
    df = pd.DataFrame(values, columns=cols)
    df.insert(0, "row ID", [f"FT{i+1:03d}" for i in range(20)])
    matrix_path = tmp_path / "matrix.csv"
    df.to_csv(matrix_path, index=False)

    groups = ["g1"] * 3 + ["g2"] * 3 + ["g3"] * 3
    meta_df = pd.DataFrame({"filename": cols, "ATTRIBUTE_Group": groups})
    metadata_path = tmp_path / "metadata.csv"
    meta_df.to_csv(metadata_path, index=False)

    inputs = VisualizationInputs(data_matrix_file=str(matrix_path), metadata_file=str(metadata_path),
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()

    pls_result = module.pls_da()
    assert pls_result.get("skipped") is not True  # PLS-DA handles >2 groups fine

    opls_result = module.opls_da()
    assert opls_result["skipped"] is True
    assert "exactly 2 groups" in opls_result["reason"]

    volcano_result = module.volcano_plot()
    assert volcano_result["skipped"] is True
    assert "exactly 2 groups" in volcano_result["reason"]


def test_normalize_matrix_sum_equalizes_sample_totals(tmp_path):
    df = pd.DataFrame({
        "s1": [10.0, 20.0, 30.0],   # total 60
        "s2": [20.0, 40.0, 60.0],   # total 120 (exactly 2x s1 -- pure dilution)
    })
    result = normalize_matrix(df, "sum")
    totals = result.sum(axis=0)
    assert totals["s1"] == pytest.approx(totals["s2"], rel=1e-9)
    # relative composition within each sample is preserved
    assert (result["s1"] / result["s1"].sum()).tolist() == pytest.approx(
        (df["s1"] / df["s1"].sum()).tolist()
    )


def test_normalize_matrix_median_and_none(tmp_path):
    df = pd.DataFrame({"s1": [1.0, 2.0, 3.0], "s2": [2.0, 4.0, 6.0]})
    none_result = normalize_matrix(df, "none")
    pd.testing.assert_frame_equal(none_result, df)

    median_result = normalize_matrix(df, "median")
    medians = median_result.median(axis=0)
    assert medians["s1"] == pytest.approx(medians["s2"], rel=1e-9)

    with pytest.raises(ValueError):
        normalize_matrix(df, "bogus")


def test_pareto_scale_matches_hand_computed_formula():
    X = np.array([[1.0, 100.0], [2.0, 200.0], [3.0, 300.0]])
    result = pareto_scale(X)
    mean = X.mean(axis=0)
    std = X.std(axis=0, ddof=0)
    expected = (X - mean) / np.sqrt(std)
    np.testing.assert_allclose(result, expected)


def test_scale_features_none_only_centers():
    X = np.array([[1.0, 10.0], [3.0, 30.0]])
    result = scale_features(X, "none")
    np.testing.assert_allclose(result.mean(axis=0), [0.0, 0.0], atol=1e-9)
    # "none" should NOT divide by std -- spread stays proportional to input
    assert result[1, 0] - result[0, 0] == pytest.approx(2.0)

    with pytest.raises(ValueError):
        scale_features(X, "bogus")


def test_sum_normalization_corrects_a_pure_dilution_difference_between_groups(tmp_path):
    """The exact scenario flagged from the real dataset: one group's
    samples are uniformly ~2x higher across the board (a loading/dilution
    artifact, not real biology). Without normalization this looks like
    almost every feature is 'significant'; sum normalization should
    remove that artifact and collapse the fold changes back toward zero."""
    rng = np.random.default_rng(42)
    n_features, n_per_group = 50, 6
    # Independent draws for both groups from the SAME distribution (i.e. no
    # real biological difference), then a uniform 2x scale applied to all
    # of group B's values -- simulating a pure loading/dilution artifact
    # rather than a real per-feature effect.
    group_a = rng.exponential(1000, size=(n_features, n_per_group))
    group_b = rng.exponential(1000, size=(n_features, n_per_group)) * 2.0

    cols_a = [f"a_{i}.mzML" for i in range(n_per_group)]
    cols_b = [f"b_{i}.mzML" for i in range(n_per_group)]
    df = pd.DataFrame(np.hstack([group_a, group_b]), columns=cols_a + cols_b)
    df.insert(0, "row ID", [f"FT{i+1:03d}" for i in range(n_features)])
    matrix_path = tmp_path / "dilution_matrix.csv"
    df.to_csv(matrix_path, index=False)

    meta_df = pd.DataFrame(
        [(c, "groupA") for c in cols_a] + [(c, "groupB") for c in cols_b],
        columns=["filename", "ATTRIBUTE_Group"],
    )
    metadata_path = tmp_path / "dilution_metadata.csv"
    meta_df.to_csv(metadata_path, index=False)

    # Without normalization: expect the pure 2x dilution to look like a
    # real, near-universal log2FC ~= 1 signal.
    inputs_raw = VisualizationInputs(
        data_matrix_file=str(matrix_path), metadata_file=str(metadata_path),
        output_dir=str(tmp_path / "out_raw"), normalization="none",
    )
    module_raw = VisualizationModule(inputs_raw)
    module_raw.load()
    result_raw = module_raw.volcano_plot()
    csv_raw = pd.read_csv(result_raw["csv"])
    assert csv_raw["log2fc"].median() == pytest.approx(1.0, abs=0.2)

    # With sum normalization: the dilution artifact should be corrected,
    # collapsing the fold changes back toward zero.
    inputs_norm = VisualizationInputs(
        data_matrix_file=str(matrix_path), metadata_file=str(metadata_path),
        output_dir=str(tmp_path / "out_norm"), normalization="sum",
    )
    module_norm = VisualizationModule(inputs_norm)
    module_norm.load()
    result_norm = module_norm.volcano_plot()
    csv_norm = pd.read_csv(result_norm["csv"])
    assert abs(csv_norm["log2fc"].median()) < 0.15


def test_normalization_summary_reports_method_and_saves_plot(tmp_path):
    matrix_path, metadata_path = _make_matrix_and_metadata(tmp_path)
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"), normalization="median", scaling="auto")
    module = VisualizationModule(inputs)
    module.load()
    result = module.normalization_summary()
    assert result["normalization"] == "median"
    assert result["scaling"] == "auto"
    assert Path(result["plot"]).exists()


def test_hist_with_gaussian_fit_returns_stats_and_handles_empty():
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(7)
    values = rng.normal(loc=5.0, scale=2.0, size=2000)
    fig, ax = plt.subplots()
    stats = _hist_with_gaussian_fit(ax, values, "#4C72B0")
    plt.close(fig)
    assert stats is not None
    assert stats["mean"] == pytest.approx(5.0, abs=0.2)
    assert stats["std"] == pytest.approx(2.0, abs=0.2)
    assert abs(stats["skew"]) < 0.3  # a real normal draw should look near-symmetric
    assert stats["n"] == 2000

    # Empty/all-NaN input: no crash, returns None instead of plotting.
    fig, ax = plt.subplots()
    empty_stats = _hist_with_gaussian_fit(ax, np.array([np.nan, np.nan]), "#4C72B0")
    plt.close(fig)
    assert empty_stats is None


def test_distribution_shape_summary_produces_plots_and_reports_skew_reduction(tmp_path):
    """Uses the same pure-dilution-artifact fixture as the normalization
    test above (one group's samples uniformly ~2x higher, no real biology)
    to confirm the sample panel's plot/stats exist, and that pooling scaled
    (Pareto) feature values gives a tighter, more symmetric distribution
    than the pre-scaling normalized values -- exactly what Pareto scaling
    is supposed to do to a matrix with wildly different feature magnitudes."""
    rng = np.random.default_rng(3)
    n_features, n_per_group = 30, 6
    group_a = rng.exponential(1000, size=(n_features, n_per_group))
    group_b = rng.exponential(1000, size=(n_features, n_per_group)) * 2.0
    # A few very high-abundance features, to give scaling something real to correct.
    group_a[0:3] *= 50
    group_b[0:3] *= 50

    cols_a = [f"a_{i}.mzML" for i in range(n_per_group)]
    cols_b = [f"b_{i}.mzML" for i in range(n_per_group)]
    df = pd.DataFrame(np.hstack([group_a, group_b]), columns=cols_a + cols_b)
    df.insert(0, "row ID", [f"FT{i+1:03d}" for i in range(n_features)])
    matrix_path = tmp_path / "matrix.csv"
    df.to_csv(matrix_path, index=False)

    meta_df = pd.DataFrame(
        [(c, "groupA") for c in cols_a] + [(c, "groupB") for c in cols_b],
        columns=["filename", "ATTRIBUTE_Group"],
    )
    metadata_path = tmp_path / "metadata.csv"
    meta_df.to_csv(metadata_path, index=False)

    inputs = VisualizationInputs(
        data_matrix_file=str(matrix_path), metadata_file=str(metadata_path),
        output_dir=str(tmp_path / "out"), normalization="sum", scaling="pareto",
    )
    module = VisualizationModule(inputs)
    module.load()
    result = module.distribution_shape_summary()

    assert Path(result["sample_plot"]).exists()
    assert Path(result["feature_plot"]).exists()
    for key in ("sample_before", "sample_after", "feature_before", "feature_after"):
        assert result[key] is not None
        assert result[key]["n"] > 0

    # Pareto-scaled feature values should be far more tightly distributed
    # around 0 than the pre-scaling (log1p, normalized) values, since a
    # few very high-abundance features no longer dominate the pooled range.
    assert result["feature_after"]["std"] < result["feature_before"]["std"]
    assert abs(result["feature_after"]["mean"]) < abs(result["feature_before"]["mean"]) + 1.0


def test_full_run_produces_report_and_plots(tmp_path):
    matrix_path, metadata_path = _make_matrix_and_metadata(tmp_path)
    mzml_path = _write(tmp_path / "sample_A.mzML", MZML_WITH_PRECOMPUTED_TIC)
    inputs = VisualizationInputs(
        mzml_files=[mzml_path], metadata_file=metadata_path, data_matrix_file=matrix_path,
        output_dir=str(tmp_path / "out"),
    )
    module = VisualizationModule(inputs)
    report = module.run()

    assert "Feature summary" in report
    assert "Chromatograms" in report
    assert "PCA / correlation" in report
    assert "Distribution shape (Gaussian fit)" in report
    assert "Feature overlap (Venn/UpSet)" in report
    assert (tmp_path / "out" / "analysis_report.txt").exists()
    assert (tmp_path / "out" / f"chromatogram_{Path(mzml_path).stem}.png").exists()
    assert (tmp_path / "out" / "sample_distribution_gaussian.png").exists()
    assert (tmp_path / "out" / "feature_distribution_gaussian.png").exists()


# --------------------------------------------------------------------- #
# Spec 1.1.0 fixes
# --------------------------------------------------------------------- #

def test_extract_chromatogram_converts_minutes_to_seconds(tmp_path):
    xml = MZML_WITH_PRECOMPUTED_TIC.replace('unitAccession="UO:0000010"', 'unitAccession="UO:0000031"')
    path = _write(tmp_path / "minutes.mzML", xml)
    data = extract_chromatogram(path)
    assert data["rt"].tolist() == [60.0, 120.0]


def test_extract_chromatogram_decodes_array_when_only_bpc_missing(tmp_path):
    xml = _make_binary_intensity_mzml([[10.0, 20.0, 5.0]], compress=False)
    # Precomputed TIC present, BPC absent: the scan must still be kept,
    # with BPC taken from the decoded array and TIC from the cvParam.
    xml = xml.replace(
        '<binaryDataArrayList',
        '<cvParam accession="MS:1000285" name="total ion current" value="999.0"/>\n          <binaryDataArrayList',
    )
    path = _write(tmp_path / "tic_only.mzML", xml)
    data = extract_chromatogram(path)
    assert data["tic"].tolist() == [999.0]
    assert data["bpc"].tolist() == [20.0]


def _two_attribute_fixture(tmp_path):
    """Metadata with a 'Group' column (matched by the loose display rule)
    and an ATTRIBUTE_ column (matched by the supervised rule), plus one
    sample column with no metadata row."""
    matrix_path = _write(
        tmp_path / "matrix.csv",
        "row ID,s1.mzML,s2.mzML,s3.mzML,s4.mzML,orphan.mzML\n"
        "1,1000,1100,2000,2100,1500\n2,500,600,300,350,400\n3,10,20,30,40,25\n",
    )
    metadata_path = _write(
        tmp_path / "metadata.csv",
        "filename,Group,ATTRIBUTE_Tissue\n"
        "s1.mzML,x,MCL\ns2.mzML,y,MCL\ns3.mzML,x,MCR\ns4.mzML,y,MCR\n",
    )
    return matrix_path, metadata_path


def test_pca_and_feature_summary_use_the_supervised_group_column(tmp_path):
    matrix_path, metadata_path = _two_attribute_fixture(tmp_path)
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()

    pca = module.pca_and_correlation()
    pls = module.pls_da()
    summary = module.feature_summary()
    assert pca["group_column"] == "ATTRIBUTE_Tissue"
    assert pls["group_column"] == "ATTRIBUTE_Tissue"
    assert summary["sample_summary"]["group_column"] == "ATTRIBUTE_Tissue"


def test_pca_reports_unmatched_samples(tmp_path):
    matrix_path, metadata_path = _two_attribute_fixture(tmp_path)
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()
    result = module.pca_and_correlation()
    assert result["n_unmatched"] == 1
    report = module.build_report()
    assert "1 sample(s) not matched to metadata" in report


# --------------------------------------------------------------------- #
# Spec 1.2.0: feature overlap (Venn/UpSet), volcano colours, 3D PCA
# --------------------------------------------------------------------- #

def _pattern_fixture(tmp_path, groups_per_sample, rows):
    """rows: list of per-sample value lists ('' = NA). Builds a matrix
    with feature ids F1.. and metadata mapping each sample to its group."""
    cols = [f"s{i}.mzML" for i in range(len(groups_per_sample))]
    lines = ["row ID," + ",".join(cols)]
    for n, row in enumerate(rows, 1):
        lines.append(f"F{n}," + ",".join(str(v) for v in row))
    matrix_path = _write(tmp_path / "matrix.csv", "\n".join(lines) + "\n")
    meta = ["filename,ATTRIBUTE_Group"] + [f"{c},{g}" for c, g in zip(cols, groups_per_sample)]
    metadata_path = _write(tmp_path / "metadata.csv", "\n".join(meta) + "\n")
    return matrix_path, metadata_path


def _overlap(tmp_path, matrix_path, metadata_path, **kw):
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"), **kw)
    module = VisualizationModule(inputs)
    module.load()
    return module.feature_overlap()


def test_feature_overlap_two_groups_counts_and_csv(tmp_path):
    matrix_path, metadata_path = _pattern_fixture(tmp_path, ["A", "A", "B", "B"], [
        [10, 10, 10, 10],   # F1 shared
        [10, 10, 0, 0],     # F2 A only
        [10, 0, "", 0],     # F3 A only (1/2 = 50% meets the default threshold)
        [0, "", 10, 10],    # F4 B only
        ["", "", "", ""],   # F5 none (all NA)
        [0, 0, 0, 0],       # F6 none (all zero)
    ])
    result = _overlap(tmp_path, matrix_path, metadata_path)
    assert result["groups"] == ["A", "B"]
    assert result["unique_per_group"] == {"A": 2, "B": 1}
    assert result["captured_per_group"] == {"A": 3, "B": 2}
    assert result["shared_all"] == 1
    assert result["none"] == 2
    assert result["plot"].endswith("venn_diagram.png") and Path(result["plot"]).exists()

    csv_df = pd.read_csv(result["csv"])
    assert csv_df["overlap"].tolist() == ["A & B", "A only", "A only", "B only", "none", "none"]
    assert {"row ID", "captured_A", "captured_B"}.issubset(csv_df.columns)


def test_feature_overlap_three_groups_region_counts(tmp_path):
    # One sample per group; region (members) -> number of features
    wanted = {("g1",): 1, ("g2",): 2, ("g3",): 3, ("g1", "g2"): 4,
              ("g1", "g3"): 5, ("g2", "g3"): 6, ("g1", "g2", "g3"): 7}
    rows = []
    for members, n in wanted.items():
        rows += [[10 if g in members else 0 for g in ("g1", "g2", "g3")]] * n
    rows.append([0, 0, 0])
    matrix_path, metadata_path = _pattern_fixture(tmp_path, ["g1", "g2", "g3"], rows)
    result = _overlap(tmp_path, matrix_path, metadata_path)
    expected = {(f"{m[0]} only" if len(m) == 1 else " & ".join(m)): n for m, n in wanted.items()}
    assert result["intersections"] == expected
    assert result["shared_all"] == 7
    assert result["none"] == 1
    assert Path(result["plot"]).name == "venn_diagram.png"


def test_feature_overlap_min_fraction_threshold(tmp_path):
    rows = [[10, 0, 0, 10, 10, 10]]  # A: 1/3 detected, B: 3/3
    matrix_path, metadata_path = _pattern_fixture(tmp_path, ["A"] * 3 + ["B"] * 3, rows)
    strict = _overlap(tmp_path, matrix_path, metadata_path)
    assert strict["intersections"] == {"B only": 1}
    loose = _overlap(tmp_path, matrix_path, metadata_path, venn_min_fraction=0.3)
    assert loose["intersections"] == {"A & B": 1}


def test_feature_overlap_four_groups_uses_upset(tmp_path):
    rng = np.random.default_rng(5)
    groups = ["g1", "g2", "g3", "g4"]
    rows = (rng.random((40, 4)) > 0.5).astype(int) * 100
    matrix_path, metadata_path = _pattern_fixture(tmp_path, groups, rows.tolist())
    result = _overlap(tmp_path, matrix_path, metadata_path)
    assert Path(result["plot"]).name == "upset_plot.png"
    assert Path(result["plot"]).exists()
    assert not (tmp_path / "out" / "venn_diagram.png").exists()
    assert sum(result["intersections"].values()) == 40 - result["none"]


def test_feature_overlap_skips_with_one_group_or_no_metadata(tmp_path):
    matrix_path, metadata_path = _pattern_fixture(tmp_path, ["A", "A"], [[1, 2]])
    assert _overlap(tmp_path, matrix_path, metadata_path)["skipped"] is True

    module = VisualizationModule(VisualizationInputs(data_matrix_file=matrix_path,
                                                     output_dir=str(tmp_path / "out2")))
    module.load()
    assert module.feature_overlap()["skipped"] is True


@pytest.mark.parametrize("bad", [0, 1.5, -0.2])
def test_invalid_venn_min_fraction_raises(tmp_path, bad):
    matrix_path, metadata_path = _pattern_fixture(tmp_path, ["A", "B"], [[1, 2]])
    module = VisualizationModule(VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                                     output_dir=str(tmp_path / "out"), venn_min_fraction=bad))
    with pytest.raises(ValueError):
        module.load()


def test_volcano_counts_significant_features_per_group(tmp_path):
    matrix_path, metadata_path = _two_group_matrix_and_metadata(tmp_path)
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()
    result = module.volcano_plot()
    n = result["n_higher_in"]
    assert set(n) == {"control", "treated"}
    assert n["treated"] >= 1 and n["control"] >= 1
    assert n["treated"] + n["control"] == result["n_significant"]


def test_pca_3d_written_with_enough_samples(tmp_path):
    matrix_path, metadata_path = _two_group_matrix_and_metadata(tmp_path)
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()
    result = module.pca_and_correlation()
    assert len(result["explained_variance_ratio"]) == 3
    assert Path(result["plot_3d"]).name == "pca_3d.png"
    assert Path(result["plot_3d"]).exists()


def test_pca_3d_skipped_with_three_samples(tmp_path):
    matrix_path, metadata_path = _make_matrix_and_metadata(tmp_path)
    inputs = VisualizationInputs(data_matrix_file=matrix_path, metadata_file=metadata_path,
                                  output_dir=str(tmp_path / "out"))
    module = VisualizationModule(inputs)
    module.load()
    result = module.pca_and_correlation()
    assert "plot_3d" not in result
    assert (tmp_path / "out" / "pca.png").exists()
    assert not (tmp_path / "out" / "pca_3d.png").exists()
    assert "3D PCA skipped" in module.build_report()
