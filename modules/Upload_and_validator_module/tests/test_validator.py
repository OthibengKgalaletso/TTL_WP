"""
Pytest suite for metabolomics_validator.py.

Run with:
    pytest -v

These tests exercise the programmatic API (InputSet / MetabolomicsValidator)
directly rather than shelling out to the CLI, and write small fixture files
on the fly with pytest's tmp_path fixture so they don't depend on the
sample_data/ folder (though sample_data/ is useful for manual/CLI testing).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from metabolomics_validator import InputSet, MetabolomicsValidator  # noqa: E402


VALID_MZML = """<?xml version="1.0" encoding="utf-8"?>
<indexedmzML>
  <mzML>
    <run>
      <spectrumList>
        <spectrum id="scan=1" index="0" defaultArrayLength="2">
          <cvParam accession="MS:1000511" name="ms level" value="1"/>
          <cvParam accession="MS:1000016" name="scan start time" value="0.1" unitAccession="UO:0000010"/>
          <binaryDataArrayList count="2">
            <binaryDataArray><cvParam accession="MS:1000514" name="m/z array"/></binaryDataArray>
            <binaryDataArray><cvParam accession="MS:1000515" name="intensity array"/></binaryDataArray>
          </binaryDataArrayList>
        </spectrum>
      </spectrumList>
    </run>
  </mzML>
</indexedmzML>
"""


def _write(path: Path, text: str) -> str:
    path.write_text(text)
    return str(path)


def _valid_required_inputs(tmp_path: Path) -> InputSet:
    mzml_path = _write(tmp_path / "sample_A.mzML", VALID_MZML)
    metadata_path = _write(
        tmp_path / "sample_metadata.csv",
        "filename,ATTRIBUTE/class\nsample_A.mzML,control\n",
    )
    matrix_path = _write(
        tmp_path / "data_matrix.csv",
        "feature_id,sample_A.mzML\nFT0001,1000\nFT0002,300\n",
    )
    return InputSet(mzml_files=[mzml_path], metadata_file=metadata_path, data_matrix_file=matrix_path)


def test_valid_required_only_passes(tmp_path):
    report = MetabolomicsValidator(_valid_required_inputs(tmp_path)).run()
    assert report.passed
    assert report.errors == []


def test_missing_required_files_fails(tmp_path):
    report = MetabolomicsValidator(InputSet()).run()
    assert not report.passed
    error_files = [i.file for i in report.errors]
    assert "mzML files" in error_files
    assert "Sample metadata" in error_files
    assert "Data matrix" in error_files


def test_optional_files_absent_are_ignored(tmp_path):
    """Omitting feature table / MGF / annotation list should never fail validation."""
    inputs = _valid_required_inputs(tmp_path)
    report = MetabolomicsValidator(inputs).run()
    assert report.passed


def test_optional_file_present_but_malformed_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.mgf_file = _write(
        tmp_path / "broken.mgf",
        "BEGIN IONS\nTITLE=x\nPEPMASS=100\nEND IONS\n",  # missing CHARGE, RTINSECONDS
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any("mgf" in i.file.lower() for i in report.errors)


def test_malformed_mzml_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    bad_mzml = _write(tmp_path / "corrupt.mzML", "<indexedmzML><mzML><run>")  # truncated
    inputs.mzml_files = [bad_mzml]
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed


def test_metadata_missing_grouping_column_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.metadata_file = _write(
        tmp_path / "no_class.csv",
        "filename,batch\nsample_A.mzML,1\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed


def test_data_matrix_non_numeric_cell_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.data_matrix_file = _write(
        tmp_path / "bad_matrix.csv",
        "feature_id,sample_A.mzML\nFT0001,not_a_number\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed


def test_mismatched_sample_names_fail_compatibility(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.data_matrix_file = _write(
        tmp_path / "mismatched_matrix.csv",
        "feature_id,sample_X.mzML\nFT0001,1000\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any(i.step.value == "Compatibility check" for i in report.errors)


def test_semicolon_delimited_metadata_and_matrix_pass(tmp_path):
    """Some regional Excel locales export CSV with ';' instead of ','."""
    mzml_path = _write(tmp_path / "sample_A.mzML", VALID_MZML)
    metadata_path = _write(
        tmp_path / "sample_metadata.csv",
        "filename;ATTRIBUTE/class\nsample_A.mzML;control\n",
    )
    matrix_path = _write(
        tmp_path / "data_matrix.csv",
        "feature_id;sample_A.mzML\nFT0001;1000\n",
    )
    inputs = InputSet(mzml_files=[mzml_path], metadata_file=metadata_path, data_matrix_file=matrix_path)
    report = MetabolomicsValidator(inputs).run()
    assert report.passed


def test_mgf_accepts_feature_id_bare_charge_and_scientific_notation(tmp_path):
    """Real MZmine/GNPS exports commonly use FEATURE_ID/SCANS instead of
    TITLE, an unsigned CHARGE, and scientific-notation intensities."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.mgf_file = _write(
        tmp_path / "gnps_style.mgf",
        "BEGIN IONS\n"
        "FEATURE_ID=1\n"
        "PEPMASS=255.233\n"
        "CHARGE=1\n"
        "SCANS=1\n"
        "RTINSECONDS=192.6\n"
        "100.5 2.34e+05 1\n"
        "150.2 5.0E+02\n"
        "END IONS\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert report.passed


def test_mgf_still_rejects_missing_tags_and_non_numeric_peaks(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.mgf_file = _write(
        tmp_path / "truly_broken.mgf",
        "BEGIN IONS\nPEPMASS=255.233\nnotanumber notanumber\nEND IONS\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed


def test_data_matrix_blank_cells_are_warning_not_error(tmp_path):
    """Blank/NA cells (non-detected features) are common in real untargeted
    metabolomics data matrices and should warn, not block validation."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.data_matrix_file = _write(
        tmp_path / "matrix_with_blanks.csv",
        "feature_id,sample_A.mzML\nFT0001,1000\nFT0002,\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert report.passed
    assert any("missing/NA" in w.message for w in report.warnings)


def test_mgf_mixed_case_tag_line_not_mistaken_for_a_peak(tmp_path):
    """Real MZmine exports include a 'Num peaks=N' line (mixed case, so
    the old ALL-CAPS-only tag check misparsed it as a malformed peak)."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.mgf_file = _write(
        tmp_path / "mzmine_style.mgf",
        "BEGIN IONS\n"
        "FEATURE_ID=2329\n"
        "RTINSECONDS=364.25\n"
        "PEPMASS=137.0952\n"
        "SCANS=2329\n"
        "Num peaks=1\n"
        "107.048422 1606\n"
        "END IONS\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert report.passed


def test_mgf_charge_optional_when_undetermined(tmp_path):
    """MZmine omits CHARGE entirely for single-scan, non-merged features —
    that alone shouldn't fail validation."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.mgf_file = _write(
        tmp_path / "no_charge.mgf",
        "BEGIN IONS\n"
        "FEATURE_ID=3841\n"
        "RTINSECONDS=472.66\n"
        "PEPMASS=163.07458\n"
        "SPECTYPE=SINGLE_SCAN\n"
        "SCANS=3841\n"
        "Num peaks=1\n"
        "103.053231 1718\n"
        "END IONS\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert report.passed


def test_annotation_list_accepts_mzmine_id_column(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.annotation_list_file = _write(
        tmp_path / "mzmine_annotations.csv",
        "id,compound_name,adduct,score\nFT0001,Glutamine,[M+H]+,0.91\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert report.passed


def test_data_matrix_accepts_mzmine_row_id_column(tmp_path):
    """MZmine's own feature-list/quant export convention names the
    identifier column 'row ID' (with a space), not feature_id/id/peak_id."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.data_matrix_file = _write(
        tmp_path / "row_id_matrix.csv",
        "row ID,sample_A.mzML\n1,1000\n2,300\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert report.passed


def test_data_matrix_full_feature_table_shape_excludes_mz_and_rt_from_sample_match(tmp_path):
    """MZmine's full feature-table-style data matrix export is shaped
    'row ID, row m/z, row retention time, <sample columns>' — the m/z and
    retention time columns describe the feature, not a sample, and must
    not be treated as unmatched sample columns during compatibility."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.data_matrix_file = _write(
        tmp_path / "full_shape_matrix.csv",
        "row ID,row m/z,row retention time,sample_A.mzML\nFT0001,255.233,192.6,1000\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert report.passed
    assert not any(
        "row m/z" in i.message or "row retention time" in i.message for i in report.issues
    )


def test_gnps_style_attribute_column_accepted(tmp_path):
    """GNPS metadata sheets use ATTRIBUTE_<anything> (e.g. ATTRIBUTE_Tissue)
    as the grouping column, not a fixed name — and are often tab-delimited
    despite the .csv extension."""
    mzml_path = _write(tmp_path / "sample_A.mzML", VALID_MZML)
    metadata_path = _write(
        tmp_path / "metadata.csv",
        "filename\tATTRIBUTE_Tissue\nsample_A.mzML\tMCL\n",
    )
    matrix_path = _write(
        tmp_path / "data_matrix.csv",
        "feature_id,sample_A.mzML\nFT0001,1000\n",
    )
    inputs = InputSet(mzml_files=[mzml_path], metadata_file=metadata_path, data_matrix_file=matrix_path)
    report = MetabolomicsValidator(inputs).run()
    assert report.passed