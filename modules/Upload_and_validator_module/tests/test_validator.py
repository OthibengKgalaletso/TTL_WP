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


# Minimal MZmine 3/4 "Export to CSV file (modular)" full feature table.
VALID_FULL_FEATURE_TABLE = (
    "id,rt,rt_range:min,rt_range:max,mz,mz_range:min,mz_range:max,height,area,"
    "datafile:sample_A.mzML:area,datafile:sample_A.mzML:height,datafile:sample_A.mzML:rt,"
    "datafile:sample_A.mzML:rt_range:min\n"  # nested field, as in real MZmine exports
    "1,3.21,3.18,3.25,255.233,255.231,255.235,1200,1000,1000,1200,3.21,3.18\n"
    "2,6.07,6.01,6.12,137.095,137.094,137.096,350,300,,,,\n"
)


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
        "feature_id,sample_A.mzML\n1,1000\n2,300\n",
    )
    feature_table_path = _write(tmp_path / "full_feature_table.csv", VALID_FULL_FEATURE_TABLE)
    return InputSet(
        mzml_files=[mzml_path], metadata_file=metadata_path,
        data_matrix_file=matrix_path, feature_table_file=feature_table_path,
    )


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
    assert "Feature table" in error_files


def test_missing_feature_table_fails(tmp_path):
    """The full MZmine feature table is a required input."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.feature_table_file = None
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any(i.file == "Feature table" and i.step.value == "Upload files" for i in report.errors)


def test_optional_files_absent_are_ignored(tmp_path):
    """Omitting the MGF / annotation list should never fail validation."""
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
    feature_table_path = _write(tmp_path / "full_feature_table.csv", VALID_FULL_FEATURE_TABLE)
    inputs = InputSet(
        mzml_files=[mzml_path], metadata_file=metadata_path,
        data_matrix_file=matrix_path, feature_table_file=feature_table_path,
    )
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
        "MSLEVEL=2\n"
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
        "FEATURE_ID=2\n"
        "MSLEVEL=2\n"
        "RTINSECONDS=364.25\n"
        "PEPMASS=137.0952\n"
        "SCANS=2\n"
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
        "FEATURE_ID=1\n"
        "MSLEVEL=2\n"
        "RTINSECONDS=472.66\n"
        "PEPMASS=163.07458\n"
        "SPECTYPE=SINGLE_SCAN\n"
        "SCANS=1\n"
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
        "id,compound_name,adduct,score\n1,Glutamine,[M+H]+,0.91\n",
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
        "row ID,row m/z,row retention time,sample_A.mzML\n1,255.233,192.6,1000\n",
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
    feature_table_path = _write(tmp_path / "full_feature_table.csv", VALID_FULL_FEATURE_TABLE)
    inputs = InputSet(
        mzml_files=[mzml_path], metadata_file=metadata_path,
        data_matrix_file=matrix_path, feature_table_file=feature_table_path,
    )
    report = MetabolomicsValidator(inputs).run()
    assert report.passed


# --------------------------------------------------------------------------- #
# MZmine-only full feature table (SPEC 2.0.0: UP-06, FQ-01..09, EX-01)
# --------------------------------------------------------------------------- #

def test_gnps_quant_table_rejected_with_reexport_message(tmp_path):
    """MZmine's reduced GNPS/IIMN quant export is not the full feature table."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.feature_table_file = _write(
        tmp_path / "iimn_gnps_quant.csv",
        "row ID,row m/z,row retention time,sample_A.mzML Peak area\n1,255.233,3.21,1000\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    ft_errors = [i for i in report.errors if i.file == inputs.feature_table_file]
    assert len(ft_errors) == 1
    assert "Export to CSV" in ft_errors[0].message


def test_full_table_missing_row_level_fields_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.feature_table_file = _write(
        tmp_path / "partial.csv",
        "id,rt,mz,datafile:sample_A.mzML:area,datafile:sample_A.mzML:height\n1,3.21,255.233,1000,1200\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any("rt_range:min" in i.message and "area" in i.message for i in report.errors)


def test_full_table_sample_missing_height_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.feature_table_file = _write(
        tmp_path / "no_height.csv",
        "id,rt,rt_range:min,rt_range:max,mz,mz_range:min,mz_range:max,height,area,"
        "datafile:sample_A.mzML:area\n"
        "1,3.21,3.18,3.25,255.233,255.231,255.235,1200,1000,1000\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any("sample_A.mzML" in i.message and "height" in i.message for i in report.errors)


def test_feature_table_duplicate_ids_fail(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    lines = VALID_FULL_FEATURE_TABLE.splitlines()
    inputs.feature_table_file = _write(
        tmp_path / "dupes.csv", "\n".join([lines[0], lines[1], lines[1]]) + "\n",
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any("Duplicate" in i.message for i in report.errors)


def test_feature_table_non_integer_id_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.feature_table_file = _write(
        tmp_path / "bad_ids.csv",
        VALID_FULL_FEATURE_TABLE.replace("\n1,", "\nFT0001,", 1),
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any("integer" in i.message for i in report.errors)


def test_feature_table_xlsx_extension_rejected(tmp_path):
    """MZmine's CSV export only writes .csv."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.feature_table_file = _write(tmp_path / "full_feature_table.xlsx", VALID_FULL_FEATURE_TABLE)
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any(i.step.value == "Validate file types" for i in report.errors)


# --------------------------------------------------------------------------- #
# MZmine-only MGF (SPEC 2.0.0: MG-04, MG-05, MG-11)
# --------------------------------------------------------------------------- #

def _mzmine_block(**overrides) -> str:
    tags = {"FEATURE_ID": "1", "MSLEVEL": "2", "RTINSECONDS": "192.6", "PEPMASS": "255.233", "SCANS": "1"}
    tags.update(overrides)
    body = "".join(f"{k}={v}\n" for k, v in tags.items() if v is not None)
    return f"BEGIN IONS\n{body}100.5 2000\nEND IONS\n"


def test_mgf_without_feature_id_fails(tmp_path):
    """TITLE/SCANS alone don't identify an MZmine feature."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.mgf_file = _write(tmp_path / "no_feature_id.mgf", _mzmine_block(FEATURE_ID=None, TITLE="x"))
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any("FEATURE_ID" in i.message for i in report.errors)


def test_mgf_without_mslevel_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.mgf_file = _write(tmp_path / "no_mslevel.mgf", _mzmine_block(MSLEVEL=None))
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any("MSLEVEL" in i.message for i in report.errors)


def test_mgf_rt_tag_instead_of_rtinseconds_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.mgf_file = _write(tmp_path / "rt_tag.mgf", _mzmine_block(RTINSECONDS=None, RT="3.21"))
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any("RTINSECONDS" in i.message for i in report.errors)


# --------------------------------------------------------------------------- #
# Compatibility with the feature table (SPEC 2.0.0: CP-08, CP-09)
# --------------------------------------------------------------------------- #

def test_metadata_sample_missing_from_feature_table_fails(tmp_path):
    inputs = _valid_required_inputs(tmp_path)
    inputs.feature_table_file = _write(
        tmp_path / "other_sample.csv",
        VALID_FULL_FEATURE_TABLE.replace("sample_A.mzML", "sample_X.mzML"),
    )
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any(
        i.step.value == "Compatibility check" and i.file == "Feature table" for i in report.errors
    )


def test_mgf_feature_id_not_in_feature_table_fails(tmp_path):
    """MGF and feature table must come from the same MZmine feature list."""
    inputs = _valid_required_inputs(tmp_path)
    inputs.mgf_file = _write(tmp_path / "other_list.mgf", _mzmine_block(FEATURE_ID="999", SCANS="999"))
    report = MetabolomicsValidator(inputs).run()
    assert not report.passed
    assert any(i.step.value == "Compatibility check" and "999" in i.message for i in report.errors)
