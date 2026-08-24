# Module: Upload & Validation

Part of the [TTL_pipeline](../../README.md) project. A standalone Python
script that validates a metabolomics dataset — mzML spectra, sample
metadata, and a data matrix (required), plus an optional feature
quantification table, MGF file, and metabolite annotation list — before it
enters downstream processing. It mirrors a four-step pipeline:

```
Upload files -> Validate file types -> Validate file integrity & format
             -> Compatibility check -> Validation passed? (YES/NO)
```

On failure it prints a structured error report grouped by step and file,
and exits with a non-zero status code.

## Module layout

```
metabolomics_validator.py   the validator (CLI + programmatic API)
run_check.sh                convenience wrapper around sample_data/valid/
tests/test_validator.py     pytest suite
sample_data/valid/          a small dataset that passes validation
sample_data/invalid/        a few files that intentionally fail, for testing
```

## Setup

Dependencies for the whole TTL_pipeline project (this module plus any
future ones) live in a single `requirements.txt` at the project root — see
the top-level [README](../../README.md) for environment setup. From the
project root, once your environment is active:

```bash
pip install -r requirements.txt
```

## Running the validator

Run these from inside `modules/upload_validation/`:

```bash
python metabolomics_validator.py \
  --mzml sample_data/valid/sample_A.mzML sample_data/valid/sample_B.mzML \
  --metadata sample_data/valid/sample_metadata.csv \
  --data-matrix sample_data/valid/data_matrix.csv \
  --feature-table sample_data/valid/feature_table.csv \
  --mgf sample_data/valid/spectra.mgf \
  --annotation-list sample_data/valid/annotation_list.csv
```

This should print `VALIDATION PASSED`. To see the failure path, swap in one
of the files from `sample_data/invalid/`, e.g.:

```bash
python metabolomics_validator.py \
  --mzml sample_data/valid/sample_A.mzML sample_data/valid/sample_B.mzML \
  --metadata sample_data/invalid/sample_metadata_missing_class.csv \
  --data-matrix sample_data/valid/data_matrix.csv
```

Only mzML files, sample metadata, and the data matrix are required; the
feature table, MGF file, and annotation list are optional and are only
validated when supplied.

## Running the tests

```bash
pytest -v
```

## Format rules enforced

See the module docstring at the top of `metabolomics_validator.py` for the
exact required columns/tags/CV accessions checked for each file type.
