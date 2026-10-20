# Module: Upload & Validation

Part of the [TTL_pipeline](../../README.md) project. A standalone Python
script that validates a metabolomics dataset — mzML spectra, sample
metadata, a data matrix and the full MZmine feature table (required), plus
an optional MZmine MGF file and metabolite annotation list — before it
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
  --feature-table sample_data/valid/mzmine_full_feature_table.csv \
  --mgf sample_data/valid/mzmine_spectra.mgf \
  --annotation-list sample_data/valid/annotation_list.csv
```

This should print `VALIDATION PASSED`. To see the failure path, swap in one
of the files from `sample_data/invalid/`, e.g.:

```bash
python metabolomics_validator.py \
  --mzml sample_data/valid/sample_A.mzML sample_data/valid/sample_B.mzML \
  --metadata sample_data/invalid/sample_metadata_missing_class.csv \
  --data-matrix sample_data/valid/data_matrix.csv \
  --feature-table sample_data/valid/mzmine_full_feature_table.csv
```

mzML files, sample metadata, the data matrix and the full feature table are
required; the MGF file and annotation list are optional and are only
validated when supplied.

### MZmine-only feature table and MGF

The feature table and MGF are accepted **only** as MZmine 3/4 output:

- **Feature table** — must be the *full* feature list, exported from MZmine
  with *Feature list → Export → Export to CSV file (modular)*, keeping all
  row types selected (`.csv` only). It must contain `id`, `mz`, `rt`,
  `rt_range:min/max`, `mz_range:min/max`, `height`, `area`, and
  `datafile:<file>:area` + `datafile:<file>:height` for every sample.
  MZmine's reduced *GNPS/IIMN quant table* (`row ID`, `row m/z`,
  `<file> Peak area`) is rejected with a re-export message.
- **MGF** — every block needs `FEATURE_ID` (integer), `MSLEVEL`,
  `RTINSECONDS` and `PEPMASS`, and every `FEATURE_ID` must exist in the
  feature table's `id` column (i.e. both exported from the same feature list).

`./run_check.sh` expects the full export at
`sample_data/valid_files/Practice_data_MeOH_POS_full_feature_table.csv`.
See `SPEC.md` (VAL-EX-01, VAL-FQ-*, VAL-MG-*, VAL-CP-08/09) for the exact rules.

## Running the tests

```bash
pytest -v
```

## Format rules enforced

See the module docstring at the top of `metabolomics_validator.py` for the
exact required columns/tags/CV accessions checked for each file type.
