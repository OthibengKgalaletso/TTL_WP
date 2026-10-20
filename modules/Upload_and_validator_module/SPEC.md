# SPEC — Upload & Validation module

| | |
|---|---|
| **Module** | `modules/Upload_and_validator_module/` |
| **Implementation** | `metabolomics_validator.py` |
| **Tests** | `tests/test_validator.py` |
| **Spec version** | 2.0.0 (MZmine-only feature table and MGF; full feature table required, 2026-10-20) |
| **Status** | Done. This module gates all downstream modules. |

This document is the **source of truth** for what the module does. Code, tests and README must agree with it. Any change in behaviour follows [§9 Change control](#9-change-control).

Requirement IDs (`VAL-xx-nn`) are stable. Never reuse or renumber one. When a requirement is removed, mark it **Retired** and leave it in place.

---

## 1. Purpose and scope

This module checks a user-supplied metabolomics dataset **before** any downstream module touches it. It answers one question: *"Is this dataset well-formed and internally consistent enough to process?"* It reports every problem it can find in one pass, so the user can fix them all at once.

**In scope:** file presence, file extensions, structural/format checks on each file, and cross-file consistency checks.

**Out of scope:** scientific quality (e.g. QC drift, poor chromatography), statistical checks, changing or fixing input files, and converting between formats. The module is **read-only** with respect to its inputs.

## 2. Pipeline position

```
Upload files -> Validate file types -> Validate file integrity & format
             -> Compatibility check -> Validation passed?
                                         YES -> hand off to downstream module
                                         NO  -> error report -> back to Upload
```

## 3. Inputs

| Input | Required | `InputSet` field | CLI flag | Accepted extensions (case-insensitive) |
|---|---|---|---|---|
| mzML files (1+) | **Yes** | `mzml_files: list[str]` | `--mzml FILE [FILE ...]` | `.mzml` |
| Sample metadata | **Yes** | `metadata_file` | `--metadata` | `.csv` `.tsv` `.txt` |
| Data matrix | **Yes** | `data_matrix_file` | `--data-matrix` | `.csv` `.tsv` `.txt` `.xlsx` `.xls` |
| Full feature table (MZmine 3/4 "Export to CSV") | **Yes** | `feature_table_file` | `--feature-table` | `.csv` |
| MGF (MS/MS spectra, MZmine export) | No | `mgf_file` | `--mgf` | `.mgf` |
| Annotation list | No | `annotation_list_file` | `--annotation-list` | `.csv` `.tsv` `.txt` `.xlsx` `.xls` |

Optional inputs are validated **only when supplied**. If an optional input is absent, that is never an issue.

`InputSet.feature_table_file` keeps its `Optional[str]` type so the dataclass signature is unchanged. Leaving it unset is reported as an Upload ERROR (VAL-UP-06).

### 3.1 Table parsing rules (all tabular inputs)

**VAL-TB-01**
- `.tsv` is read as tab-separated.
- `.csv` and `.txt` are read by sniffing the delimiter. The first line is checked for `,`, tab and `;`, and the most frequent one wins. If none of them appears, comma is used.
- A UTF-8 BOM is stripped.

**VAL-TB-02**
- `.xlsx`/`.xls` are read with `pandas.read_excel`.

**VAL-TB-03**
- All cells are read as strings.
- Only these tokens count as missing: empty string, `NA`, `N/A`, `nan`. Other spellings (e.g. `NaN`, `NULL`, `-`) are **not** treated as missing (see §8).

**VAL-TB-04**
- Column headers are matched **leniently**: names are lowercased and stripped of everything except `a–z0–9`. For example, `row ID`, `Row_ID` and `rowid` are equivalent.
- When several aliases are listed for a column, the **first column in the file** that matches any alias is used.

**VAL-TB-05**
- A table that cannot be parsed produces an ERROR `Could not parse file as a table`.

### 3.2 Filename matching rule (compatibility step)

**VAL-FN-01**
- Two names *match* if they share any variant from the set {lowercased name, lowercased name with its last extension removed}.
- So `Sample_A.mzML`, `sample_a.mzml` and `Sample_A` all match each other.
- Whitespace at the ends is trimmed. Directory paths are ignored for mzML inputs (the basename is used).

### 3.3 Accepted exporters

**VAL-EX-01**
- The feature table and the MGF are accepted **only** as MZmine 3/4 output.
- The feature table must be the **full** feature list exported with MZmine's *Feature list → Export → Export to CSV file (modular)*. It holds row-level columns (`id`, `mz`, `rt`, `rt_range:min`, …) plus per-sample `datafile:<file>:<field>` columns.
- MZmine's reduced exports, the GNPS/IIMN quant table (`row ID`, `row m/z`, `<file> Peak area`) and the legacy MZmine 2 CSV, are rejected with a targeted re-export message (VAL-FQ-08).
- The MGF must carry MZmine's per-block tags (`FEATURE_ID`, `MSLEVEL`, `RTINSECONDS`), and its feature IDs must come from the same feature list as the table (VAL-CP-09).
- XCMS and other exporters are **not** accepted for these two inputs. The data matrix, metadata and annotation list keep their exporter-agnostic aliases.

## 4. Outputs and interface contract

### 4.1 Python API (stable; the downstream GUI/Nextflow will depend on it)

```python
from metabolomics_validator import InputSet, MetabolomicsValidator
report = MetabolomicsValidator(InputSet(...)).run()   # -> ValidationReport
report.passed      # bool
report.issues      # list[Issue] in discovery order
report.errors / report.warnings
report.to_dict()   # JSON-serialisable dict, see 4.3
report.to_text()   # human-readable report, see 4.4
```

- `Issue` has these fields: `step: Step`, `severity: Severity`, `file: str`, `message: str`.
- `Severity` is `ERROR` or `WARNING`.
- `Step` values are the exact strings `"Upload files"`, `"Validate file types"`, `"Validate file integrity & format"` and `"Compatibility check"`.
- The public step methods `check_upload`, `check_file_types`, `check_integrity_and_format` and `check_compatibility` each return `bool`.

### 4.2 Pass/fail semantics

**VAL-PF-01**
- `passed` is `True` only when **every** step returns `True`.

**VAL-PF-02**
- **WARNINGs never cause failure.** ERRORs always do.

**VAL-PF-03**
- All four steps **always run**, even after an earlier step fails, so the user gets the full list of problems in one pass.
- A later step skips any file that already failed an earlier step, so it doesn't report a duplicate cascade of errors (see VAL-CP-00).

### 4.3 JSON report (`--report-format json`)

**VAL-RP-01**
- Top-level keys: `timestamp` (ISO-8601, seconds), `validation_passed`, `error_count`, `warning_count` and `issues`.
- `issues` is a list of `{step, severity, file, message}`.
- Adding keys is a MINOR change. Renaming or removing keys is MAJOR (§9).

### 4.4 Text report (default)

**VAL-RP-02**
- A header block, then `Result: VALIDATION PASSED|VALIDATION FAILED`, then the error and warning counts.
- Issues are grouped under `--- <step> ---` in pipeline order. Each issue is printed as `  [SEVERITY] <file>: <message>`.
- The footer reads either `=> Validation passed (YES). Handing off…` or `=> Validation failed (NO)… re-upload.`
- If there are no issues, it prints `No issues found across any pipeline step. Ready to proceed.`

### 4.5 CLI

**VAL-CL-01**
- Exit code `0` means passed and `1` means failed.
- Exit code `2` is an argparse usage error, raised when **none** of `--mzml`, `--metadata`, `--data-matrix` and `--feature-table` is given.
- If only some of the required flags are given, the run continues and the missing inputs are reported as Upload ERRORs (exit 1).

**VAL-CL-02**
- The report is always printed to stdout.
- `--report-out FILE` also writes it to FILE and prints `Report written to: FILE`.

## 5. Validation rules

The **Sev** column shows the issue severity. `E` = ERROR, `W` = WARNING.

### 5.1 Step 1 — Upload files

| ID | Rule | Sev |
|---|---|---|
| VAL-UP-01 | At least one mzML file is supplied. | E |
| VAL-UP-02 | A sample metadata file is supplied. | E |
| VAL-UP-03 | A data matrix file is supplied. | E |
| VAL-UP-04 | Every supplied path (required or optional) exists as a file. | E |
| VAL-UP-05 | Every supplied file is non-empty (> 0 bytes). | E |
| VAL-UP-06 | A full feature table file is supplied. | E |

### 5.2 Step 2 — Validate file types

| ID | Rule | Sev |
|---|---|---|
| VAL-FT-01 | Each existing file's extension is in the accepted set for its category (§3). Files that are missing were already reported in Step 1 and are skipped here. | E |

### 5.3 Step 3 — Integrity & format

Each file is checked independently. A file that failed Step 1 or Step 2 is skipped. For mzML files, an extension other than `.mzml` also means the file is skipped.

#### mzML (each file, streamed with `iterparse`, so large files are OK)

| ID | Rule | Sev |
|---|---|---|
| VAL-MZ-01 | The file is well-formed XML. | E |
| VAL-MZ-02 | The root element is `<mzML>` or `<indexedmzML>` (namespace ignored). | E |
| VAL-MZ-03 | The file has at least one `<spectrum>` element. | **W** |
| VAL-MZ-04 | Each `<spectrum>` has non-empty `id`, plus `index` and `defaultArrayLength` attributes. | E |
| VAL-MZ-05 | Each spectrum has a cvParam `MS:1000511` (ms level). | E |
| VAL-MZ-06 | Each spectrum has a cvParam `MS:1000016` (scan start time) **with** a `unitAccession`. | E |
| VAL-MZ-07 | Each spectrum has exactly one `MS:1000514` (m/z array) and exactly one `MS:1000515` (intensity array). Binary arrays with other accessions are ignored. | E |
| VAL-MZ-08 | Reporting cap: individual errors are reported for the first 5 bad spectra in a file, then one `...and N more spectra` ERROR. | — |

#### Sample metadata

| ID | Rule | Sev |
|---|---|---|
| VAL-MD-01 | The table has at least one row. | E |
| VAL-MD-02 | It has a filename column (aliases: `filename`, `file_name`). | E |
| VAL-MD-03 | It has a grouping column. An exact match on `Class`, `Group` or `Condition` is preferred. Otherwise the first column whose normalised name **contains** `attribute` is used (GNPS `ATTRIBUTE_<name>`). | E |
| VAL-MD-04 | There are no missing or blank values in the filename column or the grouping column. | E |
| VAL-MD-05 | Filename values are unique (exact string comparison). | E |

If the metadata fails any of these checks, it is **not** used in Step 4.

#### Data matrix

| ID | Rule | Sev |
|---|---|---|
| VAL-DM-01 | The matrix has at least 1 row and at least 2 columns. | E |
| VAL-DM-02 | It has an identifier column (aliases: `feature_id`, `id`, `peak_id`, `row ID`). | E |
| VAL-DM-03 | Every non-missing cell outside the identifier column parses as a number. This includes `row m/z` and `row retention time` when present. The message reports the count and up to 5 examples, with 1-based file row numbers that count the header. | E |
| VAL-DM-04 | Missing or NA cells (non-detected features) are allowed. Their count is reported. | **W** |
| VAL-DM-05 | Identifier values are unique. | E |
| VAL-DM-06 | The MZmine full shape is accepted. The `row m/z`/`m/z`/`mz` and `row retention time`/`retention time`/`rt` columns are feature descriptors. They are excluded from the sample columns used in Step 4. | — |

Checks run in order and stop at the first ERROR. For example, duplicate IDs (VAL-DM-05) are not reported if VAL-DM-03 already failed. A matrix that fails is not used in Step 4.

#### Full feature table (required, MZmine 3/4; see VAL-EX-01)

Per-sample columns are found by matching the **raw** header against `datafile:<file>:<field>` (case-insensitive, surrounding whitespace trimmed). `<file>` is the first `:`-separated segment after `datafile:`. `<field>` is everything after it and may itself be nested, e.g. `datafile:A.mzML:rt_range:min` has field `rt_range:min`. VAL-TB-04 normalisation is not used for these columns, because it would strip the colons.

| ID | Rule | Sev |
|---|---|---|
| VAL-FQ-08 | Wrong-export detection. Runs first. If the table has a `row ID` column **and** any column whose name contains `peak area` or `peak height`, it is MZmine's GNPS/IIMN quant export or the legacy MZmine 2 CSV. One ERROR asks the user to re-export with *Export to CSV file (modular)*, and no other feature-table checks run. | E |
| VAL-FQ-01 | It has an `id` column. *(Modified in 2.0.0: `row ID`/`feature_id` are no longer accepted.)* | E |
| VAL-FQ-02 | It has an `mz` column. *(Modified in 2.0.0: `row m/z`/`m/z` are no longer accepted.)* | E |
| VAL-FQ-03 | It has an `rt` column. *(Modified in 2.0.0: `row retention time`/`retention time` are no longer accepted.)* | E |
| VAL-FQ-04 | It has at least one per-sample `datafile:<file>:area` or `datafile:<file>:height` column. *(Modified in 2.0.0: previously any column containing `peak area`/`peak height`.)* | E |
| VAL-FQ-06 | Full-table field check. The row-level columns `rt_range:min`, `rt_range:max`, `mz_range:min`, `mz_range:max`, `height` and `area` are all present (exact header, case-insensitive). The message lists every missing one. | E |
| VAL-FQ-07 | Every sample that has any `datafile:<file>:…` column has both a `:area` and a `:height` column. The message lists the incomplete samples. | E |
| VAL-FQ-05 | ID values are unique. | E |
| VAL-FQ-09 | `id` values are non-blank integers. `mz` and `rt` values are numeric and non-blank. Per-sample `:area`/`:height` cells are numeric or blank (MZmine leaves non-detected features blank). The message reports counts and up to 5 examples, with 1-based file row numbers that count the header. | E |

Checks run in the table order above. FQ-01..04, FQ-06 and FQ-07 are all reported together. FQ-05 and FQ-09 run only if those passed. A table that fails is not used in Step 4.

#### MGF (optional, MZmine export, read line by line; see VAL-EX-01)

| ID | Rule | Sev |
|---|---|---|
| VAL-MG-01 | The file has at least one `BEGIN IONS`. | E |
| VAL-MG-02 | The `BEGIN IONS` count equals the `END IONS` count. | E |
| VAL-MG-03 | Each block has `PEPMASS=`. | E |
| VAL-MG-04 | Each block has `RTINSECONDS=`. *(Modified in 2.0.0: `RT=` is no longer accepted.)* | E |
| VAL-MG-05 | Each block has `FEATURE_ID=` with an integer value (`^\d+$`). *(Modified in 2.0.0: `TITLE=`/`SCANS=` alone no longer satisfy it. They are still allowed as extra tags.)* | E |
| VAL-MG-06 | `CHARGE=` is **optional**. When present, it must match `^\d+[+-]?$` (e.g. `1`, `1+`, `2-`). | E |
| VAL-MG-07 | Inside a block, any line that contains `=` is a tag line. Tag keys are case-insensitive. Every other non-blank line is a peak line. | — |
| VAL-MG-08 | Each peak line is `<m/z><space/tab><intensity>` followed by an optional third token. Both numbers may be signed decimals in scientific notation. | E |
| VAL-MG-09 | The file has at least one peak line in total. | E |
| VAL-MG-10 | Reporting cap: errors are reported for the first 5 bad blocks, then one `...and N more` ERROR. Bad peak lines are reported as a single count. | — |
| VAL-MG-11 | Each block has `MSLEVEL=`. Together with MG-04 and MG-05 this is the MZmine per-block fingerprint. | E |

The identifier value is taken from the first `FEATURE_ID` line in each block. These values are kept for VAL-CP-07 and VAL-CP-09, but only when the whole MGF passed.

#### Annotation list (optional)

| ID | Rule | Sev |
|---|---|---|
| VAL-AN-01 | It has an ID column (`feature_id`, `Scan`, `Query_ID`, `id`). | E |
| VAL-AN-02 | It has a name column (`Compound_Name`, `Metabolite_Name`). | E |
| VAL-AN-03 | It has an ion column (`Adduct`, `Precursor_mz`). | E |
| VAL-AN-04 | It has a score column (`Score`, `Cosine_Score`, `MSI_Level`). | E |

### 5.4 Step 4 — Compatibility check

**VAL-CP-00**
- A cross-check runs only when the files involved passed Step 3.
- If a file wasn't parsed, there's no separate compatibility issue for it, because the integrity error already reported the problem.
- Filename comparisons use VAL-FN-01.

| ID | Rule | Sev |
|---|---|---|
| VAL-CP-01 | Every supplied mzML file matches a metadata `filename` row. | E |
| VAL-CP-02 | Every metadata `filename` row matches a supplied mzML file. | **W** |
| VAL-CP-03 | Every metadata `filename` matches a data-matrix sample column. | E |
| VAL-CP-04 | Every data-matrix sample column matches a metadata `filename`. The m/z and RT descriptor columns are excluded, per VAL-DM-06. | E |
| VAL-CP-05 | The feature table and the data matrix share at least one feature ID (exact string match). | **W** |
| VAL-CP-06 | Every annotation ID exists in the data matrix IDs. The message lists up to 10 that don't. | **W** |
| VAL-CP-07 | The annotation IDs overlap with either the MGF identifiers or the data-matrix IDs. This check runs only when the annotation list, the MGF identifiers and the data matrix are all available. | **W** |
| VAL-CP-08a | Every metadata `filename` matches a feature-table sample (`<file>` from `datafile:<file>:…`). | E |
| VAL-CP-08b | Every feature-table sample matches a metadata `filename`. | **W** |
| VAL-CP-09 | Every MGF `FEATURE_ID` exists among the feature-table `id` values (exact string match). This confirms both files come from the same MZmine feature list. The message lists up to 10 that don't. | E |

## 6. Non-functional requirements

**VAL-NF-01 (Read-only)**
- The validator never modifies, moves or creates input files. The only file it writes is `--report-out` when that flag is given.

**VAL-NF-02 (Streaming)**
- mzML is parsed incrementally, with elements cleared after use. Memory must not scale with file size. The real sample files are about 55 MB each, and much larger files are expected.

**VAL-NF-03 (Independence)**
- The module runs on its own. It must not import from other `modules/*` folders, and other modules must not import it except through the public API in §4.1.
- The duplicated helpers are intentional.

**VAL-NF-04 (Dependencies)**
- The only hard runtime dependency is `pandas`. Everything else comes from the standard library.
- Reading `.xlsx`/`.xls` also needs `openpyxl` or `xlrd`, but neither is currently listed in `requirements.txt` (see §8).

**VAL-NF-05 (No interaction)**
- The validator never prompts. It must be safe to run without a terminal (CI, Nextflow).

## 7. Test traceability

Tests build their fixtures in `tmp_path` and call the Python API. They don't use `sample_data/` or the CLI.

| Test (`tests/test_validator.py`) | Covers |
|---|---|
| `test_valid_required_only_passes` | PF-01, happy path |
| `test_missing_required_files_fails` | UP-01..03 |
| `test_optional_files_absent_are_ignored` | §3 optional semantics |
| `test_optional_file_present_but_malformed_fails` | Optional input present but malformed fails; uses an MGF missing RT (MG-04) |
| `test_malformed_mzml_fails` | MZ-01 |
| `test_metadata_missing_grouping_column_fails` | MD-03 |
| `test_data_matrix_non_numeric_cell_fails` | DM-03 |
| `test_mismatched_sample_names_fail_compatibility` | CP-03 / CP-04 |
| `test_semicolon_delimited_metadata_and_matrix_pass` | TB-01 |
| `test_mgf_accepts_feature_id_bare_charge_and_scientific_notation` | MG-05, MG-06, MG-08 |
| `test_mgf_still_rejects_missing_tags_and_non_numeric_peaks` | MG-03..05, MG-08 |
| `test_data_matrix_blank_cells_are_warning_not_error` | DM-04, PF-02 |
| `test_mgf_mixed_case_tag_line_not_mistaken_for_a_peak` | MG-07 |
| `test_mgf_charge_optional_when_undetermined` | MG-06 |
| `test_annotation_list_accepts_mzmine_id_column` | AN-01 |
| `test_data_matrix_accepts_mzmine_row_id_column` | DM-02 |
| `test_data_matrix_full_feature_table_shape_excludes_mz_and_rt_from_sample_match` | DM-06, CP-04 |
| `test_gnps_style_attribute_column_accepted` | MD-03 |
| `test_missing_feature_table_fails` | UP-06 |
| `test_valid_required_only_passes` (fixture) | FQ-01..04, FQ-06, FQ-07, FQ-09 happy path, EX-01 |
| `test_gnps_quant_table_rejected_with_reexport_message` | FQ-08, EX-01 |
| `test_full_table_missing_row_level_fields_fails` | FQ-06 |
| `test_full_table_sample_missing_height_fails` | FQ-07 |
| `test_feature_table_duplicate_ids_fail` | FQ-05 |
| `test_feature_table_non_integer_id_fails` | FQ-09 |
| `test_feature_table_xlsx_extension_rejected` | FT-01 |
| `test_mgf_without_feature_id_fails` | MG-05 |
| `test_mgf_without_mslevel_fails` | MG-11 |
| `test_mgf_rt_tag_instead_of_rtinseconds_fails` | MG-04 |
| `test_metadata_sample_missing_from_feature_table_fails` | CP-08a |
| `test_mgf_feature_id_not_in_feature_table_fails` | CP-09 |

**Not yet covered by a dedicated test:**
- UP-04, UP-05
- MZ-02..08
- MD-01, MD-02, MD-04, MD-05
- DM-01, DM-05
- MG-01, MG-02, MG-09
- AN-02..04
- CP-01, CP-02, CP-05..07, CP-08b
- RP-01, RP-02, CL-01, CL-02

New or changed requirements must come with a test (§9).

## 8. Known limitations and open questions

These are recorded here as they are, so nobody mistakes them for intended behaviour. Changing any of them is a spec change (§9).

1. **Missing-value tokens are narrow (VAL-TB-03).** A cell containing `NaN`, `NULL`, `-` or `inf`-like text other than `inf`/`nan` is treated as non-numeric and fails VAL-DM-03. Example: `NaN` with a capital N fails, but `nan`, `NA` and blank pass with a warning. This was checked against the current code.
2. **Missing `.xlsx` dependency.** `.xlsx`/`.xls` are accepted extensions, but `openpyxl`/`xlrd` is not in `requirements.txt`. In a clean environment, an Excel file would fail as "Could not parse".
3. **Loose grouping-column match.** The `attribute` fallback in VAL-MD-03 matches any column *containing* `attribute`, not only the `ATTRIBUTE_` prefix.
4. **Case-sensitive uniqueness, case-insensitive matching.** VAL-MD-05 compares filenames exactly, while compatibility matching ignores case. So `A.mzML` and `a.mzml` pass the uniqueness check but both match the same file.
5. ~~Feature table isn't checked against metadata.~~ Resolved in 2.0.0 by VAL-CP-08a/b.
6. **Data-matrix ID column type.** The identifier column itself is never checked for content (e.g. blank IDs), apart from uniqueness.
7. **Repeated scan start time.** If a spectrum contains several `MS:1000016` cvParams, only the last one's `unitAccession` counts for VAL-MZ-06.
8. **`.raw` files.** The metadata error text mentions `.raw`, but only `.mzML` is accepted as spectral input.
9. **README drift.** The module README refers to `modules/upload_validation`, `sample_data/valid/` and `sample_data/invalid/`. The actual paths are `Upload_and_validator_module/`, `sample_data/valid_files/` and `sample_data/invalid_files_examples/`.
10. **FQ-06 assumes default MZmine export settings.** If a user deselects row types (e.g. RT range, m/z range) in MZmine's CSV export dialog, a genuine full table fails VAL-FQ-06. The list was confirmed against the real export `sample_data/valid_files/Practice_data_MeOH_POS_full_feature_table.csv`.
11. **Filenames containing `:`.** Per-sample columns take `<file>` as the first `:`-separated segment, so a raw file name containing `:` would be split wrongly. Such names aren't allowed on macOS or Windows. The old GNPS quant table now lives at `sample_data/invalid_files_examples/feature_table_gnps_quant_not_full.csv` as a VAL-FQ-08 example.
12. **Feature-table vs data-matrix IDs are still only a WARNING (VAL-CP-05).** The data matrix may come from another tool, so its IDs aren't required to match MZmine's.

## 9. Change control

Any change to how the validator behaves needs these steps, in this order:

1. **Spec first.** Edit this file:
   - add, modify or retire `VAL-xx-nn` rows;
   - bump the spec version;
   - add a line to the changelog (§10).
   Moving an item out of §8 counts as a change.
2. **Tests.** Add or update a test in `tests/test_validator.py` for every new or changed requirement, and update the §7 table.
3. **Code.** Change `metabolomics_validator.py`. Keep the module docstring's "Format specs" section in sync with §5.
4. **Docs.** Update the module `README.md` if usage or flags changed.
5. **Verify.** `cd modules/Upload_and_validator_module && pytest -v && ./run_check.sh`. The real sample dataset must still print `VALIDATION PASSED`, unless the change deliberately says otherwise.

**Versioning** (spec version, semver):
- **MAJOR:** anything that can break a caller or turn a previously passing dataset into a failing one. Examples:
  - a new ERROR rule;
  - a WARNING raised to ERROR;
  - a narrower set of accepted aliases or extensions;
  - a change to the `InputSet`/`ValidationReport`/`Issue`/`Step` API, the JSON keys, `Step` strings or exit codes.
- **MINOR:** changes that only make the validator more lenient, or add things. Examples:
  - new accepted aliases, extensions or optional inputs;
  - new WARNING rules;
  - an ERROR lowered to WARNING;
  - new JSON keys.
- **PATCH:** message wording, test-only changes, refactors with no behaviour change, and documentation fixes.

**Guardrails:**
- Don't weaken the required inputs (§3), the "WARNINGs never fail" rule (VAL-PF-02) or the "always run all steps" rule (VAL-PF-03) without a MAJOR bump and a stated reason.
- Real-world exporter quirks (MZmine, GNPS, XCMS, msconvert) should be handled by adding **aliases** or **leniency**, not by special-casing files.
- Record which exporter motivated each change in the changelog.

## 10. Changelog

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-10-01 | Baseline spec, written from the existing implementation and test suite. Behaviour is unchanged. |
| 2.0.0 | 2026-10-20 | **MAJOR.** Exporter: MZmine 3/4. The feature table and MGF are now accepted only as MZmine output (new VAL-EX-01). The full feature table (*Export to CSV file (modular)*) is now **required** (new VAL-UP-06; CLI usage check VAL-CL-01 includes `--feature-table`), and its extension is restricted to `.csv`. FQ-01..04 modified to MZmine 3/4 column names (`id`, `mz`, `rt`, `datafile:<file>:area|height`). Added FQ-06 (full-table row-level fields), FQ-07 (per-sample area+height), FQ-08 (GNPS/IIMN quant and legacy MZmine 2 export rejected with a re-export message) and FQ-09 (value types). MGF: MG-04 now requires `RTINSECONDS`, MG-05 requires an integer `FEATURE_ID`, and MG-11 (`MSLEVEL`) is new. Compatibility: added CP-08a/b (feature-table samples vs metadata, moved out of §8 item 5) and CP-09 (MGF FEATURE_IDs must exist in the feature table). |
