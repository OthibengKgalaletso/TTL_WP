#!/usr/bin/env python3
"""
metabolomics_validator.py

Upload & Validation Module for a metabolomics data-processing pipeline.

Mirrors this pipeline diagram:

    Upload files
        -> Validate file types
        -> Validate file integrity & format
        -> Compatibility check
        -> Validation passed? --YES--> (hand off to downstream module)
                              --NO---> Error report --> back to Upload files

Required inputs (marked * in the diagram):
    - mzML files        (one or more .mzML raw/processed spectra files)
    - Sample metadata   (.csv/.tsv/.txt table describing each sample)
    - Data matrix       (feature x sample intensity table)

Optional inputs:
    - Feature table       (e.g. MZmine/XCMS quantification export)
    - MGF file             (.mgf, MS/MS spectra for annotation)
    - Annotation list      (feature ID -> putative compound ID)

Format specs enforced by "Validate file integrity & format":

    MGF file
      - Every spectrum wrapped in BEGIN IONS / END IONS.
      - Required header tags per block: PEPMASS=, an identifier (TITLE=,
        or FEATURE_ID=/SCANS= as used by MZmine/GNPS exports), and
        RTINSECONDS= (or RT=). CHARGE= is validated for format (bare
        digit or with a '+'/'-' sign, e.g. 1 or 1+) when present, but not
        required — real exporters (e.g. MZmine) omit it for single-scan,
        non-merged features where charge can't be determined.
      - Body: numeric "[m/z] [intensity]" pairs, space- or tab-separated;
        scientific notation and an extra trailing column are tolerated.

    mzML file
      - Valid PSI-MS schema XML root (<mzML> or <indexedmzML>).
      - Each <spectrum> has id, index, defaultArrayLength attributes.
      - ms level given via cvParam accession MS:1000511.
      - scan start time given via cvParam accession MS:1000016 with a
        unit accession.
      - Exactly two binary data arrays: m/z array (MS:1000514) and
        intensity array (MS:1000515).

    Sample metadata (.csv/.tsv/.txt)
      - Required columns: filename (matches the .mzML/.raw inputs) and
        a grouping column (ATTRIBUTE/class, Class, Group, or Condition).
      - No missing values in required columns; filenames unique.

    Data matrix
      - Required columns: feature_id (or id, peak_id) plus one sample
        column per filename in the metadata.
      - Peak values must be numeric or blank; non-numeric strings are an
        error, blank/NA cells (non-detected features) are a warning.

    Feature quantification table (MZmine / XCMS style export)
      - Required columns: row ID (or feature_id), row m/z (or m/z),
        row retention time (or rt), and per-sample peak area/height
        columns.

    Metabolite annotation list
      - Required columns: feature_id (or Scan, Query_ID, id — MZmine's
        own annotation export uses plain 'id'), Compound_Name (or
        Metabolite_Name), Adduct or Precursor_mz, and a match metric
        (Score, Cosine_Score, or MSI_Level).

Usage:
    python metabolomics_validator.py \
        --mzml sample1.mzML sample2.mzML \
        --metadata sample_metadata.csv \
        --data-matrix data_matrix.csv \
        [--feature-table feature_table.csv] \
        [--mgf spectra.mgf] \
        [--annotation-list annotations.csv] \
        [--report-format text|json] [--report-out report.txt]

Exit code: 0 if validation passed, 1 otherwise (useful for CI / pipelines).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Optional

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    print("This script requires pandas. Install it with: pip install pandas", file=sys.stderr)
    raise


# --------------------------------------------------------------------------- #
# Data model
# --------------------------------------------------------------------------- #

class Severity(str, Enum):
    ERROR = "ERROR"      # blocks validation
    WARNING = "WARNING"  # flagged, does not block validation


class Step(str, Enum):
    UPLOAD = "Upload files"
    FILE_TYPE = "Validate file types"
    INTEGRITY = "Validate file integrity & format"
    COMPATIBILITY = "Compatibility check"


@dataclass
class Issue:
    step: Step
    severity: Severity
    file: str
    message: str

    def to_dict(self) -> dict:
        return {
            "step": self.step.value,
            "severity": self.severity.value,
            "file": self.file,
            "message": self.message,
        }


@dataclass
class InputSet:
    mzml_files: list[str] = field(default_factory=list)
    metadata_file: Optional[str] = None
    data_matrix_file: Optional[str] = None
    feature_table_file: Optional[str] = None
    mgf_file: Optional[str] = None
    annotation_list_file: Optional[str] = None


# Expected extensions per input category (case-insensitive)
EXPECTED_EXTENSIONS = {
    "mzML files": {".mzml"},
    "Sample metadata": {".csv", ".tsv", ".txt"},
    "Data matrix": {".csv", ".tsv", ".txt", ".xlsx", ".xls"},
    "Feature table": {".csv", ".tsv", ".txt", ".xlsx", ".xls"},
    "MGF file": {".mgf"},
    "Annotation list": {".csv", ".tsv", ".txt", ".xlsx", ".xls"},
}


# --------------------------------------------------------------------------- #
# Small text/column-matching helpers
# --------------------------------------------------------------------------- #

def _norm(s: str) -> str:
    """Lowercase and strip everything but letters/digits, for lenient header matching."""
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def _local_tag(tag: str) -> str:
    """Strip an XML namespace off an element tag, e.g. '{ns}spectrum' -> 'spectrum'."""
    return tag.rsplit("}", 1)[-1] if "}" in tag else tag


def _find_col(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    """Return the first column of df whose normalized name matches one of candidates."""
    norm_candidates = {_norm(c) for c in candidates}
    for col in df.columns:
        if _norm(col) in norm_candidates:
            return col
    return None


def _find_col_containing(df: pd.DataFrame, substrings: list[str]) -> Optional[str]:
    """Return the first column whose normalized name *contains* one of the substrings."""
    subs = [_norm(s) for s in substrings]
    for col in df.columns:
        nc = _norm(col)
        if any(s in nc for s in subs):
            return col
    return None


def _find_cols_containing(df: pd.DataFrame, substrings: list[str]) -> list[str]:
    subs = [_norm(s) for s in substrings]
    return [col for col in df.columns if any(s in _norm(col) for s in subs)]


def _filename_variants(name: str) -> set[str]:
    """Casefolded filename and its extension-stripped stem, for loose filename matching."""
    name = str(name).strip()
    variants = {name.lower()}
    variants.add(os.path.splitext(name)[0].lower())
    return variants


# --------------------------------------------------------------------------- #
# Validator
# --------------------------------------------------------------------------- #

class MetabolomicsValidator:
    """Runs the four-step validation pipeline shown in the diagram."""

    def __init__(self, inputs: InputSet):
        self.inputs = inputs
        self.issues: list[Issue] = []

        # Parsed artifacts, populated as integrity checks succeed, reused
        # in the compatibility-check step.
        self._metadata_df: Optional[pd.DataFrame] = None
        self._metadata_filename_col: Optional[str] = None
        self._matrix_df: Optional[pd.DataFrame] = None
        self._matrix_feature_col: Optional[str] = None
        self._matrix_sample_cols: list[str] = []
        self._feature_table_df: Optional[pd.DataFrame] = None
        self._feature_table_id_col: Optional[str] = None
        self._feature_table_sample_cols: list[str] = []
        self._annotation_df: Optional[pd.DataFrame] = None
        self._annotation_id_col: Optional[str] = None
        self._mgf_titles: list[str] = []

    # -- generic helpers ------------------------------------------------ #

    def _add(self, step: Step, severity: Severity, file: str, message: str) -> None:
        self.issues.append(Issue(step=step, severity=severity, file=file, message=message))

    def _read_table(self, path: str) -> Optional[pd.DataFrame]:
        """Best-effort load of a csv/tsv/txt/xlsx table, returns None on failure."""
        ext = Path(path).suffix.lower()
        try:
            if ext in (".xlsx", ".xls"):
                return pd.read_excel(path, dtype=str)
            if ext == ".tsv":
                return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["", "NA", "N/A", "nan"])
            # .csv / .txt: sniff delimiter among comma, tab, and semicolon
            # (semicolon-delimited CSV is common from some regional Excel
            # locales, even though the file extension is still .csv).
            with open(path, "r", newline="", errors="replace") as fh:
                first_line = fh.readline().lstrip("﻿")  # strip a UTF-8 BOM if present
            counts = {",": first_line.count(","), "\t": first_line.count("\t"), ";": first_line.count(";")}
            sep = max(counts, key=counts.get) if max(counts.values()) > 0 else ","
            return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, na_values=["", "NA", "N/A", "nan"], encoding="utf-8-sig")
        except Exception:
            return None

    # -- Step 1: Upload files ------------------------------------------ #

    def check_upload(self) -> bool:
        """Confirms required files were supplied and every path exists."""
        ok = True

        if not self.inputs.mzml_files:
            self._add(Step.UPLOAD, Severity.ERROR, "mzML files", "No mzML files were supplied (required).")
            ok = False
        if not self.inputs.metadata_file:
            self._add(Step.UPLOAD, Severity.ERROR, "Sample metadata", "No sample metadata file was supplied (required).")
            ok = False
        if not self.inputs.data_matrix_file:
            self._add(Step.UPLOAD, Severity.ERROR, "Data matrix", "No data matrix file was supplied (required).")
            ok = False

        all_paths = list(self.inputs.mzml_files)
        for p in (
            self.inputs.metadata_file,
            self.inputs.data_matrix_file,
            self.inputs.feature_table_file,
            self.inputs.mgf_file,
            self.inputs.annotation_list_file,
        ):
            if p:
                all_paths.append(p)

        for p in all_paths:
            if not os.path.isfile(p):
                self._add(Step.UPLOAD, Severity.ERROR, p, "File not found on disk.")
                ok = False
            elif os.path.getsize(p) == 0:
                self._add(Step.UPLOAD, Severity.ERROR, p, "File is empty (0 bytes).")
                ok = False

        return ok

    # -- Step 2: Validate file types ------------------------------------ #

    def check_file_types(self) -> bool:
        """Confirms each supplied file's extension matches what's expected."""
        ok = True

        def check_one(category: str, path: Optional[str]) -> None:
            nonlocal ok
            if not path or not os.path.isfile(path):
                return  # already reported (or optional & absent)
            ext = Path(path).suffix.lower()
            expected = EXPECTED_EXTENSIONS[category]
            if ext not in expected:
                self._add(
                    Step.FILE_TYPE, Severity.ERROR, path,
                    f"Unexpected file type '{ext or '(none)'}' for {category}; "
                    f"expected one of {sorted(expected)}.",
                )
                ok = False

        for p in self.inputs.mzml_files:
            check_one("mzML files", p)
        check_one("Sample metadata", self.inputs.metadata_file)
        check_one("Data matrix", self.inputs.data_matrix_file)
        check_one("Feature table", self.inputs.feature_table_file)
        check_one("MGF file", self.inputs.mgf_file)
        check_one("Annotation list", self.inputs.annotation_list_file)

        return ok

    # -- Step 3: Validate file integrity & format ------------------------ #

    def check_integrity_and_format(self) -> bool:
        ok = True
        ok &= self._check_mzml_files()
        ok &= self._check_metadata()
        ok &= self._check_data_matrix()
        ok &= self._check_feature_table()
        ok &= self._check_mgf()
        ok &= self._check_annotation_list()
        return ok

    # mzML: PSI-MS schema root, required spectrum attrs/CV params, binary arrays.
    def _check_mzml_files(self) -> bool:
        ok = True
        for path in self.inputs.mzml_files:
            if not os.path.isfile(path) or Path(path).suffix.lower() != ".mzml":
                continue  # already flagged upstream

            try:
                saw_root = False
                spectrum_count = 0
                bad_spectra: list[tuple[str, list[str]]] = []
                current: Optional[dict] = None

                for event, elem in ET.iterparse(path, events=("start", "end")):
                    tag = _local_tag(elem.tag)

                    if event == "start":
                        if not saw_root:
                            saw_root = True
                            if tag not in ("indexedmzML", "mzML"):
                                self._add(
                                    Step.INTEGRITY, Severity.ERROR, path,
                                    f"Root element '<{tag}>' is not a recognized mzML root "
                                    "(expected <mzML> or <indexedmzML>).",
                                )
                                ok = False
                        if tag == "spectrum":
                            current = {
                                "id": elem.get("id"),
                                "index": elem.get("index"),
                                "default_array_length": elem.get("defaultArrayLength"),
                                "has_ms_level": False,
                                "has_scan_start_time": False,
                                "binary_array_accessions": [],
                            }
                        elif tag == "cvParam" and current is not None:
                            accession = elem.get("accession")
                            if accession == "MS:1000511":
                                current["has_ms_level"] = True
                            elif accession == "MS:1000016":
                                current["has_scan_start_time"] = elem.get("unitAccession") is not None
                            elif accession in ("MS:1000514", "MS:1000515"):
                                current["binary_array_accessions"].append(accession)

                    elif event == "end":
                        if tag == "spectrum" and current is not None:
                            spectrum_count += 1
                            missing: list[str] = []
                            if not current["id"]:
                                missing.append("'id' attribute")
                            if current["index"] is None:
                                missing.append("'index' attribute")
                            if current["default_array_length"] is None:
                                missing.append("'defaultArrayLength' attribute")
                            if not current["has_ms_level"]:
                                missing.append("ms level cvParam (accession MS:1000511)")
                            if not current["has_scan_start_time"]:
                                missing.append("scan start time cvParam (accession MS:1000016) with a unit accession")
                            accessions = sorted(current["binary_array_accessions"])
                            if accessions != ["MS:1000514", "MS:1000515"]:
                                missing.append(
                                    "exactly one m/z array (MS:1000514) and one intensity array "
                                    f"(MS:1000515); found {current['binary_array_accessions'] or 'none'}"
                                )
                            if missing:
                                label = current["id"] or f"index {current['index']}"
                                bad_spectra.append((label, missing))
                            current = None
                        elem.clear()

                if spectrum_count == 0:
                    self._add(
                        Step.INTEGRITY, Severity.WARNING, path,
                        "No <spectrum> elements found; file may be empty of scan data.",
                    )
                elif bad_spectra:
                    ok = False
                    shown = bad_spectra[:5]
                    for label, missing in shown:
                        self._add(
                            Step.INTEGRITY, Severity.ERROR, path,
                            f"Spectrum '{label}' is missing: {'; '.join(missing)}.",
                        )
                    remaining = len(bad_spectra) - len(shown)
                    if remaining > 0:
                        self._add(
                            Step.INTEGRITY, Severity.ERROR, path,
                            f"...and {remaining} more spectra with missing required attributes/CV params.",
                        )

            except ET.ParseError as e:
                self._add(Step.INTEGRITY, Severity.ERROR, path, f"Malformed XML / corrupt mzML file: {e}")
                ok = False
            except Exception as e:
                self._add(Step.INTEGRITY, Severity.ERROR, path, f"Could not read file: {e}")
                ok = False
        return ok

    # Sample metadata: filename + grouping column required, no missing values, unique filenames.
    def _check_metadata(self) -> bool:
        path = self.inputs.metadata_file
        if not path or not os.path.isfile(path):
            return True
        df = self._read_table(path)
        if df is None:
            self._add(Step.INTEGRITY, Severity.ERROR, path, "Could not parse file as a table (csv/tsv/txt).")
            return False
        if df.empty:
            self._add(Step.INTEGRITY, Severity.ERROR, path, "Metadata table has no rows.")
            return False

        ok = True
        filename_col = _find_col(df, ["filename", "file_name"])
        if filename_col is None:
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                "Required column 'filename' not found (must match the input .mzML/.raw file names).",
            )
            ok = False

        # GNPS-style metadata sheets use ATTRIBUTE_<anything> as the
        # grouping column name (e.g. ATTRIBUTE_Tissue, ATTRIBUTE_Species) —
        # not a single fixed name — so any ATTRIBUTE_/ATTRIBUTE/-prefixed
        # column counts, in addition to the plain Class/Group/Condition names.
        class_col = _find_col(df, ["Class", "Group", "Condition"]) or _find_col_containing(df, ["attribute"])
        if class_col is None:
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                "Required grouping column not found (expected 'Class', 'Group', 'Condition', "
                "or a GNPS-style 'ATTRIBUTE_<name>' column, e.g. ATTRIBUTE_Tissue).",
            )
            ok = False

        if not ok:
            return False

        if df[filename_col].isna().any() or (df[filename_col].astype(str).str.strip() == "").any():
            self._add(Step.INTEGRITY, Severity.ERROR, path, "'filename' column contains missing values.")
            ok = False
        if df[class_col].isna().any() or (df[class_col].astype(str).str.strip() == "").any():
            self._add(Step.INTEGRITY, Severity.ERROR, path, f"'{class_col}' column contains missing values.")
            ok = False

        dup_mask = df[filename_col].duplicated()
        if dup_mask.any():
            dupes = df.loc[dup_mask, filename_col].tolist()
            self._add(Step.INTEGRITY, Severity.ERROR, path, f"Duplicate sample filenames found: {dupes}")
            ok = False

        if not ok:
            return False

        self._metadata_df = df
        self._metadata_filename_col = filename_col
        return True

    # Data matrix: feature_id + sample columns, strictly numeric peak values.
    def _check_data_matrix(self) -> bool:
        path = self.inputs.data_matrix_file
        if not path or not os.path.isfile(path):
            return True
        df = self._read_table(path)
        if df is None:
            self._add(Step.INTEGRITY, Severity.ERROR, path, "Could not parse file as a table (csv/tsv/txt/xlsx).")
            return False
        if df.shape[0] == 0 or df.shape[1] < 2:
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                "Data matrix must have at least one feature row and one identifier "
                "column plus one sample column.",
            )
            return False

        feature_id_col = _find_col(df, ["feature_id", "id", "peak_id", "row ID"])
        if feature_id_col is None:
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                "Required identifier column not found (expected one of: feature_id, id, peak_id, row ID).",
            )
            return False

        # Some data matrix exports use MZmine's full feature-table shape
        # ("row ID, row m/z, row retention time, <sample columns>") rather
        # than a bare "id, sample1, sample2..." matrix. When present, m/z
        # and retention time describe the feature itself, not a sample —
        # they're still validated as numeric peak-data cells below, but
        # excluded from the per-sample metadata/filename compatibility
        # check further down (see self._matrix_sample_cols).
        mz_col = _find_col(df, ["row m/z", "m/z", "mz"])
        rt_col = _find_col(df, ["row retention time", "retention time", "rt"])
        descriptive_cols = {c for c in (mz_col, rt_col) if c is not None}

        sample_cols = [c for c in df.columns if c != feature_id_col]
        numeric_part = df[sample_cols].apply(pd.to_numeric, errors="coerce")
        bad_cells = numeric_part.isna() & ~df[sample_cols].isna()
        bad_count = int(bad_cells.to_numpy().sum())
        if bad_count > 0:
            examples = []
            rows, cols = bad_cells.to_numpy().nonzero()
            for r, c in list(zip(rows, cols))[:5]:
                examples.append(f"row {r + 2}, column '{sample_cols[c]}' = {df[sample_cols].iat[r, c]!r}")
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                f"{bad_count} non-numeric value(s) found in peak data cells (must be strictly numeric); "
                f"e.g. {'; '.join(examples)}.",
            )
            return False

        # Missing/blank cells (non-detected features) are common in real
        # untargeted-metabolomics data matrices — flagged as a warning
        # (visible, non-blocking) rather than a hard failure.
        missing_count = int(df[sample_cols].isna().to_numpy().sum())
        if missing_count > 0:
            self._add(
                Step.INTEGRITY, Severity.WARNING, path,
                f"{missing_count} missing/NA value(s) found in peak data cells "
                "(commonly non-detected features — not blocking, but worth checking).",
            )

        if df[feature_id_col].duplicated().any():
            self._add(Step.INTEGRITY, Severity.ERROR, path, "Duplicate feature IDs in the identifier column.")
            return False

        self._matrix_df = df
        self._matrix_feature_col = feature_id_col
        self._matrix_sample_cols = [c for c in sample_cols if c not in descriptive_cols]
        return True

    # Feature quantification table (MZmine / XCMS style export).
    def _check_feature_table(self) -> bool:
        path = self.inputs.feature_table_file
        if not path or not os.path.isfile(path):
            return True  # optional and absent
        df = self._read_table(path)
        if df is None:
            self._add(Step.INTEGRITY, Severity.ERROR, path, "Could not parse file as a table (csv/tsv/txt/xlsx).")
            return False

        ok = True
        id_col = _find_col(df, ["row ID", "feature_id"])
        if id_col is None:
            self._add(Step.INTEGRITY, Severity.ERROR, path, "Required column not found (expected 'row ID' or 'feature_id').")
            ok = False

        mz_col = _find_col(df, ["row m/z", "m/z", "mz"])
        if mz_col is None:
            self._add(Step.INTEGRITY, Severity.ERROR, path, "Required column not found (expected 'row m/z' or 'm/z').")
            ok = False

        rt_col = _find_col(df, ["row retention time", "retention time", "rt"])
        if rt_col is None:
            self._add(Step.INTEGRITY, Severity.ERROR, path, "Required column not found (expected 'row retention time' or 'rt').")
            ok = False

        sample_cols = _find_cols_containing(df, ["peak area", "peak height"])
        if not sample_cols:
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                "No per-sample 'Peak area'/'Peak height' columns found.",
            )
            ok = False

        if not ok:
            return False

        if df[id_col].duplicated().any():
            self._add(Step.INTEGRITY, Severity.ERROR, path, f"Duplicate values in '{id_col}' column.")
            return False

        self._feature_table_df = df
        self._feature_table_id_col = id_col
        self._feature_table_sample_cols = sample_cols
        return True

    # MGF: BEGIN/END IONS blocks, required header tags, numeric peak pairs.
    def _check_mgf(self) -> bool:
        path = self.inputs.mgf_file
        if not path or not os.path.isfile(path):
            return True  # optional and absent

        # CHARGE isn't always determinable — e.g. MZmine leaves it out for
        # single-scan, non-merged features — so it's format-checked when
        # present but not required outright.
        required_tags = {"PEPMASS"}
        rt_tags = {"RTINSECONDS", "RT"}
        # A spectrum block needs *some* identifier — TITLE is the spec name,
        # but real-world exports (e.g. MZmine/GNPS) commonly use FEATURE_ID
        # or SCANS instead, so any one of these satisfies the requirement.
        id_tags = {"TITLE", "FEATURE_ID", "SCANS"}
        # Numeric pairs, allowing scientific notation (e.g. 1.23e+05) and an
        # optional trailing 3rd column (some exporters add a charge/flag
        # column per peak) — only the first two tokens (m/z, intensity) are
        # required to be numeric.
        peak_pair_re = re.compile(
            r"^-?\d+(\.\d+)?([eE][+-]?\d+)?[ \t]+-?\d+(\.\d+)?([eE][+-]?\d+)?(?:[ \t]+\S+)?$"
        )
        # Charge sign (1+ / 1-) is common but many real exports just write
        # a bare digit (e.g. CHARGE=1) — accept both.
        charge_re = re.compile(r"^\d+[+-]?$")

        try:
            begin_count = 0
            end_count = 0
            block_index = 0
            bad_blocks: list[tuple[int, list[str]]] = []
            titles: list[str] = []
            in_block = False
            tags_seen: set[str] = set()
            charge_val: Optional[str] = None
            peak_lines_total = 0
            peak_lines_bad = 0
            title_val: Optional[str] = None

            with open(path, "r", errors="replace") as fh:
                for raw_line in fh:
                    stripped = raw_line.strip()
                    if stripped == "BEGIN IONS":
                        begin_count += 1
                        block_index += 1
                        in_block = True
                        tags_seen = set()
                        charge_val = None
                        title_val = None
                        continue
                    if stripped == "END IONS":
                        end_count += 1
                        in_block = False
                        missing = sorted(required_tags - tags_seen)
                        if not (rt_tags & tags_seen):
                            missing.append("RTINSECONDS (or RT)")
                        if not (id_tags & tags_seen):
                            missing.append("TITLE (or FEATURE_ID / SCANS)")
                        if missing:
                            bad_blocks.append((block_index, missing))
                        if charge_val is not None and not charge_re.match(charge_val):
                            bad_blocks.append((block_index, [f"CHARGE value '{charge_val}' not numeric (with optional '+'/'-')"]))
                        if title_val:
                            titles.append(title_val)
                        continue
                    if not in_block or not stripped:
                        continue
                    # A numeric "[m/z] [intensity]" peak line never contains
                    # '=', so any line with one is a header/tag line — not
                    # just ones with an ALL-CAPS key (real exports include
                    # mixed-case lines like "Num peaks=1").
                    if "=" in stripped:
                        key, _, val = stripped.partition("=")
                        key = key.strip().upper()
                        tags_seen.add(key)
                        if key == "CHARGE":
                            charge_val = val.strip()
                        if key in id_tags and not title_val:
                            title_val = val.strip()
                    else:
                        peak_lines_total += 1
                        if not peak_pair_re.match(stripped):
                            peak_lines_bad += 1

            if begin_count == 0:
                self._add(Step.INTEGRITY, Severity.ERROR, path, "No 'BEGIN IONS' blocks found; not a valid MGF file.")
                return False
            if begin_count != end_count:
                self._add(
                    Step.INTEGRITY, Severity.ERROR, path,
                    f"Unbalanced MGF blocks: {begin_count} 'BEGIN IONS' vs {end_count} 'END IONS'.",
                )
                return False

            ok = True
            if bad_blocks:
                ok = False
                for idx, missing in bad_blocks[:5]:
                    self._add(Step.INTEGRITY, Severity.ERROR, path, f"Spectrum block #{idx}: missing/invalid {', '.join(missing)}.")
                if len(bad_blocks) > 5:
                    self._add(Step.INTEGRITY, Severity.ERROR, path, f"...and {len(bad_blocks) - 5} more block(s) with missing/invalid tags.")

            if peak_lines_total == 0:
                self._add(Step.INTEGRITY, Severity.ERROR, path, "No numeric '[m/z] [intensity]' peak lines found in any spectrum block.")
                ok = False
            elif peak_lines_bad > 0:
                self._add(
                    Step.INTEGRITY, Severity.ERROR, path,
                    f"{peak_lines_bad} of {peak_lines_total} peak line(s) are not valid numeric "
                    "'[m/z] [intensity]' pairs.",
                )
                ok = False

            if ok:
                self._mgf_titles = titles
            return ok

        except Exception as e:
            self._add(Step.INTEGRITY, Severity.ERROR, path, f"Could not read file: {e}")
            return False

    # Metabolite annotation list.
    def _check_annotation_list(self) -> bool:
        path = self.inputs.annotation_list_file
        if not path or not os.path.isfile(path):
            return True  # optional and absent
        df = self._read_table(path)
        if df is None:
            self._add(Step.INTEGRITY, Severity.ERROR, path, "Could not parse file as a table (csv/tsv/txt/xlsx).")
            return False

        ok = True
        # MZmine's own annotation export just calls this column 'id'.
        id_col = _find_col(df, ["feature_id", "Scan", "Query_ID", "id"])
        if id_col is None:
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                "Required column not found (expected one of: feature_id, Scan, Query_ID, id).",
            )
            ok = False

        name_col = _find_col(df, ["Compound_Name", "Metabolite_Name"])
        if name_col is None:
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                "Required column not found (expected 'Compound_Name' or 'Metabolite_Name').",
            )
            ok = False

        ion_col = _find_col(df, ["Adduct", "Precursor_mz"])
        if ion_col is None:
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                "Required column not found (expected 'Adduct' or 'Precursor_mz').",
            )
            ok = False

        score_col = _find_col(df, ["Score", "Cosine_Score", "MSI_Level"])
        if score_col is None:
            self._add(
                Step.INTEGRITY, Severity.ERROR, path,
                "Required column not found (expected one of: Score, Cosine_Score, MSI_Level).",
            )
            ok = False

        if not ok:
            return False

        self._annotation_df = df
        self._annotation_id_col = id_col
        return True

    # -- Step 4: Compatibility check ------------------------------------ #

    def check_compatibility(self) -> bool:
        """Cross-checks that the parsed files agree with each other."""
        ok = True

        metadata_variants: list[tuple[str, set[str]]] = []
        if self._metadata_df is not None and self._metadata_filename_col is not None:
            metadata_variants = [
                (str(v), _filename_variants(v))
                for v in self._metadata_df[self._metadata_filename_col].astype(str)
            ]

        mzml_variants: list[tuple[str, set[str]]] = [
            (os.path.basename(p), _filename_variants(os.path.basename(p)))
            for p in self.inputs.mzml_files
            if os.path.isfile(p)
        ]

        # metadata['filename'] <-> mzML files
        if metadata_variants and mzml_variants:
            unmatched_mzml = [
                name for name, variants in mzml_variants
                if not any(variants & mv for _, mv in metadata_variants)
            ]
            if unmatched_mzml:
                self._add(
                    Step.COMPATIBILITY, Severity.ERROR, "Sample metadata",
                    f"{len(unmatched_mzml)} mzML file(s) have no matching 'filename' row in sample "
                    f"metadata: {unmatched_mzml}",
                )
                ok = False

            unmatched_metadata = [
                name for name, mv in metadata_variants
                if not any(mv & variants for _, variants in mzml_variants)
            ]
            if unmatched_metadata:
                self._add(
                    Step.COMPATIBILITY, Severity.WARNING, "Sample metadata",
                    f"{len(unmatched_metadata)} metadata 'filename' row(s) have no matching mzML "
                    f"file supplied: {unmatched_metadata}",
                )

        # metadata['filename'] <-> data matrix sample columns
        if metadata_variants and self._matrix_sample_cols:
            matrix_variants = [(c, _filename_variants(c)) for c in self._matrix_sample_cols]

            missing_in_matrix = [
                name for name, mv in metadata_variants
                if not any(mv & cv for _, cv in matrix_variants)
            ]
            if missing_in_matrix:
                self._add(
                    Step.COMPATIBILITY, Severity.ERROR, "Data matrix",
                    f"{len(missing_in_matrix)} sample(s) in metadata are missing a matching column "
                    f"in the data matrix: {missing_in_matrix}",
                )
                ok = False

            missing_in_metadata = [
                col for col, cv in matrix_variants
                if not any(cv & mv for _, mv in metadata_variants)
            ]
            if missing_in_metadata:
                self._add(
                    Step.COMPATIBILITY, Severity.ERROR, "Data matrix",
                    f"{len(missing_in_metadata)} data-matrix sample column(s) don't match any "
                    f"metadata filename: {missing_in_metadata}",
                )
                ok = False

        # feature table <-> data matrix: feature IDs should line up.
        if self._feature_table_df is not None and self._matrix_df is not None:
            ft_ids = set(self._feature_table_df[self._feature_table_id_col].astype(str))
            matrix_ids = set(self._matrix_df[self._matrix_feature_col].astype(str))
            if not (ft_ids & matrix_ids):
                self._add(
                    Step.COMPATIBILITY, Severity.WARNING, "Feature table",
                    "No overlapping feature IDs found between the feature table and the data matrix.",
                )

        # annotation list <-> data matrix: annotated feature IDs should exist in the matrix.
        if self._annotation_df is not None and self._matrix_df is not None:
            ann_ids = set(self._annotation_df[self._annotation_id_col].astype(str))
            matrix_ids = set(self._matrix_df[self._matrix_feature_col].astype(str))
            unmatched = ann_ids - matrix_ids
            if unmatched:
                self._add(
                    Step.COMPATIBILITY, Severity.WARNING, "Annotation list",
                    f"{len(unmatched)} annotated feature ID(s) not found in the data matrix: "
                    f"{sorted(unmatched)[:10]}{'...' if len(unmatched) > 10 else ''}",
                )

        # annotation list <-> MGF: annotation IDs referencing TITLE values should exist in the MGF.
        if self._annotation_df is not None and self._mgf_titles:
            ann_ids = set(self._annotation_df[self._annotation_id_col].astype(str))
            titles = set(self._mgf_titles)
            if not (ann_ids & titles) and self._matrix_df is not None and not (ann_ids & set(self._matrix_df[self._matrix_feature_col].astype(str))):
                self._add(
                    Step.COMPATIBILITY, Severity.WARNING, "Annotation list",
                    "Annotation IDs match neither the data matrix feature IDs nor the MGF spectrum titles.",
                )

        return ok

    # -- Orchestration ---------------------------------------------------- #

    def run(self) -> "ValidationReport":
        """Runs all four pipeline steps in order and returns the report."""
        upload_ok = self.check_upload()
        types_ok = self.check_file_types()
        integrity_ok = self.check_integrity_and_format()
        compat_ok = self.check_compatibility()

        passed = upload_ok and types_ok and integrity_ok and compat_ok
        return ValidationReport(
            passed=passed,
            issues=self.issues,
            timestamp=datetime.now().isoformat(timespec="seconds"),
        )


@dataclass
class ValidationReport:
    passed: bool
    issues: list[Issue]
    timestamp: str

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.severity == Severity.ERROR]

    @property
    def warnings(self) -> list[Issue]:
        return [i for i in self.issues if i.severity == Severity.WARNING]

    def to_dict(self) -> dict:
        return {
            "timestamp": self.timestamp,
            "validation_passed": self.passed,
            "error_count": len(self.errors),
            "warning_count": len(self.warnings),
            "issues": [i.to_dict() for i in self.issues],
        }

    def to_text(self) -> str:
        lines = []
        lines.append("=" * 70)
        lines.append("METABOLOMICS DATA — UPLOAD & VALIDATION MODULE REPORT")
        lines.append("=" * 70)
        lines.append(f"Run at: {self.timestamp}")
        lines.append("")
        lines.append(f"Result: {'VALIDATION PASSED' if self.passed else 'VALIDATION FAILED'}")
        lines.append(f"Errors: {len(self.errors)}   Warnings: {len(self.warnings)}")
        lines.append("")

        if not self.issues:
            lines.append("No issues found across any pipeline step. Ready to proceed.")
        else:
            for step in (Step.UPLOAD, Step.FILE_TYPE, Step.INTEGRITY, Step.COMPATIBILITY):
                step_issues = [i for i in self.issues if i.step == step]
                if not step_issues:
                    continue
                lines.append(f"--- {step.value} ---")
                for issue in step_issues:
                    lines.append(f"  [{issue.severity.value}] {issue.file}: {issue.message}")
                lines.append("")

        lines.append("-" * 70)
        if self.passed:
            lines.append("=> Validation passed (YES). Handing off to the next pipeline module.")
        else:
            lines.append("=> Validation failed (NO). See error report above.")
            lines.append("=> Please correct the file(s) above and re-upload.")
        lines.append("=" * 70)
        return "\n".join(lines)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate a metabolomics dataset (mzML + metadata + data matrix, "
                    "plus optional feature table / MGF / annotation list) before it "
                    "enters downstream processing.",
    )
    parser.add_argument("--mzml", nargs="+", metavar="FILE", help="One or more .mzML files (required).")
    parser.add_argument("--metadata", metavar="FILE", help="Sample metadata table (required).")
    parser.add_argument("--data-matrix", metavar="FILE", help="Feature x sample intensity matrix (required).")
    parser.add_argument("--feature-table", metavar="FILE", help="Feature quantification table (optional).")
    parser.add_argument("--mgf", metavar="FILE", help="MGF file of MS/MS spectra (optional).")
    parser.add_argument("--annotation-list", metavar="FILE", help="Metabolite annotation list (optional).")
    parser.add_argument("--report-format", choices=["text", "json"], default="text", help="Report output format.")
    parser.add_argument("--report-out", metavar="FILE", help="Write the report to this file (also prints to stdout).")
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if not (args.mzml or args.metadata or args.data_matrix):
        parser.error("Provide --mzml/--metadata/--data-matrix.")

    inputs = InputSet(
        mzml_files=args.mzml or [],
        metadata_file=args.metadata,
        data_matrix_file=args.data_matrix,
        feature_table_file=args.feature_table,
        mgf_file=args.mgf,
        annotation_list_file=args.annotation_list,
    )

    validator = MetabolomicsValidator(inputs)
    report = validator.run()

    output = json.dumps(report.to_dict(), indent=2) if args.report_format == "json" else report.to_text()
    print(output)

    if args.report_out:
        Path(args.report_out).write_text(output + "\n")
        print(f"\nReport written to: {args.report_out}")

    return 0 if report.passed else 1


if __name__ == "__main__":
    sys.exit(main())