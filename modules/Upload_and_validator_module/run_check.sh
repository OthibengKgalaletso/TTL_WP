#!/usr/bin/env bash
# Convenience wrapper: validates the files in sample_data/valid/ (or a
# folder you pass as the first argument) without retyping every flag.
#
# Usage:
#   ./run_check.sh                  # uses sample_data/valid/
#   ./run_check.sh path/to/folder   # uses a different folder instead
set -euo pipefail
cd "$(dirname "$0")"

DATA_DIR="${1:-sample_data/valid_files}"
MZML_FILES=("$DATA_DIR"/*.mzML)

ARGS=(--mzml "${MZML_FILES[@]}" --metadata "$DATA_DIR/Practice_data_MeOH_POS_metadata.csv" --data-matrix "$DATA_DIR/DataMatrix_Example.csv")

[ -f "$DATA_DIR/Practice_data_MeOH_POS_iimn_gnps_quant.csv" ] && ARGS+=(--feature-table "$DATA_DIR/Practice_data_MeOH_POS_iimn_gnps_quant.csv")
[ -f "$DATA_DIR/Practice_data_MeOH_POS_iimn_gnps.mgf" ] && ARGS+=(--mgf "$DATA_DIR/Practice_data_MeOH_POS_iimn_gnps.mgf")
[ -f "$DATA_DIR/Practice_data_MeOH_POS_annotations.csv" ] && ARGS+=(--annotation-list "$DATA_DIR/Practice_data_MeOH_POS_annotations.csv")

# Prefer 'python' (always present in an active conda env) over 'python3'
# (not guaranteed in every environment), so this works regardless of which
# alias your shell happens to have.
PYBIN="python3"
command -v python >/dev/null 2>&1 && PYBIN="python"

"$PYBIN" metabolomics_validator.py "${ARGS[@]}"
