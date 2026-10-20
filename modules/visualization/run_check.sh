#!/usr/bin/env bash
# Convenience wrapper: runs the visualization module against the files in
# a folder without retyping every flag.
#
# Usage:
#   ./run_check.sh path/to/folder
#   ./run_check.sh path/to/folder ATTRIBUTE_column_name
#   ./run_check.sh path/to/folder ATTRIBUTE_column_name normalization scaling
#
#   ATTRIBUTE_column_name  which metadata column to group by for
#                          PLS-DA/OPLS-DA/volcano (skips the prompt/
#                          ambiguity). Pass "" to leave it on auto/prompt
#                          while still setting normalization/scaling.
#   normalization          none | sum (default) | median
#   scaling                none | pareto (default) | auto
#
# Like the upload_validation module's run_check.sh, this detects the
# metadata and data-matrix files by keyword in the filename rather than
# assuming fixed names, and prints what it picked before running.
set -euo pipefail
cd "$(dirname "$0")"

DATA_DIR="${1:?Usage: ./run_check.sh path/to/folder [ATTRIBUTE_column_name] [normalization] [scaling]}"
GROUP_COLUMN="${2:-}"
NORMALIZATION="${3:-}"
SCALING="${4:-}"

if [ ! -d "$DATA_DIR" ]; then
    echo "Error: folder not found: $DATA_DIR" >&2
    exit 1
fi

find_by_keyword() {
    local dir="$1" kw="$2" f base lower
    for f in "$dir"/*.csv "$dir"/*.tsv; do
        [ -f "$f" ] || continue
        base="$(basename "$f")"
        lower="$(echo "$base" | tr '[:upper:]' '[:lower:]')"
        case "$lower" in
            *"$kw"*) echo "$f"; return 0 ;;
        esac
    done
}

MZML_FILES=("$DATA_DIR"/*.mzML)

METADATA_FILE="$(find_by_keyword "$DATA_DIR" metadata)"
MATRIX_FILE="$(find_by_keyword "$DATA_DIR" matrix)"

echo "Using:"
echo "  mzML files:    ${MZML_FILES[*]}"
echo "  Metadata:      ${METADATA_FILE:-<none found — optional>}"
echo "  Data matrix:   ${MATRIX_FILE:-<none found>}"
echo "  Group column:  ${GROUP_COLUMN:-<auto-resolve / prompt>}"
echo "  Normalization: ${NORMALIZATION:-<default: sum>}"
echo "  Scaling:       ${SCALING:-<default: pareto>}"
echo

if [ -z "$MATRIX_FILE" ]; then
    echo "Error: no data matrix file found in $DATA_DIR (looked for a filename containing 'matrix')." >&2
    exit 1
fi

ARGS=(--data-matrix "$MATRIX_FILE" --output-dir "$DATA_DIR/../viz_output")
[ -n "$METADATA_FILE" ] && ARGS+=(--metadata "$METADATA_FILE")
[ -f "${MZML_FILES[0]}" ] && ARGS+=(--mzml "${MZML_FILES[@]}")
[ -n "$GROUP_COLUMN" ] && ARGS+=(--group-column "$GROUP_COLUMN")
[ -n "$NORMALIZATION" ] && ARGS+=(--normalization "$NORMALIZATION")
[ -n "$SCALING" ] && ARGS+=(--scaling "$SCALING")

PYBIN="python3"
command -v python >/dev/null 2>&1 && PYBIN="python"

"$PYBIN" visualization_module.py "${ARGS[@]}"
