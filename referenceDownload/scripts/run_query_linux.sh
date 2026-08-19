#!/usr/bin/env bash
# =============================================================================
# run_query_linux.sh
# -----------------------------------------------------------------------------
# Run an SQL query against the RNAcentral public PostgreSQL database from a
# Linux / Ubuntu shell (including WSL Ubuntu).
#
# Prerequisites:
#   - psql client  OR  Docker installed and running
#   - .env file in the parent directory (referenceDownload/) with DB credentials
#     OR export the variables manually (see below).
#
# Usage:
#   bash run_query_linux.sh --taxid 9913 --type miRNA
#   bash run_query_linux.sh --taxid 9913 --type miRNA --type tRNA --type rRNA
#   bash run_query_linux.sh --taxid 9913 --type miRNA,tRNA,rRNA
#   bash run_query_linux.sh --taxid 9913 --type miRNA --out bos_miRNA.ids.txt
#   bash run_query_linux.sh --taxid 9913 --list-types
#
# Options:
#   --taxid NUM        NCBI taxonomy ID (required)
#   --type  STRING     RNA type (e.g. miRNA, tRNA, rRNA).  Repeat the flag or
#                      supply a comma-separated list to fetch multiple types;
#                      all IDs are written to the same output file.
#                      Required unless --list-types is specified.
#   --out   FILE       Output IDs file path
#                      (default: <taxid>_<type>.ids.txt, or
#                       <taxid>_<type1>_<type2>.ids.txt for multiple types)
#   --meta-out FILE    Output metadata TSV path
#                      (default: <taxid>_<type>.metadata.tsv)
#   --no-metadata      Skip generating the sequence-metadata TSV
#   --types-out FILE   Output file for --list-types
#                      (default: <taxid>_rna_types.txt)
#   --list-types       List available RNA types for the taxid and save the
#                      summary to <taxid>_rna_types.txt
#   --use-docker       Force use of Docker even if psql is installed locally
#   --help             Show this help message
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

export PGSSLMODE=disable

# ---------------------------------------------------------------------------
# Load connection parameters from .env (if it exists)
# ---------------------------------------------------------------------------
ENV_FILE="$REPO_ROOT/.env"
if [[ -f "$ENV_FILE" ]]; then
    # shellcheck disable=SC1090
    set -a; source "$ENV_FILE"; set +a
fi

# ---------------------------------------------------------------------------
# Connection parameters – can be overridden by the .env file or by exporting
# them before calling this script.
# Get the current values from: https://rnacentral.org/help/public-database
# ---------------------------------------------------------------------------
RNACENTRAL_HOST="${RNACENTRAL_HOST:-hh-pgsql-public.ebi.ac.uk}"
RNACENTRAL_PORT="${RNACENTRAL_PORT:-5432}"
RNACENTRAL_DB="${RNACENTRAL_DB:-pfmegrnargs}"
RNACENTRAL_USER="${RNACENTRAL_USER:-reader}"
RNACENTRAL_PASSWORD="${RNACENTRAL_PASSWORD:-NWDMCE5xdipIjRrp}"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
TAXID=""
RNA_TYPES=()      # array; populated by --type (may be repeated or comma-sep)
OUT_FILE=""
META_OUT_FILE=""
TYPES_OUT_FILE=""
LIST_TYPES=false
NO_METADATA=false
USE_DOCKER=false

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
show_help() {
    sed -n '1,/^# ====/p' "$0" | grep '^#' | sed 's/^# \?//'
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --taxid)
            TAXID="$2"; shift 2 ;;
        --type)
            # Accept comma-separated list in a single --type value, e.g.
            # --type miRNA,tRNA  or  --type miRNA --type tRNA
            IFS=',' read -ra _PARTS <<< "$2"
            RNA_TYPES+=("${_PARTS[@]}")
            shift 2 ;;
        --out)
            OUT_FILE="$2"; shift 2 ;;
        --meta-out)
            META_OUT_FILE="$2"; shift 2 ;;
        --types-out)
            TYPES_OUT_FILE="$2"; shift 2 ;;
        --list-types)
            LIST_TYPES=true; shift ;;
        --no-metadata)
            NO_METADATA=true; shift ;;
        --use-docker)
            USE_DOCKER=true; shift ;;
        --help|-h)
            show_help ;;
        *)
            echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ -z "$TAXID" ]]; then
    echo "ERROR: --taxid is required." >&2
    exit 1
fi
if [[ "$LIST_TYPES" == false && "${#RNA_TYPES[@]}" -eq 0 ]]; then
    echo "ERROR: --type is required (unless --list-types is specified)." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Helper: build the SQL IN-list from RNA_TYPES array
#   e.g.  RNA_TYPES=(miRNA tRNA) -> 'miRNA','tRNA'
# ---------------------------------------------------------------------------
build_types_in() {
    local result=""
    for t in "${RNA_TYPES[@]}"; do
        [[ -n "$result" ]] && result="${result},"
        result="${result}'${t}'"
    done
    echo "$result"
}

# ---------------------------------------------------------------------------
# Build a safe underscore-joined string of all types (for file names)
# ---------------------------------------------------------------------------
types_label() {
    local IFS='_'
    echo "${RNA_TYPES[*]}"
}

# ---------------------------------------------------------------------------
# Set default output file names
# ---------------------------------------------------------------------------
if [[ -z "$OUT_FILE" && "$LIST_TYPES" == false ]]; then
    OUT_FILE="$(pwd)/${TAXID}_$(types_label).ids.txt"
fi
if [[ -z "$META_OUT_FILE" && "$LIST_TYPES" == false ]]; then
    META_OUT_FILE="$(pwd)/${TAXID}_$(types_label).metadata.tsv"
fi
if [[ -z "$TYPES_OUT_FILE" ]]; then
    TYPES_OUT_FILE="$(pwd)/${TAXID}_rna_types.txt"
fi

# ---------------------------------------------------------------------------
# Build the connection string
# ---------------------------------------------------------------------------
CONN_STR="host=${RNACENTRAL_HOST} port=${RNACENTRAL_PORT} dbname=${RNACENTRAL_DB} user=${RNACENTRAL_USER}"

# ---------------------------------------------------------------------------
# Determine whether to use local psql or Docker
# ---------------------------------------------------------------------------
if ! command -v psql &>/dev/null || [[ "$USE_DOCKER" == true ]]; then
    if command -v docker &>/dev/null; then
        USE_DOCKER=true
        echo "INFO: psql not found locally – using Docker (postgres:16 image)."
    else
        echo "ERROR: Neither 'psql' nor 'docker' found." >&2
        echo "       Install psql:  sudo apt install postgresql-client" >&2
        echo "       Or install Docker: https://docs.docker.com/get-docker/" >&2
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Helper: run a SQL file through psql (extra psql flags passed as "$@")
# ---------------------------------------------------------------------------
SQL_DIR="$REPO_ROOT/sql"

run_psql() {
    local sql_file="$1"
    shift
    if [[ "$USE_DOCKER" == true ]]; then
        # Docker – pipe the SQL file on stdin
        docker run --rm -i \
            -e PGPASSWORD="${RNACENTRAL_PASSWORD}" \
            postgres:16 psql \
            "${CONN_STR} password=${RNACENTRAL_PASSWORD}" \
            "$@" < "$sql_file"
    else
        PGPASSWORD="${RNACENTRAL_PASSWORD}" psql \
            "$CONN_STR" \
            "$@" -f "$sql_file"
    fi
}

# ===========================================================================
# --list-types: show RNA types available for the given taxid and save to file
# ===========================================================================
if [[ "$LIST_TYPES" == true ]]; then
    echo "==> Listing RNA types for taxid ${TAXID} ..."
    TMP_SQL="$(mktemp /tmp/list_types_XXXXXX.sql)"
    sed "s/:TAXID/${TAXID}/g" "$SQL_DIR/list_types_for_taxid.sql" > "$TMP_SQL"
    # Print to terminal AND save to file simultaneously
    run_psql "$TMP_SQL" | tee "$TYPES_OUT_FILE"
    rm -f "$TMP_SQL"
    echo "==> RNA type summary saved to: ${TYPES_OUT_FILE}"
    exit 0
fi

# ===========================================================================
# Fetch sequence IDs
# ===========================================================================
TYPES_IN="$(build_types_in)"
TYPES_DISPLAY="${RNA_TYPES[*]}"

echo "==> Querying IDs for taxid=${TAXID}, type(s): ${TYPES_DISPLAY} ..."
TMP_SQL="$(mktemp /tmp/ids_query_XXXXXX.sql)"
sed "s/:TAXID/${TAXID}/g; s/:RNA_TYPES/${TYPES_IN}/g" \
    "$SQL_DIR/ids_for_taxid_and_type.sql" > "$TMP_SQL"

run_psql "$TMP_SQL" -t -A > "$OUT_FILE"
rm -f "$TMP_SQL"

COUNT=$(wc -l < "$OUT_FILE")
echo "==> Done. ${COUNT} IDs written to: ${OUT_FILE}"

# ===========================================================================
# Fetch sequence metadata (TSV) – skipped when --no-metadata is passed
# Two-query approach with intelligent fallback:
#   1. Run WITH features query first (may be empty for RNA types like tRNA)
#   2. If output is empty (0 data rows), fall back to WITHOUT features query
# ===========================================================================

# Helper: run the fallback (without features) metadata query and write output.
# Uses TMP_META_ERR from the outer scope for error capture.
run_fallback_metadata_query() {
    local tmp_sql
    tmp_sql="$(mktemp /tmp/meta_query_XXXXXX.sql)"
    sed "s/:TAXID/${TAXID}/g; s/:RNA_TYPES/${TYPES_IN}/g" \
        "$SQL_DIR/sequence_metadata_for_taxid_and_type.sql" > "$tmp_sql"

    if run_psql "$tmp_sql" -A -F $'\t' -q > "$META_OUT_FILE" 2>"$TMP_META_ERR"; then
        local meta_count
        meta_count=$(wc -l < "$META_OUT_FILE")
        echo "==> Used fallback query WITHOUT features. Metadata (${meta_count} lines incl. header) written to: ${META_OUT_FILE}"
    else
        echo "WARNING: Metadata query (without features) failed – skipping metadata file." >&2
        cat "$TMP_META_ERR" >&2
        rm -f "$META_OUT_FILE"
    fi
    rm -f "$tmp_sql"
}

if [[ "$NO_METADATA" == false ]]; then
    echo "==> Querying metadata for taxid=${TAXID}, type(s): ${TYPES_DISPLAY} ..."

    TMP_META_WITH_SQL="$(mktemp /tmp/meta_with_features_XXXXXX.sql)"
    TMP_META_WITH_OUT="$(mktemp /tmp/meta_with_features_XXXXXX.tsv)"
    TMP_META_ERR="$(mktemp /tmp/meta_query_XXXXXX.err)"

    sed "s/:TAXID/${TAXID}/g; s/:RNA_TYPES/${TYPES_IN}/g" \
        "$SQL_DIR/sequence_metadata_for_taxid_and_type_with_features.sql" > "$TMP_META_WITH_SQL"

    # -A: unaligned  -F $'\t': tab separator  -q: quiet (no row-count footer)
    echo "==> Running query WITH sequence features ..."
    if run_psql "$TMP_META_WITH_SQL" -A -F $'\t' -q > "$TMP_META_WITH_OUT" 2>"$TMP_META_ERR"; then
        meta_with_count=$(wc -l < "$TMP_META_WITH_OUT")
        # A non-empty result has at least 2 lines: 1 header + 1 data row
        if [[ "$meta_with_count" -gt 1 ]]; then
            mv "$TMP_META_WITH_OUT" "$META_OUT_FILE"
            echo "==> Used query WITH features. Metadata (${meta_with_count} lines incl. header) written to: ${META_OUT_FILE}"
        else
            echo "==> WITH features query returned no data rows – falling back to query WITHOUT features ..."
            rm -f "$TMP_META_WITH_OUT"
            run_fallback_metadata_query
        fi
    else
        echo "WARNING: Metadata query (with features) failed – falling back to query WITHOUT features ..." >&2
        cat "$TMP_META_ERR" >&2
        rm -f "$TMP_META_WITH_OUT"
        run_fallback_metadata_query
    fi
    rm -f "$TMP_META_WITH_SQL" "$TMP_META_ERR"
fi
