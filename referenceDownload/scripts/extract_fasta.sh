#!/usr/bin/env bash
# =============================================================================
# extract_fasta.sh
# -----------------------------------------------------------------------------
# Extract sequences from an RNAcentral FASTA file for one or more RNA types
# using pre-generated ID list files (produced by run_query_linux.sh).
#
# Requirements:
#   - seqkit  (https://bioinf.shenwei.me/seqkit/)
#     Install: conda install -c bioconda seqkit
#          OR: sudo apt install seqkit
#          OR: download binary from https://github.com/shenwei356/seqkit/releases
#   - The RNAcentral species-specific FASTA file (downloaded from FTP)
#     See README.md for download instructions.
#
# Usage:
#   bash extract_fasta.sh --fasta rnacentral_species_specific_ids.fasta.gz \
#                         --ids   9913_miRNA.ids.txt \
#                         --out   bos_taurus_9913_miRNA.fasta
#
#   # Or process multiple types in a loop (see README for examples):
#   for TYPE in miRNA tRNA rRNA snRNA snoRNA; do
#     bash extract_fasta.sh \
#       --fasta rnacentral_species_specific_ids.fasta.gz \
#       --ids   9913_${TYPE}.ids.txt \
#       --out   bos_taurus_9913_${TYPE}.fasta
#   done
#
# Options:
#   --fasta FILE   Path to rnacentral_species_specific_ids.fasta.gz (required)
#   --ids   FILE   Path to the IDs text file (required)
#   --out   FILE   Output FASTA file path (default: output.fasta)
#   --help         Show this help
# =============================================================================

set -euo pipefail

FASTA_FILE=""
IDS_FILE=""
OUT_FILE="output.fasta"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta) FASTA_FILE="$2"; shift 2 ;;
        --ids)   IDS_FILE="$2";   shift 2 ;;
        --out)   OUT_FILE="$2";   shift 2 ;;
        --help|-h)
            sed -n '1,/^# ====/p' "$0" | grep '^#' | sed 's/^# \?//'
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
if [[ -z "$FASTA_FILE" ]]; then
    echo "ERROR: --fasta is required." >&2
    exit 1
fi
if [[ -z "$IDS_FILE" ]]; then
    echo "ERROR: --ids is required." >&2
    exit 1
fi
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "ERROR: FASTA file not found: $FASTA_FILE" >&2
    exit 1
fi
if [[ ! -f "$IDS_FILE" ]]; then
    echo "ERROR: IDs file not found: $IDS_FILE" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Check seqkit
# ---------------------------------------------------------------------------
if ! command -v seqkit &>/dev/null; then
    echo "ERROR: seqkit is not installed or not on PATH." >&2
    echo "       Install with: conda install -c bioconda seqkit" >&2
    echo "       Or download:  https://github.com/shenwei356/seqkit/releases" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Count IDs
# ---------------------------------------------------------------------------
ID_COUNT=$(grep -c . "$IDS_FILE" || true)
echo "==> Extracting ${ID_COUNT} sequences from ${FASTA_FILE} ..."
echo "    IDs file : $IDS_FILE"
echo "    Output   : $OUT_FILE"

# ---------------------------------------------------------------------------
# Extract
# ---------------------------------------------------------------------------
seqkit grep -f "$IDS_FILE" "$FASTA_FILE" > "$OUT_FILE"

OUT_COUNT=$(grep -c '^>' "$OUT_FILE" || true)
echo "==> Done. ${OUT_COUNT} sequences written to: ${OUT_FILE}"

if [[ "$ID_COUNT" -ne "$OUT_COUNT" ]]; then
    echo "WARNING: Expected ${ID_COUNT} sequences but found ${OUT_COUNT}." >&2
    echo "         Some IDs may not be present in the FASTA file." >&2
    echo "         Check that you are using the correct FASTA release." >&2
fi
