#!/usr/bin/env bash
#SBATCH --job-name=mirna_reference
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G


###############################################################################
# 5.1_mirna_reference.sh - Build unified miRNA quantification reference
#
# Description:
#   Phase V.1: Merge the known mature miRNA sequences (same FASTA supplied to
#   miRDeep2 in 4.1) with novel mature miRNA candidates discovered by miRDeep2,
#   deduplicate the combined set, and build a Bowtie1 index for per-sample
#   quantification in 5.2_mirna_quantification.sh.
#
#   Input sources (merged in this order; first occurrence of each unique
#   sequence wins):
#     1. Known mature miRNA FASTA (KNOWN_MATURE_FA) – the same unified
#        reference supplied to miRDeep2 in 4.1 (e.g. combined_mature_dedup.fasta
#        from 2.2_database_unification.sh).
#     2. Novel mature miRNA FASTA produced by get_mirdeep2_precursors.pl in
#        4.1 at the score threshold defined by NOVEL_SCORE_THRESHOLD.
#        File: ${OUTDIR}/mirdeep2/sequences${NOVEL_SCORE_THRESHOLD}/novel_mature_*.fa
#
#   USER CHOICE: set NOVEL_SCORE_THRESHOLD to one of the tier values used in
#   4.1 (SCORE_LOW / SCORE_MED / SCORE_HIGH, e.g. 1, 5, or 10).  Higher
#   thresholds include fewer but more reliable novel candidates.
#
#   Deduplication uses seqkit rmdup -s -i (by sequence, case-insensitive),
#   matching the approach used in 2.2_database_unification.sh.
#   U→T conversion is applied first with seqkit seq --rna2dna.
#
# Usage:
#   sbatch scripts/5.1_mirna_reference.sh
#
# Prerequisites:
#   - Run 4.1_mirdeep2_analysis.sh first (provides novel_mature_*.fa).
#   - KNOWN_MATURE_FA must exist (output of 2.2_database_unification.sh).
#
# Outputs (in ${REF_DIR}/mirna_quant_ref/):
#   unified_mature_all.fa          – Known + novel combined (pre-dedup, DNA)
#   unified_mature.fa              – Deduplicated mature miRNA reference (DNA)
#   unified_mature_index.*         – Bowtie1 index files
#   dup_info.txt                   – seqkit rmdup duplicate info
#   deduplication_report.txt       – Summary of deduplication statistics
#   reference_build_summary.txt    – Overall summary
###############################################################################

printf "\n\nStarted Job: 5.1-mirna_reference\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 4.1_mirdeep2_analysis.sh)
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Reference directory
REF_DIR="/storage/users/ahinsu/DD-microRNA/reference"

# Known mature miRNA FASTA – the SAME file used as KNOWN_MATURE_FA in 4.1.
# This is the unified/deduplicated reference from 2.2_database_unification.sh.
KNOWN_MATURE_FA="${REF_DIR}/db_unified/combined_mature_dedup.fasta"

# Score threshold to use when selecting novel miRNA candidates from miRDeep2.
# Must match one of the tier values used in 4.1 (SCORE_LOW / SCORE_MED / SCORE_HIGH).
# Higher values → fewer but more reliable novel candidates.
#   1  = low confidence  (predictions with score >= 1)
#   5  = medium confidence (recommended; predictions with score >= 5)
#   10 = high confidence only (predictions with score >= 10)
NOVEL_SCORE_THRESHOLD=10

# Conda environment name
CONDA_ENV="smallRNA"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

MIRDEEP2_DIR="${OUTDIR}/mirdeep2"
QUANT_REF_DIR="${REF_DIR}/mirna_quant_ref"
mkdir -p "${QUANT_REF_DIR}"

# ============================================================================
# Step 1: Locate source FASTA files
# ============================================================================

echo "Step 1: Locating source FASTA files..."

# Known mature miRNA
if [ ! -s "${KNOWN_MATURE_FA}" ]; then
    echo "ERROR: Known mature miRNA FASTA not found or empty: ${KNOWN_MATURE_FA}"
    echo "       Run 2.2_database_unification.sh first, or set KNOWN_MATURE_FA to a valid path."
    exit 1
fi
echo "  Known mature miRNAs : ${KNOWN_MATURE_FA}"

# Novel mature miRNA from miRDeep2 at the chosen score threshold
NOVEL_SEQ_DIR="${MIRDEEP2_DIR}/sequences${NOVEL_SCORE_THRESHOLD}"
if [ ! -d "${NOVEL_SEQ_DIR}" ]; then
    echo "ERROR: miRDeep2 sequences directory not found: ${NOVEL_SEQ_DIR}"
    echo "       Run 4.1_mirdeep2_analysis.sh first with SCORE_LOW/MED/HIGH"
    echo "       that covers NOVEL_SCORE_THRESHOLD=${NOVEL_SCORE_THRESHOLD}."
    exit 1
fi

# get_mirdeep2_precursors.pl writes a file named novel_mature_<timestamp>.fa
NOVEL_MATURE_FA=$(find "${NOVEL_SEQ_DIR}" -maxdepth 1 -name 'novel_mature_*.fa' \
    -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -1 | cut -d' ' -f2-)

if [ -z "${NOVEL_MATURE_FA}" ] || [ ! -s "${NOVEL_MATURE_FA}" ]; then
    echo "WARNING: No novel_mature_*.fa found in ${NOVEL_SEQ_DIR}."
    echo "         Only known mature miRNAs will be included in the reference."
    NOVEL_MATURE_FA=""
else
    echo "  Novel mature miRNAs : ${NOVEL_MATURE_FA}"
    echo "  (score threshold    : >= ${NOVEL_SCORE_THRESHOLD})"
# Normalize novel FASTA headers in-place:
    # >2_1894_x75634_coord:2:..._score:38561.9
    # -> >2_1894_coord:2:...
    TMP_NOVEL_FA="${NOVEL_MATURE_FA}.tmp"
    awk '
    /^>/ {
        h = substr($0,2)  # drop leading >
        # keep: <chr>_<id> and coord:...
        if (match(h, /^([0-9]+_[0-9]+).*(_coord:[^ ]+)/, a)) {
            # remove trailing _score:... from coord tail
            sub(/_score:.*/, "", a[2])
            print ">" a[1] a[2]
        } else {
            print $0
        }
        next
    }
    { print }
    ' "${NOVEL_MATURE_FA}" > "${TMP_NOVEL_FA}" && mv "${TMP_NOVEL_FA}" "${NOVEL_MATURE_FA}"
fi

# ============================================================================
# Step 2: Merge sources and deduplicate with seqkit rmdup
# ============================================================================

echo ""
echo "Step 2: Merging and deduplicating sequences (seqkit rmdup -s -i)..."

COMBINED_ALL="${QUANT_REF_DIR}/unified_mature_all.fa"
UNIFIED_FA="${QUANT_REF_DIR}/unified_mature.fa"
DUP_INFO="${QUANT_REF_DIR}/dup_info.txt"
DEDUP_RPT="${QUANT_REF_DIR}/deduplication_report.txt"

# Start fresh
> "${COMBINED_ALL}"

# Known mature miRNAs: convert U→T (seqkit seq --rna2dna) and tag header with source.
# seqkit rmdup keeps the first occurrence, so known sequences win when a novel
# sequence is identical to a known one.
N_KNOWN_IN=$(grep -c '^>' "${KNOWN_MATURE_FA}" 2>/dev/null || echo 0)
seqkit seq --rna2dna "${KNOWN_MATURE_FA}" \
    | awk '/^>/{print $0 "|src=known_mature"; next}{print}' \
    >> "${COMBINED_ALL}"
echo "  Known mature sequences   : ${N_KNOWN_IN}"

# Novel mature miRNAs (only if a file was found)
N_NOVEL_IN=0
if [ -n "${NOVEL_MATURE_FA}" ]; then
    N_NOVEL_IN=$(grep -c '^>' "${NOVEL_MATURE_FA}" 2>/dev/null || echo 0)
    seqkit seq --rna2dna "${NOVEL_MATURE_FA}" \
        | awk '/^>/{print $0 "|src=novel_mature"; next}{print}' \
        >> "${COMBINED_ALL}"
    echo "  Novel mature sequences   : ${N_NOVEL_IN}"
fi

N_COMBINED=$(grep -c '^>' "${COMBINED_ALL}" 2>/dev/null || echo 0)
echo "  Combined total           : ${N_COMBINED}"

echo ""
echo "Running seqkit rmdup..."
seqkit rmdup -s -i \
    -D "${DUP_INFO}" \
    -o "${UNIFIED_FA}" \
    "${COMBINED_ALL}"

N_UNIQUE=$(grep -c '^>' "${UNIFIED_FA}" 2>/dev/null || echo 0)
N_DUPS=$(( N_COMBINED - N_UNIQUE ))
echo "  ${N_COMBINED} → ${N_UNIQUE} unique sequences (${N_DUPS} duplicates removed)"
echo "  Duplicate info: ${DUP_INFO}"

cat > "${DEDUP_RPT}" << DEDUP_EOF
=================================================================
Deduplication Report – unified miRNA reference
=================================================================
Tool   : seqkit rmdup -s -i (by sequence, case-insensitive)
Source order (first occurrence kept):
  1. Known mature miRNA  (|src=known_mature)
  2. Novel mature miRNA  (|src=novel_mature)

Input sequences:
  Known mature        : ${N_KNOWN_IN}
  Novel mature        : ${N_NOVEL_IN}
  Combined total      : ${N_COMBINED}

After seqkit rmdup:
  Unique sequences    : ${N_UNIQUE}
  Duplicates removed  : ${N_DUPS}
  Duplicate info file : ${DUP_INFO}
=================================================================
DEDUP_EOF

echo ""
cat "${DEDUP_RPT}"
echo "Unified FASTA: ${UNIFIED_FA}"

# ============================================================================
# Step 3: Build Bowtie1 index
# ============================================================================

echo ""
echo "Step 3: Building Bowtie1 index for unified mature miRNA reference..."

UNIFIED_FA="${QUANT_REF_DIR}/unified_mature.fa"
INDEX_PREFIX="${QUANT_REF_DIR}/unified_mature_index"

if [ ! -s "$UNIFIED_FA" ]; then
    echo "ERROR: Unified FASTA is empty or missing. Check Step 2 output."
    exit 1
fi

bowtie-build \
    --threads ${SLURM_CPUS_PER_TASK:-4} \
    "$UNIFIED_FA" \
    "$INDEX_PREFIX"

echo "Bowtie1 index: ${INDEX_PREFIX}"

# ============================================================================
# Step 4: Summary
# ============================================================================

N_SEQS=$(grep -c '^>' "$UNIFIED_FA" 2>/dev/null || echo 0)
N_KNOWN=$(grep -c '|src=known_mature' "$UNIFIED_FA" 2>/dev/null || echo 0)
N_NOVEL=$(grep -c '|src=novel_mature' "$UNIFIED_FA" 2>/dev/null || echo 0)

cat > "${QUANT_REF_DIR}/reference_build_summary.txt" << SUMMARY_EOF
========================================================================
Phase V.1 – Unified miRNA Quantification Reference
========================================================================
Known mature miRNA sequences  : ${N_KNOWN}
Novel mature miRNA sequences  : ${N_NOVEL}  (score threshold >= ${NOVEL_SCORE_THRESHOLD})
Total (unique) sequences      : ${N_SEQS}
FASTA file                    : ${UNIFIED_FA}
Bowtie1 index prefix          : ${INDEX_PREFIX}

Next step:
  Run 5.2_mirna_quantification.sh (SLURM array) to map per-sample reads
  to this reference and compute gene counts.
========================================================================
SUMMARY_EOF

cat "${QUANT_REF_DIR}/reference_build_summary.txt"

printf "\n\nCompleted: 5.1-mirna_reference\n\n"
pwd
date
