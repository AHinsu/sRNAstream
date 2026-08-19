#!/usr/bin/env bash
#SBATCH --job-name=snrna_snorna_ref
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G



###############################################################################
# 6.1_snrna_snorna_reference.sh - Build snRNA/snoRNA/piRNA quantification references
#
# Description:
#   Phase VI.1: For each of the three RNA types (snRNA, snoRNA, piRNA),
#   deduplicate the FASTA sequences that were assembled by
#   1.5a_contamination_db_preparation.sh and stored in the contamination
#   reference directory, then build a Bowtie1 index.  These per-type indices
#   are used by 6.2_snrna_snorna_quantification.sh to quantify reads from
#   each sample independently.
#
#   Input FASTA files (written by 1.5a build_contam_index):
#     ${REF_DIR}/contamination/snrna_index_sequences.fa
#     ${REF_DIR}/contamination/snorna_index_sequences.fa
#     ${REF_DIR}/contamination/pirna_index_sequences.fa
#
#   Deduplication uses seqkit rmdup -s -i (by sequence, case-insensitive),
#   matching the approach used in 2.2_database_unification.sh.
#   U→T conversion is applied first with seqkit seq --rna2dna.
#
# Usage:
#   sbatch scripts/6.1_snrna_snorna_reference.sh
#
# Prerequisites:
#   - Run 1.5a_contamination_db_preparation.sh first (produces the FASTA files).
#
# Outputs (in ${REF_DIR}/snrna_quant_ref/):
#   snrna_dedup.fa              – Deduplicated snRNA sequences (DNA)
#   snorna_dedup.fa             – Deduplicated snoRNA sequences (DNA)
#   pirna_dedup.fa              – Deduplicated piRNA sequences (DNA)
#   snrna_dup_info.txt          – seqkit rmdup duplicate info (snRNA)
#   snorna_dup_info.txt         – seqkit rmdup duplicate info (snoRNA)
#   pirna_dup_info.txt          – seqkit rmdup duplicate info (piRNA)
#   snrna_index.*               – Bowtie1 index for snRNA
#   snorna_index.*              – Bowtie1 index for snoRNA
#   pirna_index.*               – Bowtie1 index for piRNA
#   deduplication_report.txt    – Per-type duplicate statistics
#   reference_build_summary.txt – Overall summary
###############################################################################

printf "\n\nStarted Job: 6.1-snrna_snorna_reference\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Reference directory (must match REF_DIR used in 1.5a_contamination_db_preparation.sh)
REF_DIR="/storage/users/ahinsu/DD-microRNA/reference"

# Conda environment name
CONDA_ENV="smallRNA"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

CONTAM_DIR="${REF_DIR}/contamination"
QUANT_REF_DIR="${REF_DIR}/snrna_quant_ref"
mkdir -p "${QUANT_REF_DIR}"

# ============================================================================
# Step 1: Locate source FASTA files from contamination directory
# ============================================================================

echo "Step 1: Locating contamination FASTA files..."

SNRNA_FA="${CONTAM_DIR}/snrna_index_sequences.fa"
SNORNA_FA="${CONTAM_DIR}/snorna_index_sequences.fa"
PIRNA_FA="${CONTAM_DIR}/pirna_index_sequences.fa"

echo "  snRNA  FASTA : ${SNRNA_FA}"
echo "  snoRNA FASTA : ${SNORNA_FA}"
echo "  piRNA  FASTA : ${PIRNA_FA}"

# ============================================================================
# Step 2: Deduplicate sequences for each RNA type with seqkit rmdup
# ============================================================================

echo ""
echo "Step 2: Deduplicating sequences (seqkit rmdup -s -i) for each RNA type..."

DEDUP_RPT="${QUANT_REF_DIR}/deduplication_report.txt"

{
echo "================================================================="
echo "Deduplication Report – snRNA/snoRNA/piRNA quantification references"
echo "================================================================="
echo "Tool: seqkit rmdup -s -i (by sequence, case-insensitive)"
echo ""
printf "%-12s  %8s  %8s  %12s\n" "RNA type" "Input" "Unique" "Duplicates"
echo "-----------------------------------------------------------------"
} > "${DEDUP_RPT}"

# Helper: dedup one RNA type
#   $1 – short prefix (snrna / snorna / pirna)
#   $2 – display label (snRNA / snoRNA / piRNA)
#   $3 – source FASTA path
dedup_type() {
    local prefix="$1"
    local label="$2"
    local src_fa="$3"
    local out_fa="${QUANT_REF_DIR}/${prefix}_dedup.fa"
    local dup_info="${QUANT_REF_DIR}/${prefix}_dup_info.txt"

    if [ ! -s "${src_fa}" ]; then
        echo "  WARNING: ${label} FASTA not found or empty: ${src_fa}"
        echo "           Run 1.5a_contamination_db_preparation.sh first."
        touch "${out_fa}"
        printf "%-12s  %8s  %8s  %12s\n" "${label}" "N/A" "N/A" "N/A" >> "${DEDUP_RPT}"
        return
    fi

    local n_in
    n_in=$(grep -c '^>' "${src_fa}" 2>/dev/null || echo 0)

    # Convert U→T then deduplicate
    seqkit seq --rna2dna "${src_fa}" \
        | seqkit rmdup -s -i \
            -D "${dup_info}" \
            -o "${out_fa}"

    local n_out
    n_out=$(grep -c '^>' "${out_fa}" 2>/dev/null || echo 0)
    local n_dup=$(( n_in - n_out ))

    echo "  ${label}: ${n_in} → ${n_out} unique (${n_dup} duplicates removed)"
    echo "  Duplicate info: ${dup_info}"
    printf "%-12s  %8d  %8d  %12d\n" "${label}" "${n_in}" "${n_out}" "${n_dup}" >> "${DEDUP_RPT}"
}

dedup_type "snrna"  "snRNA"  "${SNRNA_FA}"
dedup_type "snorna" "snoRNA" "${SNORNA_FA}"
dedup_type "pirna"  "piRNA"  "${PIRNA_FA}"

echo "-----------------------------------------------------------------" >> "${DEDUP_RPT}"
echo "" >> "${DEDUP_RPT}"

echo ""
echo "Deduplication report:"
cat "${DEDUP_RPT}"

# ============================================================================
# Step 3: Build Bowtie1 index for each RNA type
# ============================================================================

echo ""
echo "Step 3: Building Bowtie1 indices..."

build_index() {
    local prefix="$1"
    local label="$2"
    local fa="${QUANT_REF_DIR}/${prefix}_dedup.fa"
    local idx="${QUANT_REF_DIR}/${prefix}_index"

    if [ ! -s "${fa}" ]; then
        echo "  Skipping ${label} index (FASTA empty or missing: ${fa})"
        return
    fi

    echo ""
    echo "  Building ${label} index: ${idx}"
    bowtie-build \
        --threads ${SLURM_CPUS_PER_TASK:-4} \
        "${fa}" \
        "${idx}"
    echo "  ${label} index built: ${idx}"
}

build_index "snrna"  "snRNA"
build_index "snorna" "snoRNA"
build_index "pirna"  "piRNA"

# ============================================================================
# Step 4: Summary
# ============================================================================

N_SNRNA=$(grep  -c '^>' "${QUANT_REF_DIR}/snrna_dedup.fa"  2>/dev/null || echo 0)
N_SNORNA=$(grep -c '^>' "${QUANT_REF_DIR}/snorna_dedup.fa" 2>/dev/null || echo 0)
N_PIRNA=$(grep  -c '^>' "${QUANT_REF_DIR}/pirna_dedup.fa"  2>/dev/null || echo 0)

cat > "${QUANT_REF_DIR}/reference_build_summary.txt" << SUMMARY_EOF
========================================================================
Phase VI.1 – snRNA/snoRNA/piRNA Quantification References
========================================================================
snRNA  unique sequences : ${N_SNRNA}
snoRNA unique sequences : ${N_SNORNA}
piRNA  unique sequences : ${N_PIRNA}

FASTA files (in ${QUANT_REF_DIR}/):
  snrna_dedup.fa   → snrna_index.*
  snorna_dedup.fa  → snorna_index.*
  pirna_dedup.fa   → pirna_index.*

Contamination FASTA source directory: ${CONTAM_DIR}

Next step:
  Run 6.2_snrna_snorna_quantification.sh (SLURM array) to map per-sample
  reads to these references and compute per-sequence gene counts.
========================================================================
SUMMARY_EOF

cat "${QUANT_REF_DIR}/reference_build_summary.txt"

printf "\n\nCompleted: 6.1-snrna_snorna_reference\n\n"
pwd
date



