#!/usr/bin/env bash
#SBATCH --job-name=miranda
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --array=1-111   # Update to match the miRNA count reported by 8.1


###############################################################################
# 8.2_miranda_array.sh - Run miRanda target prediction as a SLURM array job
#
# Description:
#   Phase VIII.2 (SLURM array job): For each miRNA (one per array task), run
#   miRanda against the 3' UTR FASTA database prepared by
#   8.1_mirna_target_preparation.sh and save the raw output file.
#
#   The array task ID maps to a line in mirna_list_indexed.txt:
#     task 1 → line 1 → miRNA name 1
#     task N → line N → miRNA name N
#
#   The total number of array tasks must equal the miRNA count reported at the
#   end of 8.1_mirna_target_preparation.sh.  Update the #SBATCH --array line
#   above accordingly before submitting.
#
# Usage:
#   sbatch scripts/8.2_miranda_array.sh
#
# Prerequisites:
#   - Run 8.1_mirna_target_preparation.sh first (produces temp/ input files).
#
# Outputs (in ${OUTDIR}/miranda/raw/):
#   <miRNA_name>.txt  – Raw miRanda output for each miRNA
#
###############################################################################

printf "\n\nStarted Job: 8.2-miranda_array (task ${SLURM_ARRAY_TASK_ID})\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 8.1_mirna_target_preparation.sh)
OUTDIR="/storage/users/ahinsu/DD-microRNA/downstream/targets/miRNA"

# miRanda alignment score threshold (0–255; higher = more stringent)
MIRANDA_SCORE=140

# miRanda free-energy threshold in kcal/mol (more negative = more stringent)
MIRANDA_ENERGY=-20

# Conda environment name
CONDA_ENV="miranda"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

TEMP_DIR="${OUTDIR}/temp"
RAW_DIR="${OUTDIR}/miranda/raw"
mkdir -p "${RAW_DIR}"

# ============================================================================
# SLURM Array Job – resolve miRNA for this task
# ============================================================================

MIRNA_INDEXED_LIST="${TEMP_DIR}/mirna_list_indexed.txt"
MIRNA_FASTA="${TEMP_DIR}/mirna_selected.fa"
UTR_FASTA="${TEMP_DIR}/three_prime_utr.fa"

if [ ! -f "${MIRNA_INDEXED_LIST}" ]; then
    echo "ERROR: Indexed miRNA list not found: ${MIRNA_INDEXED_LIST}"
    echo "       Run 8.1_mirna_target_preparation.sh first."
    exit 1
fi

for f in "${MIRNA_FASTA}" "${UTR_FASTA}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: Required input file not found: ${f}"
        echo "       Run 8.1_mirna_target_preparation.sh first."
        exit 1
    fi
done

MIRNA_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MIRNA_INDEXED_LIST}")

if [ -z "${MIRNA_NAME}" ]; then
    echo "ERROR: Could not determine miRNA name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Array task   : ${SLURM_ARRAY_TASK_ID}"
echo "miRNA name   : ${MIRNA_NAME}"
echo "Score cutoff : ${MIRANDA_SCORE}"
echo "Energy cutoff: ${MIRANDA_ENERGY} kcal/mol"

# ============================================================================
# Step 1: Extract single-miRNA FASTA
# ============================================================================

echo ""
echo "Step 1: Extracting FASTA for ${MIRNA_NAME}..."

TEMP_MIRNA_FA="${RAW_DIR}/.tmp_${SLURM_ARRAY_TASK_ID}.fa"

seqkit grep -n -p "${MIRNA_NAME}" "${MIRNA_FASTA}" > "${TEMP_MIRNA_FA}"

if [ ! -s "${TEMP_MIRNA_FA}" ]; then
    echo "ERROR: No sequence found for miRNA '${MIRNA_NAME}' in ${MIRNA_FASTA}"
    rm -f "${TEMP_MIRNA_FA}"
    exit 1
fi

echo "  Extracted 1 sequence → ${TEMP_MIRNA_FA}"

# ============================================================================
# Step 2: Run miRanda
# ============================================================================

echo ""
echo "Step 2: Running miRanda for ${MIRNA_NAME}..."

OUTPUT_FILE="${RAW_DIR}/${MIRNA_NAME}.txt"

miranda \
    "${TEMP_MIRNA_FA}" \
    "${UTR_FASTA}" -strict \
    -sc "${MIRANDA_SCORE}" \
    -en "${MIRANDA_ENERGY}" \
    -out "${OUTPUT_FILE}"

rm -f "${TEMP_MIRNA_FA}"

if [ -s "${OUTPUT_FILE}" ]; then
    HIT_COUNT=$(grep -c '^>>' "${OUTPUT_FILE}" 2>/dev/null || echo 0)
    echo "  miRanda complete: ${HIT_COUNT} hits → ${OUTPUT_FILE}"
else
    echo "  miRanda complete: no hits found for ${MIRNA_NAME}"
fi

printf "\n\nCompleted: 8.2-miranda_array (task ${SLURM_ARRAY_TASK_ID}, miRNA: ${MIRNA_NAME})\n\n"
pwd
date
