#!/usr/bin/env bash
#SBATCH --job-name=fastqc_cleaned
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --array=1-72%5


###############################################################################
# 1.4_fastqc_cleaned.sh - FastQC quality control on trimmed fastq files
# 
# Description:
#   Performs quality control analysis on trimmed sequencing data using FastQC
#
# Usage: 
#   sbatch scripts/1.4_fastqc_cleaned.sh
#
# Outputs:
#   - FastQC HTML reports in ${OUTDIR}/fastqc/cleaned/
#   - FastQC zip files with detailed metrics
###############################################################################

printf "\n\nStarted Job: 1.4-fastqc_cleaned\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================
# Output directory
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Conda environment name
CONDA_ENV="smallRNA"

# The manifest file containing sample information
MANIFEST="/storage/users/ahinsu/DD-microRNA/sample_manifest.txt"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH


# Create output directory
mkdir -p "${OUTDIR}/fastqc/cleaned"

# ============================================================================
# SLURM Array Job Processing
# ============================================================================

# Get sample row (skip header) for this array task
sample_line=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n+1{print $0}' "$MANIFEST")
SAMPLE_NAME=$(echo "$sample_line" | cut -f1)

FASTQ_FILE="${OUTDIR}/fastp/${SAMPLE_NAME}_trimmed.fq.gz"


if [ ! -f "$FASTQ_FILE" ]; then
    echo "Error: FASTQ file not found: $FASTQ_FILE"
    exit 1
fi

echo "Processing sample: $SAMPLE_NAME"
echo "Input file: $FASTQ_FILE"

# Run FastQC
fastqc \
    -t ${SLURM_CPUS_PER_TASK} --noextract \
    --outdir ${OUTDIR}/fastqc/cleaned \
    $FASTQ_FILE

echo "FastQC completed for $SAMPLE_NAME"

printf "\n\nCompleted\n\n"
pwd
date
