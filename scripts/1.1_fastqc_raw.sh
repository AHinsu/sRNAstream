#!/usr/bin/env bash
#SBATCH --job-name=fastqc_raw
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=8G
#SBATCH --array=1-72%5   # 72 is number of samples


###############################################################################
# 1.1_fastqc_raw.sh - FastQC quality control on raw fastq files
# 
# Description:
#   Performs quality control analysis on raw sequencing data using FastQC
#
# Usage: 
#   sbatch scripts/1.1_fastqc_raw.sh
#
# Outputs:
#   - FastQC HTML reports in ${OUTDIR}/fastqc/raw/
#   - FastQC zip files with detailed metrics
###############################################################################

printf "\n\nStarted Job: 1.1-fastqc_raw\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Directory containing raw FASTQ files
INPUT_RAW_DIR="/storage/users/ahinsu/DD-microRNA/RawData"

# Output directory
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Conda environment name (created from environment.yml)
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
mkdir -p "${OUTDIR}/fastqc/raw"

# ============================================================================
# SLURM Array Job Processing
# ============================================================================

# Get sample row (skip header) for this array task
sample_line=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n+1{print $0}' "$MANIFEST")
SAMPLE_NAME=$(echo "$sample_line" | cut -f1)
FASTQ_FILE=$(echo "$sample_line" | cut -f2)


if [ ! -f "$FASTQ_FILE" ]; then
    echo "Error: FASTQ file not found: $FASTQ_FILE"
    exit 1
fi

echo "Processing sample: $SAMPLE_NAME"
echo "Input file: $FASTQ_FILE"

# Run FastQC
fastqc \
    -t ${SLURM_CPUS_PER_TASK} --noextract \
    --outdir ${OUTDIR}/fastqc/raw \
    $FASTQ_FILE

echo "FastQC completed for $SAMPLE_NAME"

printf "\n\nCompleted\n\n"
pwd
date
