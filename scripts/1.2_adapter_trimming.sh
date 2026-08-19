#!/usr/bin/env bash
#SBATCH --job-name=adapter_trimming
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-72


###############################################################################
# 1.2_adapter_trimming.sh - Trim adapters using cutadapt
# 
# Description:
#   Removes NEBNext Small RNA library prep kit adapters using cutadapt
#   - RNA 5' Adapter (RA5): GTTCAGAGTTCTACAGTCCGACGATC
#   - RNA 3' Adapter (RA3): AGATCGGAAGAGCACACGTCT
#   - Minimum length: 17nt, Maximum length: 30nt
#   - Quality trimming threshold: Q20 or Q30
#
# Usage: 
#   sbatch scripts/1.2_adapter_trimming.sh
#
# Outputs:
#   - Trimmed FASTQ files in ${OUTDIR}/cutadapt/
#   - Cutadapt log files with trimming statistics
###############################################################################

printf "\n\nStarted Job: 1.2-adapter_trimming\n\n"
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

# Conda environment name
CONDA_ENV="smallRNA"

# Adapter sequences (NEBNext Small RNA library prep kit)
ADAPTER_5="GTTCAGAGTTCTACAGTCCGACGATC"  # RNA 5' Adapter (RA5)
ADAPTER_3="AGATCGGAAGAGCACACGTCT"       # RNA 3' Adapter (RA3)

# Minimum length after trimming (typical for miRNA: 17-18)
MIN_LENGTH=17

# Maximum length after trimming (typical for miRNA: up to 30nt)
MAX_LENGTH=30

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
mkdir -p "${OUTDIR}/cutadapt"

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

# Run cutadapt
cutadapt \
    -g ${ADAPTER_5} \
    -a ${ADAPTER_3} \
    -m ${MIN_LENGTH} \
    -M ${MAX_LENGTH} \
    -j ${SLURM_CPUS_PER_TASK} \
    -o ${OUTDIR}/cutadapt/${SAMPLE_NAME}_cutadapt.fq.gz \
    --json ${OUTDIR}/cutadapt/${SAMPLE_NAME}_cutadapt.json \
    $FASTQ_FILE \
    > ${OUTDIR}/cutadapt/${SAMPLE_NAME}_cutadapt.log

echo "Cutadapt completed for $SAMPLE_NAME"

printf "\n\nCompleted\n\n"
pwd
date
