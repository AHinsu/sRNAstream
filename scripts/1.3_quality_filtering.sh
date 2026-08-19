#!/usr/bin/env bash
#SBATCH --job-name=quality_filtering
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-72


###############################################################################
# 1.3_quality_filtering.sh - Quality filtering with fastp
# 
# Description:
#   Performs quality filtering on cutadapt-trimmed reads using fastp
#   - Filters reads with average quality < 30
#   - Trims polyX sequences
#   - Filters by length
#
# Usage: 
#   sbatch scripts/1.3_quality_filtering.sh
#
# Outputs:
#   - Quality-filtered FASTQ files in ${OUTDIR}/fastp/
#   - HTML and JSON reports with filtering statistics
###############################################################################

printf "\n\nStarted Job: 1.3-quality_filtering\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Directory containing cutadapt-trimmed files
INPUT_DIR="/storage/users/ahinsu/DD-microRNA/results/cutadapt"

# Output directory
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Conda environment name
CONDA_ENV="smallRNA"

# Quality parameters
QUALITY=30            # Average quality threshold
MIN_LENGTH=17         # Minimum read length (typical for miRNA)
MAX_LENGTH=30        # Maximum read length
N_BASE_LIMIT=1        # Maximum number of N bases allowed

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
mkdir -p "${OUTDIR}/fastp"

# ============================================================================
# SLURM Array Job Processing
# ============================================================================

# Get sample row (skip header) for this array task
sample_line=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n+1{print $0}' "$MANIFEST")
SAMPLE_NAME=$(echo "$sample_line" | cut -f1)

FASTQ_FILE="${OUTDIR}/cutadapt/${SAMPLE_NAME}_cutadapt.fq.gz"

if [ ! -f "$FASTQ_FILE" ]; then
    echo "Error: FASTQ file not found: $FASTQ_FILE"
    exit 1
fi

echo "Processing sample: $SAMPLE_NAME"
echo "Input file: $FASTQ_FILE"

# Run fastp
fastp \
    -i $FASTQ_FILE \
    -o ${OUTDIR}/fastp/${SAMPLE_NAME}_trimmed.fq.gz \
    --trim_poly_x \
    --length_required ${MIN_LENGTH} \
    --length_limit ${MAX_LENGTH} \
    --qualified_quality_phred ${QUALITY} \
    --average_qual ${QUALITY} \
    --n_base_limit ${N_BASE_LIMIT} \
    --thread ${SLURM_CPUS_PER_TASK} \
    --html ${OUTDIR}/fastp/${SAMPLE_NAME}_fastp.html \
    --json ${OUTDIR}/fastp/${SAMPLE_NAME}_fastp.json

echo "Fastp completed for $SAMPLE_NAME"

printf "\n\nCompleted\n\n"
pwd
date
