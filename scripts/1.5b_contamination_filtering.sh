#!/usr/bin/env bash
#SBATCH --job-name=contam_filtering
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --array=1-72   # Update to match number of samples in your manifest


###############################################################################
# 1.5b_contamination_filtering.sh - Remove rRNA, tRNA, and other contamination
# 
# Description:
#   Part B of contamination filtering - runs per-sample array jobs to remove
#   contaminating RNA species from trimmed reads:
#   - Ribosomal RNA (rRNA) - most abundant contamination
#   - Transfer RNA (tRNA) - common in small RNA preps
#   - Small nuclear RNA (snRNA) - U1, U2, U4, U5, U6
#   - Small nucleolar RNA (snoRNA) - C/D box and H/ACA box
#   - Piwi-interacting RNA (piRNA) - 24-31 nt, primarily germline
#   - Other small RNAs (misc_RNA, Y_RNA, vault_RNA, ribozyme)
#   
#   Uses Bowtie with 1-2 mismatches for comprehensive contamination removal.
#   Retains only unmapped reads for downstream miRNA analysis.
#   Mapped snRNA/snoRNA/piRNA reads are saved separately for downstream
#   gene counts and DESeq2-compatible matrix building (Phase VIII).
#
#   User-provided FASTA files for each RNA type should be extracted from the
#   RNAcentral database using the scripts in the referenceDownload/ folder
#   (see referenceDownload/README.md). Bowtie indexes are built by 1.5a.
#
# Usage: 
#   sbatch scripts/1.5b_contamination_filtering.sh
#
# Requirements:
#   - Run 1.5a_contamination_db_preparation.sh first to build the Bowtie indexes
#   - Run after 1.3_quality_filtering.sh (fastp)
#
# Outputs:
#   - Cleaned FASTQ files in ${OUTDIR}/contamination_filtered/
#   - Per-sample contamination statistics
#   - snRNA/snoRNA reads in ${OUTDIR}/snrna_snorna_reads/ (for Phase VIII)
#   - piRNA reads in ${OUTDIR}/pirna_reads/ (for downstream piRNA quantification)
###############################################################################

printf "\n\nStarted Job: 1.5b-contamination_filtering\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Directory containing trimmed FASTQ files (output of 1.3_quality_filtering.sh)
INPUT_DIR="/storage/users/ahinsu/DD-microRNA/results/fastp"

# Output directory
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Reference directory containing contamination indexes (built by 1.5a)
REF_DIR="/storage/users/ahinsu/DD-microRNA/reference"

# Conda environment name
CONDA_ENV="smallRNA"

# Bowtie alignment parameters for contamination filtering
MISMATCHES=1  # Allow 1-2 mismatches to maximize contamination removal

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
mkdir -p "${OUTDIR}/contamination_filtered"

# ============================================================================
# Check for contamination reference libraries
# ============================================================================

echo "Checking for contamination reference libraries..."

# List of contamination types and their index files
CONTAM_RRNA="${REF_DIR}/contamination/rrna_index"
CONTAM_TRNA="${REF_DIR}/contamination/trna_index"
CONTAM_SNRNA="${REF_DIR}/contamination/snrna_index"
CONTAM_SNORNA="${REF_DIR}/contamination/snorna_index"
CONTAM_PIRNA="${REF_DIR}/contamination/pirna_index"
CONTAM_OTHER="${REF_DIR}/contamination/other_index"

# ============================================================================
# SLURM Array Job Processing
# ============================================================================

# Get sample row (skip header) for this array task
sample_line=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n+1{print $0}' "$MANIFEST")
SAMPLE_NAME=$(echo "$sample_line" | cut -f1)

FASTQ_FILE="${INPUT_DIR}/${SAMPLE_NAME}_trimmed.fq.gz"

if [ ! -f "$FASTQ_FILE" ]; then
    echo "Error: FASTQ file not found: $FASTQ_FILE"
    exit 1
fi

echo "Processing sample: $SAMPLE_NAME"
echo "Input file: $FASTQ_FILE"

# Initialize with input file
CURRENT_FASTQ=$FASTQ_FILE
TEMP_DIR="${OUTDIR}/contamination_filtered/temp_${SAMPLE_NAME}"
mkdir -p $TEMP_DIR

# Statistics tracking
declare -A CONTAM_STATS
TOTAL_READS=$(zcat "$FASTQ_FILE" | wc -l | awk '{print $1/4}')
CONTAM_STATS["total"]=$TOTAL_READS
REMAINING_READS=$TOTAL_READS

echo "Total input reads: $TOTAL_READS"

# ============================================================================
# Sequential Contamination Filtering
# ============================================================================

# Function to filter contamination
#
# Usage: filter_contamination <contam_type> <contam_index> <input_fq> <output_base> [save_mapped_fq]
#   contam_type    : Label for this RNA type (e.g. "snRNA")
#   contam_index   : Path to bowtie index (without .1.ebwt suffix)
#   input_fq       : Input FASTQ (gzipped)
#   output_base    : Prefix for output files
#   save_mapped_fq : (Optional) If provided, save mapped reads as uncompressed FASTQ
#                    at this path (used for Phase 8 snRNA/snoRNA analysis)
filter_contamination() {
    local contam_type=$1
    local contam_index=$2
    local input_fq=$3
    local output_base=$4
    local save_mapped_fq=${5:-""}  # Optional: path to save mapped reads as uncompressed FASTQ

    # Return path is communicated via the global NEXT_FASTQ variable so that
    # callers do not need command substitution (which would run the function in
    # a subshell and prevent CONTAM_STATS / REMAINING_READS from being updated).
    
    if [ -f "${contam_index}.1.ebwt" ]; then
        echo "" >&2
        echo "Filtering ${contam_type}..." >&2
        
        # Align to contamination index and keep unmapped reads
        bowtie \
            -v ${MISMATCHES} \
            -p ${SLURM_CPUS_PER_TASK} \
            -S \
            --un ${output_base}_unmapped.fq \
            ${contam_index} \
            -q <(zcat "$input_fq") \
            > ${output_base}_mapped.sam 2> ${output_base}_bowtie.log
        
        # Compress unmapped reads
        gzip -f ${output_base}_unmapped.fq
        
        # Count mapped (contaminating) reads
        # -F 4 excludes unmapped reads (FLAG=4) that Bowtie writes to SAM alongside the mapped reads
        local mapped_reads=$(samtools view -c -F 4 ${output_base}_mapped.sam 2>/dev/null || echo 0)
        CONTAM_STATS[$contam_type]=$mapped_reads
        
        # Update remaining reads
        REMAINING_READS=$((REMAINING_READS - mapped_reads))
        
        echo "${contam_type} reads filtered: $mapped_reads" >&2
        echo "Remaining reads: $REMAINING_READS" >&2
        
        # Save mapped reads as FASTQ if requested (for Phase 8 snRNA/snoRNA analysis)
        # -F 4 excludes unmapped reads so only truly-aligned reads are saved
        if [ -n "$save_mapped_fq" ] && [ -s "${output_base}_mapped.sam" ]; then
            samtools fastq -F 4 "${output_base}_mapped.sam" 2>/dev/null \
                > "${save_mapped_fq}"
            echo "Saved ${contam_type} reads: ${save_mapped_fq}" >&2
        fi
        
        # Clean up
        rm -f ${output_base}_mapped.sam
        
        NEXT_FASTQ="${output_base}_unmapped.fq.gz"
    else
        echo "Skipping ${contam_type} filtering (index not found)" >&2
        NEXT_FASTQ="$input_fq"
    fi
}

# Create output directory for snRNA/snoRNA reads (used by Phase 8 analysis)
mkdir -p "${OUTDIR}/snrna_snorna_reads"

# Create output directory for piRNA reads (for downstream piRNA quantification)
mkdir -p "${OUTDIR}/pirna_reads"

# Global variable used by filter_contamination to return the output path without
# relying on command substitution (which would run the function in a subshell
# and prevent CONTAM_STATS / REMAINING_READS from being updated).
NEXT_FASTQ=""

# Filter rRNA
filter_contamination "rRNA" "$CONTAM_RRNA" "$CURRENT_FASTQ" "${TEMP_DIR}/01_rrna"
CURRENT_FASTQ="$NEXT_FASTQ"

# Filter tRNA
if [ -f "$CURRENT_FASTQ" ]; then
    filter_contamination "tRNA" "$CONTAM_TRNA" "$CURRENT_FASTQ" "${TEMP_DIR}/02_trna"
    CURRENT_FASTQ="$NEXT_FASTQ"
fi

# Filter snRNA (save mapped reads for Phase 8 snRNA/snoRNA quantification)
if [ -f "$CURRENT_FASTQ" ]; then
    filter_contamination "snRNA" "$CONTAM_SNRNA" "$CURRENT_FASTQ" "${TEMP_DIR}/03_snrna" \
        "${OUTDIR}/snrna_snorna_reads/${SAMPLE_NAME}_snrna_reads.fq"
    CURRENT_FASTQ="$NEXT_FASTQ"
fi

# Filter snoRNA (save mapped reads for Phase 8 snRNA/snoRNA quantification)
if [ -f "$CURRENT_FASTQ" ]; then
    filter_contamination "snoRNA" "$CONTAM_SNORNA" "$CURRENT_FASTQ" "${TEMP_DIR}/04_snorna" \
        "${OUTDIR}/snrna_snorna_reads/${SAMPLE_NAME}_snorna_reads.fq"
    CURRENT_FASTQ="$NEXT_FASTQ"
fi

# Filter piRNA (save mapped reads for downstream piRNA gene counts and DESeq2 matrix)
if [ -f "$CURRENT_FASTQ" ]; then
    filter_contamination "piRNA" "$CONTAM_PIRNA" "$CURRENT_FASTQ" "${TEMP_DIR}/05_pirna" \
        "${OUTDIR}/pirna_reads/${SAMPLE_NAME}_pirna_reads.fq"
    CURRENT_FASTQ="$NEXT_FASTQ"
fi

# Filter other small RNAs (misc_RNA, Y_RNA, vault_RNA, ribozyme)
if [ -f "$CURRENT_FASTQ" ]; then
    filter_contamination "other" "$CONTAM_OTHER" "$CURRENT_FASTQ" "${TEMP_DIR}/06_other"
    CURRENT_FASTQ="$NEXT_FASTQ"
fi

# ============================================================================
# Finalize Output
# ============================================================================

# Move final cleaned file to output directory
if [ -f "$CURRENT_FASTQ" ]; then
    cp "$CURRENT_FASTQ" ${OUTDIR}/contamination_filtered/${SAMPLE_NAME}_cleaned.fq.gz
    echo ""
    echo "Cleaned reads saved to: ${OUTDIR}/contamination_filtered/${SAMPLE_NAME}_cleaned.fq.gz"
else
    echo "Error: No output file generated"
    exit 1
fi

# Clean up temp directory
rm -rf $TEMP_DIR

# ============================================================================
# Generate Statistics
# ============================================================================

echo ""
echo "Generating contamination statistics for $SAMPLE_NAME..."

# Export variables for Python
export SAMPLE_NAME
export OUTDIR
export TOTAL_READS
export REMAINING_READS
# Export contamination stats as JSON-like string
export CONTAM_RRNA="${CONTAM_STATS[rRNA]:-0}"
export CONTAM_TRNA="${CONTAM_STATS[tRNA]:-0}"
export CONTAM_SNRNA="${CONTAM_STATS[snRNA]:-0}"
export CONTAM_SNORNA="${CONTAM_STATS[snoRNA]:-0}"
export CONTAM_PIRNA="${CONTAM_STATS[piRNA]:-0}"
export CONTAM_OTHER="${CONTAM_STATS[other]:-0}"

python3 << 'PYTHON_STATS'
import os

sample = os.environ.get('SAMPLE_NAME')
outdir = os.environ.get('OUTDIR')

if not sample or not outdir:
    print("ERROR: Required environment variables not set")
    import sys
    sys.exit(1)

stats_file = os.path.join(outdir, 'contamination_filtered', f'{sample}_contamination_stats.txt')

# Get statistics from environment variables
total_reads = int(os.environ.get('TOTAL_READS', 0))
remaining_reads = int(os.environ.get('REMAINING_READS', 0))

# Parse contamination counts
contam_counts = {
    'rRNA': int(os.environ.get('CONTAM_RRNA', 0)),
    'tRNA': int(os.environ.get('CONTAM_TRNA', 0)),
    'snRNA': int(os.environ.get('CONTAM_SNRNA', 0)),
    'snoRNA': int(os.environ.get('CONTAM_SNORNA', 0)),
    'piRNA': int(os.environ.get('CONTAM_PIRNA', 0)),
    'other': int(os.environ.get('CONTAM_OTHER', 0))
}

total_contam = sum(contam_counts.values())

print("\n" + "="*80)
print(f"Contamination Filtering Statistics for {sample}")
print("="*80 + "\n")

with open(stats_file, 'w') as f:
    f.write(f"Contamination Filtering Statistics for {sample}\n")
    f.write("="*80 + "\n\n")
    
    f.write(f"Input Reads: {total_reads:,}\n")
    f.write(f"Output Reads (Clean): {remaining_reads:,}\n")
    f.write(f"Total Contaminating Reads Removed: {total_contam:,}\n\n")
    
    f.write("Contamination Breakdown:\n")
    for rna_type, count in sorted(contam_counts.items(), key=lambda x: x[1], reverse=True):
        if total_reads > 0:
            pct = (count / total_reads) * 100
            line = f"  {rna_type}: {count:,} ({pct:.2f}%)\n"
        else:
            line = f"  {rna_type}: {count:,}\n"
        f.write(line)
        print(line.strip())
    
    if total_reads > 0:
        clean_pct = (remaining_reads / total_reads) * 100
        contam_pct = (total_contam / total_reads) * 100
        f.write(f"\nClean Reads: {remaining_reads:,} ({clean_pct:.2f}%)\n")
        f.write(f"Contamination Rate: {contam_pct:.2f}%\n")
        print(f"\nClean Reads: {remaining_reads:,} ({clean_pct:.2f}%)")
        print(f"Contamination Rate: {contam_pct:.2f}%")

print("\n" + "="*80 + "\n")
print(f"Statistics saved to: {stats_file}")
PYTHON_STATS

echo "Contamination filtering completed for $SAMPLE_NAME"

printf "\n\nCompleted\n\n"
pwd
date
