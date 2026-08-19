#!/usr/bin/env bash
#SBATCH --job-name=genome_alignment
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-72   # Update to match number of samples in your manifest


###############################################################################
# 3.2_genome_alignment.sh - Map cleaned reads to the Ensembl reference genome
#
# Description:
#   Phase III.2: Align quality-filtered reads to the Ensembl reference genome
#   using Bowtie (v1) with two mismatches (-v 2), which is the standard
#   parameter for short small-RNA reads.
#
#   For each sample (SLURM array job):
#     1. Map contamination-filtered reads to the genome with bowtie -v 2.
#     2. Save mapped reads as a FASTQ file for downstream miRDeep2 (4.1) and
#        miRNA quantification (5.2) steps.
#     3. Parse the Bowtie log and write a per-sample mapping statistics file.
#        Each job writes only its own file so there is no concurrent-write
#        race condition on a shared output file.
#
#   After all array tasks finish, run 3.3_alignment_stats.sh to aggregate
#   the per-sample statistics into a single combined TSV.
#
# Usage:
#   sbatch scripts/3.2_genome_alignment.sh
#
# Prerequisites:
#   - Run 1.5b_contamination_filtering.sh first (produces cleaned FASTQ files).
#   - A Bowtie index for the Ensembl genome must exist (build with bowtie-build).
#     Example: bowtie-build genome.fa /path/to/index/genome (see 3.1_genome_index.sh)
#
# Outputs (per sample, in ${OUTDIR}/alignment/):
#   ${SAMPLE_NAME}_mapped.fastq         – Mapped reads as FASTQ (for Phase IV.1)
#   ${SAMPLE_NAME}_alignment.log        – Bowtie alignment log
#   ${SAMPLE_NAME}_alignment_stats.txt  – Per-sample mapping statistics summary
###############################################################################

printf "\n\nStarted Job: 3.2-genome_alignment\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Bowtie index for the Ensembl reference genome (without .1.ebwt suffix)
# Build with: bowtie-build /path/to/genome.fa /path/to/index/genome
GENOME_INDEX="/storage/users/ahinsu/DD-microRNA/reference/index/genome"

# Conda environment name
CONDA_ENV="smallRNA"

# Sample manifest (TSV with at least a 'sample_name' column in column 1)
MANIFEST="/storage/users/ahinsu/DD-microRNA/sample_manifest.txt"

# Bowtie alignment parameters
# Two mismatches (-v 2) is the standard parameter for small-RNA genome alignment
# (see miRDeep2 mapper.pl documentation and nf-core/smrnaseq pipeline).  It
# balances sensitivity to SNPs and RNA editing while remaining sufficiently
# stringent for the short (17-30 nt) reads produced by small-RNA sequencing.
MISMATCHES=2          # Allow up to 2 mismatches (standard for short small-RNA reads)
MAX_ALIGNMENTS=50     # Report up to 50 alignments per read (important for miRNAs)

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

mkdir -p "${OUTDIR}/alignment"

# ============================================================================
# SLURM Array Job – resolve sample for this task
# ============================================================================

sample_line=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n+1{print $0}' "$MANIFEST")
SAMPLE_NAME=$(echo "$sample_line" | cut -f1)

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: Could not determine sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Locate input FASTQ (prefer contamination-filtered output; fall back to fastp)
FASTQ_FILE="${OUTDIR}/contamination_filtered/${SAMPLE_NAME}_cleaned.fq.gz"
if [ ! -f "$FASTQ_FILE" ]; then
    FASTQ_FILE="${OUTDIR}/fastp/${SAMPLE_NAME}_trimmed.fq.gz"
fi
if [ ! -f "$FASTQ_FILE" ]; then
    echo "ERROR: FASTQ file not found for sample ${SAMPLE_NAME}"
    exit 1
fi

echo "Processing sample : ${SAMPLE_NAME}"
echo "Input FASTQ       : ${FASTQ_FILE}"
echo "Genome index      : ${GENOME_INDEX}"
echo "Mismatches        : ${MISMATCHES}"

# ============================================================================
# Step 1: Align reads with Bowtie (SAM output → sorted BAM)
# ============================================================================

echo ""
echo "Step 1: Running Bowtie alignment (v=2 mismatches)..."

bowtie \
    -p ${SLURM_CPUS_PER_TASK:-4} \
    -v ${MISMATCHES} \
    -q \
    -n 1 \
    -e 80 \
    -l 18 \
    -a \
    --al "${OUTDIR}/alignment/${SAMPLE_NAME}_mapped.fastq" \
    --best --strata  \
    ${GENOME_INDEX} \
     <(zcat "$FASTQ_FILE") \
      2> "${OUTDIR}/alignment/${SAMPLE_NAME}_alignment.log" > /dev/null


# ============================================================================
# Step 2: Write per-sample mapping statistics
# ============================================================================
# Each array task writes only its own file → no concurrent-write race condition.
# Run 3.3_alignment_stats.sh after all array tasks finish to aggregate these
# per-sample files into a combined all_samples_alignment_stats.tsv.

echo ""
echo "Step 2: Writing per-sample mapping statistics..."

export OUTDIR
export SAMPLE_NAME

python3 << 'PYTHON_STATS'
import os
import re

outdir    = os.environ['OUTDIR']
sample    = os.environ['SAMPLE_NAME']
align_dir = os.path.join(outdir, 'alignment')
log_file  = os.path.join(align_dir, f'{sample}_alignment.log')
stats_out = os.path.join(align_dir, f'{sample}_alignment_stats.txt')

# Extract stats from Bowtie log
total = mapped = unmapped = reported = 0

with open(log_file) as fh:
    for line in fh:
        m = re.search(r'# reads processed:\s+(\d+)', line)
        if m:
            total = int(m.group(1))

        m = re.search(r'# reads with at least one alignment:\s+(\d+)', line)
        if m:
            mapped = int(m.group(1))

        m = re.search(r'# reads that failed to align:\s+(\d+)', line)
        if m:
            unmapped = int(m.group(1))

        m = re.search(r'Reported\s+(\d+)\s+alignments', line)
        if m:
            reported = int(m.group(1))

pct_mapped   = round(mapped   / total * 100, 2) if total > 0 else 0.0
pct_unmapped = round(unmapped / total * 100, 2) if total > 0 else 0.0

# Write human-readable per-sample stats file (format parsed by 7.1_multiqc_report.sh)
with open(stats_out, 'w') as fh:
    fh.write(f"Sample                : {sample}\n")
    fh.write(f"Total reads           : {total}\n")
    fh.write(f"Mapped reads          : {mapped} ({pct_mapped}%)\n")
    fh.write(f"Unmapped reads        : {unmapped} ({pct_unmapped}%)\n")
    fh.write(f"Reported alignments   : {reported}\n")

print(f"Written per-sample stats: {stats_out}")
PYTHON_STATS


echo ""
echo "Alignment completed for ${SAMPLE_NAME}"
echo "Next: after all array tasks finish, run 3.3_alignment_stats.sh to aggregate stats."

printf "\n\nCompleted: 3.2-genome_alignment\n\n"
pwd
date
