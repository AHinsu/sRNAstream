#!/usr/bin/env bash
#SBATCH --job-name=mirna_quant
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-72   # Update to match number of samples in your manifest


###############################################################################
# 5.2_mirna_quantification.sh - Per-sample miRNA read mapping and gene counting
#
# Description:
#   Phase V.2 (SLURM array job): For each sample, map the genome-mapped reads
#   (from 3.1_genome_alignment.sh) against the unified miRNA reference built by
#   5.1_mirna_reference.sh, and compute raw read counts per miRNA sequence.
#
#   Workflow (per sample):
#     1. Map the genome-mapped FASTQ (${SAMPLE_NAME}_mapped.fastq) to the
#        unified mature miRNA Bowtie1 index using -v 1 (1 mismatch, appropriate
#        for short ~22 nt miRNA reads against a compact reference).
#     2. Compute alignment statistics from the Bowtie log and samtools flagstat.
#     3. Count aligned reads per miRNA sequence using samtools idxstats.
#     4. Write per-sample count file for aggregation by 5.3_count_matrix.sh.
#
# Usage:
#   sbatch scripts/5.2_mirna_quantification.sh
#
# Prerequisites:
#   - Run 3.2_genome_alignment.sh for all samples (produces _mapped.fastq).
#   - Run 5.1_mirna_reference.sh (produces Bowtie1 index).
#
# Outputs (per sample, in ${OUTDIR}/mirna_counts/):
#   ${SAMPLE_NAME}_mirna.bam             – Unsorted BAM (intermediate; can be deleted)
#   ${SAMPLE_NAME}_mirna_sorted.bam      – Coordinate-sorted BAM aligned to miRNA reference
#   ${SAMPLE_NAME}_mirna_sorted.bam.bai  – BAM index (written by samtools sort --write-index)
#   ${SAMPLE_NAME}_mirna_counts.txt      – Per-sequence read counts (idxstats)
#   ${SAMPLE_NAME}_mirna_alignment.log   – Bowtie alignment log
#   ${SAMPLE_NAME}_mirna_flagstat.txt    – samtools flagstat
#   ${SAMPLE_NAME}_mirna_stats.txt       – Mapping statistics summary
###############################################################################

printf "\n\nStarted Job: 5.2-mirna_quantification\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 3.1_genome_alignment.sh)
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Reference directory (must match REF_DIR used in 5.1_mirna_reference.sh)
REF_DIR="/storage/users/ahinsu/DD-microRNA/reference"

# Bowtie1 index for unified mature miRNA reference (without .1.ebwt suffix)
MIRNA_INDEX="${REF_DIR}/mirna_quant_ref/unified_mature_index"

# Conda environment name
CONDA_ENV="smallRNA"

# Sample manifest
MANIFEST="/storage/users/ahinsu/DD-microRNA/sample_manifest.txt"

# Bowtie mismatch parameter for miRNA quantification
# 1 mismatch is used here (stricter than genome alignment's 2) to reduce
# cross-mapping between closely related miRNA family members in the compact
# miRNA reference, improving per-sequence quantification accuracy.
MISMATCHES=1

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

mkdir -p "${OUTDIR}/mirna_counts"

# ============================================================================
# Resolve sample for this SLURM array task
# ============================================================================

sample_line=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n+1{print $0}' "$MANIFEST")
SAMPLE_NAME=$(echo "$sample_line" | cut -f1)

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: Could not determine sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Input: genome-mapped reads from 3.2_genome_alignment.sh
INPUT_FQ="${OUTDIR}/alignment/${SAMPLE_NAME}_mapped.fastq"

if [ ! -f "$INPUT_FQ" ]; then
    echo "ERROR: Mapped reads not found: ${INPUT_FQ}"
    echo "Run 3.2_genome_alignment.sh first."
    exit 1
fi

echo "Processing sample : ${SAMPLE_NAME}"
echo "Input FASTQ       : ${INPUT_FQ}"
echo "miRNA index       : ${MIRNA_INDEX}"

# ============================================================================
# Step 1: Map reads to unified miRNA reference
# ============================================================================

echo ""
echo "Step 1: Mapping reads to unified miRNA reference..."

bowtie \
    -p ${SLURM_CPUS_PER_TASK:-4} \
    -v ${MISMATCHES} \
    -k 1 \
    --best \
    --chunkmbs 2048 \
    -S --fullref \
    -x "${MIRNA_INDEX}" \
    -q "$INPUT_FQ" \
    2> "${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna_alignment.log" \
    | samtools view -@ ${SLURM_CPUS_PER_TASK:-4} -bS -o "${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna.bam" - 

samtools sort -@ ${SLURM_CPUS_PER_TASK:-4} --write-index \
        -o "${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna_sorted.bam" "${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna.bam"

echo "Sorted BAM: ${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna_sorted.bam"

# ============================================================================
# Step 2: Count reads per miRNA using samtools idxstats
# ============================================================================

echo ""
echo "Step 2: Counting reads per miRNA (samtools idxstats)..."

# idxstats output: ref_name, ref_length, mapped_reads, unmapped_reads
samtools idxstats \
    "${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna_sorted.bam" \
    > "${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna_counts.txt"

DETECTED=$(awk '$1 != "*" && $3 > 0' "${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna_counts.txt" | wc -l)
echo "miRNAs with ≥1 mapped read: ${DETECTED}"

# ============================================================================
# Step 3: Mapping statistics
# ============================================================================

echo ""
echo "Step 3: Calculating mapping statistics..."

samtools flagstat \
    "${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna_sorted.bam" \
    > "${OUTDIR}/mirna_counts/${SAMPLE_NAME}_mirna_flagstat.txt"

export OUTDIR
export SAMPLE_NAME

python3 << 'PYTHON_STATS'
import os
import re

outdir  = os.environ['OUTDIR']
sample  = os.environ['SAMPLE_NAME']
cnt_dir = os.path.join(outdir, 'mirna_counts')

log_file    = os.path.join(cnt_dir, f'{sample}_mirna_alignment.log')
flagstat    = os.path.join(cnt_dir, f'{sample}_mirna_flagstat.txt')
counts_file = os.path.join(cnt_dir, f'{sample}_mirna_counts.txt')
stats_out   = os.path.join(cnt_dir, f'{sample}_mirna_stats.txt')

stats = {
    'total_input_reads': 0,
    'mapped_reads': 0,
    'unmapped_reads': 0,
    'mapping_pct': 0.0,
    'mirnas_detected_1plus': 0,
    'mirnas_detected_10plus': 0,
    'mirnas_detected_100plus': 0,
}

# Parse Bowtie log
if os.path.exists(log_file):
    with open(log_file) as fh:
        for line in fh:
            m = re.search(r' reads processed: (\d+)', line)
            if m:
                stats['total_input_reads'] = int(m.group(1))
            

# Parse flagstat
if os.path.exists(flagstat):
    with open(flagstat) as fh:
        for line in fh:
            m = re.match(r'(\d+) \+ \d+ mapped', line)
            if m:
                stats['mapped_reads'] = int(m.group(1))

stats['unmapped_reads'] = stats['total_input_reads'] - stats['mapped_reads']
if stats['total_input_reads'] > 0:
    stats['mapping_pct'] = round(stats['mapped_reads'] / stats['total_input_reads'] * 100, 2)

# Parse counts (idxstats: name, length, mapped, unmapped)
if os.path.exists(counts_file):
    with open(counts_file) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            if parts[0] == '*':   # idxstats summary row – skip
                continue
            try:
                cnt = int(parts[2])
            except ValueError:
                continue
            if cnt >= 1:
                stats['mirnas_detected_1plus'] += 1
            if cnt >= 10:
                stats['mirnas_detected_10plus'] += 1
            if cnt >= 100:
                stats['mirnas_detected_100plus'] += 1

with open(stats_out, 'w') as fh:
    fh.write("=" * 60 + "\n")
    fh.write(f"miRNA Quantification Statistics – {sample}\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(f"Input reads (genome-mapped)  : {stats['total_input_reads']:>10,}\n")
    fh.write(f"Mapped to miRNA reference    : {stats['mapped_reads']:>10,}  ({stats['mapping_pct']:.2f}%)\n")
    fh.write(f"Unmapped                     : {stats['unmapped_reads']:>10,}\n\n")
    fh.write(f"Distinct miRNAs detected:\n")
    fh.write(f"  ≥1 read                    : {stats['mirnas_detected_1plus']:>10,}\n")
    fh.write(f"  ≥10 reads                  : {stats['mirnas_detected_10plus']:>10,}\n")
    fh.write(f"  ≥100 reads                 : {stats['mirnas_detected_100plus']:>10,}\n")

print(open(stats_out).read())
PYTHON_STATS


echo ""
echo "Quantification completed for ${SAMPLE_NAME}"

printf "\n\nCompleted: 5.2-mirna_quantification\n\n"
pwd
date
