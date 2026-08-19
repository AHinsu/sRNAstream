#!/usr/bin/env bash
#SBATCH --job-name=snrna_snorna_quant
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A_%a.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-72   # Update to match number of samples in your manifest


###############################################################################
# 6.2_snrna_snorna_quantification.sh - Per-sample snRNA/snoRNA/piRNA read mapping
#
# Description:
#   Phase VI.2 (SLURM array job): For each sample, map the reads that were
#   captured as snRNA, snoRNA, and piRNA during contamination filtering
#   (1.5b_contamination_filtering.sh) against the per-type Bowtie1 indices
#   built by 6.1_snrna_snorna_reference.sh, and compute raw read counts per
#   sequence using samtools idxstats.
#
#   Each RNA type is quantified independently:
#     snRNA  – reads from ${OUTDIR}/snrna_snorna_reads/${SAMPLE_NAME}_snrna_reads.fq
#     snoRNA – reads from ${OUTDIR}/snrna_snorna_reads/${SAMPLE_NAME}_snorna_reads.fq
#     piRNA  – reads from ${OUTDIR}/pirna_reads/${SAMPLE_NAME}_pirna_reads.fq
#
#   Workflow (per sample, per RNA type):
#     1. Map reads to the per-type Bowtie1 index (same parameters as 5.2).
#     2. Sort and index the BAM.
#     3. Count aligned reads per sequence with samtools idxstats.
#     4. Write per-sample count file for aggregation by 6.3_count_matrix.sh.
#
# Usage:
#   sbatch scripts/6.2_snrna_snorna_quantification.sh
#
# Prerequisites:
#   - Run 1.5b_contamination_filtering.sh for all samples (produces per-type FASTQs).
#   - Run 6.1_snrna_snorna_reference.sh (produces Bowtie1 indices).
#
# Outputs (per sample, in ${OUTDIR}/snrna_counts/, ${OUTDIR}/snorna_counts/,
#          ${OUTDIR}/pirna_counts/):
#   ${SAMPLE_NAME}_{type}.bam             – Unsorted BAM (intermediate)
#   ${SAMPLE_NAME}_{type}_sorted.bam      – Coordinate-sorted BAM
#   ${SAMPLE_NAME}_{type}_sorted.bam.bai  – BAM index
#   ${SAMPLE_NAME}_{type}_counts.txt      – Per-sequence read counts (idxstats)
#   ${SAMPLE_NAME}_{type}_alignment.log   – Bowtie alignment log
#   ${SAMPLE_NAME}_{type}_flagstat.txt    – samtools flagstat
#   ${SAMPLE_NAME}_{type}_stats.txt       – Mapping statistics summary
###############################################################################

printf "\n\nStarted Job: 6.2-snrna_snorna_quantification\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 1.5b_contamination_filtering.sh)
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Reference directory (must match REF_DIR used in 6.1_snrna_snorna_reference.sh)
REF_DIR="/storage/users/ahinsu/DD-microRNA/reference"

# Bowtie1 index prefix for each RNA type (without .1.ebwt suffix)
SNRNA_INDEX="${REF_DIR}/snrna_quant_ref/snrna_index"
SNORNA_INDEX="${REF_DIR}/snrna_quant_ref/snorna_index"
PIRNA_INDEX="${REF_DIR}/snrna_quant_ref/pirna_index"

# Conda environment name
CONDA_ENV="smallRNA"

# Sample manifest
MANIFEST="/storage/users/ahinsu/DD-microRNA/sample_manifest.txt"

# Bowtie mismatch parameter
# 1 mismatch is used (same as 5.2) for accurate per-sequence quantification
MISMATCHES=1

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

mkdir -p "${OUTDIR}/snrna_counts"
mkdir -p "${OUTDIR}/snorna_counts"
mkdir -p "${OUTDIR}/pirna_counts"

# ============================================================================
# Resolve sample for this SLURM array task
# ============================================================================

sample_line=$(awk -v n=${SLURM_ARRAY_TASK_ID} 'NR==n+1{print $0}' "$MANIFEST")
SAMPLE_NAME=$(echo "$sample_line" | cut -f1)

if [ -z "$SAMPLE_NAME" ]; then
    echo "ERROR: Could not determine sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Processing sample : ${SAMPLE_NAME}"

# Input FASTQs from 1.5b_contamination_filtering.sh
SNRNA_FQ="${OUTDIR}/snrna_snorna_reads/${SAMPLE_NAME}_snrna_reads.fq"
SNORNA_FQ="${OUTDIR}/snrna_snorna_reads/${SAMPLE_NAME}_snorna_reads.fq"
PIRNA_FQ="${OUTDIR}/pirna_reads/${SAMPLE_NAME}_pirna_reads.fq"

# ============================================================================
# Helper: map reads of one RNA type and count per sequence
#
# Usage: map_and_count <label> <index> <input_fq> <out_dir>
#   label    : RNA type label used in filenames (snrna / snorna / pirna)
#   index    : Bowtie1 index prefix (without .1.ebwt)
#   input_fq : Input FASTQ (uncompressed, from 1.5b)
#   out_dir  : Output directory for BAM, counts, and stats files
# ============================================================================

map_and_count() {
    local label="$1"
    local display="$2"
    local index="$3"
    local input_fq="$4"
    local out_dir="$5"

    echo ""
    echo "--- Processing ${display} for ${SAMPLE_NAME} ---"

    # Validate input
    if [ ! -f "${input_fq}" ]; then
        echo "  WARNING: Input FASTQ not found: ${input_fq}"
        echo "           Skipping ${display} (run 1.5b_contamination_filtering.sh first)."
        return
    fi

    if [ ! -f "${index}.1.ebwt" ]; then
        echo "  WARNING: Bowtie index not found: ${index}.1.ebwt"
        echo "           Skipping ${display} (run 6.1_snrna_snorna_reference.sh first)."
        return
    fi

    local bam="${out_dir}/${SAMPLE_NAME}_${label}.bam"
    local sorted_bam="${out_dir}/${SAMPLE_NAME}_${label}_sorted.bam"
    local counts="${out_dir}/${SAMPLE_NAME}_${label}_counts.txt"
    local log="${out_dir}/${SAMPLE_NAME}_${label}_alignment.log"
    local flagstat="${out_dir}/${SAMPLE_NAME}_${label}_flagstat.txt"
    local stats="${out_dir}/${SAMPLE_NAME}_${label}_stats.txt"

    echo "  Input FASTQ  : ${input_fq}"
    echo "  Index        : ${index}"
    echo "  Output dir   : ${out_dir}"

    # Step 1: Map reads
    echo ""
    echo "  Step 1: Mapping reads to ${display} reference..."

    bowtie \
        -p ${SLURM_CPUS_PER_TASK:-4} \
        -v ${MISMATCHES} \
        -k 1 \
        --best \
        --chunkmbs 2048 \
        -S --fullref \
        -x "${index}" \
        -q "${input_fq}" \
        2> "${log}" \
        | samtools view -@ ${SLURM_CPUS_PER_TASK:-4} -bS -o "${bam}" -

    samtools sort -@ ${SLURM_CPUS_PER_TASK:-4} --write-index \
        -o "${sorted_bam}" "${bam}"

    echo "  Sorted BAM: ${sorted_bam}"

    # Step 2: Count reads per sequence using samtools idxstats
    echo ""
    echo "  Step 2: Counting reads per sequence (samtools idxstats)..."

    samtools idxstats "${sorted_bam}" > "${counts}"

    DETECTED=$(awk '$1 != "*" && $3 > 0' "${counts}" | wc -l)
    echo "  Sequences with ≥1 mapped read: ${DETECTED}"

    # Step 3: Mapping statistics
    echo ""
    echo "  Step 3: Calculating mapping statistics..."

    samtools flagstat "${sorted_bam}" > "${flagstat}"

    # Parse stats with Python
    TYPE_LABEL="${label}" \
    TYPE_DISPLAY="${display}" \
    LOG_FILE="${log}" \
    FLAGSTAT="${flagstat}" \
    COUNTS_FILE="${counts}" \
    STATS_OUT="${stats}" \
    SAMPLE="${SAMPLE_NAME}" \
    python3 << 'PYTHON_STATS'
import os, re

sample      = os.environ['SAMPLE']
label       = os.environ['TYPE_LABEL']
display     = os.environ['TYPE_DISPLAY']
log_file    = os.environ['LOG_FILE']
flagstat    = os.environ['FLAGSTAT']
counts_file = os.environ['COUNTS_FILE']
stats_out   = os.environ['STATS_OUT']

stats = {
    'total_input_reads': 0,
    'mapped_reads': 0,
    'unmapped_reads': 0,
    'mapping_pct': 0.0,
    'detected_1plus': 0,
    'detected_10plus': 0,
    'detected_100plus': 0,
}

# Parse flagstat – this is the primary source for both total and mapped reads
# because the sorted BAM always contains all reads (mapped + unmapped) when
# Bowtie is run with -S (SAM output), making the 'in total' count reliable.
if os.path.exists(flagstat):
    with open(flagstat) as fh:
        for line in fh:
            m = re.match(r'(\d+) \+ \d+ in total', line)
            if m:
                stats['total_input_reads'] = int(m.group(1))
            m = re.match(r'(\d+) \+ \d+ mapped', line)
            if m:
                stats['mapped_reads'] = int(m.group(1))

# Fallback: parse Bowtie1 stderr log for total reads if flagstat didn't give it
if stats['total_input_reads'] == 0 and os.path.exists(log_file):
    with open(log_file) as fh:
        for line in fh:
            m = re.search(r'# reads processed: (\d+)', line)
            if m:
                stats['total_input_reads'] = int(m.group(1))

stats['unmapped_reads'] = max(0, stats['total_input_reads'] - stats['mapped_reads'])
if stats['total_input_reads'] > 0:
    stats['mapping_pct'] = round(
        stats['mapped_reads'] / stats['total_input_reads'] * 100, 2)

# Parse counts (idxstats: name, length, mapped, unmapped)
if os.path.exists(counts_file):
    with open(counts_file) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            if parts[0] == '*':  # idxstats summary row – skip
                continue
            try:
                cnt = int(parts[2])
            except ValueError:
                continue
            if cnt >= 1:
                stats['detected_1plus'] += 1
            if cnt >= 10:
                stats['detected_10plus'] += 1
            if cnt >= 100:
                stats['detected_100plus'] += 1

rna_label = display
with open(stats_out, 'w') as fh:
    fh.write("=" * 60 + "\n")
    fh.write(f"{rna_label} Quantification Statistics – {sample}\n")
    fh.write("=" * 60 + "\n\n")
    fh.write(f"Input reads                  : {stats['total_input_reads']:>10,}\n")
    fh.write(f"Mapped to {rna_label} reference : {stats['mapped_reads']:>10,}  ({stats['mapping_pct']:.2f}%)\n")
    fh.write(f"Unmapped                     : {stats['unmapped_reads']:>10,}\n\n")
    fh.write(f"Distinct sequences detected:\n")
    fh.write(f"  ≥1 read                    : {stats['detected_1plus']:>10,}\n")
    fh.write(f"  ≥10 reads                  : {stats['detected_10plus']:>10,}\n")
    fh.write(f"  ≥100 reads                 : {stats['detected_100plus']:>10,}\n")

print(open(stats_out).read())
PYTHON_STATS
}

# ============================================================================
# Map and count reads for each RNA type
# ============================================================================

map_and_count "snrna"  "snRNA"  "${SNRNA_INDEX}"  "${SNRNA_FQ}"  "${OUTDIR}/snrna_counts"
map_and_count "snorna" "snoRNA" "${SNORNA_INDEX}" "${SNORNA_FQ}" "${OUTDIR}/snorna_counts"
map_and_count "pirna"  "piRNA"  "${PIRNA_INDEX}"  "${PIRNA_FQ}"  "${OUTDIR}/pirna_counts"

echo ""
echo "Quantification completed for ${SAMPLE_NAME}"

printf "\n\nCompleted: 6.2-snrna_snorna_quantification\n\n"
pwd
date
