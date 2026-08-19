#!/usr/bin/env bash
#SBATCH --job-name=snrna_snorna_matrix
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


###############################################################################
# 6.3_count_matrix.sh - Combine per-sample snRNA/snoRNA/piRNA counts into matrices
#
# Description:
#   Phase VI.3: For each of the three RNA types (snRNA, snoRNA, piRNA),
#   aggregate the per-sample count files produced by
#   6.2_snrna_snorna_quantification.sh (samtools idxstats format) into a
#   single integer count matrix compatible with DESeq2, edgeR, or limma.
#
#   One independent set of matrices is produced per RNA type:
#     snRNA  → ${OUTDIR}/snrna_expression/
#     snoRNA → ${OUTDIR}/snorna_expression/
#     piRNA  → ${OUTDIR}/pirna_expression/
#
#   Each output directory contains:
#     - Raw integer count matrix (DESeq2 / edgeR input)
#     - CPM-normalised matrix (visualisation only)
#     - log2(CPM+1) matrix (visualisation only)
#     - Per-gene summary statistics
#     - Per-sample summary statistics
#     - Expression statistics report
#
# Usage:
#   sbatch scripts/6.3_count_matrix.sh
#
# Prerequisites:
#   - Run 6.2_snrna_snorna_quantification.sh for ALL samples first.
#
# Note: This is NOT an array job – all samples are aggregated together.
###############################################################################

printf "\n\nStarted Job: 6.3-snrna_snorna_count_matrix\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 6.2_snrna_snorna_quantification.sh)
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Sample manifest
MANIFEST="/storage/users/ahinsu/DD-microRNA/sample_manifest.txt"

# Conda environment name
CONDA_ENV="smallRNA"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

mkdir -p "${OUTDIR}/snrna_expression"
mkdir -p "${OUTDIR}/snorna_expression"
mkdir -p "${OUTDIR}/pirna_expression"

export OUTDIR
export MANIFEST

# ============================================================================
# Helper: build count matrix for one RNA type
# Called three times (snrna, snorna, pirna)
# ============================================================================

build_matrix() {
    local label="$1"        # e.g. snrna / snorna / pirna
    local display="$2"      # e.g. snRNA / snoRNA / piRNA
    local counts_dir="${OUTDIR}/${label}_counts"
    local expr_dir="${OUTDIR}/${label}_expression"

    echo ""
    echo "========================================================"
    echo "Building count matrix for ${display}"
    echo "  Counts dir : ${counts_dir}"
    echo "  Output dir : ${expr_dir}"
    echo "========================================================"

    TYPE_LABEL="${label}" \
    TYPE_DISPLAY="${display}" \
    COUNTS_DIR="${counts_dir}" \
    EXPR_DIR="${expr_dir}" \
    python3 << 'PYTHON_MATRIX'
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

outdir       = os.environ['OUTDIR']
manifest     = os.environ['MANIFEST']
label        = os.environ['TYPE_LABEL']
display      = os.environ['TYPE_DISPLAY']
counts_dir   = os.environ['COUNTS_DIR']
expr_dir     = os.environ['EXPR_DIR']

# -------------------------------------------------------------------------
# 1. Discover sample names from the manifest
# -------------------------------------------------------------------------
samples = []
if os.path.exists(manifest):
    with open(manifest) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#') or line.lower().startswith('sample'):
                continue
            samples.append(line.split('\t')[0])
print(f"Manifest samples: {len(samples)}")

# -------------------------------------------------------------------------
# 2. Load per-sample counts (samtools idxstats format)
#    Columns: ref_name, ref_length, mapped_reads, unmapped_reads
# -------------------------------------------------------------------------
all_counts = {}   # sample -> pd.Series(count, index=ref_name)

for sample in samples:
    cnt_file = os.path.join(counts_dir, f'{sample}_{label}_counts.txt')
    if not os.path.exists(cnt_file):
        print(f"  WARNING: count file missing for {sample}: {cnt_file}")
        continue
    data = {}
    with open(cnt_file) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            if parts[0] == '*':  # idxstats summary row – skip
                continue
            ref = parts[0]
            try:
                mapped = int(parts[2])
            except ValueError:
                continue
            data[ref] = mapped
    if data:
        all_counts[sample] = pd.Series(data)
        print(f"  Loaded {sample}: {len(data)} sequences")

if not all_counts:
    print(f"ERROR: No {display} count data found. Run 6.2_snrna_snorna_quantification.sh first.")
    sys.exit(1)

# -------------------------------------------------------------------------
# 3. Build count matrix
# -------------------------------------------------------------------------
count_matrix = pd.DataFrame(all_counts).fillna(0).astype(int)
count_matrix['_total'] = count_matrix.sum(axis=1)
count_matrix = count_matrix.sort_values('_total', ascending=False).drop('_total', axis=1)
print(f"\nCount matrix: {count_matrix.shape[0]} sequences x {count_matrix.shape[1]} samples")

# -------------------------------------------------------------------------
# 4. Save raw count matrix (DESeq2-ready)
# -------------------------------------------------------------------------
out_counts = os.path.join(expr_dir, f'{label}_count_matrix.tsv')
count_matrix.to_csv(out_counts, sep='\t', index_label='Genes')
print(f"Raw count matrix: {out_counts}")

# -------------------------------------------------------------------------
# 5. CPM-normalised matrix (visualisation only)
# -------------------------------------------------------------------------
cpm = count_matrix.astype(float).copy()
for col in cpm.columns:
    total = cpm[col].sum()
    if total > 0:
        cpm[col] = (cpm[col] / total) * 1e6
cpm.to_csv(os.path.join(expr_dir, f'{label}_cpm_normalized.tsv'), sep='\t', index_label='Genes', float_format='%.2f')

# -------------------------------------------------------------------------
# 6. log2(CPM+1) matrix (visualisation only)
# -------------------------------------------------------------------------
log2cpm = np.log2(cpm + 1)
log2cpm.to_csv(os.path.join(expr_dir, f'{label}_log2cpm_normalized.tsv'), sep='\t', index_label='Genes', float_format='%.3f')

# -------------------------------------------------------------------------
# 7. Per-gene summary
# -------------------------------------------------------------------------
gene_summary = pd.DataFrame({
    'Total_Counts':      count_matrix.sum(axis=1),
    'Mean_Counts':       count_matrix.mean(axis=1),
    'Median_Counts':     count_matrix.median(axis=1),
    'Max_Counts':        count_matrix.max(axis=1),
    'Samples_Detected':  (count_matrix > 0).sum(axis=1),
})
gene_summary.to_csv(
    os.path.join(expr_dir, f'{label}_gene_summary.tsv'), sep='\t', float_format='%.2f')

# -------------------------------------------------------------------------
# 8. Per-sample summary
# -------------------------------------------------------------------------
sample_summary = pd.DataFrame({
    'Total_Reads':          count_matrix.sum(axis=0),
    'Distinct_seqs_1plus':  (count_matrix > 0).sum(axis=0),
    'Distinct_seqs_10plus': (count_matrix > 9).sum(axis=0),
    'Distinct_seqs_100plus':(count_matrix > 99).sum(axis=0),
    'Mean_Counts':          count_matrix.mean(axis=0),
    'Median_Counts':        count_matrix.median(axis=0),
})
sample_summary.to_csv(
    os.path.join(expr_dir, f'{label}_sample_summary.tsv'), sep='\t', float_format='%.2f')

# -------------------------------------------------------------------------
# 9. Statistics report
# -------------------------------------------------------------------------
stats_path = os.path.join(expr_dir, f'{label}_expression_stats.txt')
with open(stats_path, 'w') as fh:
    fh.write("=" * 70 + "\n")
    fh.write(f"Phase VI.3 – {display} Count Matrix Statistics\n")
    fh.write(f"Run date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    fh.write("=" * 70 + "\n\n")
    fh.write(f"Samples processed         : {len(count_matrix.columns)}\n")
    fh.write(f"Total sequences in matrix : {len(count_matrix)}\n\n")
    fh.write("Detection thresholds (across all samples combined):\n")
    fh.write(f"  Sequences with ≥1 read  : {(count_matrix.sum(axis=1) > 0).sum()}\n")
    fh.write(f"  Sequences with ≥10 reads: {(count_matrix.sum(axis=1) >= 10).sum()}\n")
    fh.write(f"  Sequences with ≥100 reads: {(count_matrix.sum(axis=1) >= 100).sum()}\n\n")
    fh.write(f"Top 10 most abundant {display} sequences:\n")
    for idx, row in gene_summary.nlargest(10, 'Total_Counts').iterrows():
        fh.write(f"  {idx}: {int(row['Total_Counts'])} reads "
                 f"(detected in {int(row['Samples_Detected'])} samples)\n")
    fh.write(f"\nOutput files (in {expr_dir}/):\n")
    fh.write(f"  {label}_count_matrix.tsv           (DESeq2/edgeR ready)\n")
    fh.write(f"  {label}_cpm_normalized.tsv          (visualisation only)\n")
    fh.write(f"  {label}_log2cpm_normalized.tsv      (visualisation only)\n")
    fh.write(f"  {label}_gene_summary.tsv\n")
    fh.write(f"  {label}_sample_summary.tsv\n")
print(f"\n{open(stats_path).read()}")
PYTHON_MATRIX
}

# ============================================================================
# Build matrices for each RNA type
# ============================================================================

build_matrix "snrna"  "snRNA"
build_matrix "snorna" "snoRNA"
build_matrix "pirna"  "piRNA"

echo ""
echo "Count matrix preparation completed for all RNA types!"
echo ""
echo "Output directories:"
echo "  snRNA  : ${OUTDIR}/snrna_expression/"
echo "  snoRNA : ${OUTDIR}/snorna_expression/"
echo "  piRNA  : ${OUTDIR}/pirna_expression/"
echo ""

printf "\n\nCompleted: 6.3-snrna_snorna_count_matrix\n\n"
pwd
date
