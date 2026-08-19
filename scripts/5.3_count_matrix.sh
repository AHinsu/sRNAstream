#!/usr/bin/env bash
#SBATCH --job-name=mirna_count_matrix
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


###############################################################################
# 5.3_count_matrix.sh - Combine per-sample miRNA counts into a DESeq2 matrix
#
# Description:
#   Phase V.3: Aggregate the per-sample miRNA count files produced by
#   5.2_mirna_quantification.sh (samtools idxstats format) into a single
#   integer count matrix compatible with DESeq2, edgeR, or limma.
#
#   Outputs:
#     - Raw integer count matrix (DESeq2 / edgeR input)
#     - CPM-normalised matrix (visualisation only)
#     - log2(CPM+1) matrix (visualisation only)
#     - Per-gene summary statistics
#     - Per-sample summary statistics
#
# Usage:
#   sbatch scripts/5.3_count_matrix.sh
#
# Prerequisites:
#   - Run 5.2_mirna_quantification.sh for ALL samples first.
#
# Note: This is NOT an array job – all samples are aggregated together.
###############################################################################

printf "\n\nStarted Job: 5.3-mirna_count_matrix\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 5.2_mirna_quantification.sh)
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

EXPR_DIR="${OUTDIR}/mirna_expression"
mkdir -p "${EXPR_DIR}"

export OUTDIR
export MANIFEST
export EXPR_DIR

python3 << 'PYTHON_MATRIX'
import os
import sys
import re
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

outdir   = os.environ['OUTDIR']
manifest = os.environ['MANIFEST']
expr_dir = os.environ['EXPR_DIR']
counts_dir = os.path.join(outdir, 'mirna_counts')

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
#    Columns: ref_name, ref_length, mapped, unmapped
# -------------------------------------------------------------------------
all_counts = {}   # sample -> pd.Series(count, index=ref_name)

for sample in samples:
    cnt_file = os.path.join(counts_dir, f'{sample}_mirna_counts.txt')
    if not os.path.exists(cnt_file):
        print(f"  WARNING: count file missing for {sample}: {cnt_file}")
        continue
    data = {}
    with open(cnt_file) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 3 or parts[0] == '*':
                continue
            ref = parts[0]
            try:
                mapped = int(parts[2])
            except ValueError:
                continue
            data[ref] = mapped
    if data:
        all_counts[sample] = pd.Series(data)
        print(f"  Loaded {sample}: {len(data)} miRNAs")

if not all_counts:
    print("ERROR: No count data found. Run 5.2_mirna_quantification.sh first.")
    sys.exit(1)

# -------------------------------------------------------------------------
# 3. Build count matrix
# -------------------------------------------------------------------------
count_matrix = pd.DataFrame(all_counts).fillna(0).astype(int)
count_matrix['_total'] = count_matrix.sum(axis=1)
count_matrix = count_matrix.sort_values('_total', ascending=False).drop('_total', axis=1)
print(f"\nCount matrix: {count_matrix.shape[0]} miRNAs x {count_matrix.shape[1]} samples")

# -------------------------------------------------------------------------
# 4. Save raw count matrix (DESeq2-ready)
# -------------------------------------------------------------------------
out_counts = os.path.join(expr_dir, 'mirna_count_matrix.tsv')
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
cpm.to_csv(os.path.join(expr_dir, 'mirna_cpm_normalized.tsv'), sep='\t', index_label='Genes', float_format='%.2f')

# -------------------------------------------------------------------------
# 6. log2(CPM+1) matrix (visualisation only)
# -------------------------------------------------------------------------
log2cpm = np.log2(cpm + 1)
log2cpm.to_csv(os.path.join(expr_dir, 'mirna_log2cpm_normalized.tsv'), sep='\t', index_label='Genes', float_format='%.3f')

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
gene_summary.to_csv(os.path.join(expr_dir, 'mirna_gene_summary.tsv'), sep='\t', float_format='%.2f')

# -------------------------------------------------------------------------
# 8. Per-sample summary
#    Distinct miRNA counts are split into total / known / novel using the
#    |src=known_mature / |src=novel_mature tags added by 5.1_mirna_reference.sh
# -------------------------------------------------------------------------
is_known = count_matrix.index.str.contains(r'\|src=known_mature', regex=True)
is_novel = count_matrix.index.str.contains(r'\|src=novel_mature', regex=True)
known_mat = count_matrix[is_known]
novel_mat = count_matrix[is_novel]

sample_summary = pd.DataFrame({
    'Total_miRNA_Reads':              count_matrix.sum(axis=0),
    'Distinct_miRNA_1plus':           (count_matrix > 0).sum(axis=0),
    'Distinct_miRNA_1plus_known':     (known_mat > 0).sum(axis=0),
    'Distinct_miRNA_1plus_novel':     (novel_mat > 0).sum(axis=0),
    'Distinct_miRNA_10plus':          (count_matrix > 9).sum(axis=0),
    'Distinct_miRNA_10plus_known':    (known_mat > 9).sum(axis=0),
    'Distinct_miRNA_10plus_novel':    (novel_mat > 9).sum(axis=0),
    'Distinct_miRNA_100plus':         (count_matrix > 99).sum(axis=0),
    'Distinct_miRNA_100plus_known':   (known_mat > 99).sum(axis=0),
    'Distinct_miRNA_100plus_novel':   (novel_mat > 99).sum(axis=0),
    'Mean_Counts':                    count_matrix.mean(axis=0),
    'Median_Counts':                  count_matrix.median(axis=0),
})
sample_summary.to_csv(os.path.join(expr_dir, 'mirna_sample_summary.tsv'), sep='\t', float_format='%.2f')

# -------------------------------------------------------------------------
# 9. Statistics report
# -------------------------------------------------------------------------
n_ref_total = len(count_matrix)
n_ref_known = int(is_known.sum())
n_ref_novel = int(is_novel.sum())

# Row-level totals across all samples
row_totals       = count_matrix.sum(axis=1)
row_totals_known = known_mat.sum(axis=1) if not known_mat.empty else pd.Series(dtype=int)
row_totals_novel = novel_mat.sum(axis=1) if not novel_mat.empty else pd.Series(dtype=int)

def n_detected(series, threshold):
    return int((series >= threshold).sum())

stats_path = os.path.join(expr_dir, 'mirna_expression_stats.txt')
with open(stats_path, 'w') as fh:
    fh.write("=" * 70 + "\n")
    fh.write("Phase V.3 – miRNA Count Matrix Statistics\n")
    fh.write(f"Run date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    fh.write("=" * 70 + "\n\n")
    fh.write(f"Samples processed        : {len(count_matrix.columns)}\n")
    fh.write(f"Total miRNA references   : {n_ref_total}\n")
    fh.write(f"  Known miRNA references : {n_ref_known}\n")
    fh.write(f"  Novel miRNA references : {n_ref_novel}\n\n")
    fh.write("Detection thresholds (across all samples combined):\n")
    fh.write(f"{'Threshold':<20} {'Total':>8} {'Known':>8} {'Novel':>8}\n")
    fh.write("-" * 48 + "\n")
    for label, thresh in [("≥1 read", 1), ("≥10 reads", 10), ("≥100 reads", 100)]:
        t = n_detected(row_totals, thresh)
        k = n_detected(row_totals_known, thresh)
        n = n_detected(row_totals_novel, thresh)
        fh.write(f"  {label:<18} {t:>8} {k:>8} {n:>8}\n")
    fh.write("\nTop 10 most abundant miRNAs:\n")
    for idx, row in gene_summary.nlargest(10, 'Total_Counts').iterrows():
        src_tag = 'known' if '|src=known_mature' in idx else ('novel' if '|src=novel_mature' in idx else 'untagged')
        fh.write(f"  [{src_tag}] {idx}: {int(row['Total_Counts'])} reads "
                 f"(detected in {int(row['Samples_Detected'])} samples)\n")
    fh.write("\nOutput files:\n")
    fh.write(f"  mirna_count_matrix.tsv           (DESeq2/edgeR ready)\n")
    fh.write(f"  mirna_cpm_normalized.tsv          (visualisation only)\n")
    fh.write(f"  mirna_log2cpm_normalized.tsv      (visualisation only)\n")
    fh.write(f"  mirna_gene_summary.tsv\n")
    fh.write(f"  mirna_sample_summary.tsv\n")
print(f"\n{open(stats_path).read()}")

PYTHON_MATRIX



echo ""
echo "Count matrix preparation completed!"
echo "Output directory: ${EXPR_DIR}/"
echo ""

printf "\n\nCompleted: 5.3-mirna_count_matrix\n\n"
pwd
date
