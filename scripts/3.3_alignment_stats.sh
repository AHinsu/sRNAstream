#!/usr/bin/env bash
#SBATCH --job-name=alignment_stats
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G


###############################################################################
# 3.3_alignment_stats.sh - Aggregate Bowtie alignment stats from all samples
#
# Description:
#   Phase III.3: Collect the per-sample Bowtie alignment logs produced by
#   3.2_genome_alignment.sh and combine them into a single summary TSV.
#
# Usage:
#   # Submit immediately after the array job completes:
#   ARRAY_JOB_ID=$(sbatch --parsable scripts/3.2_genome_alignment.sh)
#   sbatch --dependency=afterok:${ARRAY_JOB_ID} scripts/3.3_alignment_stats.sh
#
#   # Or run manually once 3.2 has finished:
#   sbatch scripts/3.3_alignment_stats.sh
#
# Prerequisites:
#   - 3.2_genome_alignment.sh must have completed for ALL samples.
#
# Outputs (in ${OUTDIR}/alignment/):
#   all_samples_alignment_stats.tsv  – Combined per-sample mapping statistics
###############################################################################

printf "\n\nStarted Job: 3.3-alignment_stats\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 3.2_genome_alignment.sh)
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Sample manifest (TSV with at least a 'sample_name' column in column 1)
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

# ============================================================================
# Aggregate per-sample Bowtie logs into a combined TSV
# ============================================================================

echo "Aggregating alignment statistics from: ${OUTDIR}/alignment/"

export OUTDIR
export MANIFEST

python3 << 'PYTHON_AGGREGATE'
import os
import re
import glob

outdir    = os.environ['OUTDIR']
manifest  = os.environ['MANIFEST']
align_dir = os.path.join(outdir, 'alignment')
combined  = os.path.join(align_dir, 'all_samples_alignment_stats.tsv')

# -------------------------------------------------------------------------
# Load sample list from manifest (preserves order)
# -------------------------------------------------------------------------
samples = []
if os.path.exists(manifest):
    with open(manifest) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#') or line.lower().startswith('sample'):
                continue
            samples.append(line.split('\t')[0])

# Fall back: discover samples from existing Bowtie log files
if not samples:
    print("WARNING: manifest not found or empty; discovering samples from log files.")
    for log in sorted(glob.glob(os.path.join(align_dir, '*_alignment.log'))):
        s = os.path.basename(log).replace('_alignment.log', '')
        samples.append(s)

if not samples:
    print("ERROR: No samples found. Run 3.2_genome_alignment.sh first.")
    raise SystemExit(1)

print(f"Processing {len(samples)} samples...")

# -------------------------------------------------------------------------
# Parse each sample's Bowtie log
# -------------------------------------------------------------------------
rows = []
missing = []

for sample in samples:
    log_file = os.path.join(align_dir, f'{sample}_alignment.log')
    if not os.path.exists(log_file):
        print(f"  WARNING: log not found for {sample}: {log_file}")
        missing.append(sample)
        continue

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

    pct_mapped = round(mapped / total * 100, 2) if total > 0 else 0.0
    rows.append((sample, total, mapped, unmapped, pct_mapped, reported))
    print(f"  {sample}: {total} reads, {mapped} mapped ({pct_mapped}%)")

# -------------------------------------------------------------------------
# Write combined TSV (written once in a single sequential pass – no concurrent access)
# -------------------------------------------------------------------------
with open(combined, 'w') as fh:
    fh.write("sample\ttotal_reads\tmapped_reads\tunmapped_reads\tpct_mapped\treported_alignments\n")
    for row in rows:
        fh.write('\t'.join(str(x) for x in row) + '\n')

print(f"\nWritten: {combined}")
print(f"  Samples aggregated : {len(rows)}")
if missing:
    print(f"  Missing log files  : {len(missing)} ({', '.join(missing)})")

# -------------------------------------------------------------------------
# Print summary statistics
# -------------------------------------------------------------------------
if rows:
    totals   = [r[1] for r in rows]
    pcts     = [r[4] for r in rows]
    print(f"\nSummary:")
    print(f"  Mean total reads   : {sum(totals)/len(totals):,.0f}")
    print(f"  Mean mapping rate  : {sum(pcts)/len(pcts):.2f}%")
    print(f"  Min mapping rate   : {min(pcts):.2f}%")
    print(f"  Max mapping rate   : {max(pcts):.2f}%")
PYTHON_AGGREGATE

echo ""
echo "Stats aggregation completed."
echo "Combined table: ${OUTDIR}/alignment/all_samples_alignment_stats.tsv"

printf "\n\nCompleted: 3.3-alignment_stats\n\n"
pwd
date
