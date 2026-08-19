#!/usr/bin/env bash
#SBATCH --job-name=mirtrace_qc
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G


###############################################################################
# 1.6_qc_reports.sh - Comprehensive Phase I QC reports
# 
# Description:
#   Performs comprehensive miRNA quality control using miRTrace
#   and generates Phase I QC summary reports:
#   - miRNA detection and quantification
#   - RNA type distribution (miRNA, rRNA, tRNA, etc.)
#   - Complexity analysis
#   - Contamination detection
#   - Phase I Read Count Summary across all processing steps
#
# Usage: 
#   sbatch scripts/1.6_qc_reports.sh
#
# Outputs:
#   - miRTrace QC reports in ${OUTDIR}/mirtrace/
#   - HTML visualization and summary statistics
#
# Note: This is NOT an array job - processes all samples together
###############################################################################

printf "\n\nStarted Job: 1.6-qc_reports\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"
export OUTDIR

# Conda environment name
CONDA_ENV="smallRNA"

# Species for miRTrace (use 'unc' for unknown/custom species)
# Options: hsa (human), mmu (mouse), rno (rat), bta (cow), etc.
# Use 'unc' if your species is not in miRTrace's built-in database
SPECIES="bta"  # bta = Bos taurus (cattle) #Run mirtrace --list-species for list

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH


# Create output directory
mkdir -p "${OUTDIR}/mirtrace"

# ============================================================================
# Run miRTrace on all samples
# ============================================================================

echo "Running miRTrace QC on all samples"

# Create config file for miRTrace (CSV)
CONFIG_FILE="${OUTDIR}/mirtrace/mirtrace_config.csv"
> "$CONFIG_FILE"  # Clear file if exists

# Add all trimmed FASTQ files to config
for FASTQ_FILE in "${OUTDIR}"/contamination_filtered/*_cleaned.fq.gz; do
    if [ -f "$FASTQ_FILE" ]; then
        SAMPLE_NAME=$(basename "$FASTQ_FILE" _cleaned.fq.gz)
        # filename (full path),sampleName
        echo "${FASTQ_FILE},${SAMPLE_NAME}" >> "$CONFIG_FILE"
    fi
done

echo "Created config file: $CONFIG_FILE"
head $CONFIG_FILE

# Run miRTrace
mirtrace qc \
    --species ${SPECIES} \
    --protocol illumina --force \
    --output-dir ${OUTDIR}/mirtrace \
    --config $CONFIG_FILE \
    --num-threads ${SLURM_CPUS_PER_TASK}

echo "miRTrace completed"
echo "Results available in: ${OUTDIR}/mirtrace/"

# ============================================================================
# Generate Statistics
# ============================================================================

echo ""
echo "Generating miRTrace statistics..."

python3 << 'PYTHON_STATS'
import os
import json
import glob

outdir = os.environ.get('OUTDIR', '/storage/users/ahinsu/DD-microRNA/results')
mirtrace_dir = os.path.join(outdir, 'mirtrace')
stats_file = os.path.join(mirtrace_dir, 'mirtrace_stats_summary.txt')

print("\n" + "="*80)
print("miRTrace QC Statistics")
print("="*80 + "\n")

# Parse miRTrace QC summary if available
summary_file = os.path.join(mirtrace_dir, 'mirtrace-stats-contamination_basic.tsv')
results_summary = []

if os.path.exists(summary_file):
    with open(summary_file, 'r') as f:
        lines = f.readlines()
        if len(lines) > 1:
            header = lines[0].strip().split('\t')
            for line in lines[1:]:
                values = line.strip().split('\t')
                if len(values) > 0:
                    results_summary.append(dict(zip(header, values)))

# Count total samples processed
config_file = os.path.join(mirtrace_dir, 'mirtrace_config.txt')
n_samples = 0
if os.path.exists(config_file):
    with open(config_file, 'r') as f:
        n_samples = len(f.readlines())

with open(stats_file, 'w') as f:
    f.write("miRTrace Quality Control Summary\n")
    f.write("="*80 + "\n\n")
    f.write(f"Total Samples Processed: {n_samples}\n")
    
    if results_summary:
        f.write(f"\nContamination Analysis:\n")
        f.write(f"Number of samples analyzed: {len(results_summary)}\n")
        
        # If QC passed for samples
        passed = sum(1 for s in results_summary if s.get('QC', '').lower() in ['pass', 'good', 'ok'])
        f.write(f"Samples passing QC: {passed}\n")
    
    f.write(f"\nResults location: {mirtrace_dir}\n")
    f.write(f"HTML report: {os.path.join(mirtrace_dir, 'mirtrace-report.html')}\n")
    
    print(open(stats_file, 'r').read())

print("="*80 + "\n")
PYTHON_STATS

# Generate Phase I Read Count Summary
echo ""
echo "Generating Phase I comprehensive QC reports..."

python3 << 'PYTHON_REPORTS'
import os, glob

outdir = os.environ.get('OUTDIR', '')
reports_dir = os.path.join(outdir, 'reports')
os.makedirs(reports_dir, exist_ok=True)

# Read Count Summary Table
summary_file = os.path.join(reports_dir, 'phase1_read_count_summary.tsv')
contam_summary_file = os.path.join(reports_dir, 'phase1_contamination_summary.tsv')

# RNA types filtered sequentially in 1.5b (order matters for computing cumulative counts)
RNA_TYPES = ['rRNA', 'tRNA', 'snRNA', 'snoRNA', 'piRNA', 'other']

# Collect per-sample statistics
samples = {}

# Parse raw read counts from cutadapt logs (reads before any processing)
cutadapt_dir = os.path.join(outdir, 'cutadapt')
for log_file in glob.glob(os.path.join(cutadapt_dir, '*_cutadapt.log')):
    sample = os.path.basename(log_file).replace('_cutadapt.log', '')
    samples.setdefault(sample, {'contam_counts': {}})
    with open(log_file) as f:
        for line in f:
            if 'Total reads processed:' in line:
                value = line.split(':')[1].strip().replace(',', '')
                samples[sample]['raw_reads'] = int(value) if value else 0
                break

# Collect contamination stats per sample from files generated by 1.5b
contam_dir = os.path.join(outdir, 'contamination_filtered')
for stats_file in glob.glob(os.path.join(contam_dir, '*_contamination_stats.txt')):
    sample = os.path.basename(stats_file).replace('_contamination_stats.txt', '')
    samples.setdefault(sample, {'contam_counts': {}})
    with open(stats_file) as f:
        for line in f:
            line = line.strip()
            if 'Input Reads:' in line:
                value = line.split(':')[1].strip().replace(',', '')
                samples[sample]['post_trim'] = int(value) if value else 0
            elif 'Output Reads (Clean):' in line:
                value = line.split(':')[1].strip().replace(',', '')
                samples[sample]['final_clean'] = int(value) if value else 0
            else:
                # Parse contamination type lines: "  rRNA: 12345 (10.00%)"
                for rna_type in RNA_TYPES:
                    if line.startswith(rna_type + ':') and '(' in line:
                        value = line.split(':')[1].split('(')[0].strip().replace(',', '')
                        samples[sample]['contam_counts'][rna_type] = int(value) if value else 0
                        break

    # Compute cumulative post-filter read counts for each RNA type in sequential order
    running = samples[sample].get('post_trim', 0)
    for rna_type in RNA_TYPES:
        count = samples[sample]['contam_counts'].get(rna_type, 0)
        running -= count
        samples[sample][f'post_{rna_type.lower()}'] = running

# Write Read Count Summary Table (one column per filtering step)
with open(summary_file, 'w') as f:
    header = 'Sample\tRaw_Reads\tPost_Trim\tPost_rRNA\tPost_tRNA\tPost_snRNA\tPost_snoRNA\tPost_piRNA\tPost_other\tFinal_Clean\tFinal_Clean_pct\n'
    f.write(header)
    print(header.strip())
    for sample, stats in sorted(samples.items()):
        raw = stats.get('raw_reads', 'NA')
        final_clean = stats.get('final_clean', 'NA')
        if isinstance(raw, int) and raw > 0 and isinstance(final_clean, int):
            final_clean_pct = f"{(final_clean / raw * 100):.2f}"
        else:
            final_clean_pct = 'NA'
        row = '\t'.join([
            sample,
            str(raw),
            str(stats.get('post_trim', 'NA')),
            str(stats.get('post_rrna', 'NA')),
            str(stats.get('post_trna', 'NA')),
            str(stats.get('post_snrna', 'NA')),
            str(stats.get('post_snorna', 'NA')),
            str(stats.get('post_pirna', 'NA')),
            str(stats.get('post_other', 'NA')),
            str(final_clean),
            final_clean_pct,
        ]) + '\n'
        f.write(row)
        print(row.strip())

# Write Contamination Summary Table (per-type removed counts + overall rate)
with open(contam_summary_file, 'w') as f:
    header = 'Sample\t' + '\t'.join(RNA_TYPES) + '\tTotal_Removed\tFinal_Clean\tContamination_Rate_pct\n'
    f.write(header)
    for sample, stats in sorted(samples.items()):
        post_trim = stats.get('post_trim', 0)
        final_clean = stats.get('final_clean', 0)
        contam_counts = stats.get('contam_counts', {})
        total_removed = sum(contam_counts.values())
        contam_rate = f"{(total_removed / post_trim * 100):.2f}" if post_trim > 0 else 'NA'
        row = '\t'.join([
            sample,
            *[str(contam_counts.get(t, 0)) for t in RNA_TYPES],
            str(total_removed),
            str(final_clean),
            contam_rate,
        ]) + '\n'
        f.write(row)

print(f"\nPhase I QC reports written to: {reports_dir}")
PYTHON_REPORTS

printf "\n\nCompleted\n\n"
pwd
date
