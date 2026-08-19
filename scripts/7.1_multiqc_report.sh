#!/usr/bin/env bash
#SBATCH --job-name=multiqc_report
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


###############################################################################
# 7.1_multiqc_report.sh - Comprehensive MultiQC and custom pipeline report
#
# Description:
#   Phase VII: Aggregate QC and analysis metrics from all pipeline phases into
#   a comprehensive MultiQC HTML report.
#
#   The report is built using a custom MultiQC config file and seven plugin
#   modules located in ${SCRIPT_DIR}/MultiQC_customScripts/.  Standalone PNG
#   exports are generated first so they can be base64-embedded in the report.
#
#   Report section order:
#     1.  FastQC – Raw Reads              (built-in fastqc module)
#     2.  Cutadapt – Adapter Trimming     (built-in cutadapt module)
#     3.  fastp – Quality Filtering       (built-in fastp module)
#     4.  FastQC – Cleaned Reads          (built-in fastqc, 2nd pass)
#     5.  Contamination Filtering         (custom: contamination_filter)
#     6.  miRTrace                        (built-in mirtrace module)
#     7.  Mapping to Genome               (custom: genome_alignment_custom)
#     8.  miRDeep2 Analysis               (custom: mirdeep2_custom)
#     9.  miRNA Detection                 (custom: mirna_detection)
#     10. snRNA Detection                 (custom: snrna_detection)
#     11. snoRNA Detection                (custom: snorna_detection)
#     12. piRNA Detection                 (custom: pirna_detection)
#
#   Standalone PNG exports written to ${OUTDIR}/multiqc/:
#     contamination_breakdown.png
#     mirdeep2_known_novel.png
#     mirna_distinct_detected.png
#     snrna_distinct_detected.png
#     snorna_distinct_detected.png
#     pirna_distinct_detected.png
#
# Usage:
#   sbatch scripts/7.1_multiqc_report.sh
#
# Prerequisites:
#   - Complete Phase I  (1.1–1.6)
#   - Complete Phase III (3.1–3.3)
#   - Complete Phase IV (4.1)
#   - Complete Phase V  (5.1–5.3)
#   - Complete Phase VI (6.1–6.3)
#
# Outputs (in ${OUTDIR}/multiqc/):
#   multiqc_report.html          – MultiQC HTML report (all sections)
#   contamination_breakdown.png  – Contamination stacked bar (standalone)
#   mirdeep2_known_novel.png     – miRDeep2 known/novel bar chart (standalone)
#   mirna_distinct_detected.png  – Distinct miRNA bar chart (standalone)
#   snrna_distinct_detected.png  – Distinct snRNA bar chart (standalone)
#   snorna_distinct_detected.png – Distinct snoRNA bar chart (standalone)
#   pirna_distinct_detected.png  – Distinct piRNA bar chart (standalone)
#   read_tracking.tsv            – Per-sample read counts at each pipeline stage
#   pipeline_summary.txt         – Plain-text summary of all metrics
###############################################################################

printf "\n\nStarted Job: 7.1-multiqc_report\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used across all pipeline phases)
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Sample manifest
MANIFEST="/storage/users/ahinsu/DD-microRNA/sample_manifest.txt"

# Conda environment name
CONDA_ENV="smallRNA"

# Absolute path to the scripts directory that contains MultiQC_customScripts/
# (set this to the directory where your scripts live on the execution system)
SCRIPT_DIR="/storage/users/ahinsu/DD-microRNA/scripts"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

mkdir -p "${OUTDIR}/multiqc"

# ============================================================================
# Resolve path to MultiQC custom scripts and plugins
# ============================================================================
MQC_SCRIPTS="${SCRIPT_DIR}/MultiQC_customScripts"

if [ ! -d "${MQC_SCRIPTS}" ]; then
    echo "ERROR: MultiQC custom scripts directory not found: ${MQC_SCRIPTS}"
    echo "       Set SCRIPT_DIR correctly in the USER CONFIGURATION section above."
    exit 1
fi

echo "MultiQC custom scripts : ${MQC_SCRIPTS}"

# ============================================================================
# Step 1: Generate standalone Python export plots
# ============================================================================
echo ""
echo "Step 1: Generating standalone export plots..."

export OUTDIR
export MANIFEST

python3 "${MQC_SCRIPTS}/plot_contamination.py"   && echo "  [OK] contamination_breakdown.png" \
    || echo "  [WARN] plot_contamination.py failed (continuing)"

python3 "${MQC_SCRIPTS}/plot_mirdeep2.py"         && echo "  [OK] mirdeep2_known_novel.png" \
    || echo "  [WARN] plot_mirdeep2.py failed (continuing)"

python3 "${MQC_SCRIPTS}/plot_mirna_detection.py"  && echo "  [OK] mirna_distinct_detected.png" \
    || echo "  [WARN] plot_mirna_detection.py failed (continuing)"

python3 "${MQC_SCRIPTS}/plot_snrna_detection.py"  && echo "  [OK] snrna_distinct_detected.png" \
    || echo "  [WARN] plot_snrna_detection.py failed (continuing)"

python3 "${MQC_SCRIPTS}/plot_snorna_detection.py" && echo "  [OK] snorna_distinct_detected.png" \
    || echo "  [WARN] plot_snorna_detection.py failed (continuing)"

python3 "${MQC_SCRIPTS}/plot_pirna_detection.py"  && echo "  [OK] pirna_distinct_detected.png" \
    || echo "  [WARN] plot_pirna_detection.py failed (continuing)"

# ============================================================================
# Step 2: Install custom MultiQC plugin modules, then run MultiQC
# ============================================================================
echo ""
echo "Step 2: Installing custom MultiQC plugin modules..."

# Generate a sample-name replacement file so that raw FastQC reports (which are
# named after the original FASTQ file) appear in MultiQC with the sample name
# defined in the manifest (column 1 = sample name, column 2 = FASTQ path).
REPLACE_NAMES="${OUTDIR}/multiqc/replace_names.tsv"
python3 << 'PYEOF'
import os, sys

manifest = os.environ.get('MANIFEST', '')
outdir   = os.environ.get('OUTDIR', '')
outfile  = os.path.join(outdir, 'multiqc', 'replace_names.tsv')

os.makedirs(os.path.dirname(outfile), exist_ok=True)

if not manifest or not os.path.exists(manifest):
    print("WARNING: MANIFEST not set or not found – skipping replace-names generation.", file=sys.stderr)
    sys.exit(0)

written = 0
header_skipped = False
seen = set()  # avoid duplicate replacement entries
with open(manifest) as fh, open(outfile, 'w') as out:
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if not header_skipped:
            header_skipped = True
            continue  # skip the first non-comment line (header row)
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        sample_name = parts[0].strip()
        fastq_path  = parts[1].strip()
        fastq_bn    = os.path.basename(fastq_path)
        # Strip common FASTQ extensions (longest first) to get the name FastQC would use
        for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq', '.gz'):
            if fastq_bn.endswith(ext):
                fastq_bn = fastq_bn[:-len(ext)]
                break

        def _write(old, new):
            if old and old != new and old not in seen:
                out.write(f'{old}\t{new}\n')
                seen.add(old)
                return 1
            return 0

        # Map the raw FASTQ basename → sample name (for FastQC raw, which is
        # named after the input FASTQ file, not the pipeline sample name)
        written += _write(fastq_bn, sample_name)

        # Map intermediate basenames that may appear if pipeline output files
        # were created using the FASTQ basename instead of the sample name
        # (e.g. from a previous pipeline run with a different naming scheme).
        # These cover cutadapt logs, fastp JSONs, and cleaned FastQC zips
        # that could be named by the original FASTQ basename.

print(f"replace_names.tsv: {written} entries written to {outfile}")
PYEOF
echo "  [OK] replace_names.tsv generated"

# MultiQC v1.21+ removed --module-dir; custom modules must be registered as
# Python entry points (multiqc.modules.v1).  An editable install of the local
# pyproject.toml in MultiQC_customScripts/ registers all seven plugins without
# touching any other installed packages.
pip install --quiet --no-deps -e "${MQC_SCRIPTS}" \
    && echo "  [OK] multiqc-bovine-plugins installed" \
    || { echo "ERROR: pip install of custom MultiQC plugins failed"; exit 1; }

echo ""
echo "Running MultiQC..."

multiqc \
    --config   "${MQC_SCRIPTS}/multiqc_config.yaml" \
    --module fastqc \
    --module cutadapt \
    --module fastp \
    --module mirtrace \
    --module contamination_filter \
    --module genome_alignment_custom \
    --module mirdeep2_custom \
    --module mirna_detection \
    --module snrna_detection \
    --module snorna_detection \
    --module pirna_detection \
    --replace-names "${REPLACE_NAMES}" \
    --outdir   "${OUTDIR}/multiqc" \
    --filename multiqc_report_$(date +%Y%m%d_%H%M) \
    --title    "microRNA Pipeline – Comprehensive QC Report" \
    --force \
    --interactive \
    "${OUTDIR}"

echo "MultiQC report: ${OUTDIR}/multiqc/multiqc_report.html"

# ============================================================================
# Step 3: Read tracking table and plain-text pipeline summary
# ============================================================================
echo ""
echo "Step 3: Generating read tracking table and pipeline summary..."

python3 << 'PYTHON_SUMMARY'
import os
import re
import glob
import json
import sys
from datetime import datetime
from collections import defaultdict

outdir   = os.environ['OUTDIR']
manifest = os.environ['MANIFEST']
mqc_dir  = os.path.join(outdir, 'multiqc')
os.makedirs(mqc_dir, exist_ok=True)

# ── Load sample list ─────────────────────────────────────────────────────────
samples = []
if os.path.exists(manifest):
    with open(manifest) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#') or line.lower().startswith('sample'):
                continue
            samples.append(line.split('\t')[0])

if not samples:
    for f in sorted(glob.glob(
            os.path.join(outdir, 'contamination_filtered', '*_contamination_stats.txt'))):
        s = os.path.basename(f).replace('_contamination_stats.txt', '')
        if s not in samples:
            samples.append(s)

if not samples:
    print("No samples found – skipping read tracking.", file=sys.stderr)
    sys.exit(0)

print(f"Building read tracking for {len(samples)} samples...")

# ── Helper ────────────────────────────────────────────────────────────────────
def first_int(pattern, text, default=0):
    m = re.search(pattern, text)
    if m:
        try:
            return int(m.group(1).replace(',', ''))
        except ValueError:
            pass
    return default

# ── Collect per-sample counts ─────────────────────────────────────────────────
RNA_TYPES = ['rRNA', 'tRNA', 'snRNA', 'snoRNA', 'piRNA', 'other']
track = defaultdict(dict)

for s in samples:
    # Raw reads from cutadapt log
    raw = 0
    raw_log = os.path.join(outdir, 'cutadapt', f'{s}_cutadapt.log')
    if os.path.exists(raw_log):
        with open(raw_log) as fh:
            raw = first_int(r'Total reads processed:\s+([\d,]+)', fh.read())

    # QC reads from fastp JSON
    qc = 0
    fastp_json = os.path.join(outdir, 'fastp', f'{s}_fastp.json')
    if os.path.exists(fastp_json):
        try:
            with open(fastp_json) as fh:
                jdata = json.load(fh)
            qc = jdata.get('summary', {}).get('after_filtering', {}).get('total_reads', 0)
        except Exception:
            pass

    # Contamination stage counts
    cs_file = os.path.join(outdir, 'contamination_filtered', f'{s}_contamination_stats.txt')
    contam = {rt: 0 for rt in RNA_TYPES}
    clean = 0
    if os.path.exists(cs_file):
        with open(cs_file) as fh:
            content = fh.read()
        clean = first_int(r'Output Reads \(Clean\):\s+([\d,]+)', content)
        for rt in RNA_TYPES:
            contam[rt] = first_int(rf'{rt}:\s+([\d,]+)\s*\(', content)

    # miRNA-mapped reads
    mirna_stats = os.path.join(outdir, 'mirna_counts', f'{s}_mirna_stats.txt')
    mirna_mapped = 0
    if os.path.exists(mirna_stats):
        with open(mirna_stats, encoding='utf-8') as fh:
            mirna_mapped = first_int(
                r'Mapped to miRNA reference\s*:\s*([\d,]+)',
                fh.read()
            )

    track[s] = {
        'raw_reads':    raw,
        'qc_reads':     qc,
        **{rt: contam[rt] for rt in RNA_TYPES},
        'clean':        clean,
        'mirna_mapped': mirna_mapped,
    }

# ── Write TSV ─────────────────────────────────────────────────────────────────
tsv_path = os.path.join(mqc_dir, 'read_tracking.tsv')
cols = ['sample', 'raw_reads', 'qc_reads'] + RNA_TYPES + ['clean', 'mirna_mapped']
with open(tsv_path, 'w') as fh:
    fh.write('\t'.join(cols) + '\n')
    for s in samples:
        row = [s] + [str(track[s].get(c, 0)) for c in cols[1:]]
        fh.write('\t'.join(row) + '\n')
print(f"Read tracking table: {tsv_path}")

# ── Write plain-text pipeline summary ─────────────────────────────────────────
def count_seqs(fpath):
    if not os.path.exists(fpath):
        return 0
    return sum(1 for line in open(fpath) if line.startswith('>'))

summary_path = os.path.join(mqc_dir, 'pipeline_summary.txt')
with open(summary_path, 'w') as fh:
    fh.write('=' * 75 + '\n')
    fh.write('microRNA Pipeline – Comprehensive Summary\n')
    fh.write(f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
    fh.write('=' * 75 + '\n\n')
    fh.write(f'Samples processed : {len(samples)}\n\n')

    ref_dir = os.path.join(os.path.dirname(outdir), 'reference')
    for label, rpath in [
        ('Unified mature miRNA ref',
         os.path.join(ref_dir, 'mirna_quant_ref', 'unified_mature.fa')),
        ('Combined mature (dedup)',
         os.path.join(ref_dir, 'db_unified', 'combined_mature_dedup.fasta')),
    ]:
        n = count_seqs(rpath)
        fh.write(f'  {label:<40} {n:>6} sequences\n')
    fh.write('\n')

    fh.write('Read tracking (first 5 samples):\n')
    for s in samples[:5]:
        fh.write(f'  {s}: raw={track[s].get("raw_reads", 0):,}  '
                 f'clean={track[s].get("clean", 0):,}  '
                 f'mirna_mapped={track[s].get("mirna_mapped", 0):,}\n')
    if len(samples) > 5:
        fh.write(f'  ... ({len(samples) - 5} more samples in read_tracking.tsv)\n')
    fh.write('\n')

    fh.write('Output plots:\n')
    for png in sorted(glob.glob(os.path.join(mqc_dir, '*.png'))):
        fh.write(f'  {png}\n')

print(f"Pipeline summary: {summary_path}")
print("Step 3 complete.")
PYTHON_SUMMARY

echo ""
echo "============================================================"
echo "Phase VII – MultiQC report complete"
echo "============================================================"
echo ""
echo "Outputs:"
echo "  MultiQC HTML          : ${OUTDIR}/multiqc/multiqc_report.html"
echo "  Contamination plot    : ${OUTDIR}/multiqc/contamination_breakdown.png"
echo "  miRDeep2 plot         : ${OUTDIR}/multiqc/mirdeep2_known_novel.png"
echo "  miRNA distinct plot   : ${OUTDIR}/multiqc/mirna_distinct_detected.png"
echo "  snRNA distinct plot   : ${OUTDIR}/multiqc/snrna_distinct_detected.png"
echo "  snoRNA distinct plot  : ${OUTDIR}/multiqc/snorna_distinct_detected.png"
echo "  piRNA distinct plot   : ${OUTDIR}/multiqc/pirna_distinct_detected.png"
echo "  Read tracking table   : ${OUTDIR}/multiqc/read_tracking.tsv"
echo "  Pipeline summary      : ${OUTDIR}/multiqc/pipeline_summary.txt"

printf "\n\nCompleted: 7.1-multiqc_report\n\n"
pwd
date
