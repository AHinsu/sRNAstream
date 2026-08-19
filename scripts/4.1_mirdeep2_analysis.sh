#!/usr/bin/env bash
#SBATCH --job-name=mirdeep2_combined
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G


###############################################################################
# 4.1_mirdeep2_analysis.sh - Combined miRDeep2 Analysis + Result Parsing
#
# Description:
#   Phase IV: Run miRDeep2 on reads from all samples combined to identify
#   both known and novel miRNAs, then parse the results into structured tables
#   and extract per-tier FASTA files.
#
#
#   Workflow:
#     1. Collect genome-mapped FASTQ files from all samples (3.2 output) and
#        concatenate them into one combined FASTA.
#     2. Map the reads to the reference genome with mapper.pl (miRDeep2) and
#        produce the ARF alignment file.
#     3. Run miRDeep2.pl with the combined reads, genome, known mature and
#        precursor FASTAs.
#     4. Parse the miRDeep2 result CSV (Python, runs once):
#          - Correctly maps all 17 CSV columns using the canonical miRDeep2
#            column order.
#          - Classifies predictions as "known" or "novel".
#          - Applies score thresholds (SCORE_LOW / SCORE_MED / SCORE_HIGH) and
#            reports counts, Rfam alerts and randfold significance at each tier.
#          - Writes per-tier TSV tables to mirdeep2/parsed/.
#     5. Extract mature and precursor FASTA files for each score tier using
#        get_mirdeep2_precursors.pl.  Output directories:
#          mirdeep2/sequences${SCORE_LOW}/
#          mirdeep2/sequences${SCORE_MED}/
#          mirdeep2/sequences${SCORE_HIGH}/
#        Each directory contains novel_mature_*.fa and known_mature_*.fa (and
#        corresponding precursor files) that are consumed by 5.1.
#
# Usage:
#   sbatch scripts/4.1_mirdeep2_analysis.sh
#
# Prerequisites:
#   - Run 3.2_genome_alignment.sh for ALL samples first.
#   - Bowtie genome index must exist (same one used in 3.1).
#
# Key outputs (in ${OUTDIR}/mirdeep2/):
#   combined_reads.fasta             – concatenated mapped reads (FASTA)
#   collapsed_reads.fasta            – read-collapsed FASTA for miRDeep2
#   combined_mapped.arf              – miRDeep2 ARF alignment file
#   result_*.csv                     – miRDeep2 predictions (known + novel)
#   miRNAs_expressed_*.csv           – quantification table
#   pdfs/                            – PDF figures from miRDeep2
#   combined_mirdeep2.log            – miRDeep2 run log
#   parsed/mirdeep2_summary_stats.tsv– Scoring-tier summary table (21 rows × 9 cols)
#   parsed/mirdeep2_all.tsv          – All predictions (annotated)
#   parsed/mirdeep2_known.tsv        – Known miRNA subset
#   parsed/mirdeep2_novel.tsv        – Novel miRNA subset
#   parsed/mirdeep2_high_conf.tsv    – High-confidence subset (>= SCORE_HIGH)
#   parsed/mirdeep2_parsing_stats.txt– Detailed parsing statistics report
#   sequences${SCORE}/               – Per-tier mature/precursor FASTAs
###############################################################################

printf "\n\nStarted Job: 4.1-mirdeep2_combined\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 3.2_genome_alignment.sh)
OUTDIR="/storage/users/ahinsu/DD-microRNA/results"

# Reference genome FASTA file (from Ensembl)
GENOME_FASTA="/storage/users/ahinsu/DD-microRNA/reference/Bos_taurus.chromsomal.fa"

# Bowtie index for the Ensembl genome (without .1.ebwt suffix)
# Build with: bowtie-build /path/to/genome.fa /path/to/index/genome
GENOME_INDEX="/storage/users/ahinsu/DD-microRNA/reference/index/genome"

# High-confidence known mature miRNA FASTA (for miRDeep2 -d flag)
# Supply sequences from miRBase, MirGeneDB or the db_unified output from 2.2.
# Set to "" to run without a known reference (miRDeep2 will still predict novel miRNAs).
# NOTE: this same file is used by 5.1_mirna_reference.sh as the "known" source.
KNOWN_MATURE_FA="/storage/users/ahinsu/DD-microRNA/reference/db_unified/combined_mature_dedup.fasta"

# Known precursor miRNA FASTA (for miRDeep2 -p flag)
# Supply sequences from miRBase, MirGeneDB or the db_unified output from 2.2.
# Set to "" to run without precursor reference.
KNOWN_PRECURSOR_FA="/storage/users/ahinsu/DD-microRNA/reference/db_unified/combined_precursor_dedup.fasta"

# Related-species mature miRNA FASTA (optional; leave "" to skip)
# miRDeep2 uses this to boost confidence in novel predictions.
RELATED_MATURE_FA=""

# miRDeep2 score thresholds for result parsing
# Predictions with score below SCORE_LOW are excluded from output tables.
SCORE_LOW=1        # Minimum score to include a prediction at all
SCORE_MED=5        # Medium confidence
SCORE_HIGH=10      # High confidence

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

MIRDEEP2_DIR="${OUTDIR}/mirdeep2"
mkdir -p "${MIRDEEP2_DIR}"
cd "${MIRDEEP2_DIR}"

# ============================================================================
# Step 1: Concatenate mapped reads from all samples
# ============================================================================

echo "Step 1: Collecting mapped reads from all samples..."

COMBINED_FQ="${MIRDEEP2_DIR}/combined_reads.fastq"
COMBINED_FA="${MIRDEEP2_DIR}/combined_reads.fasta"

rm -f "$COMBINED_FQ" "$COMBINED_FA"

# Read sample names from manifest (skip header line)
n_samples=0
while IFS=$'\t' read -r sample_name rest; do
    [[ "$sample_name" == "sample" || "$sample_name" =~ ^# ]] && continue
    mapped_fq="${OUTDIR}/alignment/${sample_name}_mapped.fastq"
    if [ -f "$mapped_fq" ]; then
        cat "$mapped_fq" >> "$COMBINED_FQ"
        echo "  Added: ${sample_name} (${mapped_fq})"
        ((n_samples++))
    else
        echo "  WARNING: mapped reads not found for ${sample_name}: ${mapped_fq}"
    fi
done < "$MANIFEST"

if [ "$n_samples" -eq 0 ]; then
    echo "ERROR: No mapped read files found. Run 3.1_genome_alignment.sh first."
    exit 1
fi

seqkit fq2fa -j ${SLURM_CPUS_PER_TASK:-8} -o ${COMBINED_FA} "${COMBINED_FQ}"


TOTAL_READS=$(( $(wc -l < "$COMBINED_FQ") / 4 ))
echo "Combined ${n_samples} samples → ${TOTAL_READS} total reads"

rm "${COMBINED_FQ}"

awk '/^>/{sub(/ .*/, "", $0)}1' "$COMBINED_FA" > temp.fasta
mv temp.fasta "$COMBINED_FA"


# ============================================================================
# Step 2: Map collapsed reads to genome and generate ARF file
# ============================================================================

echo ""
echo "Step 2: Mapping with miRdeep2 and creating ARF format..."

COMBINED_ARF="${MIRDEEP2_DIR}/combined_mapped.arf"
COLLAPSED_FA="${MIRDEEP2_DIR}/collapsed_reads.fasta"

mapper.pl "$COMBINED_FA" -l 15 -c -q -m -v -p "${GENOME_INDEX}" \
    -s "${COLLAPSED_FA}" -t "${COMBINED_ARF}"


# ============================================================================
# Step 3: Run miRDeep2
# ============================================================================

echo ""
echo "Step 3: Running miRDeep2 on combined reads..."

###############################################
# Build argument list for miRDeep2
###############################################

# --- Known mature miRNA ---
if [[ -n "$KNOWN_MATURE_FA" && -f "$KNOWN_MATURE_FA" ]]; then
    echo "Cleaning known mature miRNA FASTA..."
    MATURE_CLEANED="matured_cleaned.fasta"
    awk '/^>/{gsub(/ /, "_")}1' "$KNOWN_MATURE_FA" > "$MATURE_CLEANED"
    MATURE_ARG="$MATURE_CLEANED"
else
    MATURE_ARG="none"
fi

# --- Known precursor miRNA ---
if [[ -n "$KNOWN_PRECURSOR_FA" && -f "$KNOWN_PRECURSOR_FA" ]]; then
    echo "Cleaning known precursor miRNA FASTA..."
    PRECURSOR_CLEANED="precursor_cleaned.fasta"
    awk '/^>/{gsub(/ /, "_")}1' "$KNOWN_PRECURSOR_FA" > "$PRECURSOR_CLEANED"
    PRECURSOR_ARG="$PRECURSOR_CLEANED"
else
    PRECURSOR_ARG="none"
fi

# --- Related mature miRNA ---
if [[ -n "$RELATED_MATURE_FA" && -f "$RELATED_MATURE_FA" ]]; then
    echo "Cleaning related mature miRNA FASTA..."
    RELATED_CLEANED="related_cleaned.fasta"
    awk '/^>/{gsub(/ /, "_")}1' "$RELATED_MATURE_FA" > "$RELATED_CLEANED"
    RELATED_ARG="$RELATED_CLEANED"
else
    RELATED_ARG="none"
fi


###############################################
# Print summary
###############################################
echo "  Genome FASTA    : ${GENOME_FASTA}"
echo "  Collapsed FASTA : ${COLLAPSED_FA}"
echo "  ARF file        : ${COMBINED_ARF}"
echo "  Known mature    : ${MATURE_ARG}"
echo "  Known precursor : ${PRECURSOR_ARG}"
echo "  Related mature  : ${RELATED_ARG}"
echo "  Score thresholds: low=${SCORE_LOW}  med=${SCORE_MED}  high=${SCORE_HIGH}"



miRDeep2.pl "${COLLAPSED_FA}" "${GENOME_FASTA}" "${COMBINED_ARF}" \
    "${MATURE_ARG}" "${RELATED_ARG}" "${PRECURSOR_ARG}" -t Cow \
    2> "${MIRDEEP2_DIR}/combined_mirdeep2.log"


echo "miRDeep2 completed. Log: ${MIRDEEP2_DIR}/combined_mirdeep2.log"

# ============================================================================
# Step 4: Parse miRDeep2 result CSV and generate statistics
# ============================================================================

echo ""
echo "Step 4: Parsing miRDeep2 results and generating statistics..."

PARSED_DIR="${MIRDEEP2_DIR}/parsed"
mkdir -p "${PARSED_DIR}"

# Count collapsed reads and ARF alignment loci for the stats report
N_COLLAPSED=$(grep -c '^>' "${COLLAPSED_FA}" 2>/dev/null || echo 0)
N_ARF_LOCI=$(wc -l < "${COMBINED_ARF}" 2>/dev/null || echo 0)

export MIRDEEP2_DIR
export PARSED_DIR
export SCORE_LOW
export SCORE_MED
export SCORE_HIGH
export N_SAMPLES="${n_samples}"
export TOTAL_COMBINED_READS="${TOTAL_READS}"
export N_COLLAPSED
export N_ARF_LOCI
export COLLAPSED_FA
export COMBINED_ARF

python3 << 'PYTHON_PARSE'
import os
import glob
from datetime import datetime

md2_dir    = os.environ['MIRDEEP2_DIR']
parsed_dir = os.environ['PARSED_DIR']
score_low  = float(os.environ.get('SCORE_LOW',  1))
score_med  = float(os.environ.get('SCORE_MED',  5))
score_high = float(os.environ.get('SCORE_HIGH', 10))
n_samples        = int(os.environ.get('N_SAMPLES', 0))
total_reads_comb = int(os.environ.get('TOTAL_COMBINED_READS', 0))
n_collapsed      = int(os.environ.get('N_COLLAPSED', 0))
n_arf_loci       = int(os.environ.get('N_ARF_LOCI', 0))

# -------------------------------------------------------------------------
# Locate the most recent miRDeep2 result CSV
# -------------------------------------------------------------------------
result_csvs = sorted(glob.glob(os.path.join(md2_dir, 'result_*.csv')))
if not result_csvs:
    print("ERROR: No miRDeep2 result CSV found in", md2_dir)
    raise SystemExit(1)
result_csv = result_csvs[-1]
print(f"Using result file: {result_csv}")

# -------------------------------------------------------------------------
# Column layout for the summary statistics table (Section 1):
#   miRDeep2 result CSVs open with a scoring-tier table that has 9 columns.
# -------------------------------------------------------------------------
SUMMARY_COLS = [
    'mirdeep2_score',                        # 0
    'novel_mirnas_reported',                 # 1
    'novel_mirnas_est_false_positives',      # 2
    'novel_mirnas_est_true_positives',       # 3
    'known_mirnas_in_species',               # 4
    'known_mirnas_in_data',                  # 5
    'known_mirnas_detected_by_mirdeep2',     # 6
    'estimated_signal_to_noise',             # 7
    'excision_gearing',                      # 8
]

# -------------------------------------------------------------------------
# Column layout for detailed predictions (Section 2, 17 columns):
#   0  provisional_id
#   1  score
#   2  estimated_probability
#   3  rfam_alert
#   4  total_reads
#   5  mature_reads
#   6  loop_reads
#   7  star_reads
#   8  significant_randfold
#   9  mirbase_miRNA
#  10  example_mirna_seed
#  11  ucsc_browser
#  12  ncbi_blastn
#  13  mature_seq
#  14  star_seq
#  15  precursor_seq
#  16  precursor_coordinate
# -------------------------------------------------------------------------
COL_NAMES = [
    'provisional_id',        # 0
    'score',                 # 1
    'estimated_probability', # 2
    'rfam_alert',            # 3
    'total_reads',           # 4
    'mature_reads',          # 5
    'loop_reads',            # 6
    'star_reads',            # 7
    'significant_randfold',  # 8
    'mirbase_miRNA',         # 9
    'example_mirna_seed',    # 10
    'ucsc_browser',          # 11
    'ncbi_blastn',           # 12
    'mature_seq',            # 13
    'star_seq',              # 14
    'precursor_seq',         # 15
    'precursor_coordinate',  # 16
]

# Tokens that indicate a header / section-marker row (not data).
# miRDeep2 uses 'provisional id' (with space) as the column header in the novel
# section and 'tag id' in the known section.  'provisional_id' / 'tag' are
# included as defensive fallbacks in case the header has been pre-processed.
PREDICTION_HEADER_TOKENS = {
    'provisional id', 'provisional_id', 'tag id', 'tag',
}

# -------------------------------------------------------------------------
# Parse the miRDeep2 result CSV in a single pass, separating the three
# sections: summary statistics, novel predictions, known predictions.
# -------------------------------------------------------------------------
with open(result_csv) as fh:
    all_lines = fh.readlines()

summary_lines = []
novel_lines   = []
known_lines   = []

# Possible parse states:
#   'summary'      – inside the scoring-tier summary table (Section 1)
#   'novel'        – inside the novel miRNA predictions block (Section 2)
#   'known'        – inside the known/mature miRNA detections block (Section 3)
#   'not_detected' – inside the "not detected" block (Section 4, not collected)
#   None           – before any recognised section
section = None

for line in all_lines:
    stripped = line.strip()
    lower    = stripped.lower()

    # ---- Detect section boundaries -----------------------------------------
    # Section 1 header: tab-delimited row whose first cell is "miRDeep2 score"
    # and which contains the expected summary column names.
    if lower.startswith('mirdeep2 score') and 'novel mirnas reported' in lower:
        section = 'summary'
        continue  # skip the header line itself

    # Section 2: novel miRNA predictions block.
    # miRDeep2 emits the marker line "novel miRNAs predicted by miRDeep2".
    elif lower == 'novel mirnas predicted by mirdeep2':
        section = 'novel'
        continue

    # Section 3: mature/known miRNA detections block.
    # miRDeep2 emits the marker line "mature miRBase miRNAs detected by miRDeep2".
    elif lower == 'mature mirbase mirnas detected by mirdeep2':
        section = 'known'
        continue

    # Section 4: miRNAs not detected – stop collecting; data has a different
    # column layout and must not be mixed with the known-miRNA rows.
    # The marker line starts with '#' and contains "not detected by miRDeep2".
    elif stripped.startswith('#') and 'not detected by mirdeep2' in lower:
        section = 'not_detected'
        continue

    # ---- Distribute lines to the appropriate bucket ------------------------
    if section == 'summary':
        summary_lines.append(line)
    elif section == 'novel':
        novel_lines.append(line)
    elif section == 'known':
        known_lines.append(line)
    # 'not_detected' lines are intentionally discarded

# -------------------------------------------------------------------------
# Parse Section 1: Summary statistics table
# -------------------------------------------------------------------------
def parse_summary(lines):
    rows = []
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        parts = [p.strip() for p in line.split('\t')]
        # Each data row starts with a numeric score value
        try:
            float(parts[0])
        except (ValueError, IndexError):
            continue  # skip header or non-numeric rows
        row = {}
        for i, col in enumerate(SUMMARY_COLS):
            row[col] = parts[i].strip() if i < len(parts) else ''
        rows.append(row)
    return rows

# -------------------------------------------------------------------------
# Parse Section 2 (novel/known): Detailed predictions
# -------------------------------------------------------------------------
def parse_predictions(lines, label):
    rows = []
    for line in lines:
        line = line.rstrip('\n')
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        parts = line.split('\t')
        if not parts:
            continue
        # Skip column-header rows (first cell matches known header tokens)
        first_cell = parts[0].strip().lower()
        if first_cell in PREDICTION_HEADER_TOKENS:
            continue
        # A valid data row must have at least the first 9 columns
        if len(parts) < 9:
            continue
        row = {'category': label}
        for i, col in enumerate(COL_NAMES):
            row[col] = parts[i].strip() if i < len(parts) else ''
        try:
            row['score'] = float(row['score'])
        except (ValueError, KeyError):
            row['score'] = 0.0
        try:
            row['total_reads'] = int(float(row['total_reads']))
        except (ValueError, KeyError):
            row['total_reads'] = 0
        try:
            row['mature_reads'] = int(float(row['mature_reads']))
        except (ValueError, KeyError):
            row['mature_reads'] = 0
        try:
            row['loop_reads'] = int(float(row['loop_reads']))
        except (ValueError, KeyError):
            row['loop_reads'] = 0
        try:
            row['star_reads'] = int(float(row['star_reads']))
        except (ValueError, KeyError):
            row['star_reads'] = 0
        rows.append(row)
    return rows

summary_rows = parse_summary(summary_lines)
novel_rows   = parse_predictions(novel_lines, 'novel')
known_rows   = parse_predictions(known_lines, 'known')

# Fallback: if no section markers detected treat everything as novel
if not novel_rows and not known_rows:
    print("Warning: could not detect novel/known sections; treating all rows as novel.")
    novel_rows = parse_predictions(all_lines, 'novel')

all_rows = novel_rows + known_rows

print(f"Parsed {len(summary_rows)} summary statistics rows")
print(f"Parsed {len(novel_rows)} novel predictions and {len(known_rows)} known miRNAs")

# -------------------------------------------------------------------------
# Write TSV tables
# -------------------------------------------------------------------------
def write_tsv(rows, cols, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r.get(c, '')) for c in cols) + '\n')

PREDICTION_COLS = COL_NAMES + ['category']

write_tsv(summary_rows, SUMMARY_COLS,
          os.path.join(parsed_dir, 'mirdeep2_summary_stats.tsv'))
write_tsv(all_rows,     PREDICTION_COLS,
          os.path.join(parsed_dir, 'mirdeep2_all.tsv'))
write_tsv(known_rows,   PREDICTION_COLS,
          os.path.join(parsed_dir, 'mirdeep2_known.tsv'))
write_tsv(novel_rows,   PREDICTION_COLS,
          os.path.join(parsed_dir, 'mirdeep2_novel.tsv'))

high_conf_rows = [r for r in all_rows if r['score'] >= score_high]
write_tsv(high_conf_rows, PREDICTION_COLS,
          os.path.join(parsed_dir, 'mirdeep2_high_conf.tsv'))

print(f"Written TSV tables to {parsed_dir}/")

# -------------------------------------------------------------------------
# Statistics
# -------------------------------------------------------------------------
def tier_stats(rows, threshold):
    above  = [r for r in rows if r['score'] >= threshold]
    known  = [r for r in above if r['category'] == 'known']
    novel  = [r for r in above if r['category'] == 'novel']
    # rfam_alert: miRDeep2 writes 'yes' when an Rfam match is found
    rfam   = [r for r in above
              if r.get('rfam_alert', '').strip().lower() == 'yes']
    # significant_randfold: miRDeep2 writes 'yes' for significant p-value
    randf  = [r for r in above
              if r.get('significant_randfold', '').strip().lower() == 'yes']
    reads  = sum(r.get('total_reads', 0) for r in above)
    return {
        'total': len(above), 'known': len(known), 'novel': len(novel),
        'rfam_alert': len(rfam), 'significant_randfold': len(randf),
        'total_reads': reads,
    }

stats_path = os.path.join(parsed_dir, 'mirdeep2_parsing_stats.txt')
with open(stats_path, 'w') as fh:
    fh.write("=" * 70 + "\n")
    fh.write("Phase IV – miRDeep2 Analysis & Parsing Statistics\n")
    fh.write(f"Run date   : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    fh.write(f"Source file: {result_csv}\n")
    fh.write("=" * 70 + "\n\n")

    fh.write("Input summary:\n")
    fh.write(f"  Samples combined              : {n_samples}\n")
    fh.write(f"  Total raw reads (all samples) : {total_reads_comb:,}\n")
    fh.write(f"  Collapsed (unique) reads      : {n_collapsed:,}\n")
    fh.write(f"  Genome alignment loci (ARF)   : {n_arf_loci:,}\n\n")

    fh.write("Summary statistics table:\n")
    fh.write(f"  Score tiers extracted         : {len(summary_rows)}\n\n")

    fh.write("All predictions (no score filter):\n")
    fh.write(f"  Known miRNAs  : {len(known_rows)}\n")
    fh.write(f"  Novel miRNAs  : {len(novel_rows)}\n")
    fh.write(f"  Total         : {len(all_rows)}\n\n")

    for label, thresh in [('Low', score_low), ('Medium', score_med), ('High', score_high)]:
        s = tier_stats(all_rows, thresh)
        fh.write(f"Score >= {thresh:.0f} ({label} confidence):\n")
        fh.write(f"  Total predictions         : {s['total']}\n")
        fh.write(f"    Known miRNAs            : {s['known']}\n")
        fh.write(f"    Novel miRNAs            : {s['novel']}\n")
        fh.write(f"  Total reads at loci       : {s['total_reads']:,}\n")
        fh.write(f"  Rfam alerts               : {s['rfam_alert']}\n")
        fh.write(f"  Significant randfold      : {s['significant_randfold']}\n\n")

    fh.write("Output files:\n")
    fh.write(f"  {parsed_dir}/mirdeep2_summary_stats.tsv\n")
    fh.write(f"  {parsed_dir}/mirdeep2_all.tsv\n")
    fh.write(f"  {parsed_dir}/mirdeep2_known.tsv\n")
    fh.write(f"  {parsed_dir}/mirdeep2_novel.tsv\n")
    fh.write(f"  {parsed_dir}/mirdeep2_high_conf.tsv\n")
    fh.write(f"  {parsed_dir}/mirdeep2_parsing_stats.txt\n\n")

    fh.write("Next step:\n")
    fh.write("  Run 5.1_mirna_reference.sh to build the unified miRNA quantification\n")
    fh.write("  reference by merging these known mature miRNAs with the selected\n")
    fh.write("  novel candidates.\n")

print("\n" + open(stats_path).read())
PYTHON_PARSE

# ============================================================================
# Step 5: Extract FASTA sequences using get_mirdeep2_precursors.pl
# ============================================================================

echo ""
echo "Step 5: Extracting mature and precursor sequences per score tier..."

# Find newest miRDeep2 result CSV (timestamped filename)
RESULT_CSV=$(find "${MIRDEEP2_DIR}" -maxdepth 1 -name 'result_*.csv' -printf '%T@ %p\n' \
    2>/dev/null | sort -rn | head -1 | cut -d' ' -f2-)

if [ -z "${RESULT_CSV}" ]; then
    echo "ERROR: No miRDeep2 result CSV found in ${MIRDEEP2_DIR}"
    exit 1
fi

echo "  Result CSV: ${RESULT_CSV}"

# get_mirdeep2_precursors.pl requires the basename and must be run from the
# directory that contains the result file.
RESULT_BASENAME=$(basename "${RESULT_CSV}")
cd "${MIRDEEP2_DIR}" || exit 1

mkdir -p "sequences${SCORE_LOW}" "sequences${SCORE_MED}" "sequences${SCORE_HIGH}"

# -m: extract mature sequences  -d: extract precursor sequences
get_mirdeep2_precursors.pl -r "${RESULT_BASENAME}" -s "${SCORE_LOW}"  -m -d -o "sequences${SCORE_LOW}"
get_mirdeep2_precursors.pl -r "${RESULT_BASENAME}" -s "${SCORE_MED}"  -m -d -o "sequences${SCORE_MED}"
get_mirdeep2_precursors.pl -r "${RESULT_BASENAME}" -s "${SCORE_HIGH}" -m -d -o "sequences${SCORE_HIGH}"

echo ""
echo "FASTA extraction complete."
echo "  sequences${SCORE_LOW}/  – score >= ${SCORE_LOW}  (low confidence)"
echo "  sequences${SCORE_MED}/  – score >= ${SCORE_MED}  (medium confidence)"
echo "  sequences${SCORE_HIGH}/ – score >= ${SCORE_HIGH} (high confidence)"
echo ""
echo "Each directory contains:"
echo "  novel_mature_*.fa    – novel mature miRNA sequences"
echo "  known_mature_*.fa    – known mature miRNA sequences (if applicable)"
echo "  novel_precursor_*.fa – novel precursor sequences"
echo "  known_precursor_*.fa – known precursor sequences (if applicable)"
echo ""
echo "Use these files in 5.1_mirna_reference.sh to build the quantification index."

printf "\n\nCompleted: 4.1-mirdeep2_combined\n\n"
pwd
date

