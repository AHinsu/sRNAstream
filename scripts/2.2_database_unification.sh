#!/usr/bin/env bash
#SBATCH --job-name=db_unification
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


###############################################################################
# 2.2_database_unification.sh - Database Unification from User-Supplied FASTAs
#
# Description:
#   Phase II.2: Accept one or more user-supplied FASTA files for mature miRNA
#   and precursor/hairpin miRNA sequences, then:
#
#   1. Count sequences in each input file and display a summary table.
#
#   2. Convert RNA to DNA (U → T) for every input file using prinseq-lite.pl
#      (falling back to Python if not available).
#
#   3. Plot a combined sequence-length density distribution showing all input
#      files overlaid (equivalent to the density plot in script 2.1).
#
#   4. Concatenate all mature inputs into one file, and all precursor inputs
#      into another.
#
#   5. Run seqkit rmdup -s -i (by-sequence, case-insensitive) on each combined
#      FASTA.  Save the duplicate info file (-D) alongside the deduped output.
#
#   6. Plot the length distribution and sequence counts of the final deduped
#      files (before/after dedup bar chart + per-file input counts).
#
#   7. Write a summary report.
#
#   Output is written to DB_UNIFIED_OUTDIR, which is separate from the
#   DB_OUTDIR used by 2.1_database_preparation.sh.
#
# Usage:
#   sbatch scripts/2.2_database_unification.sh
#
# Prerequisites:
#   - seqkit on PATH (provides rmdup)
#   - Python 3 with matplotlib and pandas (provided by the smallRNA conda env)
#   - prinseq-lite.pl on PATH (optional; Python fallback used if absent)
#
# User configuration:
#   Edit MATURE_FASTAS, PRECURSOR_FASTAS, and DB_UNIFIED_OUTDIR below.
#
# Outputs (in ${DB_UNIFIED_OUTDIR}/):
#   mature_dna_<Label>.fasta            – Per-source mature FASTA after U→T
#   precursor_dna_<Label>.fasta         – Per-source precursor FASTA after U→T
#   length_density_<Label>.png          – Per-source length density plots
#   length_density_combined.png         – Overlay density plot (all sources)
#   combined_mature_all.fasta           – All mature sequences concatenated
#   combined_mature_dedup.fasta         – Mature sequences after seqkit rmdup
#   combined_mature_dup_info.txt        – seqkit rmdup duplicate info (mature)
#   combined_precursor_all.fasta        – All precursor sequences concatenated
#   combined_precursor_dedup.fasta      – Precursor sequences after seqkit rmdup
#   combined_precursor_dup_info.txt     – seqkit rmdup duplicate info (precursor)
#   per_file_sequence_counts.png        – Bar chart: sequences per input file
#   length_density_final.png            – Length distribution of deduped files
#   dedup_stats.png                     – Before/after dedup sequence counts
#   db_unification_summary.txt          – Overall summary
###############################################################################

printf "\n\nStarted Job: 2.2-database_unification\n\n"
pwd
date
printf "\n\n"

SCRIPT_DIR="/storage/users/ahinsu/DD-microRNA/scripts"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory for unified database files (different from 2.1's DB_OUTDIR)
DB_UNIFIED_OUTDIR="/storage/users/ahinsu/DD-microRNA/reference/db_unified"

# Conda environment name
CONDA_ENV="smallRNA"

# -------------------------------------------------------------------------
# Mature miRNA FASTA files (1 or more entries).
# Format: each entry is "Label:/absolute/path/to/file.fa"
#   • Use filesystem-safe characters in Label (letters, digits, _, -).
#   • U→T conversion is applied automatically.
#   • Tip: run db_mirbase_filter.py first if using unfiltered miRBase files.
MATURE_FASTAS=(
    "miRBase_mature:/storage/users/ahinsu/DD-microRNA/reference/mirbase/bta_mature.fa"
    "MirGeneDB_mature:/storage/users/ahinsu/DD-microRNA/reference/mirgenedb/mirgenedb_mature.fa"
    )

# Precursor/hairpin miRNA FASTA files (1 or more entries).
# Same format as MATURE_FASTAS.
PRECURSOR_FASTAS=(
    "miRBase_hairpin:/storage/users/ahinsu/DD-microRNA/reference/mirbase/bta_hairpin.fa"
    "MirGeneDB_hairpin:/storage/users/ahinsu/DD-microRNA/reference/mirgenedb/mirgenedb_hairpin.fa"
)

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

mkdir -p "${DB_UNIFIED_OUTDIR}"

# ============================================================================
# Helper: sanitise a label string for use in filenames
# ============================================================================
sanitize_label() {
    echo "${1//[^a-zA-Z0-9_-]/_}"
}

# ============================================================================
# Helper: convert RNA to DNA (U → T)
# Accepts a full output path (including .fasta extension).
# Uses prinseq-lite.pl if available; falls back to Python.
# ============================================================================
rna_to_dna() {
    local input_fa="$1"
    local output_fa="$2"   # caller must supply full path with .fasta extension

    if [ ! -s "$input_fa" ]; then
        echo "  Skipping U→T conversion: file not found or empty: $input_fa"
        return
    fi

    if command -v prinseq-lite.pl &> /dev/null; then
        echo "  Using prinseq-lite.pl for U→T conversion..."
        local out_prefix="${output_fa%.fasta}"
        prinseq-lite.pl -fasta "$input_fa" -dna_rna dna -out_format 1 \
            -out_good "$out_prefix" -out_bad null
        if [ $? -ne 0 ] || [ ! -s "$output_fa" ]; then
            echo "  ERROR: prinseq-lite.pl conversion failed for $input_fa" >&2
            return 1
        fi
    else
        echo "  prinseq-lite.pl not found – using Python for U→T conversion..."
        python3 "${SCRIPT_DIR}/db_rna_to_dna.py" "$input_fa" "$output_fa"
        if [ $? -ne 0 ] || [ ! -s "$output_fa" ]; then
            echo "  ERROR: Python U→T conversion failed for $input_fa" >&2
            return 1
        fi
    fi

    local n_seqs
    n_seqs=$(grep -c '^>' "$output_fa" 2>/dev/null || echo 0)
    echo "  Written ${n_seqs} sequences to ${output_fa}"
}

# ============================================================================
# Step 1: Count sequences per input file and display summary table
# ============================================================================

echo "============================================================"
echo "Step 1: Counting sequences per input file"
echo "============================================================"
echo ""

printf "  %-30s  %-12s  %s\n" "Label" "Sequences" "File"
printf "  %-30s  %-12s  %s\n" "-----" "---------" "----"

for entry in "${MATURE_FASTAS[@]}"; do
    label="${entry%%:*}"
    raw_fa="${entry#*:}"
    if [ -s "$raw_fa" ]; then
        n=$(grep -c '^>' "$raw_fa" 2>/dev/null || echo 0)
    else
        n=0
    fi
    printf "  %-30s  %-12s  %s\n" "[mature] ${label}" "$n" "$raw_fa"
done

for entry in "${PRECURSOR_FASTAS[@]}"; do
    label="${entry%%:*}"
    raw_fa="${entry#*:}"
    if [ -s "$raw_fa" ]; then
        n=$(grep -c '^>' "$raw_fa" 2>/dev/null || echo 0)
    else
        n=0
    fi
    printf "  %-30s  %-12s  %s\n" "[precursor] ${label}" "$n" "$raw_fa"
done
echo ""

# ============================================================================
# Step 2: Convert each input FASTA from RNA to DNA (U → T)
# ============================================================================

echo "============================================================"
echo "Step 2: Converting input FASTA files (U → T)"
echo "============================================================"

MATURE_DNA_FILES=()
MATURE_SOURCE_PAIRS=()
PRECURSOR_DNA_FILES=()
PRECURSOR_SOURCE_PAIRS=()

for entry in "${MATURE_FASTAS[@]}"; do
    label="${entry%%:*}"
    raw_fa="${entry#*:}"
    safe_label="$(sanitize_label "${label}")"
    dna_fa="${DB_UNIFIED_OUTDIR}/mature_dna_${safe_label}.fasta"

    echo "Processing mature – ${label}..."
    rna_to_dna "$raw_fa" "$dna_fa"

    if [ -s "$dna_fa" ]; then
        MATURE_DNA_FILES+=("$dna_fa")
        MATURE_SOURCE_PAIRS+=("${label}:${dna_fa}")
    fi
done

for entry in "${PRECURSOR_FASTAS[@]}"; do
    label="${entry%%:*}"
    raw_fa="${entry#*:}"
    safe_label="$(sanitize_label "${label}")"
    dna_fa="${DB_UNIFIED_OUTDIR}/precursor_dna_${safe_label}.fasta"

    echo "Processing precursor – ${label}..."
    rna_to_dna "$raw_fa" "$dna_fa"

    if [ -s "$dna_fa" ]; then
        PRECURSOR_DNA_FILES+=("$dna_fa")
        PRECURSOR_SOURCE_PAIRS+=("${label}:${dna_fa}")
    fi
done

# Build pipe-separated strings consumed by Python helper scripts
MATURE_SOURCES="$(IFS='|'; echo "${MATURE_SOURCE_PAIRS[*]}")"
PRECURSOR_SOURCES="$(IFS='|'; echo "${PRECURSOR_SOURCE_PAIRS[*]}")"
export MATURE_SOURCES PRECURSOR_SOURCES

echo ""
echo "Mature sources   : ${MATURE_SOURCES:-<none>}"
echo "Precursor sources: ${PRECURSOR_SOURCES:-<none>}"

# ============================================================================
# Step 3: Plot combined sequence-length density distribution (all input files)
# ============================================================================

echo ""
echo "============================================================"
echo "Step 3: Plotting combined sequence-length density distribution"
echo "============================================================"

export DB_OUTDIR="${DB_UNIFIED_OUTDIR}"

python3 "${SCRIPT_DIR}/db_density_plot.py"

# ============================================================================
# Step 4: Concatenate sequences per category
# ============================================================================

echo ""
echo "============================================================"
echo "Step 4: Combining sequences per category"
echo "============================================================"

COMBINED_MATURE="${DB_UNIFIED_OUTDIR}/combined_mature_all.fasta"
COMBINED_PRECURSOR="${DB_UNIFIED_OUTDIR}/combined_precursor_all.fasta"

> "${COMBINED_MATURE}"
for fa in "${MATURE_DNA_FILES[@]}"; do
    [ -s "$fa" ] && cat "$fa" >> "${COMBINED_MATURE}"
done
COMBINED_MATURE_CNT=$(grep -c '^>' "${COMBINED_MATURE}" 2>/dev/null || echo 0)
echo "Combined mature sequences   : ${COMBINED_MATURE_CNT}"

> "${COMBINED_PRECURSOR}"
for fa in "${PRECURSOR_DNA_FILES[@]}"; do
    [ -s "$fa" ] && cat "$fa" >> "${COMBINED_PRECURSOR}"
done
COMBINED_PRECURSOR_CNT=$(grep -c '^>' "${COMBINED_PRECURSOR}" 2>/dev/null || echo 0)
echo "Combined precursor sequences: ${COMBINED_PRECURSOR_CNT}"

# ============================================================================
# Step 5: Deduplicate with seqkit rmdup (by sequence, case-insensitive)
#         Save duplicate info file with -D flag
# ============================================================================

echo ""
echo "============================================================"
echo "Step 5: Running seqkit rmdup"
echo "============================================================"

COMBINED_MATURE_DEDUP="${DB_UNIFIED_OUTDIR}/combined_mature_dedup.fasta"
MATURE_DUP_INFO="${DB_UNIFIED_OUTDIR}/combined_mature_dup_info.txt"
COMBINED_PRECURSOR_DEDUP="${DB_UNIFIED_OUTDIR}/combined_precursor_dedup.fasta"
PRECURSOR_DUP_INFO="${DB_UNIFIED_OUTDIR}/combined_precursor_dup_info.txt"

echo "Running seqkit rmdup on combined mature..."
if [ -s "${COMBINED_MATURE}" ]; then
    seqkit rmdup -s -i \
        -D "${MATURE_DUP_INFO}" \
        -o "${COMBINED_MATURE_DEDUP}" \
        "${COMBINED_MATURE}"
    COMBINED_MATURE_DEDUP_CNT=$(grep -c '^>' "${COMBINED_MATURE_DEDUP}" 2>/dev/null || echo 0)
    echo "  Mature : ${COMBINED_MATURE_CNT} → ${COMBINED_MATURE_DEDUP_CNT} after dedup"
    echo "  Duplicate info: ${MATURE_DUP_INFO}"
else
    echo "  WARNING: Combined mature FASTA is empty – skipping seqkit rmdup"
fi

echo "Running seqkit rmdup on combined precursor..."
if [ -s "${COMBINED_PRECURSOR}" ]; then
    seqkit rmdup -s -i \
        -D "${PRECURSOR_DUP_INFO}" \
        -o "${COMBINED_PRECURSOR_DEDUP}" \
        "${COMBINED_PRECURSOR}"
    COMBINED_PRECURSOR_DEDUP_CNT=$(grep -c '^>' "${COMBINED_PRECURSOR_DEDUP}" 2>/dev/null || echo 0)
    echo "  Precursor: ${COMBINED_PRECURSOR_CNT} → ${COMBINED_PRECURSOR_DEDUP_CNT} after dedup"
    echo "  Duplicate info: ${PRECURSOR_DUP_INFO}"
else
    echo "  WARNING: Combined precursor FASTA is empty – skipping seqkit rmdup"
fi

# ============================================================================
# Step 6: Plot length distributions and sequence counts of the final files
# ============================================================================

echo ""
echo "============================================================"
echo "Step 6: Plotting final length distributions and dedup stats"
echo "============================================================"

export DB_OUTDIR="${DB_UNIFIED_OUTDIR}"
export MATURE_SOURCES PRECURSOR_SOURCES
export COMBINED_MATURE COMBINED_MATURE_DEDUP
export COMBINED_PRECURSOR COMBINED_PRECURSOR_DEDUP

python3 "${SCRIPT_DIR}/db_unification_plot.py"

# ============================================================================
# Step 7: Summary report
# ============================================================================

echo ""
echo "============================================================"
echo "Step 7: Writing summary report"
echo "============================================================"

SUMMARY="${DB_UNIFIED_OUTDIR}/db_unification_summary.txt"

{
    echo "========================================================================"
    echo "Phase II.2 – Database Unification Summary"
    echo "Run date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================================================"
    echo ""
    echo "Output directory: ${DB_UNIFIED_OUTDIR}"
    echo ""
    echo "Mature input files:"
    for entry in "${MATURE_FASTAS[@]}"; do
        label="${entry%%:*}"
        raw_fa="${entry#*:}"
        n=$(grep -c '^>' "$raw_fa" 2>/dev/null || echo 0)
        printf "  %-30s  %6d sequences   %s\n" "$label" "$n" "$raw_fa"
    done
    echo ""
    echo "Precursor input files:"
    for entry in "${PRECURSOR_FASTAS[@]}"; do
        label="${entry%%:*}"
        raw_fa="${entry#*:}"
        n=$(grep -c '^>' "$raw_fa" 2>/dev/null || echo 0)
        printf "  %-30s  %6d sequences   %s\n" "$label" "$n" "$raw_fa"
    done
    echo ""
    echo "Combined (pre-dedup):"
    printf "  %-30s  %6d sequences\n" "Mature combined" \
        "$(grep -c '^>' "${COMBINED_MATURE}" 2>/dev/null || echo 0)"
    printf "  %-30s  %6d sequences\n" "Precursor combined" \
        "$(grep -c '^>' "${COMBINED_PRECURSOR}" 2>/dev/null || echo 0)"
    echo ""
    echo "After seqkit rmdup (final files):"
    printf "  %-30s  %6d sequences   %s\n" "Mature deduped" \
        "$(grep -c '^>' "${COMBINED_MATURE_DEDUP}" 2>/dev/null || echo 0)" \
        "${COMBINED_MATURE_DEDUP}"
    printf "  %-30s  %6d sequences   %s\n" "Precursor deduped" \
        "$(grep -c '^>' "${COMBINED_PRECURSOR_DEDUP}" 2>/dev/null || echo 0)" \
        "${COMBINED_PRECURSOR_DEDUP}"
    echo ""
    echo "Duplicate info files:"
    echo "  ${MATURE_DUP_INFO}"
    echo "  ${PRECURSOR_DUP_INFO}"
    echo ""
    echo "Plots:"
    for png in "${DB_UNIFIED_OUTDIR}"/length_density_*.png \
               "${DB_UNIFIED_OUTDIR}"/per_file_sequence_counts.png \
               "${DB_UNIFIED_OUTDIR}"/length_density_final.png \
               "${DB_UNIFIED_OUTDIR}"/dedup_stats.png; do
        [ -f "$png" ] && echo "  $png"
    done
} | tee "${SUMMARY}"

printf "\n\nCompleted: 2.2-database_unification\n\n"
pwd
date

