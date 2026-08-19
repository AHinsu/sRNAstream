#!/usr/bin/env bash
#SBATCH --job-name=db_preparation
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G


###############################################################################
# 2.1_database_preparation.sh - Database Files Preparation
#
# Description:
#   Phase II.1: Prepare and QC all miRNA reference database files.
#
#   1. (Optional) Extract miRNA gene sequences (gene_biotype "miRNA" and
#      "miRNA_primary_transcript") from the user-supplied Ensembl genome FASTA
#      and GTF using bedtools getfasta.  Leave GENOME_FASTA or GTF_FILE empty
#      to skip this step.
#
#   2. Convert RNA sequences (U → T) for all external FASTA files supplied by
#      the user via the MATURE_FASTAS and PRECURSOR_FASTAS arrays (see USER
#      CONFIGURATION below).  Conversion uses the prinseq-lite.pl script
#      bundled with the conda environment, falling back to Python if it is not
#      found on PATH.
#
#      Note: miRBase files must be pre-filtered by species before being listed
#      here.  Use scripts/db_mirbase_filter.py for that step.
#
#   3. Plot a density plot showing the number of sequences per length for every
#      processed FASTA file (Ensembl + all user-supplied files after U→T
#      conversion).
#
#   4. Run duplication statistics separately for:
#        - mature sequences  (exact duplicate sequences)
#        - precursor sequences (exact duplicate sequences)
#      and report the number of unique vs duplicate sequences per source.
#
#   4b. Combine all sequences per category and run seqkit rmdup on the merged
#       FASTA, then run prinseq-lite.pl (-stats_info, -stats_len, -stats_dupl)
#       on both the combined and deduped FASTAs, and plot dedup statistics.
#
# Usage:
#   sbatch scripts/2.1_database_preparation.sh
#
# Prerequisites:
#   - bedtools on PATH (or loaded via conda/module)  [only needed for Ensembl]
#   - Python 3 with matplotlib and pandas (provided by the smallRNA conda env)
#   - seqkit on PATH
#   - prinseq-lite.pl on PATH (optional; Python fallback used if absent)
#
# Outputs (in ${DB_OUTDIR}/):
#   ensembl_mirna.fasta                        – Ensembl miRNA sequences (DNA)
#   mature_dna_<Label>.fasta                   – Mature FASTA after U→T per source
#   precursor_dna_<Label>.fasta                – Precursor FASTA after U→T per source
#   length_density_<Label>.png                 – Density plots per FASTA file
#   length_density_combined.png                – Overlay density plot (all sources)
#   duplication_stats_mature.txt               – Duplication report for mature sequences
#   duplication_stats_precursor.txt            – Duplication report for precursor sequences
#   combined_mature_all.fasta                  – All mature sequences merged
#   combined_mature_dedup.fasta                – Merged mature sequences after seqkit rmdup
#   prinseq_dedup_mature.fasta                 – Merged mature sequences after prinseq dedup
#   combined_hairpin_all.fasta                 – All hairpin/precursor sequences merged
#   combined_hairpin_dedup.fasta               – Merged hairpin sequences after seqkit rmdup
#   prinseq_dedup_hairpin.fasta                – Merged hairpin sequences after prinseq dedup
#   prinseq_combined_mature_all.txt            – prinseq-lite.pl stats (combined mature)
#   prinseq_combined_mature_dedup.txt          – prinseq-lite.pl stats (seqkit deduped mature)
#   prinseq_combined_mature_prinseqDedup.txt   – prinseq-lite.pl stats (prinseq deduped mature)
#   prinseq_combined_hairpin_all.txt           – prinseq-lite.pl stats (combined hairpin)
#   prinseq_combined_hairpin_dedup.txt         – prinseq-lite.pl stats (seqkit deduped hairpin)
#   prinseq_combined_hairpin_prinseqDedup.txt  – prinseq-lite.pl stats (prinseq deduped hairpin)
#   combined_dedup_stats.png                   – Bar chart: before/after seqkit rmdup/prinseq dedup
#   db_preparation_summary.txt                 – Overall summary statistics
###############################################################################

printf "\n\nStarted Job: 2.1-database_preparation\n\n"
pwd
date
printf "\n\n"

# Resolve the directory that contains this script so Python helper scripts
# can be located relative to it, regardless of the working directory.
SCRIPT_DIR="/storage/users/ahinsu/DD-microRNA/scripts"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Reference genome FASTA (Ensembl) – leave empty ("") to skip Ensembl extraction
GENOME_FASTA="/storage/users/ahinsu/DD-microRNA/reference/Bos_taurus.chromsomal.fa"

# Ensembl GTF annotation – leave empty ("") to skip Ensembl extraction
GTF_FILE="/storage/users/ahinsu/DD-microRNA/reference/Bos_taurus.cleaned.gtf"

# Output directory for all database files
DB_OUTDIR="/storage/users/ahinsu/DD-microRNA/reference/db_prep"

# Conda environment name
CONDA_ENV="smallRNA"

# -------------------------------------------------------------------------
# Additional mature miRNA FASTA files (0 or more entries).
# Format: each entry is "Label:/absolute/path/to/file.fa"
#   • U→T conversion is applied automatically.
#   • Use filesystem-safe characters in Label (letters, digits, _, -).
#   • Tip: pre-filter miRBase files by species with db_mirbase_filter.py
#     before listing them here.
# Leave the array empty () to provide no additional mature files.
MATURE_FASTAS=(
    "miRBase_mature:/storage/users/ahinsu/DD-microRNA/reference/mirbase/bta_mature.fa"
    "MirGeneDB_mature:/storage/users/ahinsu/DD-microRNA/reference/mirgenedb/mirgenedb_mature.fa"
    "RNAcentral:/storage/users/ahinsu/DD-microRNA/reference/rnacentral/rnacentral_mirna.fa"
)

# Additional precursor/hairpin miRNA FASTA files (0 or more entries).
# Same format as MATURE_FASTAS.  Ensembl sequences (if extracted above) are
# automatically included and do not need to be listed here.
# Leave the array empty () to provide no additional precursor files.
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

# Create output directory
mkdir -p "${DB_OUTDIR}"

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
        # prinseq-lite.pl appends .fasta to -out_good prefix
        local out_prefix="${output_fa%.fasta}"
        prinseq-lite.pl -fasta "$input_fa" -dna_rna dna -out_format 1 \
            -out_good "$out_prefix" -out_bad null
        if [ $? -ne 0 ] || [ ! -s "$output_fa" ]; then
            echo "  ERROR: prinseq-lite.pl conversion failed for $input_fa" >&2
            return 1
        fi
        # result is now at ${out_prefix}.fasta = $output_fa
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
# Step 1: Extract Ensembl miRNA sequences
# ============================================================================

echo "============================================================"
echo "Step 1: Extracting Ensembl miRNA sequences from GTF + FASTA"
echo "============================================================"

ENSEMBL_GTF="${DB_OUTDIR}/ensembl_mirna.gtf"
ENSEMBL_FA="${DB_OUTDIR}/ensembl_mirna.fasta"
ENSEMBL_BED="${DB_OUTDIR}/ensembl_mirna.bed"

# Convert GTF → BED6, then extract sequences with bedtools getfasta
gtf_to_bed_and_extract() {
    local gtf="$1"
    local bed="$2"
    local fasta_out="$3"
    local label="$4"

    if [ ! -s "$gtf" ]; then
        echo "  Skipping ${label}: GTF file is empty."
        touch "$fasta_out"
        return
    fi

    # Build BED6 (0-based start)
    awk 'BEGIN{OFS="\t"} {
        chr=$1; start=$4-1; end=$5; score=0; strand=$7;
        name="";
        match($9, /gene_name "([^"]+)"/, a); if (a[1]!="") name=a[1];
        else { match($9, /gene_id "([^"]+)"/, b); name=b[1]; }
        print chr, start, end, name, score, strand
    }' "$gtf" > "$bed"

    if command -v bedtools &> /dev/null; then
        bedtools getfasta \
            -fi "$GENOME_FASTA" \
            -bed "$bed" \
            -name \
            -s \
            -fo "$fasta_out" 2>/dev/null
    else
        echo "  WARNING: bedtools not found – cannot extract ${label} sequences"
        touch "$fasta_out"
        return
    fi

    local n_seqs
    n_seqs=$(grep -c '^>' "$fasta_out" 2>/dev/null || echo 0)
    echo "  Extracted ${n_seqs} ${label} sequences → ${fasta_out}"
}


# Skip this block if either file is missing
if [ ! -f "$GTF_FILE" ] || [ ! -f "$GENOME_FASTA" ]; then
    echo "Skipping block because required file(s) are missing:"
    [ ! -f "$GTF_FILE" ] && echo "  - Missing GTF: $GTF_FILE"
    [ ! -f "$GENOME_FASTA" ] && echo "  - Missing FASTA: $GENOME_FASTA"
else
    # Both files exist — run this block
# Extract gene entries where gene_biotype is "miRNA"
grep -P '\tgene\t' "$GTF_FILE" \
  | grep -P 'gene_biotype "miRNA"' \
  > "$ENSEMBL_GTF" || touch "$ENSEMBL_GTF"

ENSEMBL_CNT=$(grep -c '^' "$ENSEMBL_GTF" 2>/dev/null || echo 0)
echo "Found ${ENSEMBL_CNT} miRNA gene entries in Ensembl GTF"

gtf_to_bed_and_extract \
    "$ENSEMBL_GTF" \
    "$ENSEMBL_BED" \
    "$ENSEMBL_FA" \
    "Ensembl miRNA"
    
fi



# ============================================================================
# Step 2: Convert external database FASTA files from RNA to DNA (U → T)
# ============================================================================

echo ""
echo "============================================================"
echo "Step 2: Converting user-supplied FASTA files (U → T)"
echo "============================================================"

# Arrays to collect the processed file paths and label:path pairs for Python.
MATURE_DNA_FILES=()
MATURE_SOURCE_PAIRS=()
PRECURSOR_DNA_FILES=()
PRECURSOR_SOURCE_PAIRS=()

# Process each entry in MATURE_FASTAS
for entry in "${MATURE_FASTAS[@]}"; do
    label="${entry%%:*}"
    raw_fa="${entry#*:}"
    # Sanitise label for use in filenames
    safe_label="$(sanitize_label "${label}")"
    dna_fa="${DB_OUTDIR}/mature_dna_${safe_label}.fasta"

    echo "Processing mature – ${label}..."
    rna_to_dna "$raw_fa" "$dna_fa"

    if [ -s "$dna_fa" ]; then
        MATURE_DNA_FILES+=("$dna_fa")
        MATURE_SOURCE_PAIRS+=("${label}:${dna_fa}")
    fi
done

# Process each entry in PRECURSOR_FASTAS
for entry in "${PRECURSOR_FASTAS[@]}"; do
    label="${entry%%:*}"
    raw_fa="${entry#*:}"
    safe_label="$(sanitize_label "${label}")"
    dna_fa="${DB_OUTDIR}/precursor_dna_${safe_label}.fasta"

    echo "Processing precursor – ${label}..."
    rna_to_dna "$raw_fa" "$dna_fa"

    if [ -s "$dna_fa" ]; then
        PRECURSOR_DNA_FILES+=("$dna_fa")
        PRECURSOR_SOURCE_PAIRS+=("${label}:${dna_fa}")
    fi
done

# Prepend Ensembl to precursor list if it was extracted successfully
if [ -s "${ENSEMBL_FA}" ]; then
    PRECURSOR_DNA_FILES=("${ENSEMBL_FA}" "${PRECURSOR_DNA_FILES[@]}")
    PRECURSOR_SOURCE_PAIRS=("Ensembl:${ENSEMBL_FA}" "${PRECURSOR_SOURCE_PAIRS[@]}")
fi

# Build pipe-separated strings consumed by Python helper scripts
MATURE_SOURCES="$(IFS='|'; echo "${MATURE_SOURCE_PAIRS[*]}")"
PRECURSOR_SOURCES="$(IFS='|'; echo "${PRECURSOR_SOURCE_PAIRS[*]}")"
export MATURE_SOURCES PRECURSOR_SOURCES

echo ""
echo "Mature sources   : ${MATURE_SOURCES:-<none>}"
echo "Precursor sources: ${PRECURSOR_SOURCES:-<none>}"

# ============================================================================
# Step 3: Plot density plots (sequences per length)
# ============================================================================

echo ""
echo "============================================================"
echo "Step 3: Plotting sequence-length density plots"
echo "============================================================"

export DB_OUTDIR ENSEMBL_FA MATURE_SOURCES PRECURSOR_SOURCES

python3 "${SCRIPT_DIR}/db_density_plot.py"

# ============================================================================
# Step 4: Duplication statistics for mature and precursor sequences
# ============================================================================

echo ""
echo "============================================================"
echo "Step 4: Duplication statistics"
echo "============================================================"

export DB_OUTDIR ENSEMBL_FA MATURE_SOURCES PRECURSOR_SOURCES

python3 "${SCRIPT_DIR}/db_duplication_stats.py"

# ============================================================================
# Step 4b: Combined deduplication – merge all sources and run seqkit rmdup
#          then run prinseq-lite.pl stats and plot results
# ============================================================================

echo ""
echo "============================================================"
echo "Step 4b: Combined deduplication (seqkit rmdup + prinseq-lite.pl)"
echo "============================================================"

COMBINED_MATURE="${DB_OUTDIR}/combined_mature_all.fasta"
COMBINED_MATURE_DEDUP="${DB_OUTDIR}/combined_mature_dedup.fasta"
COMBINED_HAIRPIN="${DB_OUTDIR}/combined_hairpin_all.fasta"
COMBINED_HAIRPIN_DEDUP="${DB_OUTDIR}/combined_hairpin_dedup.fasta"

# Combine all processed mature FASTA files
> "${COMBINED_MATURE}"
for fa in "${MATURE_DNA_FILES[@]}"; do
    [ -s "$fa" ] && cat "$fa" >> "${COMBINED_MATURE}"
done
COMBINED_MATURE_CNT=$(grep -c '^>' "${COMBINED_MATURE}" 2>/dev/null || echo 0)
echo "Combined mature sequences (all sources): ${COMBINED_MATURE_CNT}"

# Combine all processed precursor FASTA files (Ensembl + user-supplied)
> "${COMBINED_HAIRPIN}"
for fa in "${PRECURSOR_DNA_FILES[@]}"; do
    [ -s "$fa" ] && cat "$fa" >> "${COMBINED_HAIRPIN}"
done
COMBINED_HAIRPIN_CNT=$(grep -c '^>' "${COMBINED_HAIRPIN}" 2>/dev/null || echo 0)
echo "Combined hairpin sequences (all sources): ${COMBINED_HAIRPIN_CNT}"

# Run seqtk rmdup on combined FASTAs
echo "Running seqtk rmdup on combined mature..."
seqkit rmdup -s -i -o "${COMBINED_MATURE_DEDUP}" "${COMBINED_MATURE}"
COMBINED_MATURE_DEDUP_CNT=$(grep -c '^>' "${COMBINED_MATURE_DEDUP}" 2>/dev/null || echo 0)
echo "  Mature after dedup: ${COMBINED_MATURE_DEDUP_CNT}"

echo "Running seqtk rmdup on combined hairpin..."
seqkit rmdup -s -i -o "${COMBINED_HAIRPIN_DEDUP}" "${COMBINED_HAIRPIN}" 
COMBINED_HAIRPIN_DEDUP_CNT=$(grep -c '^>' "${COMBINED_HAIRPIN_DEDUP}" 2>/dev/null || echo 0)
echo "  Hairpin after dedup: ${COMBINED_HAIRPIN_DEDUP_CNT}"

# Run prinseq-lite.pl stats on combined and deduped FASTAs
echo "Running prinseq-lite.pl stats for combined mature (all sources)..."
if [ -s "${COMBINED_MATURE}" ]; then
    prinseq-lite.pl -fasta "${COMBINED_MATURE}" \
        -stats_info -stats_len -stats_dupl \
        -out_good null -out_bad null \
        > "${DB_OUTDIR}/prinseq_combined_mature_all.txt" 2>&1
    echo "  Written: ${DB_OUTDIR}/prinseq_combined_mature_all.txt"
fi

echo "Running prinseq-lite.pl dedup for combined mature..."
if [ -s "${COMBINED_MATURE}" ]; then
    # -derep 14 removes exact duplicates (1), 5'-duplicates (2), 3'-duplicates (4),
    # and reverse-complement duplicates (8): 1+2+4+8 = 15; 14 omits exact-only (1+2+4+8-1=14)
    # NOTE: -out_good takes a filename PREFIX; prinseq-lite.pl appends .fasta automatically,
    # so the output will be written to ${DB_OUTDIR}/prinseq_dedup_mature.fasta
    prinseq-lite.pl -fasta "${COMBINED_MATURE}" \
        -derep 14 -out_format 1 -out_good "${DB_OUTDIR}/prinseq_dedup_mature" -out_bad null
    prinseq-lite.pl -fasta "${DB_OUTDIR}/prinseq_dedup_mature.fasta" \
        -stats_info -stats_len -stats_dupl \
        -out_good null -out_bad null \
        > "${DB_OUTDIR}/prinseq_combined_mature_prinseqDedup.txt" 2>&1
    echo "  Written: ${DB_OUTDIR}/prinseq_combined_mature_prinseqDedup.txt"
fi


echo "Running prinseq-lite.pl stats for combined mature (after dedup)..."
if [ -s "${COMBINED_MATURE_DEDUP}" ]; then
    prinseq-lite.pl -fasta "${COMBINED_MATURE_DEDUP}" \
        -stats_info -stats_len -stats_dupl \
        -out_good null -out_bad null \
        > "${DB_OUTDIR}/prinseq_combined_mature_dedup.txt" 2>&1
    echo "  Written: ${DB_OUTDIR}/prinseq_combined_mature_dedup.txt"
fi


# for Hairpin
echo "Running prinseq-lite.pl stats for combined hairpin (all sources)..."
if [ -s "${COMBINED_HAIRPIN}" ]; then
    prinseq-lite.pl -fasta "${COMBINED_HAIRPIN}" \
        -stats_info -stats_len -stats_dupl \
        -out_good null -out_bad null \
        > "${DB_OUTDIR}/prinseq_combined_hairpin_all.txt" 2>&1
    echo "  Written: ${DB_OUTDIR}/prinseq_combined_hairpin_all.txt"
fi

echo "Running prinseq-lite.pl dedup for combined hairpin..."
if [ -s "${COMBINED_HAIRPIN}" ]; then
    # -derep 14: removes 5'-duplicates (2), 3'-duplicates (4), and reverse-complement duplicates (8)
    # NOTE: -out_good takes a filename PREFIX; prinseq-lite.pl appends .fasta automatically,
    # so the output will be written to ${DB_OUTDIR}/prinseq_dedup_hairpin.fasta
    prinseq-lite.pl -fasta "${COMBINED_HAIRPIN}" \
        -derep 14 -out_format 1 -out_good "${DB_OUTDIR}/prinseq_dedup_hairpin" -out_bad null
    prinseq-lite.pl -fasta  "${DB_OUTDIR}/prinseq_dedup_hairpin.fasta" \
        -stats_info -stats_len -stats_dupl \
        -out_good null -out_bad null \
        > "${DB_OUTDIR}/prinseq_combined_hairpin_prinseqDedup.txt" 2>&1
    echo "  Written: ${DB_OUTDIR}/prinseq_combined_hairpin_prinseqDedup.txt"
fi

echo "Running prinseq-lite.pl stats for combined hairpin (after dedup)..."
if [ -s "${COMBINED_HAIRPIN_DEDUP}" ]; then
    prinseq-lite.pl -fasta "${COMBINED_HAIRPIN_DEDUP}" \
        -stats_info -stats_len -stats_dupl \
        -out_good null -out_bad null \
        > "${DB_OUTDIR}/prinseq_combined_hairpin_dedup.txt" 2>&1
    echo "  Written: ${DB_OUTDIR}/prinseq_combined_hairpin_dedup.txt"
fi

# Plot combined deduplication statistics
PRINSEQ_DEDUP_MATURE="${DB_OUTDIR}/prinseq_dedup_mature.fasta"
PRINSEQ_DEDUP_HAIRPIN="${DB_OUTDIR}/prinseq_dedup_hairpin.fasta"

export COMBINED_MATURE COMBINED_MATURE_DEDUP PRINSEQ_DEDUP_MATURE
export COMBINED_HAIRPIN COMBINED_HAIRPIN_DEDUP PRINSEQ_DEDUP_HAIRPIN

python3 "${SCRIPT_DIR}/db_combined_dedup_plot.py"

# ============================================================================
# Step 5: Summary report
# ============================================================================

echo ""
echo "============================================================"
echo "Step 5: Writing summary report"
echo "============================================================"

export DB_OUTDIR GENOME_FASTA GTF_FILE ENSEMBL_FA MATURE_SOURCES PRECURSOR_SOURCES

python3 "${SCRIPT_DIR}/db_summary.py"

printf "\n\nCompleted: 2.1-database_preparation\n\n"
pwd
date
