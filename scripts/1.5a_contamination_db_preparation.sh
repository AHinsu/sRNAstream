#!/usr/bin/env bash
#SBATCH --job-name=contam_db_prep
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G


###############################################################################
# 1.5a_contamination_db_preparation.sh - Build contamination reference indexes
# 
# Description:
#   Part A of contamination filtering - prepares Bowtie indexes for all RNA
#   contamination types (rRNA, tRNA, snRNA, snoRNA, piRNA, and other small RNAs).
#   This is a one-time setup step that must be run before 1.5b.
#
#   Sequences are sourced from two places:
#     1. User-provided FASTA files extracted from the RNAcentral database using
#        the scripts in the referenceDownload/ folder (see referenceDownload/README.md)
#     2. Ensembl GTF annotation (gene_biotype-based extraction via bedtools)
#
#   NOTE on Ensembl GTF and RNA biotypes:
#   Ensembl GTF files store RNA biotypes (rRNA, tRNA, snRNA, snoRNA, etc.) as
#   attributes in column 9 (e.g., gene_biotype "rRNA"), NOT as the feature type
#   in column 3 (which only contains: gene, transcript, exon, CDS, UTR, etc.).
#   This script extracts sequences from the GTF using gene_biotype matching and
#   supplements them with user-provided FASTA files from RNAcentral
#   (recommended for more comprehensive coverage).
#
# Usage: 
#   sbatch scripts/1.5a_contamination_db_preparation.sh
#
# Run before:
#   1.5b_contamination_filtering.sh
#
# Outputs:
#   - Bowtie indexes in ${REF_DIR}/contamination/
###############################################################################

printf "\n\nStarted Job: 1.5a-contamination_db_preparation\n\n"
pwd
date

printf "\n\n"

# ============================================================================
# USER CONFIGURATION - Edit these variables for your project
# ============================================================================

# Reference genome FASTA file (from Ensembl)
GENOME_FASTA="/storage/users/ahinsu/DD-microRNA/reference/Bos_taurus.chromsomal.fa"

# GTF annotation file (from Ensembl)
GTF_FILE="/storage/users/ahinsu/DD-microRNA/reference/Bos_taurus.cleaned.gtf"

# Reference directory (contamination indexes will be written here)
REF_DIR="/storage/users/ahinsu/DD-microRNA/reference"

# Conda environment name
CONDA_ENV="smallRNA"

# ============================================================================
# EXTERNAL DATABASE FASTA FILES (optional but recommended)
# ============================================================================
# NOTE: Ensembl GTF files encode RNA biotypes in the attributes field (column 9)
# as gene_biotype "rRNA", NOT as the feature type in column 3 (which only
# contains: gene, transcript, exon, CDS, UTR, etc.). As a result, Ensembl GTF
# annotations for some RNA types (particularly tRNA and rRNA) may be incomplete.
# For comprehensive contamination filtering, extract species-specific FASTA
# files from the RNAcentral database using the scripts in the referenceDownload/
# folder (see referenceDownload/README.md for full instructions).
#
# All user-provided FASTA files below should be extracted from RNAcentral using
# the SQL queries and helper scripts in referenceDownload/. The RNA type names
# used in the RNAcentral queries correspond to those below (e.g. rRNA, tRNA,
# snRNA, snoRNA, piRNA, misc_RNA, Y_RNA, vault_RNA).
#
# Leave a path as empty string ("") to use the Ensembl GTF alone for that type.
#
# rRNA sequences (RECOMMENDED - Ensembl rRNA annotations are often sparse):
#   Extract from RNAcentral with RNA type "rRNA" for your species taxid.
#   See referenceDownload/README.md - Step 2 and Step 4.
#
# tRNA sequences (RECOMMENDED - Ensembl tRNA annotations are often sparse):
#   Extract from RNAcentral with RNA type "tRNA" for your species taxid.
#   See referenceDownload/README.md - Step 2 and Step 4.
#
# snRNA sequences:
#   Extract from RNAcentral with RNA type "snRNA" for your species taxid.
#   See referenceDownload/README.md - Step 2 and Step 4.
#
# snoRNA sequences:
#   Extract from RNAcentral with RNA type "snoRNA" for your species taxid.
#   See referenceDownload/README.md - Step 2 and Step 4.
#
# piRNA sequences:
#   Extract from RNAcentral with RNA type "piRNA" for your species taxid.
#   See referenceDownload/README.md - Step 2 and Step 4.
#   piRNA reads are saved separately for downstream gene counts and DESeq2
#   analysis (similar to snRNA/snoRNA in Phase VIII).
#
# Other small RNAs (misc_RNA, Y_RNA, vault_RNA, ribozyme):
#   Extract from RNAcentral with the relevant RNA types for your species taxid.
#   See referenceDownload/README.md - Step 2 and Step 4.
#
#
#
# Paths to user-provided species-specific FASTA files (set to "" to skip)
# Replace filenames and path in the example paths below with your path and filenames.
USER_PIRNA_FA="${REF_DIR}/external/bos_taurus_9913_piRNA.fasta"      # e.g., "${REF_DIR}/external/bos_taurus_9913_piRNA.fasta"
USER_RRNA_FA="${REF_DIR}/external/bos_taurus_9913_rRNA.fasta"       # e.g., "${REF_DIR}/external/bos_taurus_9913_rRNA.fasta"
USER_TRNA_FA="${REF_DIR}/external/bos_taurus_9913_tRNA.fasta"       # e.g., "${REF_DIR}/external/bos_taurus_9913_tRNA.fasta"
USER_SNRNA_FA="${REF_DIR}/external/bos_taurus_9913_snRNA.fasta"      # e.g., "${REF_DIR}/external/bos_taurus_9913_snRNA.fasta"
USER_SNORNA_FA="${REF_DIR}/external/bos_taurus_9913_snoRNA.fasta"     # e.g., "${REF_DIR}/external/bos_taurus_9913_snoRNA.fasta"
USER_OTHER_FA="${REF_DIR}/external/bos_taurus_9913_otherRNA.fasta"      # e.g., "${REF_DIR}/external/bos_taurus_9913_otherRNA.fasta"


# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH


# ============================================================================
# Check required files
# ============================================================================

if [ ! -f "$GENOME_FASTA" ]; then
    echo "Error: Genome FASTA file not found: $GENOME_FASTA"
    echo "Please set GENOME_FASTA to the path of your reference genome FASTA file."
    exit 1
fi

if [ ! -f "$GTF_FILE" ]; then
    echo "Warning: GTF file not found: $GTF_FILE"
    echo "GTF-based extraction will be skipped; only user-provided FASTA files will be used."
fi

echo "Reference genome FASTA: $GENOME_FASTA"
echo "GTF annotation file:    $GTF_FILE"

# Create output directory for contamination indexes
mkdir -p "${REF_DIR}/contamination"

# Contamination index paths
CONTAM_RRNA="${REF_DIR}/contamination/rrna_index"
CONTAM_TRNA="${REF_DIR}/contamination/trna_index"
CONTAM_SNRNA="${REF_DIR}/contamination/snrna_index"
CONTAM_SNORNA="${REF_DIR}/contamination/snorna_index"
CONTAM_PIRNA="${REF_DIR}/contamination/pirna_index"
CONTAM_OTHER="${REF_DIR}/contamination/other_index"

# ---------------------------------------------------------------------------
# Helper: extract transcript features from Ensembl GTF by gene_biotype
#
# In Ensembl GTF, RNA biotypes are stored as attributes in column 9
# (e.g., gene_biotype "rRNA"), not as the feature type in column 3.
# This function filters transcript records by gene_biotype and converts
# them to BED format, then uses bedtools to extract FASTA sequences.
#
# Uses the script-level GTF_FILE and GENOME_FASTA variables.
#
# Usage: extract_from_ensembl_gtf <biotype_pattern> <out_gtf> <out_bed> <out_fa>
#   biotype_pattern : ERE pattern matching gene_biotype value(s), e.g. "rRNA" or "snRNA|snoRNA"
# ---------------------------------------------------------------------------
extract_from_ensembl_gtf() {
    local biotype_pattern=$1
    local out_gtf=$2
    local out_bed=$3
    local out_fa=$4

    if [ ! -f "$GTF_FILE" ]; then
        echo "  Ensembl GTF not found at $GTF_FILE - skipping GTF extraction"
        touch "$out_fa"
        return
    fi
    if ! command -v bedtools &> /dev/null; then
        echo "  bedtools not found - skipping GTF extraction"
        touch "$out_fa"
        return
    fi
    if [ ! -f "$GENOME_FASTA" ]; then
        echo "  Genome FASTA not found at $GENOME_FASTA - skipping GTF extraction"
        touch "$out_fa"
        return
    fi

    echo "  Extracting from Ensembl GTF (gene_biotype matching: ${biotype_pattern})..."
    # Filter transcript features by gene_biotype attribute value.
    # Parentheses around pat are required for ERE alternation (e.g. "snRNA|snoRNA").
    awk -v pat="$biotype_pattern" \
        '$3 == "transcript" && $0 ~ ("gene_biotype \"(" pat ")\"")' \
        "$GTF_FILE" > "$out_gtf" || touch "$out_gtf"

    if [ -s "$out_gtf" ]; then
        # Convert GTF to BED; extract transcript_id from attributes for the name field.
        # Uses POSIX awk: match() sets RSTART/RLENGTH; "transcript_id \"" is 15 chars.
        awk 'BEGIN{OFS="\t"} {
            tid="unknown"
            if (match($0, /transcript_id "[^"]+"/) ) tid=substr($0, RSTART+15, RLENGTH-16)
            print $1,$4-1,$5,tid,$6,$7
        }' "$out_gtf" > "$out_bed"
        bedtools getfasta -fi "$GENOME_FASTA" -bed "$out_bed" -name -s -fo "$out_fa" || touch "$out_fa"
        local n
        n=$(grep -c '^>' "$out_fa" 2>/dev/null || echo 0)
        echo "  Extracted ${n} sequences from Ensembl GTF"
    else
        echo "  No matching transcripts found in Ensembl GTF for gene_biotype: ${biotype_pattern}"
        touch "$out_fa"
    fi
}

# ---------------------------------------------------------------------------
# Helper: merge user-provided FASTA with Ensembl-extracted FASTA and build
# a bowtie index. At least one source must be non-empty to build the index.
#
# Usage: build_contam_index <rna_type> <user_fa> <ensembl_fa> <contam_index>
# ---------------------------------------------------------------------------
build_contam_index() {
    local rna_type=$1
    local user_fa=$2
    local ensembl_fa=$3
    local contam_index=$4
    local merged_fa="${contam_index}_sequences.fa"

    > "$merged_fa"  # create/empty the merged file

    if [ -n "$user_fa" ] && [ -f "$user_fa" ] && [ -s "$user_fa" ]; then
        echo "  Adding user-provided ${rna_type} sequences: $user_fa"
        cat "$user_fa" >> "$merged_fa"
    else
        if [ -n "$user_fa" ]; then
            echo "  User-provided ${rna_type} FASTA not found or empty: $user_fa"
        fi
    fi

    if [ -f "$ensembl_fa" ] && [ -s "$ensembl_fa" ]; then
        echo "  Adding Ensembl-derived ${rna_type} sequences"
        cat "$ensembl_fa" >> "$merged_fa"
    fi

    if [ -s "$merged_fa" ]; then
        echo "  Building bowtie index for ${rna_type}..."
        bowtie-build "$merged_fa" "$contam_index"
    else
        echo "  Warning: No ${rna_type} sequences available - index will not be built."
        local upper_type
        upper_type=$(echo "$rna_type" | tr '[:lower:]' '[:upper:]')
        echo "           To enable ${rna_type} filtering, provide a FASTA via USER_${upper_type}_FA"
        echo "           or ensure GTF_FILE contains ${rna_type} gene_biotype annotations."
    fi
}

# ---------------------------------------------------------------------------
# Build contamination indexes (skip if index already exists)
# ---------------------------------------------------------------------------

# rRNA
if [ ! -f "${CONTAM_RRNA}.1.ebwt" ]; then
    echo ""
    echo "rRNA index not found - building..."
    extract_from_ensembl_gtf "rRNA|Mt_rRNA" \
        "${REF_DIR}/contamination/rrna_ensembl.gtf" \
        "${REF_DIR}/contamination/rrna_ensembl.bed" \
        "${REF_DIR}/contamination/rrna_ensembl.fa"
    build_contam_index "rRNA" "$USER_RRNA_FA" "${REF_DIR}/contamination/rrna_ensembl.fa" "$CONTAM_RRNA"
fi

# tRNA
if [ ! -f "${CONTAM_TRNA}.1.ebwt" ]; then
    echo ""
    echo "tRNA index not found - building..."
    extract_from_ensembl_gtf "tRNA|Mt_tRNA" \
        "${REF_DIR}/contamination/trna_ensembl.gtf" \
        "${REF_DIR}/contamination/trna_ensembl.bed" \
        "${REF_DIR}/contamination/trna_ensembl.fa"
    build_contam_index "tRNA" "$USER_TRNA_FA" "${REF_DIR}/contamination/trna_ensembl.fa" "$CONTAM_TRNA"
fi

# snRNA
if [ ! -f "${CONTAM_SNRNA}.1.ebwt" ]; then
    echo ""
    echo "snRNA index not found - building..."
    extract_from_ensembl_gtf "snRNA" \
        "${REF_DIR}/contamination/snrna_ensembl.gtf" \
        "${REF_DIR}/contamination/snrna_ensembl.bed" \
        "${REF_DIR}/contamination/snrna_ensembl.fa"
    build_contam_index "snRNA" "$USER_SNRNA_FA" "${REF_DIR}/contamination/snrna_ensembl.fa" "$CONTAM_SNRNA"
fi

# snoRNA
if [ ! -f "${CONTAM_SNORNA}.1.ebwt" ]; then
    echo ""
    echo "snoRNA index not found - building..."
    extract_from_ensembl_gtf "snoRNA" \
        "${REF_DIR}/contamination/snorna_ensembl.gtf" \
        "${REF_DIR}/contamination/snorna_ensembl.bed" \
        "${REF_DIR}/contamination/snorna_ensembl.fa"
    build_contam_index "snoRNA" "$USER_SNORNA_FA" "${REF_DIR}/contamination/snorna_ensembl.fa" "$CONTAM_SNORNA"
fi

# piRNA (reads saved separately for downstream gene counts and DESeq2 matrix)
if [ ! -f "${CONTAM_PIRNA}.1.ebwt" ]; then
    echo ""
    echo "piRNA index not found - building..."
    extract_from_ensembl_gtf "piRNA" \
        "${REF_DIR}/contamination/pirna_ensembl.gtf" \
        "${REF_DIR}/contamination/pirna_ensembl.bed" \
        "${REF_DIR}/contamination/pirna_ensembl.fa"
    build_contam_index "piRNA" "$USER_PIRNA_FA" "${REF_DIR}/contamination/pirna_ensembl.fa" "$CONTAM_PIRNA"
fi

# Other small RNAs (misc_RNA, Y_RNA, vault_RNA, ribozyme)
if [ ! -f "${CONTAM_OTHER}.1.ebwt" ]; then
    echo ""
    echo "Other small RNA index not found - building..."
    extract_from_ensembl_gtf "misc_RNA|Y_RNA|vault_RNA|ribozyme|scaRNA" \
        "${REF_DIR}/contamination/other_ensembl.gtf" \
        "${REF_DIR}/contamination/other_ensembl.bed" \
        "${REF_DIR}/contamination/other_ensembl.fa"
    build_contam_index "other" "$USER_OTHER_FA" "${REF_DIR}/contamination/other_ensembl.fa" "$CONTAM_OTHER"
fi

echo ""
echo "Contamination database preparation complete."
echo "Indexes written to: ${REF_DIR}/contamination/"
echo "Run 1.5b_contamination_filtering.sh to filter samples."

printf "\n\nCompleted\n\n"
pwd
date
