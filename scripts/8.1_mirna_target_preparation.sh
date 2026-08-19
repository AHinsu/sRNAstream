#!/usr/bin/env bash
#SBATCH --job-name=miranda_prep
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


###############################################################################
# 8.1_mirna_target_preparation.sh - Build input files for miRNA target prediction
#
# Description:
#   Phase VIII.1: Prepare the miRNA sequences and 3' UTR reference files
#   required by the miRanda target prediction array job (8.2_miranda_array.sh).
#   Run this script once before submitting 8.2.
#
#   Workflow:
#     1. Extract miRNA sequences from the unified mature miRNA FASTA for each
#        name listed in MIRNA_LIST; write a combined FASTA and an indexed name
#        list (one miRNA per line) for array-job indexing in 8.2.
#     2. Extract 3' UTR coordinates from the Ensembl GTF annotation and write
#        them as a BED file.
#     3. Retrieve the corresponding FASTA sequences with bedtools getfasta and
#        convert to uppercase (required by miRanda).
#
#   At the end the script prints the number of unique miRNAs that will be
#   tested – use this value to set the --array limit in 8.2_miranda_array.sh.
#
# Usage:
#   sbatch scripts/8.1_mirna_target_preparation.sh
#
# Prerequisites:
#   - Unified mature miRNA FASTA (from 5.1_mirna_reference.sh).
#   - Ensembl genome FASTA and GTF annotation.
#
# Outputs (in ${OUTDIR}/temp/):
#   mirna_selected.fa      – Extracted miRNA sequences (FASTA)
#   mirna_list_indexed.txt – One miRNA name per line; line N = array task N in 8.2
#   three_prime_utr.bed    – 3' UTR regions in BED format
#   three_prime_utr.fa     – 3' UTR sequences in FASTA format (uppercase DNA)
#
# Next step:
#   Update #SBATCH --array in 8.2_miranda_array.sh to 1-<miRNA count> and submit:
#   sbatch scripts/8.2_miranda_array.sh
###############################################################################

printf "\n\nStarted Job: 8.1-mirna_target_preparation\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory
OUTDIR="/storage/users/ahinsu/DD-microRNA/downstream/targets/miRNA"

# File containing miRNA names to predict targets for (one name per line)
MIRNA_LIST="/storage/users/ahinsu/DD-microRNA/downstream/mirna_list.txt"

# Unified mature miRNA FASTA (output of 5.1_mirna_reference.sh)
MIRNA_FASTA="/storage/users/ahinsu/DD-microRNA/reference/mirna_quant_ref/unified_mature.fa"

# Ensembl genome FASTA
ENSEMBL_FASTA="/storage/users/ahinsu/DD-microRNA/reference/Bos_taurus.chromsomal.fa"

# Ensembl GTF annotation
ENSEMBL_GTF="/storage/users/ahinsu/DD-microRNA/reference/Bos_taurus.cleaned.gtf"

# Conda environment name
CONDA_ENV="smallRNA"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

TEMP_DIR="${OUTDIR}/temp"
mkdir -p "${TEMP_DIR}"

# ============================================================================
# Step 1: Extract miRNA sequences from the unified mature FASTA
# ============================================================================

echo "Step 1: Extracting miRNA sequences from FASTA..."

MIRNA_FASTA_OUT="${TEMP_DIR}/mirna_selected.fa"
MIRNA_INDEXED_LIST="${TEMP_DIR}/mirna_list_indexed.txt"

# Validate input files
for f in "${MIRNA_LIST}" "${MIRNA_FASTA}" "${ENSEMBL_FASTA}" "${ENSEMBL_GTF}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: Input file not found: ${f}"
        exit 1
    fi
done

> "${MIRNA_FASTA_OUT}"
> "${MIRNA_INDEXED_LIST}"

sort -u "${MIRNA_LIST}" -o "${MIRNA_LIST}"


while IFS= read -r mirna_name; do
    # Skip blank lines and comments
    [[ -z "${mirna_name}" || "${mirna_name}" =~ ^# ]] && continue

    mirna_tmp="${TEMP_DIR}/mirna_selected.${mirna_name//[^A-Za-z0-9_.-]/_}.fa"
    > "${mirna_tmp}"

    seqkit grep -n -p "${mirna_name}" "${MIRNA_FASTA}" > "${mirna_tmp}" 2>/dev/null || true
    if [ -s "${mirna_tmp}" ]; then
        cat "${mirna_tmp}" >> "${MIRNA_FASTA_OUT}"
        echo "${mirna_name}" >> "${MIRNA_INDEXED_LIST}"
    else
        echo "  WARNING: miRNA '${mirna_name}' not found in ${MIRNA_FASTA}"
    fi

    rm -f "${mirna_tmp}"
done < "${MIRNA_LIST}"

MIRNA_COUNT=$(grep -c '^>' "${MIRNA_FASTA_OUT}" 2>/dev/null || echo 0)
echo "  Extracted ${MIRNA_COUNT} miRNA sequences  → ${MIRNA_FASTA_OUT}"
echo "  Indexed name list               → ${MIRNA_INDEXED_LIST}"

if [ "${MIRNA_COUNT}" -eq 0 ]; then
    echo "ERROR: No miRNA sequences were extracted. Check MIRNA_LIST and MIRNA_FASTA."
    exit 1
fi

# ============================================================================
# Step 2: Extract 3' UTR coordinates from the Ensembl GTF
# ============================================================================

echo ""
echo "Step 2: Extracting 3' UTR coordinates from GTF..."

UTR_BED="${TEMP_DIR}/three_prime_utr.bed"

awk '
    BEGIN { OFS="\t" }
    $3 == "three_prime_utr" {
        gene_id = ""; transcript_id = ""
        for (i = 9; i <= NF; i++) {
            if ($i ~ /gene_id/) {
                match($(i+1), /"([^"]+)"/, arr); gene_id = arr[1]
            }
            if ($i ~ /transcript_id/) {
                match($(i+1), /"([^"]+)"/, arr); transcript_id = arr[1]
            }
        }
        # BED format: chrom, chromStart (0-based), chromEnd, name, score, strand
        print $1, $4-1, $5, gene_id ":" transcript_id, "0", $7
    }
' "${ENSEMBL_GTF}" | sort -k1,1 -k2,2n > "${UTR_BED}"

UTR_COUNT=$(wc -l < "${UTR_BED}")
echo "  Extracted ${UTR_COUNT} 3' UTR regions → ${UTR_BED}"

if [ "${UTR_COUNT}" -eq 0 ]; then
    echo "WARNING: No 3' UTR regions found in GTF. Check that the GTF contains three_prime_UTR features."
fi

# ============================================================================
# Step 3: Extract 3' UTR FASTA sequences with bedtools
# ============================================================================

echo ""
echo "Step 3: Extracting 3' UTR FASTA sequences with bedtools getfasta..."

UTR_FASTA="${TEMP_DIR}/three_prime_utr.fa"

bedtools getfasta \
    -fi "${ENSEMBL_FASTA}" \
    -bed "${UTR_BED}" \
    -name \
    -s \
    > "${TEMP_DIR}/utr_blocks.fa" 2>/dev/null || {
        echo "ERROR: bedtools getfasta failed. Check FASTA/GTF chromosome name compatibility."
        exit 1
    }


#Merge sequence from subsequent blocks
awk '
    /^>/ {
        # Header looks like: >GENE:TX:chr:start-end(+)
        split($0, a, ":")
        key = a[1] ":" a[2]   # keep only gene:transcript
        current = key
        next
    }
    {
        seq[current] = seq[current] $0
    }
    END {
        for (k in seq) {
            print ">" k
            print seq[k]
        }
    }
' "${TEMP_DIR}/utr_blocks.fa" \
> "${TEMP_DIR}/three_prime_utr.fa"

rm "${TEMP_DIR}/utr_blocks.fa"

# miRanda requires uppercase DNA sequences
seqkit seq -u "${UTR_FASTA}" > "${UTR_FASTA}.upper" && mv "${UTR_FASTA}.upper" "${UTR_FASTA}"

echo "  3' UTR FASTA sequences → ${UTR_FASTA}"

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "========================================================================"
echo "Phase VIII.1 – miRNA Target Prediction Preparation"
echo "========================================================================"
echo "miRNA sequences extracted  : ${MIRNA_COUNT}"
echo "3' UTR regions             : ${UTR_COUNT}"
echo ""
echo "Output files:"
echo "  ${MIRNA_FASTA_OUT}"
echo "  ${MIRNA_INDEXED_LIST}"
echo "  ${UTR_BED}"
echo "  ${UTR_FASTA}"
echo ""
echo "Next step:"
echo "  Number of array jobs for 8.2_miranda_array.sh : ${MIRNA_COUNT}"
echo "  Update #SBATCH --array=1-${MIRNA_COUNT} in 8.2_miranda_array.sh, then submit:"
echo "  sbatch scripts/8.2_miranda_array.sh"
echo "========================================================================"

printf "\n\nCompleted: 8.1-mirna_target_preparation\n\n"
pwd
date
