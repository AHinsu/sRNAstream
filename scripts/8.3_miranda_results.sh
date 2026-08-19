#!/usr/bin/env bash
#SBATCH --job-name=miranda_results
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


###############################################################################
# 8.3_miranda_results.sh - Filter miRanda predictions and build summary tables
#
# Description:
#   Phase VIII.3: Apply user-defined score and energy thresholds to the
#   all_predictions.tsv file produced by 8.2a_miranda_results_parsed.sh,
#   then build two summary tables from the filtered predictions:
#
#   Workflow:
#     1. Filter all_predictions.tsv by SCORE_THRESHOLD and ENERGY_THRESHOLD,
#        writing the passing rows to parsed/filtered_predictions.tsv.
#     2. Build a per-miRNA target summary (mirna_target_summary.tsv):
#        miRNA | Num_Targets | Target_mRNA_List (targets joined by "; ").
#     3. Build a per-miRNA gene summary (mirna_gene_summary.tsv):
#        miRNA | Num_Genes | Gene_ID_List (unique gene IDs joined by "; ",
#        where gene ID = substring of target ID before the first ":").
#
# Usage:
#   # Submit after 8.2a has finished:
#   sbatch scripts/8.3_miranda_results.sh
#
# Prerequisites:
#   - Run 8.2a_miranda_results_parsed.sh first to produce all_predictions.tsv.
#
# Outputs (in ${OUTDIR}/miranda/):
#   parsed/filtered_predictions.tsv  – Predictions passing score/energy filters
#   mirna_target_summary.tsv         – miRNA → Num_Targets + "; "-joined target list
#   mirna_gene_summary.tsv           – miRNA → Num_Genes  + "; "-joined gene ID list
###############################################################################

printf "\n\nStarted Job: 8.3-miranda_results\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 8.1, 8.2 and 8.2a)
OUTDIR="/storage/users/ahinsu/DD-microRNA/downstream/targets/miRNA"

# Minimum alignment score threshold (inclusive; predictions >= this value pass)
SCORE_THRESHOLD=150

# Maximum binding energy threshold in kcal/mol (inclusive; predictions <= this
# value pass; use a negative number, e.g. -20 means energy must be <= -20)
ENERGY_THRESHOLD=-25

# ============================================================================
# Resolve paths
# ============================================================================

ALL_PREDICTIONS="${OUTDIR}/miranda/all_predictions.tsv"
PARSED_DIR="${OUTDIR}/miranda/parsed"
FILTERED_TSV="${PARSED_DIR}/filtered_predictions.tsv"
# Add thresholds to summary filenames
SUMMARY_FILE="${PARSED_DIR}/target_summary_st${SCORE_THRESHOLD}_et${ENERGY_THRESHOLD}.tsv"
GENE_SUMMARY_FILE="${PARSED_DIR}/gene_summary_st${SCORE_THRESHOLD}_et${ENERGY_THRESHOLD}.tsv"
mkdir -p "${PARSED_DIR}"

if [ ! -f "${ALL_PREDICTIONS}" ]; then
    echo "ERROR: Input file not found: ${ALL_PREDICTIONS}"
    echo "       Run 8.2a_miranda_results_parsed.sh first."
    exit 1
fi

# ============================================================================
# Step 1: Filter all_predictions.tsv by score and energy thresholds
# ============================================================================

echo ""
echo "Step 1: Filtering predictions (score >= ${SCORE_THRESHOLD}, energy <= ${ENERGY_THRESHOLD} kcal/mol)..."

awk -v FS="\t" -v OFS="\t" -v score_thr="${SCORE_THRESHOLD}" -v energy_thr="${ENERGY_THRESHOLD}" '
    BEGIN { OFS="\t" }
    NR == 1 { print; next }
    $3 >= score_thr && $4 <= energy_thr { print }
' "${ALL_PREDICTIONS}" > "${FILTERED_TSV}"

FILTERED_COUNT=$(tail -n +2 "${FILTERED_TSV}" | wc -l)
echo "  ${FILTERED_COUNT} predictions passed filters → ${FILTERED_TSV}"

# ============================================================================
# Step 2: Build per-miRNA target summary table
# ============================================================================

echo ""
echo "Step 2: Building per-miRNA target summary table..."

echo -e "miRNA_ID\tNum_Targets\tTarget_mRNA_List" > "${SUMMARY_FILE}"

awk -v FS="\t" -v OFS="\t" '
    NR > 1 {
        mir    = $1
        target = $2
        if (mir in targets) {
            targets[mir] = targets[mir] "; " target
            counts[mir]++
        } else {
            targets[mir] = target
            counts[mir]  = 1
        }
    }
    END {
        for (mir in targets) {
            print mir "\t" counts[mir] "\t" targets[mir]
        }
    }
' "${FILTERED_TSV}" | sort >> "${SUMMARY_FILE}"

MIRNA_WITH_TARGETS=$(tail -n +2 "${SUMMARY_FILE}" | wc -l)
echo "  ${MIRNA_WITH_TARGETS} miRNAs with target predictions → ${SUMMARY_FILE}"

# ============================================================================
# Step 3: Build per-miRNA gene summary table
# ============================================================================

echo ""
echo "Step 3: Building per-miRNA gene summary table..."

echo -e "miRNA_ID\tNum_Genes\tGene_ID_List" > "${GENE_SUMMARY_FILE}"

# Gene ID = portion of target ID before the first ":".
# Gene IDs are deduplicated per miRNA and joined with "; ".
awk -v FS="\t" -v OFS="\t" '
    NR > 1 {
        mir    = $1
        target = $2
        gene   = (match(target, /:/) ? substr(target, 1, RSTART - 1) : target)
        key    = mir SUBSEP gene
        if (!(key in seen)) {
            seen[key] = 1
            if (mir in genes) {
                genes[mir]      = genes[mir] "; " gene
                gene_counts[mir]++
            } else {
                genes[mir]      = gene
                gene_counts[mir] = 1
            }
        }
    }
    END {
        for (mir in genes) {
            print mir "\t" gene_counts[mir] "\t" genes[mir]
        }
    }
' "${FILTERED_TSV}" | sort >> "${GENE_SUMMARY_FILE}"

MIRNA_WITH_GENES=$(tail -n +2 "${GENE_SUMMARY_FILE}" | wc -l)
echo "  ${MIRNA_WITH_GENES} miRNAs with gene predictions → ${GENE_SUMMARY_FILE}"

# ============================================================================
# Final summary
# ============================================================================

echo ""
echo "========================================================================"
echo "Phase VIII.3 – miRNA Target Prediction Filtering and Summarisation"
echo "========================================================================"
echo "Score threshold            : >= ${SCORE_THRESHOLD}"
echo "Energy threshold           : <= ${ENERGY_THRESHOLD} kcal/mol"
echo "Predictions after filter   : ${FILTERED_COUNT}"
echo "miRNAs with targets        : ${MIRNA_WITH_TARGETS}"
echo "miRNAs with genes          : ${MIRNA_WITH_GENES}"
echo ""
echo "Output files:"
echo "  ${FILTERED_TSV}"
echo "  ${SUMMARY_FILE}"
echo "  ${GENE_SUMMARY_FILE}"
echo "========================================================================"

printf "\n\nCompleted: 8.3-miranda_results\n\n"
pwd
date
