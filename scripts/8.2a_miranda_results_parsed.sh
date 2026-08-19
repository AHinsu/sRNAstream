#!/usr/bin/env bash
#SBATCH --job-name=miranda_results_parsed
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%A.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%A.%x.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G


###############################################################################
# 8.2a_miranda_results_parsed.sh - Parse, merge and filter miRanda target predictions
#
# Description:
#   Phase VIII.2a: Aggregate the per-miRNA raw miRanda output files produced by
#   8.2_miranda_array.sh into consolidated result tables and summary statistics.
#   Run this script after all 8.2 array tasks have finished.
#
#   This script produces the same outputs as 8.3_miranda_results.sh with the
#   following corrections:
#     - Target IDs have the leading ">" (from FASTA headers) stripped.
#     - Both Max_Score and Max_Energy are extracted from miRanda hit lines
#       and reported with correct column headers.
#
#   Workflow:
#     1. Parse each raw miRanda output file and extract miRNA name, target ID,
#        max score and max binding energy into a single flat TSV (all_predictions.tsv).
#     2. Build a per-miRNA summary table listing the number of predicted targets
#        and a semicolon-separated target list (mirna_target_summary.tsv).
#     3. Generate statistics: total predictions, energy distribution, top 10
#        miRNAs by number of targets (target_prediction_statistics.txt).
#     4. Write per-miRNA detailed TSV files in the parsed/ sub-directory.
#     5. Remove temporary files created by 8.1_mirna_target_preparation.sh.
#
#
# Prerequisites:
#   - Run 8.2_miranda_array.sh for ALL miRNAs first.
#
# Outputs (in ${OUTDIR}/miranda/):
#   all_predictions.tsv             – All miRNA-target predictions (miRNA, target, max score, max energy)
#   mirna_target_summary.tsv        – miRNA → number of targets + semicolon-separated list
#   target_prediction_statistics.txt – Summary statistics and energy distribution
#   parsed/<miRNA>.detailed.tsv     – Per-miRNA predictions sorted by max energy
###############################################################################

printf "\n\nStarted Job: 8.2a-miranda_results_parsed\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Output directory (must match OUTDIR used in 8.1 and 8.2)
OUTDIR="/storage/users/ahinsu/DD-microRNA/downstream/targets/miRNA"

# miRanda parameters used in 8.2 (for the statistics report only)
MIRANDA_SCORE=140
MIRANDA_ENERGY=-20

# Conda environment name
CONDA_ENV="miranda"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

RAW_DIR="${OUTDIR}/miranda/raw"
PARSED_DIR="${OUTDIR}/miranda/parsed"
mkdir -p "${PARSED_DIR}"

if [ ! -d "${RAW_DIR}" ]; then
    echo "ERROR: Raw miRanda output directory not found: ${RAW_DIR}"
    echo "       Run 8.2_miranda_array.sh first."
    exit 1
fi

RAW_COUNT=$(find "${RAW_DIR}" -maxdepth 1 -name '*.txt' -type f | wc -l)
if [ "${RAW_COUNT}" -eq 0 ]; then
    echo "ERROR: No miRanda output files found in ${RAW_DIR}"
    echo "       Run 8.2_miranda_array.sh first."
    exit 1
fi
echo "Found ${RAW_COUNT} miRanda output files in ${RAW_DIR}"

# ============================================================================
# Step 1: Parse raw miRanda output files
# ============================================================================

echo ""
echo "Step 1: Parsing miRanda output files..."

PARSED_TSV="${OUTDIR}/miranda/all_predictions.tsv"
echo -e "miRNA\tTarget_ID\tMax_Score\tMax_Energy_kcal_mol" > "${PARSED_TSV}"

for miranda_file in "${RAW_DIR}"/*.txt; do
    [ -f "${miranda_file}" ] || continue
    mirna_name=$(basename "${miranda_file}" .txt)

    # Extract target hits from miRanda output.
    # miRanda hit lines start with ">>" and contain tab-separated fields:
    #   >>miRNA  Target  TotalScore  TotalEnergy  MaxScore  MaxEnergy  Strand  ...
    # The target ID may carry a leading ">" from the FASTA header – strip it.
    awk -v mir="${mirna_name}" '
        BEGIN { OFS="\t" }
        /^>>/ {
            target = $2
            sub(/^>/, "", target)
            score  = $5
            energy = $6
            if (target != "" && score != "" && energy != "") {
                print mir, target, score, energy
            }
        }
    ' "${miranda_file}" >> "${PARSED_TSV}"
done

PREDICTION_COUNT=$(tail -n +2 "${PARSED_TSV}" | wc -l)
MIRNA_WITH_TARGETS=$(cut -f1 "${PARSED_TSV}" | sort -u | wc -l)
echo "  Parsed ${PREDICTION_COUNT} miRNA-target predictions → ${PARSED_TSV}"

# ============================================================================
# Step 2: Build per-miRNA summary table
# ============================================================================

#echo ""
#echo "Step 2: Building per-miRNA summary table..."
#
#SUMMARY_FILE="${OUTDIR}/miranda/mirna_target_summary.tsv"
#echo -e "miRNA_ID\tNum_Targets\tTarget_mRNA_List" > "${SUMMARY_FILE}"
#
#awk '
#    NR > 1 {
#        mir    = $1
#        target = $2
#        if (mir in targets) {
#            targets[mir] = targets[mir] ";" target
#            counts[mir]++
#        } else {
#            targets[mir] = target
#            counts[mir]  = 1
#        }
#    }
#    END {
#        for (mir in targets) {
#            print mir "\t" counts[mir] "\t" targets[mir]
#        }
#    }
#' "${PARSED_TSV}" | sort >> "${SUMMARY_FILE}"
#
#MIRNA_WITH_TARGETS=$(tail -n +2 "${SUMMARY_FILE}" | wc -l)
#echo "  ${MIRNA_WITH_TARGETS} miRNAs with target predictions → ${SUMMARY_FILE}"

# ============================================================================
# Step 3: Generate statistics report
# ============================================================================

echo ""
echo "Step 3: Generating statistics report..."

STATS_FILE="${OUTDIR}/miranda/target_prediction_statistics.txt"

AVG_TARGETS="N/A"
if [ "${MIRNA_WITH_TARGETS}" -gt 0 ]; then
    AVG_TARGETS=$(echo "scale=2; ${PREDICTION_COUNT} / ${MIRNA_WITH_TARGETS}" | bc -l 2>/dev/null || echo "N/A")
fi

cat > "${STATS_FILE}" << STATS_EOF
================================================================================
miRNA Target Prediction Statistics
================================================================================

Pipeline Parameters:
  miRanda Score Threshold:    ${MIRANDA_SCORE}
  miRanda Energy Threshold:   ${MIRANDA_ENERGY} kcal/mol

Prediction Results:
  Total Predictions:          ${PREDICTION_COUNT}
  miRNAs with Targets:        ${MIRNA_WITH_TARGETS}
  Average Targets per miRNA:  ${AVG_TARGETS}

Top 10 miRNAs by Number of Targets:
STATS_EOF

tail -n +2 "${SUMMARY_FILE}" | sort -t$'\t' -k2 -rn | head -10 | \
    awk '{printf "  %d. %s\t%s targets\n", NR, $1, $2}' >> "${STATS_FILE}"

cat >> "${STATS_FILE}" << STATS_EOF

Energy Distribution of Predictions:
STATS_EOF

awk 'NR > 1 { print $4 }' "${PARSED_TSV}" | sort -n | \
    awk '{
        if      ($1 <= -25) vs++
        else if ($1 <= -20) s++
        else if ($1 <= -15) m++
        else if ($1 <= -10) w++
        else               vw++
    }
    END {
        printf "  Very Strong (<=-25 kcal/mol)    : %d\n", vs ? vs : 0
        printf "  Strong (-24 to -20 kcal/mol)    : %d\n", s  ? s  : 0
        printf "  Moderate (-20 to -15 kcal/mol)  : %d\n", m  ? m  : 0
        printf "  Weak (-15 to -10 kcal/mol)      : %d\n", w  ? w  : 0
        printf "  Very Weak (> -10 kcal/mol)      : %d\n", vw ? vw : 0
    }' >> "${STATS_FILE}"



cat >> "${STATS_FILE}" << STATS_EOF

Output Files:
  - all_predictions.tsv          : All miRNA-target predictions (miRNA, target, max score, max energy)
  - mirna_target_summary.tsv     : miRNA → number of targets + semicolon-separated list
  - parsed/<miRNA>.detailed.tsv  : Per-miRNA predictions sorted by max binding energy

================================================================================

Stats and numbers from the R file
--------------------------------------------------------------------------------

STATS_EOF

cat "${STATS_FILE}"

# ============================================================================
# Step 4: Write per-miRNA detailed TSV files
# ============================================================================

echo ""
echo "Step 4: Writing per-miRNA detailed files..."

while IFS=$'\t' read -r mirna_id num_targets target_list; do
    [ "${mirna_id}" = "miRNA_ID" ] && continue
    detailed_file="${PARSED_DIR}/${mirna_id}.detailed.tsv"
    echo -e "Target_ID\tMax_Score\tMax_Energy_kcal_mol" > "${detailed_file}"
    awk -v mir="${mirna_id}" '$1 == mir { print $2 "\t" $3 "\t" $4 }' "${PARSED_TSV}" \
        | sort -k3 -n >> "${detailed_file}"
done < "${SUMMARY_FILE}"

echo "  Per-miRNA detailed files written → ${PARSED_DIR}/"

# ============================================================================
# Step 5: Remove temporary preparation files
# ============================================================================

echo ""
echo "Step 5: Removing temporary files..."

TEMP_DIR="${OUTDIR}/temp"
if [ -d "${TEMP_DIR}" ]; then
    rm -rf "${TEMP_DIR}"
    echo "  Removed ${TEMP_DIR}"
else
    echo "  Temporary directory not found (already removed?): ${TEMP_DIR}"
fi

# ============================================================================
# Final summary
# ============================================================================

echo ""
echo "========================================================================"
echo "Phase VIII.2a – miRNA Target Prediction Results (Parsed)"
echo "========================================================================"
echo "Total predictions          : ${PREDICTION_COUNT}"
echo "miRNAs with targets        : ${MIRNA_WITH_TARGETS}"
echo "Average targets per miRNA  : ${AVG_TARGETS}"
echo ""
echo ""
echo "Running R script for plotting and summarising"
echo ""
echo ""
module load apps/R-4.2.1.tcl

Rscript /storage/users/ahinsu/DD-microRNA/scripts/miranda_summary.R "${PARSED_TSV}" >> "${STATS_FILE}" #The last part is if you want to pass arguments



printf "\n\nCompleted: 8.2a-miranda_results_parsed\n\n"
pwd
date
