#!/usr/bin/env bash
#SBATCH --job-name=genome_index
#SBATCH --output=/storage/users/ahinsu/sbatch_logs/job-%j.%x.out
#SBATCH --error=/storage/users/ahinsu/sbatch_logs/job-%j.%x.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G


###############################################################################
# 3.1_genome_index.sh - Map cleaned reads to the Ensembl reference genome
#
# Description:
#   Phase III.1: Prepares genome index with Bowtie for mapping
#
# Usage:
#   sbatch scripts/3.1_genome_index.sh
#
################################################################################

printf "\n\nStarted Job: 3.2-genome_alignment\n\n"
pwd
date
printf "\n\n"

# ============================================================================
# USER CONFIGURATION – edit these variables for your project
# ============================================================================

# Reference genome FASTA file (from Ensembl)
GENOME_FASTA="/storage/users/ahinsu/DD-microRNA/reference/Bos_taurus.chromsomal.fa"

# Output directory
OUTDIR="/storage/users/ahinsu/DD-microRNA/reference/index"

# Conda environment name
CONDA_ENV="smallRNA"

# ============================================================================
# Load required modules and activate conda environment
# ============================================================================
module load apps/anaconda-4.7.12.tcl
eval "$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
unset PYTHONPATH

mkdir -p ${OUTDIR}

#Build index

bowtie-build --threads ${SLURM_CPUS_PER_TASK:-4}  \
    ${GENOME_FASTA} ${OUTDIR}/genome


printf "\n\nCompleted: 3.1-genome_index\n\n"
pwd
date
