#!/usr/bin/env python3
"""Write a summary report for Phase II.1 database preparation.

Reads DB_OUTDIR, GENOME_FASTA, GTF_FILE, MATURE_SOURCES, PRECURSOR_SOURCES,
and ENSEMBL_FA from environment variables set by the calling shell script.

Environment variables:
    DB_OUTDIR         – directory where the summary is written
    GENOME_FASTA      – path to Ensembl genome FASTA (optional)
    GTF_FILE          – path to Ensembl GTF (optional)
    MATURE_SOURCES    – pipe-separated 'Label:path' pairs for mature FASTAs
    PRECURSOR_SOURCES – pipe-separated 'Label:path' pairs for precursor FASTAs
    ENSEMBL_FA        – optional path to Ensembl extracted FASTA

Usage:
    DB_OUTDIR=/path/to/db \\
        MATURE_SOURCES="miRBase_mature:/path1.fa" \\
        PRECURSOR_SOURCES="miRBase_hairpin:/path2.fa" \\
        python3 db_summary.py
"""
import os
import glob
from datetime import datetime

db_outdir    = os.environ['DB_OUTDIR']
genome_fasta = os.environ.get('GENOME_FASTA', '(not provided)')
gtf_file     = os.environ.get('GTF_FILE',     '(not provided)')


def parse_sources(env_var):
    """Parse pipe-separated 'Label:path' pairs from an environment variable."""
    raw = os.environ.get(env_var, '').strip()
    sources = []
    if not raw:
        return sources
    for entry in raw.split('|'):
        entry = entry.strip()
        if ':' in entry:
            label, path = entry.split(':', 1)
            sources.append((label.strip(), path.strip()))
    return sources

def count_seqs(fasta_path):
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        return 0
    with open(fasta_path) as fh:
        return sum(1 for line in fh if line.startswith('>'))

summary_path = os.path.join(db_outdir, 'db_preparation_summary.txt')

with open(summary_path, 'w') as fh:
    fh.write("=" * 80 + "\n")
    fh.write("Phase II.1 – Database Preparation Summary\n")
    fh.write(f"Run date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    fh.write("=" * 80 + "\n\n")

    fh.write("Input files:\n")
    fh.write(f"  Genome FASTA : {genome_fasta}\n")
    fh.write(f"  GTF          : {gtf_file}\n\n")

    # Collect all sources from env vars
    ensembl_fa = os.environ.get('ENSEMBL_FA', '').strip()
    mature_sources    = parse_sources('MATURE_SOURCES')
    precursor_sources = parse_sources('PRECURSOR_SOURCES')
    if ensembl_fa:
        # Ensure Ensembl is included exactly once: prefer ENSEMBL_FA over any
        # Ensembl entry that may already be present in PRECURSOR_SOURCES.
        precursor_sources = [
            (label, path) for (label, path) in precursor_sources
            if label != 'Ensembl'
        ]
        precursor_sources = [('Ensembl', ensembl_fa)] + precursor_sources
    fh.write("Mature sequence counts per source:\n")
    for label, path in mature_sources:
        n = count_seqs(path)
        fh.write(f"  {label:<35} {n:>6} sequences   {path}\n")

    fh.write("\nPrecursor sequence counts per source:\n")
    for label, path in precursor_sources:
        n = count_seqs(path)
        fh.write(f"  {label:<35} {n:>6} sequences   {path}\n")

    fh.write("\nPlots:\n")
    for png in sorted(glob.glob(os.path.join(db_outdir, 'length_density_*.png'))):
        fh.write(f"  {png}\n")

    fh.write("\nDuplication statistics:\n")
    for fname in ['duplication_stats_mature.txt', 'duplication_stats_precursor.txt']:
        fh.write(f"  {os.path.join(db_outdir, fname)}\n")

    fh.write("\nCombined deduplication:\n")
    for fname in [
        'combined_mature_all.fasta', 'combined_mature_dedup.fasta',
        'prinseq_dedup_mature.fasta',
        'combined_hairpin_all.fasta', 'combined_hairpin_dedup.fasta',
        'prinseq_dedup_hairpin.fasta',
        'combined_dedup_stats.png',
    ]:
        fpath = os.path.join(db_outdir, fname)
        if os.path.exists(fpath):
            fh.write(f"  {fpath}\n")

    fh.write("\nprinseq-lite.pl statistics:\n")
    for fname in [
        'prinseq_combined_mature_all.txt', 'prinseq_combined_mature_dedup.txt',
        'prinseq_combined_mature_prinseqDedup.txt',
        'prinseq_combined_hairpin_all.txt', 'prinseq_combined_hairpin_dedup.txt',
        'prinseq_combined_hairpin_prinseqDedup.txt',
    ]:
        fpath = os.path.join(db_outdir, fname)
        if os.path.exists(fpath):
            fh.write(f"  {fpath}\n")

    fh.write("\nNext steps:\n")
    fh.write("  Run 2.2_database_unification.sh to further combine and deduplicate database files.\n")

print(open(summary_path).read())
