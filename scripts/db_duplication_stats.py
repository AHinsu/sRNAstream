#!/usr/bin/env python3
"""Compute duplication statistics for mature and precursor miRNA FASTA files.

Reads source file lists from environment variables set by the calling shell
script, along with the output directory.

Environment variables:
    DB_OUTDIR         – directory where reports are written
    MATURE_SOURCES    – pipe-separated 'Label:path' pairs for mature FASTAs
    PRECURSOR_SOURCES – pipe-separated 'Label:path' pairs for precursor FASTAs
    ENSEMBL_FA        – optional path to Ensembl FASTA (prepended to precursor)

Usage:
    DB_OUTDIR=/path/to/db \\
        MATURE_SOURCES="miRBase_mature:/path1.fa|MirGeneDB_mature:/path2.fa" \\
        PRECURSOR_SOURCES="miRBase_hairpin:/path3.fa" \\
        python3 db_duplication_stats.py
"""
import os
from collections import Counter

db_outdir = os.environ['DB_OUTDIR']


def parse_sources(env_var):
    """Parse pipe-separated 'Label:path' pairs from an environment variable.

    Returns a list of (label, path) tuples.
    """
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

def parse_fasta(fasta_path):
    """Return dict of seq_id -> sequence (uppercase, T-based)."""
    seqs = {}
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        return seqs
    current_id, current_seq = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(current_seq).upper()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = ''.join(current_seq).upper()
    return seqs

def duplication_report(label, fasta_path):
    seqs = parse_fasta(fasta_path)
    if not seqs:
        return None

    total = len(seqs)
    seq_counts = Counter(seqs.values())
    unique = sum(1 for v in seq_counts.values() if v == 1)
    duplicated_seqs = sum(1 for v in seq_counts.values() if v > 1)
    total_dup_entries = sum(v for v in seq_counts.values() if v > 1)
    dup_entries = total - unique  # entries that are duplicates of something

    report = {
        'label': label,
        'path': fasta_path,
        'total_entries': total,
        'unique_sequences': len(seq_counts),
        'fully_unique': unique,
        'duplicated_sequence_groups': duplicated_seqs,
        'duplicate_entries': dup_entries,
        'duplication_rate_pct': round((dup_entries / total * 100) if total > 0 else 0, 2),
    }
    return report

# Build source lists from environment variables.
# Mature sources come from MATURE_SOURCES; precursor sources from PRECURSOR_SOURCES
# with Ensembl prepended if ENSEMBL_FA is set.
mature_sources = parse_sources('MATURE_SOURCES')

precursor_sources = parse_sources('PRECURSOR_SOURCES')
ensembl_fa = os.environ.get('ENSEMBL_FA', '').strip()
if ensembl_fa:
    # Avoid adding a duplicate Ensembl entry if it is already present
    has_ensembl = any(
        label == 'Ensembl' and path == ensembl_fa
        for label, path in precursor_sources
    )
    if not has_ensembl:
        precursor_sources = [('Ensembl', ensembl_fa)] + precursor_sources

def write_duplication_stats(sources, out_path):
    reports = [duplication_report(label, path) for label, path in sources]
    reports = [r for r in reports if r is not None]

    with open(out_path, 'w') as fh:
        fh.write("=" * 70 + "\n")
        fh.write("Duplication Statistics\n")
        fh.write("=" * 70 + "\n\n")
        for r in reports:
            fh.write(f"Source: {r['label']}\n")
            fh.write(f"  File:                       {r['path']}\n")
            fh.write(f"  Total entries:              {r['total_entries']}\n")
            fh.write(f"  Unique sequences:           {r['unique_sequences']}\n")
            fh.write(f"  Fully unique (no dups):     {r['fully_unique']}\n")
            fh.write(f"  Duplicate sequence groups:  {r['duplicated_sequence_groups']}\n")
            fh.write(f"  Duplicate entries:          {r['duplicate_entries']}\n")
            fh.write(f"  Duplication rate:           {r['duplication_rate_pct']} %\n")
            fh.write("\n")

        # Cross-database deduplication summary
        if len(reports) > 1:
            fh.write("-" * 70 + "\n")
            fh.write("Cross-database deduplication (combined unique sequences):\n")
            all_seqs = set()
            per_source = {}
            for label, path in sources:
                parsed = parse_fasta(path)
                seq_set = set(parsed.values())
                per_source[label] = seq_set
                all_seqs.update(seq_set)
            fh.write(f"  Total unique sequences across all sources: {len(all_seqs)}\n")
            for label, seq_set in per_source.items():
                fh.write(f"  {label}: {len(seq_set)} unique sequences\n")

    print(f"Written: {out_path}")
    for r in reports:
        print(f"  {r['label']}: {r['total_entries']} entries, "
              f"{r['unique_sequences']} unique, "
              f"{r['duplication_rate_pct']}% duplication rate")

print("\nMature sequence duplication statistics:")
write_duplication_stats(
    mature_sources,
    os.path.join(db_outdir, 'duplication_stats_mature.txt')
)

print("\nPrecursor sequence duplication statistics:")
write_duplication_stats(
    precursor_sources,
    os.path.join(db_outdir, 'duplication_stats_precursor.txt')
)
