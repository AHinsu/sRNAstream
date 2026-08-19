#!/usr/bin/env python3
"""Plot sequence-length density distributions for miRNA database FASTA files.

Reads source file lists from environment variables set by the calling shell
script, along with the output directory.

Environment variables:
    DB_OUTDIR         – directory where plots are written
    MATURE_SOURCES    – pipe-separated 'Label:path' pairs for mature FASTAs
    PRECURSOR_SOURCES – pipe-separated 'Label:path' pairs for precursor FASTAs
    ENSEMBL_FA        – optional path to the Ensembl extracted FASTA

Usage:
    DB_OUTDIR=/path/to/db \\
        MATURE_SOURCES="miRBase_mature:/path1.fa|MirGeneDB_mature:/path2.fa" \\
        PRECURSOR_SOURCES="miRBase_hairpin:/path3.fa" \\
        python3 db_density_plot.py
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd

db_outdir = os.environ['DB_OUTDIR']


def parse_sources(env_var):
    """Parse pipe-separated 'Label:path' pairs from an environment variable.

    Returns an ordered dict of {label: path}.
    """
    raw = os.environ.get(env_var, '').strip()
    sources = {}
    if not raw:
        return sources
    for entry in raw.split('|'):
        entry = entry.strip()
        if ':' in entry:
            label, path = entry.split(':', 1)
            sources[label.strip()] = path.strip()
    return sources


# Build ordered dict of all FASTA files to plot.
# Order: Ensembl (if provided) → mature sources → precursor sources.
fasta_files = {}

ensembl_fa = os.environ.get('ENSEMBL_FA', '').strip()
if ensembl_fa:
    fasta_files['Ensembl'] = ensembl_fa

fasta_files.update(parse_sources('MATURE_SOURCES'))
fasta_files.update(parse_sources('PRECURSOR_SOURCES'))

def read_lengths(fasta_path):
    """Return list of sequence lengths from a FASTA file."""
    lengths = []
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        return lengths
    seq = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if seq:
                    lengths.append(len(''.join(seq)))
                seq = []
            else:
                seq.append(line)
        if seq:
            lengths.append(len(''.join(seq)))
    return lengths

all_lengths = {}
for label, fpath in fasta_files.items():
    lengths = read_lengths(fpath)
    if lengths:
        all_lengths[label] = lengths
        print(f"  {label}: {len(lengths)} sequences, lengths {min(lengths)}–{max(lengths)} nt")
    else:
        print(f"  {label}: no sequences found – skipping plot")

# Individual density plot per source
for label, lengths in all_lengths.items():
    fig, ax = plt.subplots(figsize=(8, 4))
    df = pd.DataFrame({'Length (nt)': lengths})
    length_counts = df['Length (nt)'].value_counts().sort_index()
    ax.bar(length_counts.index, length_counts.values, color='steelblue', edgecolor='white', linewidth=0.5)
    ax.set_xlabel('Sequence length (nt)', fontsize=12)
    ax.set_ylabel('Number of sequences', fontsize=12)
    ax.set_title(f'Sequence length distribution – {label.replace("_", " ")}', fontsize=13)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    plt.tight_layout()
    out_png = os.path.join(db_outdir, f'length_density_{label}.png')
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_png}")

# Combined overlay density plot for all sources
if all_lengths:
    fig, ax = plt.subplots(figsize=(10, 5))
    colors = plt.get_cmap('tab10').colors
    for i, (label, lengths) in enumerate(all_lengths.items()):
        color = colors[i % len(colors)]
        df = pd.DataFrame({'length': lengths})
        counts = df['length'].value_counts().sort_index()
        ax.plot(counts.index, counts.values, label=label.replace('_', ' '),
                color=color, linewidth=1.5, marker='o', markersize=3)
    ax.set_xlabel('Sequence length (nt)', fontsize=12)
    ax.set_ylabel('Number of sequences', fontsize=12)
    ax.set_title('Sequence length distribution – all databases', fontsize=13)
    ax.legend(fontsize=8, framealpha=0.7)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    plt.tight_layout()
    out_png = os.path.join(db_outdir, 'length_density_combined.png')
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"  Saved combined plot: {out_png}")

print("Density plots completed.")
