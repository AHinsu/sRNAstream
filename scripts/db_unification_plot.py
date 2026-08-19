#!/usr/bin/env python3
"""Plot deduplication statistics and length distributions for script 2.2.

Generates three figures:
  1. Per-file sequence counts bar chart (mature and precursor, side-by-side).
  2. Length distributions of the final deduped mature and precursor FASTAs.
  3. Before/after seqkit rmdup sequence counts (mature and precursor).

Environment variables:
    DB_OUTDIR              – directory where plots are written
    MATURE_SOURCES         – pipe-separated 'Label:path' for raw mature FASTAs
    PRECURSOR_SOURCES      – pipe-separated 'Label:path' for raw precursor FASTAs
    COMBINED_MATURE        – path to combined (pre-dedup) mature FASTA
    COMBINED_MATURE_DEDUP  – path to deduped mature FASTA
    COMBINED_PRECURSOR     – path to combined (pre-dedup) precursor FASTA
    COMBINED_PRECURSOR_DEDUP – path to deduped precursor FASTA

Usage:
    DB_OUTDIR=/path/to/db \\
        MATURE_SOURCES="label1:/path1.fa|label2:/path2.fa" \\
        PRECURSOR_SOURCES="label3:/path3.fa" \\
        COMBINED_MATURE=/path/combined_mature_all.fasta \\
        COMBINED_MATURE_DEDUP=/path/combined_mature_dedup.fasta \\
        COMBINED_PRECURSOR=/path/combined_precursor_all.fasta \\
        COMBINED_PRECURSOR_DEDUP=/path/combined_precursor_dedup.fasta \\
        python3 db_unification_plot.py
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd

db_outdir = os.environ['DB_OUTDIR']


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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


def count_seqs(fasta_path):
    """Count FASTA header lines (sequences) in a file."""
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        return 0
    with open(fasta_path) as fh:
        return sum(1 for line in fh if line.startswith('>'))


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


# ---------------------------------------------------------------------------
# Read env vars
# ---------------------------------------------------------------------------

mature_sources    = parse_sources('MATURE_SOURCES')
precursor_sources = parse_sources('PRECURSOR_SOURCES')

combined_mature        = os.environ.get('COMBINED_MATURE', '')
combined_mature_dedup  = os.environ.get('COMBINED_MATURE_DEDUP', '')
combined_precursor     = os.environ.get('COMBINED_PRECURSOR', '')
combined_precursor_dedup = os.environ.get('COMBINED_PRECURSOR_DEDUP', '')


# ---------------------------------------------------------------------------
# Figure 1: Per-file sequence counts (bar chart)
# ---------------------------------------------------------------------------

mat_labels = [lbl for lbl, _ in mature_sources]
mat_counts = [count_seqs(path) for _, path in mature_sources]
pre_labels = [lbl for lbl, _ in precursor_sources]
pre_counts = [count_seqs(path) for _, path in precursor_sources]

print("\nPer-file sequence counts:")
for lbl, cnt in zip(mat_labels, mat_counts):
    print(f"  Mature    – {lbl}: {cnt}")
for lbl, cnt in zip(pre_labels, pre_counts):
    print(f"  Precursor – {lbl}: {cnt}")

if mat_labels or pre_labels:
    n_bars = len(mat_labels) + len(pre_labels)
    fig, axes = plt.subplots(1, 2, figsize=(max(10, n_bars * 1.5 + 4), 5),
                             squeeze=False)
    colors = plt.get_cmap('tab10').colors

    for ax, labels, counts, title in [
        (axes[0, 0], mat_labels, mat_counts, 'Mature miRNA – sequences per input file'),
        (axes[0, 1], pre_labels, pre_counts, 'Precursor miRNA – sequences per input file'),
    ]:
        if labels:
            bar_colors = [colors[i % len(colors)] for i in range(len(labels))]
            bars = ax.bar(labels, counts, color=bar_colors, edgecolor='white', linewidth=0.5)
            for bar, count in zip(bars, counts):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() * 1.01,
                        str(count), ha='center', va='bottom', fontsize=9)
            ax.set_ylabel('Number of sequences', fontsize=11)
            ax.tick_params(axis='x', rotation=30)
            ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        else:
            ax.text(0.5, 0.5, 'No files provided', ha='center', va='center',
                    transform=ax.transAxes, fontsize=12, color='grey')
        ax.set_title(title, fontsize=12)

    plt.tight_layout()
    out_png = os.path.join(db_outdir, 'per_file_sequence_counts.png')
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_png}")


# ---------------------------------------------------------------------------
# Figure 2: Length distributions of the final deduped FASTA files
# ---------------------------------------------------------------------------

dedup_files = {
    'Mature (deduped)':    combined_mature_dedup,
    'Precursor (deduped)': combined_precursor_dedup,
}

dedup_lengths = {}
for lbl, path in dedup_files.items():
    lengths = read_lengths(path)
    if lengths:
        dedup_lengths[lbl] = lengths

if dedup_lengths:
    n_panels = len(dedup_lengths)
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 5), squeeze=False)
    colors = plt.get_cmap('tab10').colors

    for idx, (lbl, lengths) in enumerate(dedup_lengths.items()):
        ax = axes[0, idx]
        df = pd.DataFrame({'Length (nt)': lengths})
        length_counts = df['Length (nt)'].value_counts().sort_index()
        ax.bar(length_counts.index, length_counts.values,
               color=colors[idx % len(colors)], edgecolor='white', linewidth=0.5)
        ax.set_xlabel('Sequence length (nt)', fontsize=12)
        ax.set_ylabel('Number of sequences', fontsize=12)
        ax.set_title(f'Length distribution – {lbl}', fontsize=13)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        print(f"  {lbl}: {len(lengths)} sequences, "
              f"lengths {min(lengths)}–{max(lengths)} nt")

    plt.tight_layout()
    out_png = os.path.join(db_outdir, 'length_density_final.png')
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_png}")


# ---------------------------------------------------------------------------
# Figure 3: Before/after seqkit rmdup sequence counts
# ---------------------------------------------------------------------------

n_mature_all   = count_seqs(combined_mature)
n_mature_dedup = count_seqs(combined_mature_dedup)
n_pre_all      = count_seqs(combined_precursor)
n_pre_dedup    = count_seqs(combined_precursor_dedup)

print(f"\nMature:    {n_mature_all} combined → {n_mature_dedup} after seqkit rmdup")
print(f"Precursor: {n_pre_all} combined → {n_pre_dedup} after seqkit rmdup")

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

for ax, label, n_all, n_dedup in [
    (axes[0], 'Mature',    n_mature_all, n_mature_dedup),
    (axes[1], 'Precursor', n_pre_all,    n_pre_dedup),
]:
    categories = ['Before dedup\n(combined)', 'After dedup\n(seqkit rmdup)']
    counts = [n_all, n_dedup]
    bars = ax.bar(categories, counts,
                  color=['steelblue', 'darkorange'],
                  edgecolor='white', linewidth=0.5)
    ax.set_title(f'{label} miRNA – Deduplication', fontsize=13)
    ax.set_ylabel('Number of sequences', fontsize=11)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() * 1.01,
                str(count), ha='center', va='bottom', fontsize=11)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

plt.tight_layout()
out_png = os.path.join(db_outdir, 'dedup_stats.png')
fig.savefig(out_png, dpi=150)
plt.close(fig)
print(f"Saved: {out_png}")

print("\nUnification plots completed.")
