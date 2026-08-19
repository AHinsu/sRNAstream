#!/usr/bin/env python3
"""Plot deduplication statistics from combined mature and hairpin FASTA files.

Reads COMBINED_MATURE, COMBINED_MATURE_DEDUP, PRINSEQ_DEDUP_MATURE,
COMBINED_HAIRPIN, COMBINED_HAIRPIN_DEDUP, and PRINSEQ_DEDUP_HAIRPIN from
environment variables set by the calling shell script.

Usage:
    COMBINED_MATURE=... COMBINED_MATURE_DEDUP=... PRINSEQ_DEDUP_MATURE=... \\
    COMBINED_HAIRPIN=... COMBINED_HAIRPIN_DEDUP=... PRINSEQ_DEDUP_HAIRPIN=... \\
    python3 db_combined_dedup_plot.py
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def count_seqs(fasta_path):
    """Count sequences (FASTA headers) in a FASTA file."""
    if not fasta_path or not os.path.exists(fasta_path) or \
            os.path.getsize(fasta_path) == 0:
        return 0
    with open(fasta_path) as fh:
        return sum(1 for line in fh if line.startswith('>'))


combined_mature        = os.environ.get('COMBINED_MATURE', '')
combined_mature_dedup  = os.environ.get('COMBINED_MATURE_DEDUP', '')
prinseq_dedup_mature   = os.environ.get('PRINSEQ_DEDUP_MATURE', '')
combined_hairpin       = os.environ.get('COMBINED_HAIRPIN', '')
combined_hairpin_dedup = os.environ.get('COMBINED_HAIRPIN_DEDUP', '')
prinseq_dedup_hairpin  = os.environ.get('PRINSEQ_DEDUP_HAIRPIN', '')

# Determine output directory from the path of the first provided file
db_outdir = os.path.dirname(combined_mature or combined_hairpin) or '.'

# Count sequences before and after deduplication
n_mature_all            = count_seqs(combined_mature)
n_mature_seqkit         = count_seqs(combined_mature_dedup)
n_mature_prinseq        = count_seqs(prinseq_dedup_mature)
n_hairpin_all           = count_seqs(combined_hairpin)
n_hairpin_seqkit        = count_seqs(combined_hairpin_dedup)
n_hairpin_prinseq       = count_seqs(prinseq_dedup_hairpin)

print(f"Mature sequences:  {n_mature_all} total → {n_mature_seqkit} after seqkit rmdup"
      f" → {n_mature_prinseq} after prinseq dedup")
print(f"Hairpin sequences: {n_hairpin_all} total → {n_hairpin_seqkit} after seqkit rmdup"
      f" → {n_hairpin_prinseq} after prinseq dedup")

# Bar chart: before dedup / after seqkit rmdup / after prinseq dedup
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, label, n_all, n_seqkit, n_prinseq in [
    (axes[0], 'Mature',  n_mature_all,  n_mature_seqkit,  n_mature_prinseq),
    (axes[1], 'Hairpin', n_hairpin_all, n_hairpin_seqkit, n_hairpin_prinseq),
]:
    categories = [
        'Before dedup\n(combined)',
        'After dedup\n(seqkit rmdup)',
        'After dedup\n(prinseq)',
    ]
    counts = [n_all, n_seqkit, n_prinseq]
    bars = ax.bar(categories, counts,
                  color=['steelblue', 'darkorange', 'seagreen'],
                  edgecolor='white', linewidth=0.5)
    ax.set_title(f'{label} miRNA – Deduplication', fontsize=13)
    ax.set_ylabel('Number of sequences', fontsize=11)
    for bar, count in zip(bars, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() * 1.01,
            str(count),
            ha='center', va='bottom', fontsize=10,
        )
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

plt.tight_layout()
out_png = os.path.join(db_outdir, 'combined_dedup_stats.png')
fig.savefig(out_png, dpi=150)
plt.close(fig)
print(f"Saved: {out_png}")
print("Combined deduplication plot completed.")
