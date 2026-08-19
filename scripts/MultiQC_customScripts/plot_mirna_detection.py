#!/usr/bin/env python3
"""
plot_mirna_detection.py – Standalone export: distinct miRNA detected per sample

Reads mirna_counts/*_mirna_stats.txt produced by 5.2_mirna_quantification.sh
and exports a grouped bar chart showing the number of distinct miRNAs detected
per sample at read-count thresholds ≥1, ≥10, and ≥100.

Environment variables
---------------------
OUTDIR   (required) – pipeline results directory
MANIFEST (optional) – used to sort samples in manifest order
"""

import os
import sys
import re
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Reuse shared helper from plot_contamination (same OUTDIR env convention)
THRESHOLDS = [1, 10, 100]
TH_LABELS = ["\u22651", "\u226510", "\u2265100"]
TH_COLORS = ["#4e79a7", "#f28e2b", "#e15759"]


def load_manifest_order(manifest_path):
    samples = []
    if not manifest_path or not os.path.exists(manifest_path):
        return samples
    with open(manifest_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.lower().startswith("sample"):
                continue
            samples.append(line.split("\t")[0])
    return samples


def parse_distinct_counts(fpath):
    """Return (d1, d10, d100) distinct detected counts from a *_stats.txt file."""
    d1 = d10 = d100 = 0
    with open(fpath, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            # Handles both "≥1 read" and "≥1 reads" (miRNA uses "read", others "reads")
            m = re.search(r"\u22651\s+read", line)
            if m:
                val = re.search(r":\s*([\d,]+)", line)
                if val:
                    d1 = int(val.group(1).replace(",", ""))
            m = re.search(r"\u226510\s+reads", line)
            if m:
                val = re.search(r":\s*([\d,]+)", line)
                if val:
                    d10 = int(val.group(1).replace(",", ""))
            m = re.search(r"\u2265100\s+reads", line)
            if m:
                val = re.search(r":\s*([\d,]+)", line)
                if val:
                    d100 = int(val.group(1).replace(",", ""))
    return d1, d10, d100


def plot_distinct_detected(
    outdir, subdirectory, filename_suffix, sample_label, rna_label, out_filename
):
    """Core plotting function reused by all RNA-type detection scripts."""
    mqc_dir = os.path.join(outdir, "multiqc")
    os.makedirs(mqc_dir, exist_ok=True)

    stats_files = sorted(
        glob.glob(os.path.join(outdir, subdirectory, f"*{filename_suffix}"))
    )
    if not stats_files:
        print(
            f"WARNING: No {sample_label} stats files found in {subdirectory}/ – skipping plot.",
            file=sys.stderr,
        )
        return

    manifest = os.environ.get("MANIFEST", "")
    data = {}
    for fpath in stats_files:
        bn = os.path.basename(fpath)
        sample = re.sub(re.escape(filename_suffix) + "$", "", bn)
        d1, d10, d100 = parse_distinct_counts(fpath)
        data[sample] = (d1, d10, d100)

    if not data:
        return

    manifest_order = load_manifest_order(manifest)
    if manifest_order:
        samples = [s for s in manifest_order if s in data]
        samples += [s for s in sorted(data) if s not in samples]
    else:
        samples = sorted(data.keys())

    d1_vals   = [data[s][0] for s in samples]
    d10_vals  = [data[s][1] for s in samples]
    d100_vals = [data[s][2] for s in samples]

    xs = np.arange(len(samples))
    bar_w = 0.25
    fig_width = max(10, len(samples) * 0.6)
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    ax.bar(xs - bar_w, d1_vals,   bar_w, label=f"{TH_LABELS[0]} read",  color=TH_COLORS[0], edgecolor="white")
    ax.bar(xs,         d10_vals,  bar_w, label=f"{TH_LABELS[1]} reads", color=TH_COLORS[1], edgecolor="white")
    ax.bar(xs + bar_w, d100_vals, bar_w, label=f"{TH_LABELS[2]} reads", color=TH_COLORS[2], edgecolor="white")

    ax.set_xticks(xs)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"Distinct {rna_label} detected", fontsize=11)
    ax.set_title(
        f"Distinct {rna_label} Detected per Sample at Read-Count Thresholds",
        fontsize=12,
    )
    ax.legend(title="Min reads", fontsize=9)
    plt.tight_layout()

    out_path = os.path.join(mqc_dir, out_filename)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def main():
    outdir = os.environ.get("OUTDIR")
    if not outdir:
        print("ERROR: OUTDIR environment variable is not set.", file=sys.stderr)
        sys.exit(1)

    plot_distinct_detected(
        outdir=outdir,
        subdirectory="mirna_counts",
        filename_suffix="_mirna_stats.txt",
        sample_label="miRNA",
        rna_label="miRNAs",
        out_filename="mirna_distinct_detected.png",
    )


if __name__ == "__main__":
    main()
