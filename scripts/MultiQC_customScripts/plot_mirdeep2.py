#!/usr/bin/env python3
"""
plot_mirdeep2.py – Standalone export: miRDeep2 known/novel bar chart

Reads mirdeep2/parsed/mirdeep2_all.tsv produced by 4.1_mirdeep2_analysis.sh
and exports a grouped bar chart showing known and novel miRNA counts at each
score threshold ≥1, ≥2, …, ≥10.

Environment variables
---------------------
OUTDIR (required) – pipeline results directory
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCORE_THRESHOLDS = list(range(1, 11))  # 1 through 10
KNOWN_COLOR = "#4e79a7"
NOVEL_COLOR = "#f28e2b"
BAR_WIDTH = 0.35


def main():
    outdir = os.environ.get("OUTDIR")
    if not outdir:
        print("ERROR: OUTDIR environment variable is not set.", file=sys.stderr)
        sys.exit(1)

    mqc_dir = os.path.join(outdir, "multiqc")
    os.makedirs(mqc_dir, exist_ok=True)

    tsv_path = os.path.join(outdir, "mirdeep2", "parsed", "mirdeep2_all.tsv")
    if not os.path.exists(tsv_path):
        print(f"WARNING: {tsv_path} not found – skipping miRDeep2 plot.", file=sys.stderr)
        return

    # Parse TSV: need 'score' and 'category' columns
    header = None
    rows = []
    with open(tsv_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if header is None:
                header = line.split("\t")
                continue
            parts = line.split("\t")
            row = dict(zip(header, parts))
            try:
                row["score"] = float(row.get("score", 0))
            except (ValueError, TypeError):
                row["score"] = 0.0
            rows.append(row)

    if not rows:
        print("WARNING: No rows parsed from mirdeep2_all.tsv – skipping plot.", file=sys.stderr)
        return

    # Count known and novel for each threshold
    known_counts = []
    novel_counts = []
    for thresh in SCORE_THRESHOLDS:
        above = [r for r in rows if r["score"] >= thresh]
        known_counts.append(sum(1 for r in above if r.get("category", "") == "known"))
        novel_counts.append(sum(1 for r in above if r.get("category", "") == "novel"))

    # Plot
    xs = np.arange(len(SCORE_THRESHOLDS))
    fig, ax = plt.subplots(figsize=(12, 5))

    bars_known = ax.bar(xs - BAR_WIDTH / 2, known_counts, BAR_WIDTH,
                        label="Known miRNA", color=KNOWN_COLOR, edgecolor="white")
    bars_novel = ax.bar(xs + BAR_WIDTH / 2, novel_counts, BAR_WIDTH,
                        label="Novel miRNA", color=NOVEL_COLOR, edgecolor="white")

    # Value labels above bars
    for bar in list(bars_known) + list(bars_novel):
        h = bar.get_height()
        if h > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                h + max(known_counts + novel_counts) * 0.01,
                str(int(h)),
                ha="center", va="bottom", fontsize=8,
            )

    ax.set_xticks(xs)
    ax.set_xticklabels([f"\u2265{t}" for t in SCORE_THRESHOLDS], fontsize=10)
    ax.set_xlabel("miRDeep2 score threshold", fontsize=11)
    ax.set_ylabel("Number of miRNA predictions", fontsize=11)
    ax.set_title(
        "miRDeep2: Known and Novel miRNA Predictions by Score Threshold (Phase IV)",
        fontsize=12,
    )
    ax.legend(fontsize=10)
    plt.tight_layout()

    out_path = os.path.join(mqc_dir, "mirdeep2_known_novel.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
