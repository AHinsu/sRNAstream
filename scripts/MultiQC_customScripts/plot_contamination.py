#!/usr/bin/env python3
"""
plot_contamination.py – Standalone export: contamination filtering stacked bar chart

Reads *_contamination_stats.txt files produced by 1.5b_contamination_filtering.sh
and writes a publication-ready stacked bar chart to ${OUTDIR}/multiqc/.

Environment variables
---------------------
OUTDIR   (required) – pipeline results directory
MANIFEST (optional) – used only to sort samples in manifest order
"""

import os
import re
import glob
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
RNA_TYPES = ["rRNA", "tRNA", "snRNA", "snoRNA", "piRNA", "other", "clean"]

COLORS = {
    "rRNA":   "#e15759",
    "tRNA":   "#f28e2b",
    "snRNA":  "#76b7b2",
    "snoRNA": "#edc948",
    "piRNA":  "#b07aa1",
    "other":  "#bab0ac",
    "clean":  "#59a14f",
}

LABELS = {
    "rRNA":   "rRNA",
    "tRNA":   "tRNA",
    "snRNA":  "snRNA",
    "snoRNA": "snoRNA",
    "piRNA":  "piRNA",
    "other":  "Other RNA",
    "clean":  "Clean reads",
}


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
def parse_contamination_stats(fpath):
    """Return dict of {rna_type: count} from a *_contamination_stats.txt file."""
    counts = {}
    rna_keys = {"rRNA", "tRNA", "snRNA", "snoRNA", "piRNA", "other"}
    with open(fpath, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("Output Reads (Clean):"):
                val = line.split(":", 1)[1].strip().replace(",", "")
                try:
                    counts["clean"] = int(val)
                except ValueError:
                    pass
            else:
                m = re.match(r"^(\w+):\s+([\d,]+)\s*\(", line)
                if m and m.group(1) in rna_keys:
                    try:
                        counts[m.group(1)] = int(m.group(2).replace(",", ""))
                    except ValueError:
                        pass
    return counts


def load_manifest_order(manifest_path):
    """Return ordered list of sample names from the manifest (skip header)."""
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    outdir = os.environ.get("OUTDIR")
    if not outdir:
        print("ERROR: OUTDIR environment variable is not set.", file=sys.stderr)
        sys.exit(1)

    manifest = os.environ.get("MANIFEST", "")
    mqc_dir = os.path.join(outdir, "multiqc")
    os.makedirs(mqc_dir, exist_ok=True)

    # Discover stats files
    stats_files = sorted(
        glob.glob(os.path.join(outdir, "contamination_filtered", "*_contamination_stats.txt"))
    )
    if not stats_files:
        print("WARNING: No contamination stats files found – skipping plot.", file=sys.stderr)
        return

    # Parse each file
    data = {}
    for fpath in stats_files:
        bn = os.path.basename(fpath)
        sample = re.sub(r"_contamination_stats\.txt$", "", bn)
        parsed = parse_contamination_stats(fpath)
        if parsed:
            data[sample] = parsed

    if not data:
        print("WARNING: Could not parse any contamination stats – skipping plot.", file=sys.stderr)
        return

    # Order samples
    manifest_order = load_manifest_order(manifest)
    if manifest_order:
        samples = [s for s in manifest_order if s in data]
        samples += [s for s in sorted(data) if s not in samples]
    else:
        samples = sorted(data.keys())

    # Build matrix (samples × rna_types)
    matrix = np.array(
        [[data[s].get(rt, 0) for rt in RNA_TYPES] for s in samples],
        dtype=float,
    )

    # ---------------------------------------------------------------------------
    # Plot
    # ---------------------------------------------------------------------------
    fig_width = max(10, len(samples) * 0.55)
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    xs = np.arange(len(samples))
    bar_width = 0.7
    bottom = np.zeros(len(samples))

    for i, rt in enumerate(RNA_TYPES):
        vals = matrix[:, i]
        ax.bar(xs, vals, bar_width, bottom=bottom, color=COLORS[rt], label=LABELS[rt])
        bottom += vals

    ax.set_xticks(xs)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Read count", fontsize=11)
    ax.set_title("Contamination Filtering: Read Fate per Sample (Phase I.5b)", fontsize=12)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.85, ncol=2)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
    plt.tight_layout()

    out_path = os.path.join(mqc_dir, "contamination_breakdown.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
