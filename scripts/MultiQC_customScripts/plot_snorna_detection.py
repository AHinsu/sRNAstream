#!/usr/bin/env python3
"""
plot_snorna_detection.py – Standalone export: distinct snoRNA detected per sample

Reads snorna_counts/*_snorna_stats.txt produced by 6.2_snrna_snorna_quantification.sh
and exports a grouped bar chart showing distinct snoRNA sequences detected per
sample at ≥1, ≥10, and ≥100 read thresholds.

Environment variables
---------------------
OUTDIR   (required) – pipeline results directory
MANIFEST (optional) – used to sort samples in manifest order
"""

import os
import sys

# Insert the MultiQC_customScripts directory on sys.path so that
# plot_mirna_detection can be imported when this script is run standalone
# (i.e. not as part of an installed package).
sys.path.insert(0, os.path.dirname(__file__))
from plot_mirna_detection import plot_distinct_detected


def main():
    outdir = os.environ.get("OUTDIR")
    if not outdir:
        print("ERROR: OUTDIR environment variable is not set.", file=sys.stderr)
        sys.exit(1)

    plot_distinct_detected(
        outdir=outdir,
        subdirectory="snorna_counts",
        filename_suffix="_snorna_stats.txt",
        sample_label="snoRNA",
        rna_label="snoRNAs",
        out_filename="snorna_distinct_detected.png",
    )


if __name__ == "__main__":
    main()
