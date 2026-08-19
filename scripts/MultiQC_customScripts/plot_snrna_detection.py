#!/usr/bin/env python3
"""
plot_snrna_detection.py – Standalone export: distinct snRNA detected per sample

Reads snrna_counts/*_snrna_stats.txt produced by 6.2_snrna_snorna_quantification.sh
and exports a grouped bar chart showing distinct snRNA sequences detected per
sample at ≥1, ≥10, and ≥100 read thresholds.

Environment variables
---------------------
OUTDIR   (required) – pipeline results directory
MANIFEST (optional) – used to sort samples in manifest order
"""

import os
import sys

# Import the shared plotting helper from the co-located module
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
        subdirectory="snrna_counts",
        filename_suffix="_snrna_stats.txt",
        sample_label="snRNA",
        rna_label="snRNAs",
        out_filename="snrna_distinct_detected.png",
    )


if __name__ == "__main__":
    main()
