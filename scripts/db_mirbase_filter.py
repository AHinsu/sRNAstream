#!/usr/bin/env python3
"""Filter a miRBase FASTA file to keep only sequences for a given species code.

Usage:
    python3 db_mirbase_filter.py <input.fa> <output.fa> <species_code>
"""
import sys

in_fa, out_fa, code = sys.argv[1], sys.argv[2], sys.argv[3].lower()
keep = False
with open(in_fa) as fh_in, open(out_fa, 'w') as fh_out:
    for line in fh_in:
        if line.startswith('>'):
            keep = line.lower().startswith(f'>{code}-') or f' {code}-' in line.lower()
        if keep:
            fh_out.write(line)
