#!/usr/bin/env python3
"""Convert RNA sequences to DNA (U → T) in a FASTA file.

Usage:
    python3 db_rna_to_dna.py <input.fa> <output.fa>
"""
import sys

in_fa, out_fa = sys.argv[1], sys.argv[2]
with open(in_fa) as fh_in, open(out_fa, 'w') as fh_out:
    for line in fh_in:
        if line.startswith('>'):
            fh_out.write(line)
        else:
            fh_out.write(line.upper().replace('U', 'T'))
