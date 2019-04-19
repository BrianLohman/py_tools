#!/usr/bin/env python

# with bed file, check to see if variants are in possible CpG context

# USAGE: cpg_lookup.py [bed file] [reference genome] 

from pyfaidx import Fasta
import sys

# define I/O
bed_file = sys.argv[1]
ref_file = sys.argv[2]
cpg_bed = open('cpg_only_' + bed_file, 'w')
fa = Fasta(ref_file, read_ahead = 100000, as_raw = True, sequence_always_upper = True)

# loop over variants in bed file
k = 0
for chrom, start, end, ref, alt, sample,  in (line.strip().split('\t') for line in open(bed_file)):
    k += 1
    if k % 10000 == 0:
        print(str(k) + ' variants screened', file = sys.stderr)
    
    if ref == 'C':
        three_prime_base = fa[chrom][int(end)]
        if three_prime_base == 'G':
            print(chrom, start, end, ref, alt, sample, sep = '\t', file = cpg_bed)
        else:
            continue

    if ref == 'G':
        three_prime_base = fa[chrom][int(start) - 1]
        if three_prime_base == 'C':
            print(chrom, start, end, ref, alt, sample, sep = '\t', file = cpg_bed)
        else:
            continue
    
    else:
        continue

