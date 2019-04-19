#!/usr/bin/env python

# with bed file, check to see if variants are in possible CpG context

# USAGE: cpg_lookup.py [variant file] [reference genome] 

from pyfaidx import Fasta
import sys

# define I/O
bed_file = sys.argv[1]
ref_file = sys.argv[2]
cpg_bed = open('cpg_only_' + bed_file, 'w')
fa = Fasta(ref_file, read_ahead = 100000, as_raw = True, sequence_always_upper = True)

# loop over variants in variant file
k = 0
for snp, sample, batch, allele, consequence, impact, symbol  in (line.strip().split('\t') for line in open(bed_file)):
    k += 1
    if k % 100 == 0:
        print(str(k) + ' variants screened', file = sys.stderr)

    if snp == 'ID':
        continue
    
    chrom, end, ref, alt = snp.strip().split(':')
    start = int(end) - 1

    if ref == 'C':
        assert(ref == fa[chrom][start])
        three_prime_base = fa[chrom][int(end)]
        if three_prime_base == 'G':
            print(snp, sample, batch, allele, consequence, impact, symbol,  sep = '\t', file = cpg_bed)
        else:
            continue

    if ref == 'G':
        assert(ref == fa[chrom][start])
        three_prime_base = fa[chrom][int(start) - 1]
        if three_prime_base == 'C':
            print(snp, sample, batch, allele, consequence, impact, symbol, sep = '\t', file = cpg_bed)
        else:
            continue
    
    else:
        continue

