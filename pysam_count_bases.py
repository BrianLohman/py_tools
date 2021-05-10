#!/usr/bin/env python

import pysam
import sys

# name output based on input and print heades
outfile = open(sys.argv[1][:-4]+"_bases.txt", 'w')
print('\t'.join(["A", "T", "G", "C", "N", "Length"]), file = outfile)

# read in BAM
bamfile = pysam.AlignmentFile(sys.argv[1], "rb")

# iterate over BAM
for r in bamfile.fetch():
    if r.is_secondary or r.is_supplementary:
        continue
    else:
        read = str(r.get_forward_sequence())
        print('\t'.join(str(x) for x in [read.count('A'), read.count('T'), read.count('G'), read.count('C'), read.count('N'), len(read)]), file = outfile)
