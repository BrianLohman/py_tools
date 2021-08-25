#!/usr/bin/env python

# Filter de novos from readback phasing

# Usage: filter_readback_phasing.py [input bed file] [output filename]

import sys

out = open(sys.argv[2], 'w')

with open(sys.argv[1], 'r') as bed:
    for chrom, start, end, ref, alt, sample_id, candidate_sites, dad_evidence, mom_evidence in (x.strip().split("\t") for x in bed):
        if chrom == 'chrom':
            print(chrom, start, end, ref, alt, sample_id, candidate_sites, dad_evidence, mom_evidence, 'origin', sep = "\t", file = out)
            continue

        if candidate_sites == 'NA':
            continue

        if int(mom_evidence) >= 2 and int(dad_evidence) >= 2:
            continue

        if int(dad_evidence) == 0 and int(mom_evidence) == 0:
            continue

        if int(dad_evidence) == 1 and int(mom_evidence) == 1:
            continue

        if int(dad_evidence) + int(mom_evidence) < 2:
            continue

        if int(mom_evidence) > int(dad_evidence):
            print(chrom, start, end, ref, alt, sample_id, candidate_sites, dad_evidence, mom_evidence, 'mom', sep = "\t", file = out)
            continue

        if int(dad_evidence) > int(mom_evidence):
            print(chrom, start, end, ref, alt, sample_id, candidate_sites, dad_evidence, mom_evidence, 'dad', sep = "\t", file = out)
            continue

        else:
            print("fail: variant not account for by logic " + chrom, start, end, file = sys.stderr)
            continue
