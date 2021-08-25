#!/usr/bin/env python

# convert simons format list of de novos to bed

import sys

# dict to convert simons id to IID scheme
d = {}
with open(sys.argv[2], "r") as id_dict:
    for IID, simons_id in (line.strip().split("\t") for line in id_dict):
        d[simons_id] = IID

problem_samples = []

out = open("all_simons_deNovos.bed", "w")

with open(sys.argv[1], "r") as simons:
    for ID, SampleID, Batch, Allele, Consequence, Impact, Symbol in (
        line.strip("\r\n").split("\t") for line in simons
    ):
        if ID == "ID":
            continue
        if Batch != "P231":
            continue
        chrom, end, ref, alt = ID.split(":")
        start = int(end) - len(ref)

        if SampleID in d:
            print(chrom, start, end, ref, alt, SampleID, sep="\t", file=out)
        else:
            problem_samples.append(SampleID)
            continue
            # print(chrom, start, end, ref, alt, SampleID,  sep = '\t', file = out)

print(set(problem_samples))
