import sys

# convert from bim to bed format
# no chr prefix

# USAGE: python bim-to-bed.py [.bim]

bim = sys.argv[1]

out = open(sys.argv[1][:-3] + "bed", "w")

for toks in (x.strip().split("\t") for x in open(bim)):
    chrom = 0
    start = 3
    end = int(toks[start]) + 1
    out.write("\t".join([str(chrom), str(toks[start]), str(end), "\n"]))
