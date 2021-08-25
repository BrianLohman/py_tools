from __future__ import print_function
import sys

pos_to_rs = {}

log = open(sys.argv[1][:-4] + "_add_rsID.log", "w")

for i, line in enumerate(
    open("/scratch/ucgd/lustre/work/u0806040/data/hg37_all_snps150_chr_pos_rsID.txt")
):
    fields = line.strip().split("\t")
    chrom = fields[0]
    pos = fields[2]
    rsid = fields[3]
    pos_to_rs[(chrom, pos)] = rsid

for i, line in enumerate(open(sys.argv[1])):
    fields = line.strip().split("\t")
    chrom = fields[0]
    pos = fields[3]
    ref = fields[4]
    alt = fields[5]
    try:
        rsid = pos_to_rs[(chrom, pos)]
    except KeyError:
        print("%s not found" % (",".join([chrom, pos])), file=log)
        continue
    print("\t".join([chrom, rsid, "0", pos, ref, alt]))
