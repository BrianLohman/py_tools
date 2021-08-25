from __future__ import print_function
import sys
import glob

# convert the rsID to chr:BP in summary stats (.assoc files)

# USAGE: python gwas_rsID_2_chrBP.py

# read in gwas files
gwas = glob.glob("*.rsID.hg38.assoc")

# loop through gwas files swapping rsID for chr:BP
# print results to .chrBP.hg38.assoc

# conversion
for g in gwas:
    out = open(g[:-16] + ".chrBP.hg38.assoc", "w")
    out.write("\t".join(["SNP", "CHR", "BP", "A1", "A2", "BETA", "P\n"]))
    for i, fields in enumerate(x.strip().split("\t") for x in open(g)):
        if i == 0 and "SNP" in fields:
            continue
        chrom = fields[1]
        pos = fields[2]
        ref = fields[3]
        alt = fields[4]
        beta = fields[5]
        p = fields[6]
        snp_ID = ":".join([str(chrom), str(pos)])

        out.write("\t".join([snp_ID, chrom, pos, ref, alt, beta, p + "\n"]))
    out.close()
