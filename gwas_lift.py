from __future__ import print_function
import sys
import glob

# USAGE: python gwas_lift.py [target (hg38) dictionary]

# make empty dictionary
rs_to_pos = {}

# populate dictionary
for toks in (x.strip().split("\t") for x in open(sys.argv[1])):
    chrom = toks[0]
    hg38_pos = toks[2]
    rsid = toks[3]
    rs_to_pos[rsid] = [chrom, hg38_pos]

# read in gwas files
gwas = glob.glob("*.assoc")

# loop through gwas files, looking up rsIDs in dictionary
# print unmapped rsIDs to .unMapped file
# print results to .lifted.hg38 file

for g in gwas:
    out = open(g[:-6] + ".hg38.assoc", "w")
    log = open(g + ".log", "w")
    out.write("\t".join(["SNP", "CHR", "BP", "A1", "A2", "BETA", "P\n"]))
    for i, toks in enumerate(x.strip().split("\t") for x in open(g)):
        if i == 0 and "SNP" in toks:
            continue
        A1 = toks[3]
        A2 = toks[4]
        beta = toks[5]
        p = toks[6]
        rsid = toks[0]
        try:
            snp = rs_to_pos[rsid]
        except KeyError:
            print("%s not found" % rsid, file=log)
            continue
        chrom = snp[0]
        hg38_pos = snp[1]
        out.write("\t".join([rsid, chrom, hg38_pos, A1, A2, beta, p + "\n"]))
    out.close()
    log.close()
