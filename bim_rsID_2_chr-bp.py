import sys

assoc = sys.argv[1]
out = open(sys.argv[1] + ".renamed", "w")

for toks in (x.strip().split("\t") for x in open(assoc)):
    chrom = 0
    bp = 3
    A2 = 5
    A1 = 4
    snp_ID = ":".join([toks[chrom],toks[bp]])

    out.write("\t".join([toks[chrom],snp_ID,'0',toks[bp],toks[A1], toks[A2],"\n"]))
