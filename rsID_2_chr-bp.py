import sys

assoc = sys.argv[1]
out = open(sys.argv[1] + ".renamed", "w")

for i, toks in enumerate(x.strip().split("\t") for x in open(assoc)):
    if i == 0 and "SNP" in toks:
        out.write("\t".join(["SNP", "A1", "A2", "BETA", "P", "CHR", "BP\n"]))
        continue
    chrom = 5
    bp = 6
    P = 4
    beta = 3
    A2 = 2
    A1 = 1

    snp_ID = ":".join([toks[chrom], toks[bp]])

    out.write(
        "\t".join(
            [
                snp_ID,
                toks[A1],
                toks[A2],
                toks[beta],
                toks[P],
                toks[chrom],
                toks[bp],
                "\n",
            ]
        )
    )
