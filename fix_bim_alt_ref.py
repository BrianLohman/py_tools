import sys
from pyfaidx import Fasta

# USAGE
# python fix_bim_ref_alt.py [input bim file] [reference fasta] [output prefix]

bim = sys.argv[1]
fa = Fasta(sys.argv[2], as_raw=True)
prefix = sys.argv[3]

# params for .bim
chrom_index = 0
pos_index = 3
a_index = 5
b_index = 4

assert a_index > b_index, "specify a > b"

bad = open(prefix + ".bad.bim", "w")
nbad = 0

out = open(prefix + ".fixed.bim", "w")

seq = ""
last_chrom = None
for i, toks in enumerate(x.rstrip().split("\t") for x in open(bim)):
    if i == 0 and "CHR" in toks:
        continue

    if "chr" + toks[chrom_index] != last_chrom:
        last_chrom = "chr" + toks[chrom_index]
        seq = fa[last_chrom]

    pos = int(toks[pos_index]) - 1
    ref = seq[pos]  # ref sequence for the given chrom, pos

    # error checking
    if not toks[a_index] in "ACTG":
        nbad += 1
        bad.write("\t".join(toks) + "\t" + str(i) + "\n")
        continue
    if not toks[b_index] in "ACTG":
        nbad += 1
        bad.write("\t".join(toks) + "\t" + str(i) + "\n")
        continue

    if not ref in (toks[a_index], toks[b_index]):
        nbad += 1
        bad.write("\t".join(toks) + "\t" + str(i) + "\n")
        continue

    # if the a_index doesn't hold the reference allele, then we flip a and b.
    if toks[a_index] != ref:
        toks[a_index], toks[b_index] = toks[b_index], toks[a_index]
    out.write("\t".join(toks) + "\n")

print("n_bad: %d out of %d" % (nbad, i))
