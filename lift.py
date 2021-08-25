from __future__ import print_function
import sys
from pyfaidx import Fasta
import gzip
import toolshed as ts

# USAGE: python lift.py [hg38 rsIDs and positions] [hg19/37 .bim with rsIDs] [hg37 ref] [hg38 ref]


def read_lookup(fname):
    d = {}
    k = 0
    for l in open(fname):
        chrom, start, end, rs = l.rstrip().split("\t")
        if chrom.startswith("chr"):
            chrom = chrom[3:]
        d[rs] = (chrom, end)
        k += 1
        if k % 10000000 == 0:
            print(k, file=sys.stderr)
    return d


L = read_lookup(sys.argv[1])
fa = Fasta(sys.argv[3], read_ahead=100000, as_raw=True, sequence_always_upper=True)
fa38 = Fasta(sys.argv[4], read_ahead=100000, as_raw=True, sequence_always_upper=True)

print("read %d locations into lookup keyed by rsid" % len(L), file=sys.stderr)


def fix(fa, chrom, pos, ref, alt, prefix="", rs=None):
    pos = int(pos)
    faref = fa[prefix + chrom][pos - 1]
    if faref != ref:
        # assert faref == alt, (chrom, pos, ref, alt, faref)
        if faref != alt:
            print(rs, chrom, pos, ref, alt, "setting ref to ->", faref, file=log)
            return faref, alt

        return alt, ref
    return ref, alt


log = open(sys.argv[2] + ".log", "w")

for (chrom, rs, value, pos, ref, alt) in (
    l.rstrip().split("\t") for l in open(sys.argv[2])
):

    # fix ref/alt order on original reference
    ref, alt = fix(fa, chrom, pos, ref, alt)

    try:
        oc, oe = L[rs]
    except KeyError:
        print("%s not found" % rs, file=log)
        continue

    chrom = oc
    pos = oe

    # update ref and alt to hg38. in some cases these are flipped
    ref, alt = fix(fa38, chrom, pos, ref, alt, prefix="chr", rs=rs)

    print("\t".join(("chr" + chrom, rs, value, pos, ref, alt)))
