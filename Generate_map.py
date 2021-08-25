pos_to_rs = {}
for line in open("rs2chr_bp.tab"):
    fields = line.strip().split("\t")
    chrom = fields[0]
    pos = fields[2]
    rsid = fields[1]
    pos_to_rs[(chrom, pos)] = rsid

for line in open("simons.1kgsites.bim"):
    fields = line.strip().split("\t")
    chrom = fields[0]
    pos = fields[3]
    ref = fields[4]
    alt = fields[5]
    rsid = pos_to_rs[(chrom, pos)]

    print("\t".join([chrom, rsid, "0", pos, ref, alt]))
