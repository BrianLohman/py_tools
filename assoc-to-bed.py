import sys

assoc_sorted = sys.argv[1]

out = open(sys.argv[1] + ".bed", "w")

for i, toks in enumerate(x.strip().split("\t") for x in open(assoc_sorted)):
	if i == 0 and "CHR" in toks: continue
	chrom = 5
	end = 6

	start = int(toks[end]) - 1

	chrom_full = "".join(['chr', str(toks[chrom])])

	out.write("\t".join([chrom_full, str(start), toks[end], "\n"]))
