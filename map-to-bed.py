import sys

assoc_sorted = sys.argv[1]

out = open(sys.argv[1] + ".bed", "w")

for i, toks in enumerate(x.strip().split("\t") for x in open(assoc_sorted)):
	chrom = 0
	start = 3

	#print(toks[chrom])
	#print(toks[start])

	end = int(toks[start]) + 1
	#print(end)

	out.write("\t".join([toks[chrom], toks[start], str(end), "\n"]))
