# build regions file for extracting gene models in pathways of interest from vcf file

import sys

# USAGE: python build_regions.py [thinned gene sorted gene models]

# assign input file:
gene_models = sys.argv[1]

# build dictionary of gene name and associated pathway
gene_to_pathway = {}

for toks in (x.strip().split("\t") for x in open("genes_by_pathway.txt")):
    gene_to_pathway[toks[0]] = toks[1]

# make space for the output to go
out_fh = "regions_with_name" + gene_models[-9:]
print(out_fh)
out = open(out_fh, "w")

# read in the lines from the gene models
for toks in (x.strip().split("\t") for x in open(gene_models)):
    chrom = toks[0]
    start = toks[2]
    end = toks[3]
    gene = toks[4]

    # check that the gene is in the above dictionary
    try:
        pathway = gene_to_pathway[gene]
    except KeyError:
        continue

    # check that the smaller of the start/end comes first
    if int(start) > int(end):
        start = toks[3]
        end = toks[2]

    # print the chrom, start, end, name, and pathway
    out.write("\t".join([chrom, start, end, gene, pathway + "\n"]))
