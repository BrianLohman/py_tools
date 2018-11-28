# reduce the number of entries in the gene models file to 1 per gene
# For each gene, list the smallest and largest position in bp

# USAGE: python reduce_gene_models.py [gene_models] 

# ASSUMES: 1) the list of gene models is sorted based on gene name, 2) the same gene does not occur on more than one chromosome

import sys

gene_models = sys.argv[1]

out = open("thinned_" + gene_models, "w")

current_gene = None
current_chrom = None
start, end = [], []

for i, toks in enumerate(x.strip().split('\t') for x in open(gene_models)):
    if i == 0:
        current_gene = toks[4]
	current_chrom = toks[0]
    if toks[4] == current_gene:
        start.append(int(toks[2]))
        end.append(int(toks[3]))
    else:
        out.write('\t'.join([current_chrom, toks[1], str(min(start)), str(max(end)), current_gene+'\n'])) 
        start, end = [], []
        current_gene = toks[4]
	current_chrom = toks[0]
        start.append(int(toks[2]))
        end.append(int(toks[3]))
