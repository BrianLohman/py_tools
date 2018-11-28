# intersect the list of genes in a pathway and the list of SFARI genes, given a gene score filter

import sys

# USAGE: python intersect_pathway_SFARI.py [pathway bed file] 

# read in the SFARI gene list as a dictionary
sfari = {}

for toks in (x.strip().split('\t') for x in open('SFARI-Gene_genes_export18-09-2018.txt')):
	sfari[toks[1]] = toks[5]

# dictioinary look up of genes in the pathway file
for toks in (x.strip().split('\t') for x in open(sys.argv[1])):
	gene = toks[3]

	try:
		score = sfari[gene]
	except KeyError:
		continue

	print gene,score
