#!/usr/bin/env python

# read in a data frame of genes and counts of variants in that gene
# read in a data frame samples and counts of variants in said genes
# return a list of samples with the highest Spearman Rho (sample variant counts vs specified list)

# USAGE:
# python replicate_by_variant_count.py [table: rows = genes, columns = count of variants, gene neame] [table: rows = genes, columns = individuals]

import sys
import pandas
import scipy.stats

# read in the pre-specified list of genes and coutns of variatns therein
gene_importance = pandas.read_table(sys.argv[1], header = None)
gene_importance.columns = ['count', 'gene']

# read in table of samples and their variant counts
samples = pandas.read_table(sys.argv[2])

# reorder the rows to match the standard
samples = samples.set_index('gene')
samples = samples.reindex(index = gene_importance['gene'])
samples = samples.reset_index()

# check gene order
if list(samples['gene']) != list(gene_importance['gene']):
	print("sample order error!")
	sys.exit(1)

# calculate spearman rho between each sample and the standard
iid = list(samples)[1:]
result = {}
for i in iid:
	tmp  = scipy.stats.spearmanr(a = gene_importance['count'], b = samples[i])
	result[i] = tmp[0]

# convert dict to data frame
x = pandas.DataFrame.from_dict(result, orient = 'index')
x.columns = ['cor']
#print(x.head())

# sort data frame base on rho
y = x.sort_values(['cor'], ascending = False)
#print(y.head())

# print to file
y.to_csv(sys.argv[2][:-4]+"_sample_rho.txt", sep = '\t', header = False)
