# read in variant tables and filter down to genes in given list
# write the number of variants that each individual has in these genes to file

# USAGE: python fast_filter_variant_table.py -v [variant table] -g [genes of interest] -o [output name] 

import dask.dataframe as ddf
import sys
import argparse

# define arguments
parser = argparse.ArgumentParser(description='sum the number of variants per gene in an individual')
parser.add_argument('-v', '--variants', dest = 'variants', help = 'variant table')
parser.add_argument('-g', '--genes', dest = 'genes', help = 'genes of interest')
parser.add_argument('-o', '--out', dest = 'out', help = 'output file name')
args = parser.parse_args()

# read in table of variants
variants = ddf.read_table(args.variants, blocksize = 500e6) # 1 GB blocks
#print variants.head()

# define out file
out = open(args.out, "w")

# define list of genes of interest
with open(args.genes) as g:
	genes_of_interest = g.read().splitlines()

#print genes

# make list of EUR proband IDs
master = ddf.read_table("/scratch/ucgd/lustre/work/u0806040/data/15_Jan_19_Simons_master_ancestry_corrected_PRS.txt", blocksize = 25e6, dtype={'other_dx_axis_i': 'object', 'other_dx_axis_ii': 'object', 'other_dx_icd': 'object'}) # 25 MB blocks

# text to numeric for speed
#master['family_member'] = master['family_member'].astype('category')
#master['ancestry.prediction'] = master['ancestry.prediction'].astype('category')

# sample filters
# probands
probands = master.loc[master['family_member'] == 'p1']

# of EUR ancestry
eur_probands = probands.loc[probands['ancestry.prediction'] == 'EUR']

# collect IIDs
proband_ids = eur_probands['IID']

print proband_ids.head()

# variant filters
# medium and and high impact
med_high = variants[variants['impact'].isin(['MED', 'HIGH'])]

# variants to those in gens of interest
voi = med_high[med_high ['gene'].isin(genes_of_interest)]

# reorganize data frame so that rows are genes of interest, columns are IIDs and value are coutns of variants
# drop metadata columns that are not needed in output
voi = voi.drop(voi.columns[0:4], axis = 'columns')
voi = voi.drop(voi.columns[1:4], axis = 'columns')

# convert back to pandas now that the data frame is small
voi = voi.compute()

# replace -1 (missing) with 0 in preparation for summing columns
voi = voi.replace(to_replace = -1, value = 0)

# make list of samples
samples = list(voi)[1:]

# empty table for gene and sums by individusl
table = {}

# loop through genes, summing the genotypes for each indivudal
for g in genes_of_interest:
	tmp = voi.loc[voi['gene'] == g ]
	tmp = tmp.drop(['gene'], axis = 'columns')
	table[g] = tmp.sum()

# write to file
out.write('\t'.join(['gene', '\t'.join([str(i) for i in samples])])+'\n')
for g in genes_of_interest:
	out.write('\t'.join([g, '\t'.join([str(i) for i in list(table[g])])+'\n']))
