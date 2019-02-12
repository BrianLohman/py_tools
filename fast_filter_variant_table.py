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
variants = ddf.read_table(args.variants, blocksize = 10e8) # 1 GB blocks
#print variants.head()

# define out file
out = open(args.out, "w")

# define list of genes of interest
with open(args.genes) as g:
	genes = g.read().splitlines()

#print genes

# make list of EUR proband IDs
master = ddf.read_table("/scratch/ucgd/lustre/work/u0806040/data/15_Jan_19_Simons_master_ancestry_corrected_PRS.txt", blocksize = 25e6, dtype={'other_dx_axis_i': 'object', 'other_dx_axis_ii': 'object', 'other_dx_icd': 'object'}) # 25 MB blocks

# text to numeric for speed
#master['family_member'] = master['family_member'].astype('category').compute()
#master['ancestry.prediction'] = master['ancestry.prediction'].astype('category').compute()

# filters
probands = master.loc[master['family_member'] == 'p1'] # probands
probands.compute()

eur_probands = probands.loc[probands['ancestry.prediction'] == 'EUR'] # of EUR ancestry
eur_probands.compute()

proband_ids = probands['IID'] # collect IDs
proband_ids.compute()

print proband_ids
