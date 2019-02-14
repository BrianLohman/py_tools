# read in variant tables and filter down to genes in given list
# write the number of variants that each individual has in these genes to file

# USAGE: python fast_filter_variant_table.py -v [variant table] -g [genes of interest] -o [output name] 

# note that you should consider upgrading to python 3.x from 2.x 
# stopping development of dask for python 2.x is an ongoing discussion
# https://github.com/dask/dask/issues/4047

import dask.dataframe as ddf
from dask.distributed import Client
import sys
import argparse

# define arguments
parser = argparse.ArgumentParser(description='sum the number of variants per gene in an individual')
parser.add_argument('-v', '--variants', dest = 'variants', help = 'variant table')
parser.add_argument('-g', '--genes', dest = 'genes', help = 'genes of interest')
parser.add_argument('-o', '--out', dest = 'out', help = 'output file name')
args = parser.parse_args()

# if this isn't working for you, see https://distributed.dask.org/en/latest/setup.html
# or talk to your friendly IT professional
client = Client()
client.restart()

# read in table of variants
variants = ddf.read_table(args.variants, blocksize = 50e6) # 50 MB blocks
#print variants.head()
print 'variants read in'

# define list of genes of interest
with open(args.genes) as g:
	genes_of_interest = g.read().splitlines()

# make list of EUR proband IDs
# I know nothing about this size, but consider passing in an index, which will make operations later on faster
# if it makes sense for what your trying to do
# e.g. I imagine indexing by gene might be useful
# see http://docs.dask.org/en/latest/dataframe-performance.html#use-the-index
master = ddf.read_table("/scratch/ucgd/lustre/work/u0806040/data/15_Jan_19_Simons_master_ancestry_corrected_PRS.txt",
            blocksize = 25e6, # I would try leaving this off...dask can figure out the optimal block sizes on its own
            dtype={'other_dx_axis_i': 'object', 'other_dx_axis_ii': 'object', 'other_dx_icd': 'object'}
            ) # 25 MB blocks

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

# variant filters
# medium and and high impact
med_high = variants[variants['impact'].isin(['MED', 'HIGH'])]

# variants to those in gens of interest
voi = med_high[med_high ['gene'].isin(genes_of_interest)]

# reorganize data frame so that rows are genes of interest, columns are IIDs and value are coutns of variants
# drop metadata columns that are not needed in output
voi = voi.drop(voi.columns[0:4], axis = 'columns')
voi = voi.drop(voi.columns[1:4], axis = 'columns')

print 'filtered variants'

# convert back to pandas now that the data frame is small
voi = voi.compute()
print 'converted back to pandas'
print voi.shape 

# replace -1 (missing) with 0 in preparation for summing columns
voi = voi.replace(to_replace = -1, value = 0)

# make list of samples
samples = list(voi)[1:]

# empty table for gene and sums by individusl
table = {}

print 'summing variants per gene for each sample'

# loop through genes, summing the genotypes for each indivudal
for g in genes_of_interest:
	tmp = voi.loc[voi['gene'] == g ]
	tmp = tmp.drop(['gene'], axis = 'columns')
	table[g] = tmp.sum()

print 'witing to file'

# write to file
with open(args.out, "w") as out:
    out = open(args.out, "w")
    out.write('\t'.join(['gene', '\t'.join([str(i) for i in samples])])+'\n')
    for g in genes_of_interest:
        out.write('\t'.join([g, '\t'.join([str(i) for i in list(table[g])])+'\n']))