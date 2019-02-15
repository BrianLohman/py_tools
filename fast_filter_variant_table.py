# read in variant tables and filter down to genes in given list
# write the number of variants that each individual has in these genes to file

# USAGE: python fast_filter_variant_table.py -v [variant table] -g [genes of interest] -o [output name] 

import dask.dataframe as ddf
from dask.distributed import Client
import sys
import argparse
import pandas

# set up resources for dask: 32 workers, each 1 thread
# this works on kingspeak31 which has many more threads than cores, so we need to account for that
client = Client(processes = False, n_workers = 32, threads_per_worker = 1, memory_limit = '4.8GB')

# define arguments
parser = argparse.ArgumentParser(description='sum the number of variants per gene in an individual')
parser.add_argument('-v', '--variants', dest = 'variants', help = 'variant table')
parser.add_argument('-g', '--genes', dest = 'genes', help = 'genes of interest')
parser.add_argument('-o', '--out', dest = 'out', help = 'output file name')
args = parser.parse_args()

# define list of genes of interest
with open(args.genes) as g:
	genes_of_interest = g.read().splitlines()

# make list of EUR proband IDs
master = pandas.read_table("/scratch/ucgd/lustre/work/u0806040/data/15_Jan_19_Simons_master_ancestry_corrected_PRS.txt", dtype={'other_dx_axis_i': 'object', 'other_dx_axis_ii': 'object', 'other_dx_icd': 'object'})
probands = master.loc[master['family_member'] == 'p1']
eur_probands = probands.loc[probands['ancestry.prediction'] == 'EUR']
proband_ids = eur_probands['IID']

# read in table of variants
print 'reading in variants'
variants = ddf.read_table(args.variants)
#print variants.head()

# filter variants
print 'setting up variant filters'
# medium and and high impact
variants[variants.impact.isin(['MED', 'HIGH'])]

# in gens of interest
variants[variants.gene.isin(genes_of_interest)]

# reorganize data frame so that rows are genes of interest, columns are IIDs and value are coutns of variants
# drop metadata columns that are not needed in output
variants.drop(variants.columns[0:4], axis = 'columns')
variants.drop(variants.columns[1:4], axis = 'columns')

# convert back to pandas now that the data frame is small
print 'computing and returing pandas data frame'
voi = variants.compute()

## ERRORS GALORE
#distributed.worker - WARNING - gc.collect() took 1.651s. This is usually a sign that the some tasks handle too many Python objects at the same time. Rechunking the work into smaller tasks might help.
#distributed.worker - WARNING - Worker is at 103% memory usage. Pausing worker.  Process memory: 5.05 GB -- Worker memory limit: 4.80 GB

# MEMORY USAGE CLIMBS OVER TIME

# replace -1 (missing) with 0 in preparation for summing columns
voi.replace(to_replace = -1, value = 0, inplace = True)

# make list of samples
samples = list(voi)[1:]

# empty table for gene and sums by individusl
table = {}

print 'summing variants per gene for each sample'

# loop through genes, summing the genotypes for each indivudal
for g in genes_of_interest:
	tmp = voi.loc[voi.gene == g ]
	tmp = tmp.drop(['gene'], axis = 'columns')
	table[g] = tmp.sum()

print 'witing to file'

# write to file
with open(args.out, "w") as out:
	out.write('\t'.join(['gene', '\t'.join([str(i) for i in samples])])+'\n')
	for g in genes_of_interest:
		out.write('\t'.join([g, '\t'.join([str(i) for i in list(table[g])])+'\n']))
