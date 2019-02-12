# read in variant tables and filter down to genes in given list
# write the number of variants that each individual has in these genes to file

import dask
import sys
import argparse

parser = argparse.ArgumentParser(description='sum the number of variants per gene in an individual')

parser.add_argument('-v', '--variants', dest = 'variants', help = 'variant table')
parser.add_argument('-g', '--genes', dest = 'genes', help = 'genes of interest')
parser.add_argument('-o', '--out', dest = 'out', help = 'output file name')
args = parser.parse_args()

# read in table of variants
variants = pandas.read_table(args.variants, low_memory = False)

# define out file
out = open(args.out, "w")

# define list of genes of interest
with open(args.genes) as g:
	genes = g.read().splitlines()

master = pandas.read_table("/scratch/ucgd/lustre/work/u0806040/data/15_Jan_19_Simons_master_ancestry_corrected_PRS.txt", low_memory = False)
probands = master.loc[master['family_member'] == 'p1']
probands = probands.loc[probands['ancestry.prediction'] == 'EUR']
probands = probands['IID']


