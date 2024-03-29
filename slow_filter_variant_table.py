#!/usr/bin/env python
# read in a variant table, filter down to med and high impact variants in certain genes
import sys

import numpy
import pandas

# read in table of variants
variants = pandas.read_table(sys.argv[1], low_memory=False)

# define out file
out = open(sys.argv[1][:-4] + "_variants_in_genes_by_individual.txt", "w")

# define list of genes of interest
genes = [
    "PLXNA3",
    "CACNA1A",
    "PLXNB1",
    "EP300",
    "CREBBP",
    "PLXNA4",
    "ROBO1",
    "OPRM1",
    "GNAS",
    "SLIT3",
    "ITPR1",
    "EPHB6",
    "NOS1",
    "NTRK1",
    "CACNA1C",
    "ROBO2",
    "GRID1",
    "GRM8",
    "PLCB1",
    "SRGAP3",
    "SEMA5A",
    "PTGER3",
]

# get list of EUR proband IDs
master = pandas.read_table(
    "/scratch/ucgd/lustre/work/u0806040/data/15_Jan_19_Simons_master_ancestry_corrected_PRS.txt",
    low_memory=False,
)
probands = master.loc[master["family_member"] == "p1"]
probands = probands.loc[probands["ancestry.prediction"] == "EUR"]
probands = probands["IID"]

# filter down variants to mediuml and high impact
voi = variants[variants["impact"].isin(["MED", "HIGH"])]

# filter down variants to those in gens of interest
voi = voi[voi["gene"].isin(genes)]

# reorganize data frame so that rows are genes of interest, columns are IIDs and value are coutns of variants
voi = voi.drop(voi.columns[0:4], axis="columns")
voi = voi.drop(voi.columns[1:4], axis="columns")

# replace -1 (missing) with 0 in preparation for summing columns
voi = voi.replace(to_replace=-1, value=0)

# make list of samples
samples = list(voi)[1:]

# empty table for gene and sums by individusl
table = {}

# loop through genes, summing the genotypes for each indivudal
for g in genes:
    tmp = voi.loc[voi["gene"] == g]
    tmp = tmp.drop(["gene"], axis="columns")
    table[g] = tmp.sum()

# write to file
out.write("\t".join(["gene", "\t".join([str(i) for i in samples])]) + "\n")
for g in genes:
    out.write("\t".join([g, "\t".join([str(i) for i in list(table[g])]) + "\n"]))
