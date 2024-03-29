#!/usr/bin/env python

# assign individuals to a bin by their mutation load
# uses the input of variant_table.nim
# currently only considering probands of EUR ancestry, as defined by peddy

import argparse
import pandas

parser = argparse.ArgumentParser(
    description="assign individuals to bins based on mutation load, after applying filters"
)

parser.add_argument("-v", "--variants", dest="variants", help="variant table")
parser.add_argument(
    "-g", "--gnomad", dest="gnomad", default=1, type=float, help="gnomad AF threshold"
)
parser.add_argument(
    "-s", "--sfari", dest="sfari", default=5, type=int, help="sfari score threshold"
)
parser.add_argument("-o", "--out", dest="outfile", help="name of results file")
args = parser.parse_args()

# read in variant table
variant_table = pandas.read_table(args.variants, low_memory=False)

# drop any duplicates and write to file
variant_table = variant_table.drop_duplicates()
variant_table.to_csv("dedup_" + args.variants, sep="\t", index=False)

# sample down to variants that pass filters
# predicated impact
variant_table = variant_table.loc[variant_table["impact"].isin(["MED", "HIGH"])]

# sfari score
variant_table = variant_table.loc[variant_table["sfari_score"] <= args.sfari]

# gnomad allele frequency
variant_table = variant_table.loc[variant_table["gnomAD_AF"] <= args.gnomad]

# sum accros individuals into new data frame
# swap all occurances of -1 for 0
samples = variant_table.iloc[:, 8:]

samples = samples.replace(to_replace=-1, value=0)
ids = list(samples)
counts = pandas.DataFrame(samples.sum(axis=0))
counts["IID"] = ids
counts.columns = ["n_variants", "IID"]

# add the PRS by ID
simons = pandas.read_table(
    "/scratch/ucgd/lustre/work/u0806040/data/15_Jan_19_Simons_master_ancestry_corrected_PRS.txt",
    low_memory=False,
)
master = simons.merge(counts, on="IID")
probands = master.loc[master["family_member"] == "p1"]
probands = probands.loc[probands["ancestry.prediction"] == "EUR"]
probands = probands.drop(probands.columns[0:122], axis="columns")

# assign bins
summary = probands["n_variants"].describe()


def which_quartile(i):
    if i <= summary[4]:
        return "1"
    if i >= summary[6]:
        return "3"
    return "2"


probands["group"] = probands["n_variants"].apply(which_quartile)

# drop the n_variants columnn
# probands = probands.drop(['n_variants'], axis = 1)

# write to file
probands.to_csv(args.outfile, sep="\t", index=False)
