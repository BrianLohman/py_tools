#!/usr/bin/env python
# Given 1)a table with PRS and 2) variants as input:
# 	randomly split half
# 	write files for clinco to test each half: differences in mean PRS by group (defined by n-variants)
# USAGE: python shuffle_split_test.py --input [variant table] --prefix [out file prefix]
import argparse
import sys

import numpy as np
import pandas

parser = argparse.ArgumentParser(
    description="join two variant tables, randomly split half, write clinco input to file"
)
parser.add_argument("--file_a", dest="file_a", help="first variant table")
parser.add_argument("--prefix", dest="prefix", help="prefix for output file names")
args = parser.parse_args()

file_a = pandas.read_table(args.file_a, low_memory=False)

# filters: only high and medium impact mutations in genes with decent sfari score (<=5)
file_a = file_a.loc[file_a["impact"].isin(["MED", "HIGH"])]
file_a = file_a.loc[file_a["sfari_score"] <= 5]

# get rid of meta data
samples = file_a.iloc[:, 8:]

# replace missing with 0
samples = samples.replace(to_replace=-1, value=0)

# make list of samples
ids = list(samples)

# count variatns per individual
counts = pandas.DataFrame(samples_a.sum(axis=0))

# attach sample ids
counts["IID"] = ids_a

# name new columns
counts.columns = ["n_variants", "IID"]

# drop duplicates
counts = counts.drop_duplicates()

# add in prs
simons = pandas.read_table(
    "/scratch/ucgd/lustre/work/u0806040/data/15_Jan_19_Simons_master_ancestry_corrected_PRS.txt",
    low_memory=False,
)
master = simons.merge(counts, on="IID")
probands = master.loc[master["family_member"] == "p1"]
probands = probands.loc[probands["ancestry.prediction"] == "EUR"]
probands = probands.drop(probands.columns[0:122], axis="columns")

print(str(probands.shape))

# randomly split half
x = np.random.rand(len(probands)) < 0.5

half_1 = probands[x]
print(str(len(half_1)) + " probands in first half")

half_2 = probands[~x]
print(str(len(half_2)) + " probands in second half")

# FIRST HALF
# assign bins
summary = half_1["n_variants"].describe()


def which_quartile(i):
    if i <= summary[4]:
        return "1"
    if i >= summary[6]:
        return "3"
    return "2"


half_1["group"] = half_1["n_variants"].apply(which_quartile)

# drop the n_variants column
half_1 = half_1.drop(["n_variants"], axis=1)

# drop the middle group
half_1 = half_1.loc[half_1["group"].isin(["1", "3"])]

# write to file
half_1.to_csv(
    sys.argv[3] + "_half_1_EUR_probands_for_clinco.txt", sep="\t", index=False
)

# SECOND HALF
half_2["group"] = half_2["n_variants"].apply(which_quartile)

half_2 = half_2.drop(["n_variants"], axis=1)

half_2 = half_2.loc[half_2["group"].isin(["1", "3"])]

half_2.to_csv(
    args.prefix + "_half_2_EUR_probands_for_clinco.txt", sep="\t", index=False
)
