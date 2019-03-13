# Given two table with PRS and n-variants as input:
# 	combine
# 	randomly split half
#	test each half: differences in mean PRS by group (defined by n-variants)

# USAGE: python shuffle_split_test.py --file_a [variant table] --file_b [variant table] --prefix [out file prefix]

import pandas
import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = 'join two variant tables, randomly split half, write clinco input to file')
parser.add_argument('--file_a', dest = 'file_a', help = 'first variant table')
parser.add_argument('--file_b', dest = 'file_b', help = 'second variant table')
parser.add_argument('--prefix', dest = 'prefix', help = 'prefix for output file names')
args = parser.parse_args()

file_a = pandas.read_table(args.file_a, low_memory = False)
file_b = pandas.read_table(args.file_b, low_memory = False)

# filters: only high and medium impact mutations in genes with decent sfari score (<=5)
file_a = file_a.loc[file_a['impact'].isin(['MED', 'HIGH'])]
file_a = file_a.loc[file_a['sfari_score'] <= 5 ]

file_b = file_b.loc[file_b['impact'].isin(['MED', 'HIGH'])]
file_b = file_b.loc[file_b['sfari_score'] <= 5 ]

# get rid of meta data
samples_a = file_a.iloc[:,8:]
samples_b = file_b.iloc[:,8:]

# replace missing with 0
samples_a = samples_a.replace(to_replace = -1, value = 0)
samples_b = samples_b.replace(to_replace = -1, value = 0)

# make list of samples
ids_a = list(samples_a)
ids_b = list(samples_b)

# count variatns per individual
counts_a = pandas.DataFrame(samples_a.sum(axis = 0))
counts_b = pandas.DataFrame(samples_b.sum(axis = 0))

# attach sample ids
counts_a['IID'] = ids_a
counts_b['IID'] = ids_b

# name new columns
counts_a.columns = ['n_variants', 'IID']
counts_b.columns = ['n_variants', 'IID']

# merge
joint = pandas.concat([counts_a, counts_b], ignore_index = True)
joint = joint.drop_duplicates()

# add in prs
simons = pandas.read_table("/scratch/ucgd/lustre/work/u0806040/data/15_Jan_19_Simons_master_ancestry_corrected_PRS.txt", low_memory = False)
master = simons.merge(joint, on = 'IID')
probands = master.loc[master['family_member'] == 'p1']
probands = probands.loc[probands['ancestry.prediction'] == 'EUR']
probands = probands.drop(probands.columns[0:122], axis = 'columns')

print(str(probands.shape))

# randomly split half
x = np.random.rand(len(probands)) < 0.5

half_1 = probands[x]
print(str(len(half_1)) + " probands in first half")

half_2 = probands[~x]
print(str(len(half_2)) + " probands in second half")

# FIRST HALF
# assign bins
summary = half_1['n_variants'].describe()

def which_quartile(i):
	if i <= summary[4]:
		return "1"
	if i >= summary[6]:
		return "3"
	return "2"

half_1['group'] = half_1['n_variants'].apply(which_quartile)

# drop the n_variants column
half_1 = half_1.drop(['n_variants'], axis = 1)

# drop the middle group
half_1 = half_1.loc[half_1['group'].isin(["1","3"])]

# write to file
half_1.to_csv(sys.argv[3]+"_half_1_EUR_probands_for_clinco.txt", sep = '\t', index = False)

# SECOND HALF
half_2['group'] = half_2['n_variants'].apply(which_quartile)

half_2 = half_2.drop(['n_variants'], axis = 1)

half_2 = half_2.loc[half_2['group'].isin(["1","3"])]

half_2.to_csv(args.prefix+"_half_2_EUR_probands_for_clinco.txt", sep = '\t', index = False)
