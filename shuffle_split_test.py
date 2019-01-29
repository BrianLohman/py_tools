# Given two table with PRS and n-variants as input:
# 	combine
# 	randomly split half
#	test each half: differences in mean PRS by group (defined by n-variants)

import pandas
import sys
import numpy as np

file_a = pandas.read_table(sys.argv[1])
file_b = pandas.read_table(sys.argv[2])

# merge
joint = file_a.merge(file_b, on = 'IID')

# randomly split half
x = np.random.rand(len(joint)) < 0.5

half_1 = joint[x]
half_2 = joint[~x]


# FIRST HALF
# assign bins
summary = half_1['n_variants'].describe()

def which_quartile:
	if i <= summary[4]:
		return "1"
	if i >= summary[6]:
		return "3"
	return "2"

half_1['group'] = half_1['n_variants'].apply(which_quartile)

# drop the n_variants column
half_1 = half_1.drop(['n_variants'], axis = 1)

# drop the middle group
half_1 = half_1[half_1.group != 2]

# write to file
half_1.to_csv("half_1_EUR_probands_for_clinco.txt", sep = '\t', index = False)

# SECOND HALF
half_2['group'] = half_2['n_variants'].apply(which_quartile)

half_2 = half_2.drop(['n_variants'], axis = 1)

half_2 = half_2[half_2.group != 1]

half_2.to_csv("half_2_EUR_probands_for_clinco.txt", sep = '\t', index = False)
