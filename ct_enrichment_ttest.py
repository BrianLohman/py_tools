#!/usr/bin/env python

import pandas as pd
from scipy import stats
import numpy as np
import sys

# phenotype data
m = pd.read_table("June2019_EUR_master_PRS_pheno.txt", sep = '\t')

# variant info
t = pd.read_table("SFARI.FullBreakdown.txt")

# dictionary of probands and siblings in each cagegory of the mutation spectrum
groups = {}

for c in cols:
    tmp = t.loc[ :,t.columns.isin(["Proband", "Sibling", str(c)])]
    tmp = tmp[tmp.iloc[:,2] > 0]
    n = tmp.columns[2]
    k = n.split(".")[0]

    if k == "P1":
        groups[n] = tmp.iloc[:,0]
    if k == "S1":
        groups[n] = tmp.iloc[:,1]

# pick traits to use
traits = ["bmi", "head_circumference", "phrase_delay", "word_delay", "srs_total", "scq_life_total", "father_age_birth"]
probands = []
siblings = []

for key in groups.keys():
    if key.split(".")[0] == "P1":
        probands.append(key)
    if key.split(".")[0] == "S1":
        siblings.append(key)

# empty dict for results
result = {}

# run stats
for group in groups.keys():
    result[group] = {}
    for trait in traits:
        if group.split(".")[0] == "P1":
            p = m.loc[m.family_member == "p1" ,]
            p = m.loc[m.family_member == "p1" ,]
            case = p.loc[p.IID.isin(groups[group]), trait].values.astype(int)
            control = p.loc[~p.IID.isin(groups[group]), trait].values.astype(int)
        elif group.split(".")[0] == "S1":
            p = m.loc[m.family_member == "s1" ,]
            p = m.loc[m.family_member == "s1" ,]
            case = p.loc[p.IID.isin(groups[group]), trait].values.astype(int)
            control = p.loc[~p.IID.isin(groups[group]), trait].values.astype(int)
        else:
            sys.exit("Invalid group")

        result[group][trait] = {"case_mean": np.mean(case), "case_sd": np.std(case), "control_mean": np.mean(control), "control_sd": np.std(control), "pval": stats.ttest_ind(case, control, equal_var = False)[1]}

o = open("CT_deNovo_ttest_results.txt", "w")
print("\t".join(["group", "trait", "pval"]), file = o)

for group in groups.keys():
    for trait in traits:
    if result[group][trait]['pval'] < 0.05:
        print(group, trait, result[group][trait]['pval'], sep = "\t", file = o)

o.close()
