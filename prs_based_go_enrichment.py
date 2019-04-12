#!/usr/bin/env python

# REQUIRES 15_Jan_19_Simons_master_ancestry_corrected_PRS.txt AND med_high_GATK_variants.txt IN WORKING DIRECTORY

import pandas
import random
import requests
import argparse
from collections import defaultdict
from collections import Counter

parser = argparse.ArgumentParser(description='Test for GO enrichment with Panther')
parser.add_argument('-n', '--sample_size', dest = 'n', help = 'number of individals in each group')
parser.add_argument('-o', '--result', dest = 'o', help = 'name of file contaning results')
parser.add_arugment('-p', '--prs', dest = 'prs', help = 'PRS or other column to generate groups')
parser.add_argument('-q', '--quartile', dest = 'q', help = 'quartile (or other such -ile), either "lower" or "upper"')
args = parser.parse_args()

# Read in Simons data
master = pandas.read_table('15_Jan_19_Simons_master_ancestry_corrected_PRS.txt')
probands = master.loc[master.phenotype == 2]
proband_ids = list(probands.simons_id)

# Read in table of variants (med and high impact from GATK)
variants = pandas.read_table('med_high_GATK_variants.txt')

# select probands based on some criteria
prs = probands.sort_values(by=args.prs+"_ancestry_resid")

# select quartile of interest and get IDs
if args.q == "lower":
	sample = prs.head(int(args.n)).simons_id
if args.q == "upper":
	sample = prs.tail(int(args.n)).simons_id

# Select variants from those probands
variants_of_interest = variants.loc[variants.SampleID.isin(sample)]

# Select genes with at least 1 variant
gene_symbols_of_interest = list(set(variants_of_interest.SYMBOL))
print(str(len(gene_symbols_of_interest))+' unique genes of interest')

# Write genes of interest to file
with open('./test_genes.txt', 'w') as f:
    for item in gene_symbols_of_interest:
        f.write("%s\n" % item)

# Submit a request to the Panther database for GO enrichment
# Assign values to form fields
form_data = {'organism':'Homo sapiens',
       'geneList': ('test_genes.txt', open('./test_genes.txt', 'rb')),
       'enrichmentType': 'fullGO_process',
        'correction': 'FDR',
       'type':'enrichment'}

# Send request to Panther
response = requests.post('http://pantherdb.org/webservices/garuda/tools/enrichment/VER_2/enrichment.jsp?', files = form_data)

# Convert Panther response(bytes) to pandas data farme
text = response.content.decode("utf-8")

header_keys = ['Id','Name','GeneId','P-value', 'FDR']
header2idx = dict(zip(range(len(header_keys)), header_keys))

result = defaultdict(list)

for i, line in enumerate(text.split('\n')):
    if line in ('\r', ''):
        continue

    for idx, toks in enumerate(line.split('\t')):
        key = header2idx[idx]
        result[key].append(toks)

#print([(k, len(result[k])) for k in result])

# Convert to pandas data frame
df = pandas.DataFrame.from_dict(result)

# Remove duplicate column headers
df = df.iloc[1:]

# Drop gene column and subsequent duplicate rows
df = df.drop(['GeneId'], axis = 1)
df = df.drop_duplicates()

# Convert FDR to float
df.FDR = df.FDR.astype(float)

# Only significant after FDR correction
df = df[df['FDR'] < 0.05]

# check to see if anything is left                                                                         .
if df.shape[0] == 0:
        print("No significant enrichment")

else:
	# Sort by FDR and print to file
	df = df.sort_values(by=['FDR'])
	df.to_csv(args.o, index = False, sep = '\t')
