#!/usr/bin/env python

# Simulate the expected functional enrichment (GO BP) of the average proband in SSC

import pandas
import random
import requests
from collections import defaultdict
from collections import Counter
import statistics as stats

# Read in Simons data
master = pandas.read_table('15_Jan_19_Simons_master_ancestry_corrected_PRS.txt')
probands = master.loc[master.phenotype == 2]
proband_ids = list(probands.simons_id)

# Read in table of variants (med and high impact from GATK)
variants = pandas.read_table('med_high_GATK_variants.txt')

# Make space to store GO terms during the loop
d = {}

for i in range(1,1001):
    sample = random.sample(proband_ids, 413)

    # Select randomly sampled probands
    variants_of_interest = variants.loc[variants.SampleID.isin(sample)]

    # Select genes with at least 1 variant
    gene_symbols_of_interest = list(set(variants_of_interest.SYMBOL))
    print(str(len(gene_symbols_of_interest))+' unique genes of interest on iteration '+str(i))

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
    
    # convert FDR to float
    df.FDR = df.FDR.astype(float)
   
    # drop gene name and remove duplicates
    df = df.drop('GeneId', axis = 1)
    df = df.drop_duplicates()
 
    # list of uniq GO groups
    uniq_go = list(set(df.Name))

    # build dictionary of go name and FDR
    for g in uniq_go:
        if g in d:
            d[g] = d[g].append(df[df.Name == g].FDR)
        else:
            d[g] = df[df.Name == g].FDR

# gather results
out = open('1000_simulations_GO_mean_sd.txt', 'w')

for key, values in d.items():
    if len(values) <= 1:
        print('\t'.join([key, str(stats.mean(values)), 'NA']), file = out)
    else:
        print('\t'.join([key, str(stats.mean(values)), str(stats.stdev(values))]), file = out)
