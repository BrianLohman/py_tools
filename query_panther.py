#!/usr/bin/env python

# Submit a request to the Panther database for GO enrichment
import pandas
import requests
from collections import defaultdict

# assign values to form fields
form_data = {'organism':'Homo sapiens',
       'geneList': ('test_genes.txt', open('./test_genes.txt', 'rb')),
       'enrichmentType': 'fullGO_process',
        'correction': 'FDR',
       'type':'enrichment'}

# send request to Panther
response = requests.post('http://pantherdb.org/webservices/garuda/tools/enrichment/VER_2/enrichment.jsp?', files = form_data)

# convert Panther response(bytes) to pandas data farme
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

# convert to pandas data frame
df = pandas.DataFrame.from_dict(result)

# remove duplicate column headers
df = df.iloc[1:]

#print(df.head())

# sort by FDR
#print(df.sort_values(by=['FDR']))

# print results to file
df.sort_values(by=['FDR']).to_csv('test_out.txt', sep = '\t', index = False)
