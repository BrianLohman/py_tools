#!/usr/bin/env python

# Gather kraken2 reports from a directory and assemble them into 1 large table
# No filters applied
# Rows are observations from the kraken2 report
# Columns are individuals
# 0 is used where a taxon was not observed in a given sample

# USAGE: collect_kraken2_reports.py --dir [directory to look for reports in] --out [results file name] 

import argparse
import pandas
import glob

# build arg parser
parser = argparse.ArgumentParser(description = 'parse and assemble kraken2 reports')
parser.add_argument('--dir', help = 'directory to look for kraken2 reports in')
parser.add_argument('--out', help = 'output filename')
args = parser.parse_args()

# glob in files
report_list = glob.glob(args.dir+'/*_kraken2_report.txt')

samples = []
# loop through kraken2 reports
for report in report_list:
	print(report)
	
	# get sample ID from filename
	sample_name = report.split('_')[-3]
	sample_name = sample_name.split('/')[-1]
	
	# read in as pandas dataframe
	sample = pandas.read_csv(report, sep="\t")
	sample.columns = ["percent_total", "reads", "taxon_reads", "taxon", "NCBI_taxon_ID","name"]
	
	# strip white space from name columns
	sample.name = sample.name.str.strip()
	
	# drop undesired columns
	sample = sample[['name','reads']]	
	sample.columns = ['name', sample_name]

	# append current sample to list of samples
	samples.append(sample)

# concat
master = pandas.concat(samples).groupby('name', as_index = False, sort = False).first().set_index('name').fillna('0')

# print to file
master.to_csv(args.out, sep = "\t")
