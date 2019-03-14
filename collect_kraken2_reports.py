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

# empty data frame for results
master = pandas.DataFrame()

# loop through kraken2 reports
for report in report_list:
	# get sample ID from filename
	sample = report.split('_')[0]
	sample = sample.split('/')[-1]
	
	# read in as pandas dataframe
	df = pandas.read_csv(report, sep="\t")
	df.columns = ["percent_total", "reads", "taxon_reads", "taxon", "NCBI_taxon_ID","name"]
	
	# strip white space from name columns
	df.name = df.name.str.strip()
	
	# set name as index
	df = df.set_index(df.name)
	
	# keep only the reads as smaller data frame
	df = df['reads']
	df = pandas.DataFrame(df)
	df.columns = [sample]

	# merge each report with the prior	
	master = pandas.merge(master, df, how = 'outer', right_index = True, left_index = True)

	# replace NaN with 0
	master = master.fillna(value = 0)

# print to file
master.to_csv(args.out, sep = "\t")
