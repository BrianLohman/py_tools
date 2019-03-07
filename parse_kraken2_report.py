#!/usr/bin/env python

# parse the output of a kraken2 report
# apply various filters

import argparse

parser = argparse.ArgumentParser(description = 'parse a kraken2 report, applying filters, and summarizing')
parser.add_argument('-r', '--min-reads', dest = 'min_reads', type = int, default = 0, help = 'minimum number of reads required to include')
parser.add_argument('-t', '--taxon', dest = 'taxon', type = list, default = ['U', 'R', 'D', 'k', 'P', 'C', 'O', 'F', 'G', 'S'],  help = 'taxon levels to include')
parser.add_argument('-k', '--kraken2-report', dest = 'report', help = 'kraken2 report')
parser.add_argument('-o', '--out', dest = 'out', help = 'filename of output')
args = parser.parse_args()

# open outfile
with open(args.out, "w") as out, open(args.report) as report:

	# read in report
	for (percent_total, reads, taxon_reads, taxon, NCBI_taxon_ID, name) in (l.strip().split("\t") for l in report):
		if int(reads) < args.min_reads:
			continue
		if taxon not in args.taxon:
			continue

		print(percent_total, reads, taxon_reads, taxon, NCBI_taxon_ID, name.strip(), sep = '\t', file = out)
