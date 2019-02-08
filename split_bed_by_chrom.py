# split a bed file by chrom
# output a new bed file for each chrom, keeping the input file prefix

import sys
import pandas

bedfile = sys.argv[1]
prefix = sys.argv[2]

bed = pandas.read_table(bedfile, header = None)
bed.columns = ['chrom', 'start', 'end']

# list of chroms
chroms = bed.chrom.unique()

# iterate over list, subsetting and writing to file
for chrom in chroms:
	tmp = bed.loc[bed['chrom'] == chrom ]
	tmp.to_csv(prefix+'_'+chrom+'.bed', index = False, sep = '\t')
