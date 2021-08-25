import sys

import numpy as np
import pandas

# USAGE: python list_build_regions.py sorted_hapmap_combined_chroms_hg37.out no_hap_sorted_regions_with_name_hg37.tab 2 > pathway_haplotype_regions_hg37.txt

# ASSUMES: the list of genes and hapmap markers have been sorted by chromsome

# define input
hapmap = sys.argv[1]
genes = sys.argv[2]

# set rate of recombination threshold
threshold = int(sys.argv[3])

# read in hapmap data and filter by threshold
hapmap = pandas.read_table(hapmap, header=None)
hapmap.columns = ["chrom", "bp", "rate"]
hapmap["bp"] = pandas.to_numeric(hapmap["bp"])
hapmap["rate"] = pandas.to_numeric(hapmap["rate"])

hapmap = hapmap[hapmap["rate"] > threshold]
# print hapmap.head(10)

# read in gene list
genes = pandas.read_table(genes, header=None)
genes.columns = ["chrom", "start", "end", "gene", "pathway"]
# print genes.head(10)

# loop through chromosomes
for chrom in set(genes.chrom):
    if chrom == "chrY" or chrom == "chr":  # there is no Y chrom in the hapmap data
        continue
    chrom_hapmap = hapmap[hapmap.chrom == chrom]
    chrom_hapmap = chrom_hapmap.sort_values("bp")
    # print chrom_hapmap.head(10)
    chrom_genes = genes[genes.chrom == chrom]

    # loop through genes within chromosome, subtracting location of start and end from every hapmap marker
    for i, row in chrom_genes.iterrows():
        diffs_start = row.start - chrom_hapmap.bp
        diffs_end = row.end - chrom_hapmap.bp

        # check for the lowest value that is > 0, becuase column is sorted
        min_val, min_idx = None, None
        for i, d in enumerate(diffs_start):
            if i == 0:
                min_val = d
                continue

            if d > min_val:
                continue
            if d < min_val:
                min_val = d
                min_idx = i

            if d < 0:
                break

        # check for the first value that is < 0, becuase column is sorted
        max_val, max_idx = None, None
        for i, d in enumerate(diffs_end):
            if d > 0:
                continue
            else:
                max_val = d
                max_idx = i
                break

        # for some cases the end position is off the hapmap number line
        # if true, print the end position
        if max_idx is None:
            # print min_idx, len(chrom_hapmap) - 1
            # 1/0
            print(
                join(
                    [
                        row.chrom,
                        str(min_position["bp"]),
                        str(row.end),
                        row.gene,
                        row.pathway,
                    ]
                ),
                sep="\t",
            )
            continue

        # lookup chosen index in array
        # print "min_val = ", min_val, "min_idx = ", min_idx
        # print "max_val = ", max_val, "max_idx = ", max_idx
        min_position = chrom_hapmap.iloc[min_idx].to_dict()
        max_position = chrom_hapmap.iloc[max_idx].to_dict()

        # print results to file along with unchanged tokens
        print(
            join(
                [
                    row.chrom,
                    str(min_position["bp"]),
                    str(max_position["bp"]),
                    row.gene,
                    row.pathway,
                ]
            ),
            sep="\t",
        )
