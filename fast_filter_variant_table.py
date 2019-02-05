#!/usr/bin/pyspark

# read in variant tables and filter down to genes in given list
# write the number of variants that each individual has in these genes to file


# build RDD
master = sc.textFile("./neural_table.txt", 5).map(lambda line: line.split("\t"))


