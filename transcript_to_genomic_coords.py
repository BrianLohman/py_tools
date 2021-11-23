#!/usr/bin/env python
"""
Convert transcript coordinates to genomic coordinates
USAGE: python transcript_to_genomic_coords.py --transcripts [transcripts.txt] --queries [query.txt] --output [output.txt]

Algorithm: genomic coordinate = genomic_start + query_start + sum(deletions up to relevant matching segment)

Asumptions:
1. The transcript names are unique
2. Transcript alignments have no leading hard or soft clips
3. Transcripts are mapped in the forward orientation
"""
import argparse
import re
import sys
from os import sep

parser = argparse.ArgumentParser(
    description="Convert transcript coordinates to genomic corrdinates"
)
parser.add_argument(
    "-t",
    "--transcripts",
    dest="transcripts",
    help="A four column (tab-separated) file containing the transcripts",
)
parser.add_argument(
    "-q",
    "--queries",
    dest="query",
    help="A two column (tab-separated) file indicating a set of queries",
)
parser.add_argument(
    "-o",
    "--output",
    dest="output",
    help="File name of output. Four column tab separated file with one row for each of the input queries",
)
args = parser.parse_args()

# Dictionary of transcript chr, start and cigar, keyed by transcript name
transcript_dict = {}

for transcript_id, chr, transcript_start, cigar in (
    x.strip().split("\t") for x in open(args.transcripts)
):
    transcript_dict[(transcript_id)] = [chr, transcript_start, cigar]
# Make empty list for results
results = []

# Loop through queries
for query_id, query_start in (x.strip().split("\t") for x in open(args.query)):
    # check that query_id is in transcript_dict
    if query_id in transcript_dict.keys():
        # get chr, transcript_start, and cigar
        chr, transcript_start, cigar = transcript_dict[query_id]
    else:
        sys.exit("query transcript ID not in transcript dictionary: " + query_id)

    # check that query_start is less than length of transcript
    if int(query_start) > sum(map(int, re.findall(r"(\d+)[MI]", cigar))):
        sys.exit("query start is greater than transcript length: " + query_start)

    # cigars that lead with a match
    if cigar.split("M")[0].isdecimal():
        # if the query start is in the first matching segment
        if int(query_start) < int(cigar.split("M")[0]):
            genomic_coord = int(transcript_start) + int(query_start)
            results.append("\t".join([query_id, query_start, chr, str(genomic_coord)]))
        # if the query starts after the first matching segment
        else:
            # how many matches or insertions are needed to get to query start?

            # break cigar into segments
            l = list(re.findall(r"(\d+)([A-Z])", cigar))

            # set the length of the first matching segment
            c = int(l[0][0])

            # empty space to sum deletions
            d = 0

            # loop through the cigar segments to find how many matches or insertions are needed to get to query start
            # use entire cigar when necessary
            for x, i in enumerate(l):
                if x < len(l) - 1:
                    if "D" in l[x + 1]:
                        d += int(l[x + 1][0])
                        continue
                    else:
                        y = int(l[x][0])
                else:
                    y = sum(map(int, re.findall(r"(\d+)[D]", cigar)))
                # test to see if the query start is less than the sum of first match and next match
                if int(query_start) < int(c) + int(y):
                    genomic_coord = int(transcript_start) + int(query_start) + d
                    results.append(
                        "\t".join([query_id, query_start, chr, str(genomic_coord)])
                    )
                    break
                else:
                    continue
    else:
        sys.exit(
            "cigar does not start with match. Are soft/hard clips included? "
            + transcript_dict[query_id]
        )
# Write results to output file
print("\n".join(results), file=open(args.output, "w"))
