#!/usr/bin/env python

import sys

# counters
error_count = 0
unmapped_count = 0

# name output based on input
outfile = open(sys.argv[1][:-3] + "_bases.txt", "w")

for toks in (x.strip().split("\t") for x in open(sys.argv[1])):
    # skip headers
    if str(toks[0]).startswith("@"):
        continue

    # skip reads that don't map to reference
    if toks[2] == "*":
        continue

    # skip reads with unmapped flags
    if toks[2] in '"{}"'.format([77, 73, 69, 117, 153, 141, 137, 133]):
        unmapped_count += 1
        continue

    if toks[1] in '"{}"'.format([65, 129]):
        print(
            "\t".join(
                [
                    str(toks[9].count("A")),
                    str(toks[9].count("T")),
                    str(toks[9].count("G")),
                    str(toks[9].count("C")),
                    str(toks[9].count("N")),
                    str(len(toks[9])),
                ]
            ),
            file=outfile,
        )
    elif toks[1] in '"{}"'.format([99, 97, 89, 83, 81, 163, 161, 147, 145, 113]):
        print(
            "\t".join(
                [
                    str(toks[9].count("T")),
                    str(toks[9].count("A")),
                    str(toks[9].count("C")),
                    str(toks[9].count("G")),
                    str(toks[9].count("N")),
                    str(len(toks[9])),
                ]
            ),
            file=outfile,
        )
    else:
        error_count += 1

# print error counts to screen
print("Reads unmapped " + str(unmapped_count))
print("Reads with errors " + str(error_count))
