#!/usr/bin/env python
import sys

# USAGE: python subset_snps_bim_UCSC.py [.bim file with snps of interest] [master list of snps] > [snps in .bim]

# to speed lookup times, the snps of interst from a master dictionary from UCSC
# for use before lift.py, a dictionary based method of liftover based on rsIDs

rsids = set()
for line in open(sys.argv[1]):
    fields = line.split("\t", 3)
    rsids.add(fields[1])

print >>sys.stderr, len(rsids)
found = set()

# write errors to log
log = open(sys.argv[1][:-4] + ".subset.log", "w")

keep = {}
for line in open(sys.argv[2]):
    fields = line.rstrip().split("\t")
    if not fields[3] in rsids:
        continue
    if fields[3] in keep:
        if "alt" in fields[0]:
            continue
        try:
            assert "alt" in keep[fields[3]], (line, keep[fields[3]])
        except AssertionError as err:
            log.write(str(err) + "\n")
            continue
    keep[fields[3]] = line
    found.add(fields[3])

for k, l in keep.iteritems():
    print(l),
