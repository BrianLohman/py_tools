import sys
from collections import defaultdict

import toolshed as ts
from quicksect import IntervalTree


# STATUS: INCOMPLETE, not currently developing

# USAGE: python tree_build_regions.py hapmap_combined_chroms.txt regions_with_name_hg37.tab
hapmap = sys.argv[1]
regions = sys.argv[2]
# construct the tree
tree = defaultdict(IntervalTree)

# loop over file and add intervals to tree
for toks in ts.reader(hapmap):
    tree[toks["chr"]].add(
        int(toks["position"]), int(toks["position"]), other=toks["rate_cm_mb"]
    )

# search for an interval in the tree
for toks in (x.strip().split("\t") for x in open(regions)):
    chrom = toks[0]
    start = int(toks[1])
    end = int(toks[2])
    gene = toks[3]
    pathway = toks[4]

    my_search = tree[chrom].search(start, end)
    print(my_search)
