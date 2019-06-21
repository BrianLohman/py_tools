#!/usr/bin/env python

from __future__ import print_function
import sys
import cyvcf2
print(cyvcf2, cyvcf2.__version__, file=sys.stderr)
from cyvcf2 import VCF, Writer
import time

sites = []
for toks in (x.rstrip().split("\t") for x in open("1000G.mac5eur.hm3.bim")):
    sites.append((toks[0], int(toks[3]), toks[4], toks[5]))


def match_sites(args):
    vcf_path, sites = args
    vcf = VCF(vcf_path)
    matches = []
    t0 = time.time()
    for s in sites:
        for v in vcf("%s:%d-%d" % (s[0], s[1] - 1, s[1])):
            if len(v.ALT) > 1: continue
            ra = (v.REF, v.ALT[0])
            if ra != (s[2], s[3]) and ra != (s[3], s[2]):
                continue
            if v.start != s[1] - 1: continue
            matches.append(v)

    print("found %d out of %d sites in %.1f seconds" % (len(matches),
        len(sites), time.time() - t0), file=sys.stderr)
    return matches


print("n-sites: %d" % len(sites), file=sys.stderr)
vcf_path = "/scratch/ucgd/lustre/work/u6000771/Projects/2016/UtahAutism/utah-519.vcf.gz"

#import multiprocessing as mp
from multiprocessing.dummy import Pool

p = Pool(16)

wtr = Writer("-", VCF(vcf_path))

step = 500
for res in p.imap(match_sites, ((vcf_path, sites[idx:idx+step]) for idx in range(0, len(sites), step))):

    for v in res:
        wtr.write_record(v)

wtr.close()
