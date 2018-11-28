from __future__ import print_function
import sys
import glob

# convert the chr:bp to rsID in summary stats (.assoc files)

# USAGE: python gwas_chrBP_2_rsID.py [source (hg37) dictionary]

# empty dictionary
pos_to_rs = {}

# build dictionary
for i, line in enumerate(open(sys.argv[1])):
    fields=line.strip().split('\t')
    chrom=fields[0]
    pos=fields[2]
    rsid=fields[3]
    pos_to_rs[(chrom,pos)]=rsid

# read in gwas files
gwas = glob.glob('*.assoc')

# loop through gwas files, looking up chr and BP in dictionary
# print unmapped rsIDs to _add_rsID.log file
# print results to .hg38.assoc

# lookup
for g in gwas:
    out = open(g[:-6] + ".rsID.assoc", "w")
    log = open(g[:-6] + "_add_rsID.log", "w")
    out.write("\t".join(["SNP", "CHR", "BP", "A1", "A2", "BETA", "P\n"]))
    for i, fields in enumerate(x.strip().split('\t') for x in open(g)):
        if i == 0 and "SNP" in fields:
            continue
        chrom=fields[1]
        pos=fields[2]
        ref=fields[3]
        alt=fields[4]
        beta=fields[5]
        p=fields[6]
        try:
            rsid=pos_to_rs[(chrom,pos)]
        except KeyError:
            print("%s not found" % (','.join([chrom,pos])), file=log)
            continue
        out.write('\t'.join([rsid,chrom,pos,ref,alt,beta,p+'\n']))
    out.close()
    log.close()
