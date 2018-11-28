import glob

# read in files
gene_list = glob.glob('*.txt')

# blank space for final result
out = open('genes_by_pathway.tab', "w")

# loop through the files, extracting gene names as column 1 and file name as column 2
for g in gene_list:
    tmp = []
    for i, toks in enumerate(x.strip().split() for x in open(g)):
        if i == 0: continue
        gene = toks[1].split(';')[0]
        print gene
        pathway = g[:-13]
        print pathway

        tmp = [gene,pathway]
        out.write('\t'.join(tmp)+'\n') 
