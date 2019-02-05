# read in a variant table, filter down to med and high impact variants in certain genes

import pandas

variants = pandas.read_table(sys.argv[1], low_memory = False)
genes = ['PLXNA3','CACNA1A','PLXNB1','EP300','CREBBP','PLXNA4','ROBO1','OPRM1','GNAS','SLIT3','ITPR1','EPHB6','NOS1','NTRK1','CACNA1C','ROBO2','GRID1','GRM8','PLCB1','SRGAP3','SEMA5A', 'PTGER3']

voi = variants[variants['impact'].isin(['MED', 'HIGH'])]

voi = voi[voi['gene'].isin(genes)]


