import sys

# USAGE: python split_by_pathway.py [regions file with pathway]

regions = sys.argv[1]

signaling = open("cell_signaling_" + regions[-8:] , "w")

metabolism = open("metabolism_" + regions[-8:] , "w")

neural = open("neural_" + regions[-8:] , "w")

immune = open("immune_" + regions[-8:] , "w")

structure = open("cell_structure_" + regions[-8:] , "w")

signaling_sub = ['Calcium_Signaling_Pathway', 'MAPK_signaling_pathway', 'cell_adhesion_molecules', 'wnt_signaling_pathway', 'mTOR_signaling_pathway', 'vascular_smooth_muscle_contratction', 'ErbB_signaling', 'ECM-receptor_interaction', 'GnRH_signaling_pathway']
neural_sub = ['neuroactive_ligand-reeptor_interaction', 'long-term_potentiation', 'neurotrophin_signaling_pathway', 'long-term_depression', 'axon_guidance']
metabolism_sub = ['melanogenesis', 'lysine_degradation', 'tyrosine_metabolism', 'Phenylalanine_metabolism', 'Aldosterone-regulated_sodium_reabsorption', 'Tryptophan_metabolism', 'Arginine_and_proline_metabolism', 'Ubiquitin_mediated_proteolysi']
immune_sub = ['Leukocyte_transendothelial_migration', 'chemokine_signaling_pathway', 'B_cell_receptor_signaling_pathway']
structure_sub = ['focal_adhesion', 'gap_juction', 'regulation_of_actin_cytoskeleton', 'adherens_junction', 'Endocytosis']

for (chrom, start, end, gene, pathway) in (l.strip().split('\t') for l in open(regions)):
    #print chrom, start, end, gene, pathway
    if pathway in signaling_sub:
        signaling.write('\t'.join([chrom, start, end, gene, pathway+'\n']))

    if pathway in neural_sub:
        neural.write('\t'.join([chrom, start, end, gene, pathway+'\n']))

    if pathway in metabolism_sub:
        metabolism.write('\t'.join([chrom, start, end, gene, pathway+'\n']))
    
    if pathway in immune_sub:
        immune.write('\t'.join([chrom, start, end, gene, pathway+'\n']))

    if pathway in structure_sub:
        structure.write('\t'.join([chrom, start, end, gene, pathway+'\n']))
