#!/usr/bin/env python3
import anndata as ad
import matplotlib.pyplot as plt
import metacells as mc
import numpy as np
import os
import sys
import pandas as pd
import scipy.sparse as sp
import seaborn as sb
from math import hypot
from matplotlib.collections import LineCollection
from IPython.display import set_matplotlib_formats
import re

def main():
    os.chdir('./Metacells_TotalLP/')
    ## Load the AnnData merged of all samples
    adall = ad.read_h5ad("adalltot_v2_s_n_dec22.h5ad")
    adall.X=adall.X.astype("float64", casting='same_kind', copy=True)
    adall.X=adall.X.astype("float32", casting='same_kind', copy=True)
    # AnnData object with n_obs x n_vars = 489621 x 21932
    mc.ut.set_name(adall, 'CDsmall')
    print(adall.shape)

    ## 1st sort: remove dying cells + red cells
    hb=['HBB','HBA1','HBA2']
    mts=['MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB','MTRNR2L11','MTRNR2L12','MTRNR2L13','MTRNR2L6','MTRNR2L10','MTRNR2L8','MTRNR2L7','MTRNR2L5','MTRNR2L4','MTRNR2L1','MTRNR2L3']

    excluded_gene_names=mts+hb
    excluded_gene_patterns = ['MT-.*','MT1.*']
    mc.pl.analyze_clean_genes(adall,excluded_gene_names=excluded_gene_names,
                            excluded_gene_patterns=excluded_gene_patterns,random_seed=42)
    mc.pl.pick_clean_genes(adall)
    #set CDsmall.var[clean_gene]: 20630 true (94.06%) out of 21932 bools


    full = adall
    mc.pl.analyze_clean_cells(full,
                            properly_sampled_min_cell_total=800,
                            properly_sampled_max_cell_total=25000,
                            properly_sampled_max_excluded_genes_fraction=0.2)    # Exclude all cells with > 20% of the genes hb, mts
    # set CDsmall.obs[clean_cell]:  240819 true (49.18%) out of 489621 bools
    mc.pl.pick_clean_cells(full)
    adall = mc.pl.extract_clean_data(full)

    ## 2nd sort: remove epithelial cells 
    epith=['PLA2G2A','CLCA1','REG4','S100A14','ITLN1','ELF3','PIGR','EPCAM','REG1B','REG1A','REG3A','FABP1','RBP2','SST','FABP2','SPINK1','FABP6','AGR2','AGR3','CLDN3','CLDN4','DEFA6','DEFA5','SPINK4','ALDOB','LCN2','MUC2','KRT8','KRT18','TSPAN8','OLFM4','GPX2','IFI27','PHGR1','MT1G','CLDN7','KRT19','FXYD3','LGALS4','FCGBP','TFF3','TFF1']
    excluded_gene_names=epith
    excluded_gene_patterns = ["None"]
    mc.pl.analyze_clean_genes(adall,excluded_gene_names=excluded_gene_names,excluded_gene_patterns=excluded_gene_patterns,
                            properly_sampled_min_gene_total=0,
                            noisy_lonely_min_gene_normalized_variance=100,  #### NEW
                            random_seed=42)
    mc.pl.pick_clean_genes(adall)

    full = adall
    mc.pl.analyze_clean_cells(full,properly_sampled_min_cell_total=800,
                            properly_sampled_max_cell_total=25000, 
                            properly_sampled_max_excluded_genes_fraction=0.1)   
    mc.pl.pick_clean_cells(full)
    # 232996 true (96.75%) out of 240819 bools

    adall = mc.pl.extract_clean_data(full)

    ##3rd sort: Selection of genes variables between individuals to remove from the metacell calculations:
    list(filter(lambda x: re.search(r'IGHD', x), list(adall.var_names)))
    trois=["MALAT1","XIST","JCHAIN"]
    hla=['HLA-F','HLA-G','HLA-A','HLA-E','HLA-C','HLA-B','HLA-DRB5','HLA-DRB1','HLA-DQA1','HLA-DQB1','HLA-DQB1-AS1','HLA-DQA2','HLA-DQB2','HLA-DOB','HLA-DMB','HLA-DMA','HLA-DOA','HLA-DPA1','HLA-DPB1']
    # Using excluded genes, Table S5 of the Leader et al papaer: Single-cell analysis of human non-small cell lung cancer lesions refines tumor classification and patient stratification
    text_file = open("Leader_genes.txt", "r")

    Leaderexcl = text_file.read().split(',')
    excluded_gene_names=trois+hla+Leaderexcl
    excluded_gene_patterns = ['RPS.*','RPL.*','RPP.*','IGLJ.+','IGLV.+','IGKJ.+','IGHV.+','IGHJ.+'] # ADAPTATION NEW : NO MORE EXCLUSION OF IGHD.+
    
    mc.pl.analyze_clean_genes(adall, excluded_gene_names=excluded_gene_names,excluded_gene_patterns=excluded_gene_patterns,
                            properly_sampled_min_gene_total=0,
                            noisy_lonely_min_gene_normalized_variance=100, #### NEW
                            random_seed=42)

    mc.pl.pick_clean_genes(adall)
    #set CDsmall.clean.clean.var[clean_gene]: 20129 true (97.77%) out of 20589 bools

    full = adall
    mc.pl.analyze_clean_cells(full,properly_sampled_min_cell_total=800,
                            properly_sampled_max_cell_total=25000,
                            properly_sampled_max_excluded_genes_fraction=0.95)   
    mc.pl.pick_clean_cells(full)
    # set CDsmall.clean.clean.obs[clean_cell]:  232842 true (99.93%) out of 232996 bools     
    clean = mc.pl.extract_clean_data(full)
    list(filter(lambda x: re.search(r'IGHD', x), list(clean.var_names)))
    #  text_file2 = open("genes1156_topPC1_noBexclu.txt", "r")
    #   suspect_gene_names= text_file2.read().split(',')    # V3 : remove it

    suspect_gene_names=['None']    
    suspect_gene_patterns = ['None']
    suspect_genes_mask = mc.tl.find_named_genes(clean, names=suspect_gene_names,
                                                patterns=suspect_gene_patterns)
    suspect_gene_names = sorted(clean.var_names[suspect_genes_mask])
    mc.pl.relate_genes(clean, random_seed=42)

    module_of_genes = clean.var['related_genes_module']
    suspect_gene_modules = np.unique(module_of_genes[suspect_genes_mask])          
    suspect_gene_modules = suspect_gene_modules[suspect_gene_modules >= 0]
    print(suspect_gene_modules)

    similarity_of_genes = mc.ut.get_vv_frame(clean, 'related_genes_similarity')
    forbidden_genes_mask = suspect_genes_mask
    forbidden_gene_names=sorted(clean.var_names[forbidden_genes_mask])
    print(' '.join(forbidden_gene_names))
    len(' '.join(forbidden_gene_names))
    #0

    # Calculation of metacells
    max_parallel_piles = mc.pl.guess_max_parallel_piles(clean)
    print(max_parallel_piles)
    #       288
    # max_parallel_piles = 200
    mc.pl.set_max_parallel_piles(max_parallel_piles)  
    mc.pl.divide_and_conquer_pipeline(clean,forbidden_gene_names=forbidden_gene_names,  
                                    random_seed=123456, target_metacell_size=40000)     

    metacells = mc.pl.collect_metacells(clean, name='CDsmall.metacells')
    mc.pl.compute_umap_by_features(metacells, max_top_feature_genes=1000,
                                min_dist=2.0, random_seed=42)

    #9091 floats
    umap_x = mc.ut.get_o_numpy(metacells, 'umap_x')
    umap_y = mc.ut.get_o_numpy(metacells, 'umap_y')
    plot = sb.scatterplot(x=umap_x, y=umap_y)

    umap_edges = sp.coo_matrix(mc.ut.get_oo_proper(metacells, 'obs_outgoing_weights'))
    min_long_edge_size = 4
    sb.set()
    plot = sb.scatterplot(x=umap_x, y=umap_y)
    plot.get_figure().savefig("cleanalltot_v2_s_n_dec22.png") 

    clean.write('all_clean_all_alltot_v2_s_n_dec22.h5ad') 
    metacells.write('all__alltot_v2_s_n_dec22.h5ad')



if __name__ == "__main__":
    main()