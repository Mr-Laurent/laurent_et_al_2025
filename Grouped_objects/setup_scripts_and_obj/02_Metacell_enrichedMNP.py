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

def main():
    os.chdir('./Metacells_EnrichMNP/')
    ## Load the AnnData merged of the 7 MNP enriched samples
    adall = ad.read_h5ad("ad16nov22_reAug24.h5ad")
    adall.X=adall.X.astype("float32", casting='same_kind', copy=True)

    mc.ut.set_name(adall, 'CDsmall')
    print(adall.shape)

    ## 1st sort: remove dying cells + red cells
    hb=['HBB','HBA1','HBA2']
    mts=['MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB','MTRNR2L11','MTRNR2L12','MTRNR2L13','MTRNR2L6','MTRNR2L10','MTRNR2L8','MTRNR2L7','MTRNR2L5','MTRNR2L4','MTRNR2L1','MTRNR2L3']

    excluded_gene_names=mts+hb
    excluded_gene_patterns = ['MT-.*','MT1.*']

    mc.pl.analyze_clean_genes(adall,
        excluded_gene_names=excluded_gene_names,
        excluded_gene_patterns=excluded_gene_patterns,
        random_seed=42)

    mc.pl.pick_clean_genes(adall)
    #set CDsmall.var[clean_gene]: 30351 true (82.95%) out of 36591 bools

    full = adall
    mc.pl.analyze_clean_cells(full,
        properly_sampled_min_cell_total=800,
        properly_sampled_max_cell_total=25000,
        properly_sampled_max_excluded_genes_fraction=0.2)     # Exclude all cells with > 20% of the genes hb, mts
    # set CDsmall.obs[clean_cell]: 31934 true (26.73%) out of 119472 bools
   
    mc.pl.pick_clean_cells(full)
    adall = mc.pl.extract_clean_data(full)

    ## 2nd sort: remove epithelial cells 
    epith=['PLA2G2A','CLCA1','REG4','S100A14','ITLN1','ELF3','PIGR','EPCAM','REG1B','REG1A','REG3A','FABP1','RBP2','SST','FABP2','SPINK1','FABP6','AGR2','AGR3','CLDN3','CLDN4','DEFA6','DEFA5','SPINK4','ALDOB','LCN2','MUC2','KRT8','KRT18','TSPAN8','OLFM4','GPX2','IFI27','PHGR1','MT1G','CLDN7','KRT19','FXYD3','LGALS4','FCGBP','TFF3','TFF1']
    excluded_gene_names=epith
    excluded_gene_patterns = ["None"]

    mc.pl.analyze_clean_genes(adall,
        excluded_gene_names=excluded_gene_names,
        excluded_gene_patterns=excluded_gene_patterns,
        properly_sampled_min_gene_total=0,   
        noisy_lonely_min_gene_normalized_variance=100,  
        random_seed=42)

    mc.pl.pick_clean_genes(adall)

    full = adall
    mc.pl.analyze_clean_cells(
        full,
        properly_sampled_min_cell_total=800,
        properly_sampled_max_cell_total=25000,
        properly_sampled_max_excluded_genes_fraction=0.1)   

    mc.pl.pick_clean_cells(full)
    # 31425 true (98.41%) out of 31934 bools

    adall = mc.pl.extract_clean_data(full)

    ##3rd sort: Selection of genes variables between individuals to remove from the metacell calculations:
    trois=["MALAT1","XIST","JCHAIN"]
    hla=['HLA-F','HLA-G','HLA-A','HLA-E','HLA-C','HLA-B','HLA-DRB5','HLA-DRB1','HLA-DQA1','HLA-DQB1','HLA-DQB1-AS1','HLA-DQA2','HLA-DQB2','HLA-DOB','HLA-DMB','HLA-DMA','HLA-DOA','HLA-DPA1','HLA-DPB1']
    # Using excluded genes, Table S5 of the Leader et al paper: "Single-cell analysis of human non-small cell lung cancer lesions refines tumor classification and patient stratification"
    text_file = open("Leader_genes.txt", "r")
    #  /!\ Need to remove the last newline for it to work
    Leaderexcl = text_file.read().split(',')
    excluded_gene_names=trois+hla+Leaderexcl
    excluded_gene_patterns = ['RPS.*','RPL.*','RPP.*','IGLJ.+','IGLV.+','IGKJ.+','IGHV.+','IGHJ.+'] # We keep IGHD.+ genes

    mc.pl.analyze_clean_genes(adall,
        excluded_gene_names=excluded_gene_names,
        excluded_gene_patterns=excluded_gene_patterns,
        properly_sampled_min_gene_total=0,    
        noisy_lonely_min_gene_normalized_variance=100, 
        random_seed=42)

    mc.pl.pick_clean_genes(adall)
    #set CDsmall.clean.clean.var[clean_gene]: 30043 true (99.12%) out of 30310 bools

    full = adall
    mc.pl.analyze_clean_cells(full,
        properly_sampled_min_cell_total=800,
        properly_sampled_max_cell_total=25000,
        properly_sampled_max_excluded_genes_fraction=0.95)    # 0.95 as we have a lot of genes and that I want to keep plasmablasts

    mc.pl.pick_clean_cells(full)
    # set CDsmall.clean.clean.obs[clean_cell]: 31419 true (99.98%) out of 31425 bools    

    clean = mc.pl.extract_clean_data(full)

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
    max_parallel_piles = 200
    mc.pl.set_max_parallel_piles(max_parallel_piles)  
    mc.pl.divide_and_conquer_pipeline(clean,
                                    forbidden_gene_names=forbidden_gene_names,  #  This step is really long
                                    random_seed=123456,
                                    target_metacell_size=5000)     

    metacells = mc.pl.collect_metacells(clean, name='CDsmall.metacells')
    mc.pl.compute_umap_by_features(metacells, max_top_feature_genes=1000,
                                min_dist=2.0, random_seed=42)
    #1545 floats

    clean.write('all_clean_all_16nov22.h5ad') 
    metacells.write('all__16nov22.h5ad')

if __name__ == "__main__":
    main()