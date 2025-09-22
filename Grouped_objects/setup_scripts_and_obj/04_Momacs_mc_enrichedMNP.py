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

os.chdir('./Metacells_EnrichMNP/')
clean = ad.read_h5ad("all_clean_all_16nov22_delStro.h5ad")
excl = pd.read_csv("cells_to_keep_Macs_inMNPenrich.csv")
npex=excl.to_numpy()[:,1]
clean=clean[npex]
max_parallel_piles = 200
mc.pl.set_max_parallel_piles(max_parallel_piles)  
mc.pl.divide_and_conquer_pipeline(clean,forbidden_gene_names=[''],  
                                  random_seed=123456,
                                  target_metacell_size=5000)      
metacells = mc.pl.collect_metacells(clean, name='CDsmall.metacells')
mc.pl.compute_umap_by_features(metacells, max_top_feature_genes=1000,
                               min_dist=2.0, random_seed=42)   # 802 floats
umap_x = mc.ut.get_o_numpy(metacells, 'umap_x')
umap_y = mc.ut.get_o_numpy(metacells, 'umap_y')
plot = sb.scatterplot(x=umap_x, y=umap_y)
umap_edges = sp.coo_matrix(mc.ut.get_oo_proper(metacells, 'obs_outgoing_weights'))
min_long_edge_size = 4
sb.set()
plot = sb.scatterplot(x=umap_x, y=umap_y)

plot.get_figure().savefig("clean16nov22Macs.png") 
clean.write('all_clean_all_16nov22_Macs.h5ad') 
metacells.write('all__16nov22_Macs.h5ad')
