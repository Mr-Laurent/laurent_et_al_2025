#!/usr/bin/python

import anndata as ad
import os
import sys
import pandas as pd
import metacells as mc
import argparse

os.chdir('~/Grouped_objects/setup_scripts_and_obj/')

def main():
    parser = argparse.ArgumentParser(description="Pipeline metacells paramétrable")
    parser.add_argument("--clean", required=True,
                        help=".h5ad entry file (clean)")
    parser.add_argument("--tokeep", required=True,
                        help="csv file of cells to subset")
    parser.add_argument("--target_size", type=int, default=5000,
                        help="target metacell size (default: 5000)")
    parser.add_argument("--clean_out", required=True,
                        help="clean (cells) h5ad output name")
    parser.add_argument("--metacells_out", required=True,
                        help="metacells h5ad output name")

    args = parser.parse_args()

    # Readand load files
    clean = ad.read_h5ad(args.clean)
    tokeep = pd.read_csv(args.excl)
    npex=tokeep.to_numpy()[:,1]
    clean=clean[npex]

    # Make the metacell model
    max_parallel_piles = 200
    mc.pl.set_max_parallel_piles(max_parallel_piles)  
    mc.pl.divide_and_conquer_pipeline(clean,
                                    forbidden_gene_names=[''],  
                                    random_seed=123456,
                                    target_metacell_size=5000)      
    metacells = mc.pl.collect_metacells(clean, name='CDsmall.metacells')
    mc.pl.compute_umap_by_features(metacells, max_top_feature_genes=1000,
                                min_dist=2.0, random_seed=42)  
    
    # Save the cell + metacell objects
    clean.write(args.clean_out)
    metacells.write(args.metacells_out)
    
if __name__ == "__main__":
    main()
