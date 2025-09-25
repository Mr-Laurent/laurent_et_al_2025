#!/usr/bin/python
# python 01_Merge_totalLP_h5.py alltot_s_n_dec22

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pandas as pd
import scipy.sparse as sp
import seaborn as sb
import re
from math import hypot
from matplotlib.collections import LineCollection
from IPython.display import set_matplotlib_formats

os.chdir('~/Grouped_objects/setup_scripts_and_obj/')
directory = './h5ad_totalLP/'

def main():
  # get the sample names from the list of all the .h5ad files 
  grp=[]
  for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".h5ad"): 
      
        base_name, file_ext = os.path.splitext(filename)
        grp.append(base_name)

        # Old file names: split by the point, take part 0, then split by the _ and take the last part (-1)
        # grp.append(filename.split('.')[0].split('_')[-1])

  # Import .h5ad files in python, with a variable corresponding name)    
  for k in grp:
    exec(f"ad{k}t=ad.AnnData.transpose(ad.read_h5ad('{directory}{k}.h5ad'))")

  # Get list of anndata variables
  patvar=[a for a in dir() if re.match(r'ad.+',a)]
  # transform it as the list of anndata matrices so we can concatenate them
  pat2=[]
  for i in range(0,len(patvar)): pat2.append(eval(patvar[i]))
  ad_all = pat2[0].concatenate(*pat2[1:], join='inner',batch_categories=patvar, uns_merge="first", index_unique='-')
  # AnnData object with n_obs x n_vars = 420025 x 36591    obs: 'batch'
  #  and ad_all.obs   gives the name of the object the cell is coming from: adXXXt

  ad.AnnData.write(ad_all,filename="ad%s.h5ad" % sys.argv[1])


if __name__ == "__main__":
    main()
