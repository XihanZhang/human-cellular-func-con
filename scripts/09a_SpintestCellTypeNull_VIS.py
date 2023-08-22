#!/bin/python

from netneurotools import datasets, freesurfer
import pandas as pd
import numpy as np
from nibabel.freesurfer import read_annot, read_geometry
def _decode_list(vals):
    """ List decoder
    """
    return [val.decode() if hasattr(val, 'decode') else val for val in vals]

# this script will use native space AHBA sample coordinates (x,y,z) to map onto
# a 32k midthickness file that is spatially aligned with native cortical
# geometry (vertices are aligned across individuals)



### Load data and atalas
## 1. set dirs
base_dir = '/gpfs/milgram/project/holmes/xz555/gradient_shift'
gradient_dir = base_dir+'/data/Gradients_Margulies2016'
figure_dir = base_dir+'/figures/surface_plots'
out_dir = base_dir+'/data/NetworkEnrichment'

## 2. load the 18 cell-types
LakeDFC_schaefer400 = pd.read_csv(gradient_dir+'/schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3.csv',header=None,index_col=0)
LakeVIS_schaefer400 = pd.read_csv(gradient_dir+'/schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3.csv',header=0,index_col=0)

## 3. load the annot files of Schaefer400
annot = datasets.fetch_schaefer2018(version='fsaverage', data_dir=None, url=None, resume=True, verbose=1)['400Parcels7Networks']



### Prepare cell type fractions parceled in Schaefer 400
## 1. make 0 to be a small none-zero value
LakeDFC_schaefer400.replace(0,0.0001,inplace=True)
LakeVIS_schaefer400.replace(0,0.0001,inplace=True)

## 2. rename the column names according to annot file
vn, _, names_l = read_annot(annot.lh)
names_l = _decode_list(names_l)
vn, _, names_r = read_annot(annot.rh)
names_r = _decode_list(names_r)
LakeVIS_schaefer400.columns = names_l[1:] + names_r[1:]

## 3. make name list of NaN column (parcel index and name)
# index list of NaN column
NaN_Lake_schaefer400_index = [0] + LakeDFC_schaefer400.columns[LakeDFC_schaefer400.isna().any()].tolist()
# name list of NaN column
NaN_Lake_schaefer400_name = [names_l[0]] + LakeVIS_schaefer400.columns[LakeVIS_schaefer400.isna().any()].tolist()

## 4. make name list of row names (cell-type)
cell_types_DFC = list(LakeDFC_schaefer400.index.values)
cell_types_VIS = list(LakeVIS_schaefer400.index.values)

## 5. drop NaN columns
LakeDFC_schaefer400_drop = LakeDFC_schaefer400.dropna(axis=1)
LakeVIS_schaefer400_drop = LakeVIS_schaefer400.dropna(axis=1)
# save it
file_name = '/Empirical_LakeDFC_Schaefer400_NaNDropped.csv'
file_path = out_dir + file_name
LakeDFC_schaefer400_drop.to_csv(file_path)
file_name = '/Empirical_LakeVIS_Schaefer400_NaNDropped.csv'
file_path = out_dir + file_name
LakeVIS_schaefer400_drop.to_csv(file_path)


### Spin them~
# Cornblath method, more here: https://markello-spatialnulls.netlify.app/spatial_nulls.html
# loop through the 18 cell types
for this_cell in cell_types_VIS:
    print(this_cell)
    brain = LakeVIS_schaefer400_drop.loc[this_cell].to_numpy()
    nulls = freesurfer.spin_data(brain, lhannot=annot.lh, rhannot=annot.rh,n_rotate=5000,
                                 version='fsaverage',drop=NaN_Lake_schaefer400_name,verbose=True)
    file_name = '/Null_' + this_cell + '_Schaefer400_NaNDropped.csv'
    file_path = out_dir + file_name
    pd.DataFrame(nulls).to_csv(file_path)
    
