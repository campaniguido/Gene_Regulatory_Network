import scanpy as sc
import numpy as np
import pandas as pd
import os
from GENIE3 import GENIE3
import matplotlib.pyplot as plt
import sys
import dask.dataframe as dd

path_h5ad = '/Users/ieo6943/Documents/Guido/Albi/clustered.h5ad/MDA_chemo/data/SCT/clustered.h5ad'
##
path_genes = '/Users/ieo6943/Downloads/genes.csv'

path_h5ad = '/hpcnfs/data/PGP/acossa/archive_Cellula/MDA_chemo/data/SCT/clustered.h5ad'
path_results = '/hpcnfs/scratch/PGP/gcampani/ALBI/GRN/'
path_genes = '/hpcnfs/scratch/PGP/acossa/breast_albi/MDA/single_cell/scripts/GRNs/genes.csv'

origin = sys.argv[1]
treatment = sys.argv[2]
nthreads = int(sys.argv[3])
type_GRN = sys.argv[4]
size_subset = int(sys.argv[5])
n_trees = int(sys.argv[6])

#origin = 'PT'
#treatment = 'treated'
#nthreads = int(20)
#size_subset = 30

os.makedirs(path_results, exist_ok=True)

cluster_h5ad = sc.read_h5ad(path_h5ad)

genes = list(pd.read_csv(path_genes, index_col=0).index)

adata_NT = cluster_h5ad[cluster_h5ad.obs['condition']==f'{origin}, {treatment}']
adata_NT = adata_NT[:, adata_NT.var_names.isin(genes)]
expression_data = pd.DataFrame(adata_NT.X.toarray(), index=adata_NT.obs_names, columns=adata_NT.var_names)
expression_data_sample = expression_data.sample(n=size_subset, random_state=2)
expression_data_sample = expression_data

zero_columns = list((expression_data_sample == 0).all())
non_zero_columns = list((expression_data_sample != 0).any())
zero_genes = expression_data_sample.iloc[:,zero_columns].columns
filtered_data = expression_data_sample.iloc[:,non_zero_columns]
filtered_data.to_csv(os.path.join(path_results, f'cell_{origin}_{treatment}.csv'), sep=',', index=True, header=True)



#Corr matrix
if type_GRN=='corr':
    corr_expression_data = filtered_data.corr()
    ddf = dd.from_pandas(filtered_data, npartitions=nthreads)
    correlation_matrix = ddf.corr().compute()
    np.savetxt(os.path.join(path_results, f'corr_{origin}_{treatment}.csv'), correlation_matrix, delimiter=',')

#  GENIE3
if type_GRN=='GRN':
    grn_matrix = GENIE3(filtered_data.values, nthreads=nthreads,ntrees=n_trees)
    df_grn = pd.DataFrame(grn_matrix,index=filtered_data.columns, columns=filtered_data.columns)
    df_grn.to_csv(os.path.join(path_results, f'GRN_{origin}_{treatment}_{size_subset}.csv'), sep=',', index=True, header=True)

