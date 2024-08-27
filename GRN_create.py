import scanpy as sc
import numpy as np
import pandas as pd
import os
from GENIE3 import GENIE3
import umap
import matplotlib.pyplot as plt
import sys
#import dask.dataframe as dd


#path_h5ad = '/Users/ieo6943/Documents/Guido/Albi/clustered.h5ad/MDA_chemo/data/SCT/clustered.h5ad'
#path_results = '/Users/ieo6943/Downloads/'
#path_genes = '/Users/ieo6943/Downloads/genes.csv'
#origin = 'PT'
#treatment = 'treated'
#nthreads = 5

path_h5ad = '/hpcnfs/data/PGP/acossa/archive_Cellula/MDA_chemo/data/SCT/clustered.h5ad'
path_results = '/hpcnfs/scratch/PGP/gcampani/ALBI/GRN/'
path_genes = '/hpcnfs/scratch/PGP/acossa/breast_albi/MDA/single_cell/scripts/GRNs/genes.csv'

origin = sys.argv[1]
treatment = sys.argv[2]
nthreads = sys.argv[3]
type_GRN = sys.argv[4]

os.makedirs(path_results, exist_ok=True)

cluster_h5ad = sc.read_h5ad(path_h5ad)
genes = list(pd.read_csv(path_genes, index_col=0).index)

adata_NT = cluster_h5ad[cluster_h5ad.obs['condition']==f'{origin}, {treatment}']
adata_NT = cluster_h5ad[:, cluster_h5ad.var_names.isin(genes)]
expression_data = pd.DataFrame(adata_NT.X.toarray(), index=adata_NT.obs_names, columns=adata_NT.var_names)

#Corr matrix
#if type_GRN=='corr':
#    #corr_expression_data = expression_data.corr()
#    # Converti il DataFrame Pandas in un DataFrame Dask
#    ddf = dd.from_pandas(expression_data, npartitions=nthreads)
#    correlation_matrix = ddf.corr().compute()
#    np.savetxt(os.path.join(path_results, f'corr_{origin}_{treatment}.csv'), correlation_matrix, delimiter=',')

#  GENIE3 
if type_GRN=='GRN':
    grn_matrix = GENIE3(expression_data.values, nthreads=nthreads)
    np.savetxt(os.path.join(path_results, f'GRN_{origin}_{treatment}.csv'), grn_matrix, delimiter=',')




#UMAP
#conditions = list(adata_NT.obs['condition'].unique())
#colors = ['red', 'blue', 'green', 'orange', 'purple']
#condition_color_dict = dict(zip(conditions, colors))
#
#colors = adata_NT.obs['condition'].map(condition_color_dict)
#
#
#
## Esegui UMAP
#reducer = umap.UMAP(n_components=2)
#embedding = reducer.fit_transform(expression_data)
#
## Visualizza i risultati
#plt.figure(figsize=(10, 8))
#plt.scatter(embedding[:, 0], embedding[:, 1], c=colors, s=0.1)
#plt.title("UMAP Projection")
#plt.xlabel("UMAP1")
#plt.ylabel("UMAP2")
#plt.show()
#