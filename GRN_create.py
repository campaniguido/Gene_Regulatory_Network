import scanpy as sc
import numpy as np
import pandas as pd
import os
from GENIE3 import GENIE3
import umap
import matplotlib.pyplot as plt
import sys
import psutil
#import dask.dataframe as dd

###############################################
#Managing resources
###############################################
#def log_resources(stage):
#    memory_used = psutil.Process().memory_info().rss / (1024 ** 3)
#    cpu_percent = psutil.cpu_percent(interval=1)
#    with open("resource_usage.log", "a") as log_file:
#        log_file.write(f"{stage}: Memoria usata: {memory_used:.2f} GB, "
#                       f"CPU: {cpu_percent}%\n")
#################################################

#log_resources("begining Script")

#path_h5ad = '/Users/ieo6943/Documents/Guido/Albi/clustered.h5ad/MDA_chemo/data/SCT/clustered.h5ad'
#
#path_genes = '/Users/ieo6943/Downloads/genes.csv'

path_h5ad = '/hpcnfs/data/PGP/acossa/archive_Cellula/MDA_chemo/data/SCT/clustered.h5ad'
path_results = '/hpcnfs/scratch/PGP/gcampani/ALBI/GRN/'
path_genes = '/hpcnfs/scratch/PGP/acossa/breast_albi/MDA/single_cell/scripts/GRNs/genes.csv'

origin = sys.argv[1]
treatment = sys.argv[2]
nthreads = int(sys.argv[3])
type_GRN = sys.argv[4]

origin = 'PT'
treatment = 'untreated'
nthreads = int(20)
size_subset = 1000

#path_results = f'/Users/ieo6943/Documents/Guido/Albi/data_my_GRN/{size_subset}'
os.makedirs(path_results, exist_ok=True)

#log_resources("before ANNdata")
cluster_h5ad = sc.read_h5ad(path_h5ad)
#log_resources("after ANNdata")

genes = list(pd.read_csv(path_genes, index_col=0).index)

adata_NT = cluster_h5ad[cluster_h5ad.obs['condition']==f'{origin}, {treatment}']
adata_NT = adata_NT[:, adata_NT.var_names.isin(genes)]
expression_data = pd.DataFrame(adata_NT.X.toarray(), index=adata_NT.obs_names, columns=adata_NT.var_names)
expression_data_sample = expression_data.sample(n=size_subset, random_state=2)

#Corr matrix
#if type_GRN=='corr':
#    #corr_expression_data = expression_data.corr()
#    # Converti il DataFrame Pandas in un DataFrame Dask
#    ddf = dd.from_pandas(expression_data, npartitions=nthreads)
#    correlation_matrix = ddf.corr().compute()
#    np.savetxt(os.path.join(path_results, f'corr_{origin}_{treatment}.csv'), correlation_matrix, delimiter=',')

#  GENIE3 

zero_columns = list((expression_data_sample == 0).all())
non_zero_columns = list((expression_data_sample != 0).any())

zero_genes = expression_data_sample.iloc[:,zero_columns].columns
filtered_data = expression_data_sample.iloc[:,non_zero_columns]
filtered_data.to_csv(os.path.join(path_results, f'cell_{origin}_{treatment}.csv'), sep=',', index=True, header=True)

#log_resources("before Network")
if type_GRN=='GRN':
    grn_matrix = GENIE3(filtered_data.values, nthreads=nthreads,ntrees=500)
    df_grn = pd.DataFrame(grn_matrix,index=filtered_data.columns, columns=filtered_data.columns)
    df_grn.to_csv(os.path.join(path_results, f'GRN_{origin}_{treatment}_1000_sample.csv'), sep=',', index=True, header=True)

#log_resources("After Network")




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