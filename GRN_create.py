import scanpy as sc
import pandas as pd
import os
from GENIE3 import GENIE3
import sys
from scipy.spatial.distance import pdist, squareform
import scipy.sparse as sp

#path_h5ad = '/Users/ieo6943/Documents/Guido/Albi/clustered.h5ad/MDA_chemo/data/SCT/clustered.h5ad'
#path_genes = '/Users/ieo6943/Downloads/genes.csv'
#path_results = '/Users/ieo6943/Documents/Guido/Albi/data_my_GRN/'

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
#treatment = 'untreated'
#nthreads = 8
#size_subset = 'all_cell'
#n_trees = 500
genes_sub = 'all_genes'

path_results_size = os.path.join(path_results, str(size_subset))
path_results_subset = os.path.join(path_results_size, genes_sub)
os.makedirs(path_results_subset, exist_ok=True)

cluster_h5ad = sc.read_h5ad(path_h5ad)
sc.pp.normalize_total(cluster_h5ad, target_sum=1e4)
sc.pp.log1p(cluster_h5ad)
sc.pp.scale(cluster_h5ad)
cluster_h5ad.X = sp.csr_matrix(cluster_h5ad.X)
cluster_h5ad.layers["norm"] = cluster_h5ad.X.copy()
#aggiungere riga per salvare

genes = list(pd.read_csv(path_genes, index_col=0).index)

#aggiungere riga in caso di subset o all genes
all_data = pd.DataFrame(cluster_h5ad.X.toarray(), index=cluster_h5ad.obs_names, columns=cluster_h5ad.var_names)
top_genes = list(all_data.sum(axis=0).sort_values(ascending=False).head(100).index)
genes_of_interest = ['MT1X', 'LOX', 'APOE', 'SOD2', 'ENO2', 'PGK1', 'NNMT', 'ARF5','SERPINF1', 'PAEP', 'LCN1', 'CST1', 'CST4', 'FABP5', 'TMEM158', 'PLAUR','LDHA', 'TAGLN', 'FXYD5', 'EIF4EBP1']
#genes_subset = top_genes + genes_of_interest
genes_subset = genes

adata_NT = cluster_h5ad[cluster_h5ad.obs['condition']==f'{origin}, {treatment}']
adata_NT = adata_NT[:, adata_NT.var_names.isin(genes_subset)]
expression_data = pd.DataFrame(adata_NT.X.toarray(), index=adata_NT.obs_names, columns=adata_NT.var_names)
#expression_data_sample = expression_data.sample(n=size_subset, random_state=2)
expression_data_sample = expression_data

zero_columns = list((expression_data_sample == 0).all())
non_zero_columns = list((expression_data_sample != 0).any())
#zero_genes = expression_data_sample.iloc[:,zero_columns].columns
filtered_data = expression_data_sample.iloc[:,non_zero_columns]
filtered_data.to_csv(os.path.join(path_results_subset, f'cell_subset_{origin}_{treatment}.csv'), sep=',', index=True, header=True)


#Corr matrix
if type_GRN=='corr':
    data = filtered_data.values
    correlation_matrix = 1 - squareform(pdist(data.T, metric='correlation'))
    correlation_df = pd.DataFrame(correlation_matrix, index=filtered_data.columns, columns=filtered_data.columns)
    correlation_df.to_csv(os.path.join(path_results_subset, f'corr_{origin}_{treatment}.csv'), sep=',', index=True, header=True)
    

#  GENIE3
if type_GRN=='GENIE':
    grn_matrix = GENIE3(filtered_data.values, nthreads=nthreads,ntrees=n_trees)
    df_grn = pd.DataFrame(grn_matrix,index=filtered_data.columns, columns=filtered_data.columns)
    df_grn.to_csv(os.path.join(path_results_subset, f'GENIE_{origin}_{treatment}.csv'), sep=',', index=True, header=True)









