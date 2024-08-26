import scanpy as sc
import numpy as np
import pandas as pd
import os
from GENIE3 import GENIE3
import sys


#path_h5ad = '/Users/ieo6943/Documents/Guido/Albi/clustered.h5ad/MDA_chemo/data/SCT/clustered.h5ad'
#path_results = '/Users/ieo6943/Downloads/'
#path_genes = '/Users/ieo6943/Downloads/genes.csv'

path_h5ad = '/hpcnfs/data/PGP/acossa/archive_Cellula/MDA_chemo/data/SCT/clustered.h5ad'
path_results = '/hpcnfs/scratch/PGP/gcampani/ALBI/GRN/'
path_genes = '/hpcnfs/scratch/PGP/acossa/breast_albi/MDA/single_cell/scripts/GRNs/genes.csv'

origin = sys.argv[1]
treatment = sys.argv[2]

os.makedirs(path_results, exist_ok=True)

cluster_h5ad = sc.read_h5ad(path_h5ad)
genes = list(pd.read_csv(path_genes, index_col=0).index)

adata_NT = cluster_h5ad[cluster_h5ad.obs['condition']==f'{origin}, {treatment}']
adata_NT = adata_NT[:, adata_NT.var_names.isin(genes)]

expression_data = pd.DataFrame(adata_NT.X.toarray(), index=adata_NT.obs_names, columns=adata_NT.var_names)

# Esegui GENIE3 per inferire la GRN
grn_matrix = GENIE3(expression_data.values)
np.savetxt(os.path.join(path_results, f'GRN_{origin}_{treatment}.csv'), grn_matrix, delimiter=',')