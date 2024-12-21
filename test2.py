# Import necessary packages
import os
import copy
import anndata as ad
import pandas as pd
import numpy as np
import scimap as sm
import scanpy as sc
import matplotlib.pyplot as plt
import bbknn
# Set the working directory
np.random.seed(42)

# Import data
D0 = pd.read_csv("data/D0_batch_2_mesmer.seg.whole.cell.mask_raw_signal_value.csv", header=0, index_col=0)
D10 = pd.read_csv("data/D10_batch_2_mesmer.seg.whole.cell.mask_raw_signal_value.csv", header=0, index_col=0)
D20 = pd.read_csv("data/D20_batch_2_mesmer.seg.whole.cell.mask_raw_signal_value.csv", header=0, index_col=0)
metadata_D0 = pd.read_csv("data/D0_batch_2_mesmer.seg.whole.cell.mask_metadata.csv", header=0, index_col=0)
metadata_D10 = pd.read_csv("data/D10_batch_2_mesmer.seg.whole.cell.mask_metadata.csv", header=0, index_col=0)
metadata_D20 = pd.read_csv("data/D20_batch_2_mesmer.seg.whole.cell.mask_metadata.csv", header=0, index_col=0)

# Create AnnData Object
adata_D0 = sc.AnnData(X=D0, obs=metadata_D0)
adata_D10 = sc.AnnData(X=D10, obs=metadata_D10)
adata_D20 = sc.AnnData(X=D20, obs=metadata_D20)
adata_D0 = adata_D0[adata_D0.obs.total_counts < 6000, :]
adata_D10 = adata_D10[adata_D10.obs.total_counts < 15000, :]
adata_D20 = adata_D20[adata_D20.obs.total_counts < 75000, :]
sc.pp.normalize_total(adata_D0, target_sum=1e4) 
sc.pp.normalize_total(adata_D10, target_sum=1e4)
sc.pp.normalize_total(adata_D20, target_sum=1e4)
sc.pp.log1p(adata_D0) 
sc.pp.log1p(adata_D10)
sc.pp.log1p(adata_D20)
adata_D20.X[0,:].sum()
adata_merged = ad.concat([adata_D0, adata_D10, adata_D20], join='outer')
adata_merged.obs_names_make_unique()
adata_merged.raw = adata_merged
sc.pp.highly_variable_genes(adata_merged, n_top_genes=2000) # 預設
sc.pp.regress_out(adata_merged, ['total_counts']) 
sc.pp.scale(adata_merged, max_value=10)
sc.tl.pca(adata_merged, svd_solver='arpack')
sc.pl.pca(adata_merged, color=['PD1','pH2AX'])
sc.pl.pca_variance_ratio(adata_merged, log=True) # PCs to the total variance in the data
# This gives us information about how many PCs we should consider in order to compute 
# the neighborhood relations of cells, e.g. clustering function sc.tl.louvain() or tSNE sc.tl.tsne().

variance_ratio = adata_merged.uns['pca']['variance_ratio']
print('variance_ratio:', variance_ratio)
cumulative_variance_ratio = np.cumsum(variance_ratio)
print('cumulative_variance_ratio:', cumulative_variance_ratio)
n_pcs = np.where(cumulative_variance_ratio >= 0.95)[0][0] + 1
print('Selected number of PCs:', n_pcs)
adata_merged.write('merged_before_harmony.h5ad')
sc.external.pp.bbknn(adata_merged, batch_key='Sample')
sc.tl.umap(adata_merged, random_state=42)
sc.tl.leiden(adata_merged, resolution=0.3, random_state=42)  
sc.pl.umap(adata_merged, color=['leiden','Sample'])