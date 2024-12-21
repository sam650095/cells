import scanpy as sc
import bbknn
import logging

logging.basicConfig(level=logging.INFO)

# 加載示例數據集
adata = sc.datasets.pbmc3k()  # 這是一個 PBMC (外周血單核細胞) 的小型數據集

# 篩選和標準化數據
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 計算高變基因
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var['highly_variable']].copy()  # 顯式轉換為副本

# 標準化高變基因
sc.pp.scale(adata, max_value=10)

# PCA 降維
sc.tl.pca(adata, svd_solver='arpack')

# 使用 bbknn 替代標準最近鄰搜尋（適用於批次效應校正）
adata.obs['batch'] = ['batch1'] * (adata.n_obs // 2) + ['batch2'] * (adata.n_obs // 2)
print("Batch labels set successfully.")
print(f"adata shape: {adata.shape}")
# 確保 BBKNN 正常運行
try:
    bbknn.bbknn(adata, batch_key='batch')
    print("BBKNN completed successfully.")
except Exception as e:
    print("Error during BBKNN execution:", e)

# UMAP 可視化
sc.tl.umap(adata)
sc.pl.umap(adata, color='batch')
