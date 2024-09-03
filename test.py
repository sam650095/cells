import matplotlib.pyplot as plt
import scimap as sm 
import matplotlib
import anndata as ad
import time
start = time.time()
# 設置
chosen_column = 'phenotype'
matplotlib.use('Agg')
adata = ad.read_h5ad('adata_phenotyping.h5ad')
plt.rcParams['figure.figsize'] = [15, 10]

# 清除任何現有的圖形
plt.clf()

print('plot start')
# 創建 Voronoi 圖
sm.pl.voronoi(adata, 
              color_by=chosen_column, 
              voronoi_edge_color='black', 
              voronoi_line_width=0.3, 
              voronoi_alpha=0.8, 
              size_max=5000, 
              overlay_points=None, 
              plot_legend=True, 
              legend_size=6)

# 保存圖形
plt.savefig(f'interactions_voronoi_{chosen_column}.png')
plt.close()  # 關閉圖形以釋放內存
end = time.time()

# 檢查數據
print(f"Data shape: {adata.shape}")
print(f"Available columns: {adata.obs.columns}")
print(f"Unique values in {chosen_column}: {adata.obs[chosen_column].unique()}")
print(format(end - start)+'秒')
# # 如果還是不行，嘗試使用 scatter plot
# plt.figure(figsize=(15, 10))
# sc = plt.scatter(adata.obsm['spatial'][:, 0], 
#                  adata.obsm['spatial'][:, 1], 
#                  c=adata.obs[chosen_column].astype('category').cat.codes, 
#                  s=1, 
#                  alpha=0.8)
# plt.colorbar(sc)
# plt.title(f'Scatter plot of {chosen_column}')
# plt.savefig(f'scatter_{chosen_column}.png')
# plt.close()