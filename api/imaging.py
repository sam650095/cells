from django.conf import settings

import os
import matplotlib.pyplot as plt
import scanpy as sc
import matplotlib
matplotlib.use('Agg')

# preprocessing
def produce_and_save_img(adata, saveplace, name):
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)
    save_dir = os.path.join(settings.MEDIA_ROOT, 'tempimage', saveplace)
    os.makedirs(save_dir, exist_ok=True)
    save_name = f'{adata.uns["prefix"]}_{name}.png'
    save_path = os.path.join(save_dir, save_name)
    plt.savefig(save_path)
    plt.close()
    return adata, save_name

def replaceimage():
    origin_dir = os.path.join(settings.MEDIA_ROOT, 'tempimage', 'origin')
    preview_dir = os.path.join(settings.MEDIA_ROOT, 'tempimage', 'preview')
    replaced = []
    for preview_file in os.listdir(preview_dir):
        prefix = preview_file.split('_')[0]
        preview_path = os.path.join(preview_dir, preview_file)
        
        for origin_file in os.listdir(origin_dir):
            if origin_file.startswith(prefix + '_'):
                origin_path = os.path.join(origin_dir, origin_file)
                os.replace(preview_path, origin_path)
                replaced.append(f"Replaced {origin_file}")
                break
    image_files = sorted(os.listdir(origin_dir))
    image_names = [file for file in image_files if file.lower().endswith(('.png', '.jpg', '.jpeg'))]
    return image_names

# clustering
def plot_clustering_pca(adata):
    pca_dir = os.path.join(settings.MEDIA_ROOT, 'pca_img')
    os.makedirs(pca_dir, exist_ok=True)

    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.regress_out(adata, ['total_counts'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, show=False)
    pca_plot_path = os.path.join(pca_dir, 'clustering_pca.jpg')
    plt.savefig(pca_plot_path, bbox_inches='tight')
    plt.close() 
    return 'clustering_pca.jpg'

def clustering_result(adata, npcs):
    if 'leiden_R' not in adata.obs.columns:
        adata.obs['leiden_R'] = adata.obs['leiden']
    
    # 創建儲存目錄
    save_dir = os.path.join(settings.MEDIA_ROOT, 'cluster_result')
    os.makedirs(save_dir, exist_ok=True)
    
    # Summary表格
    # summary_df_cluster = summary_cluster(adata)
    # dataframe_to_image(summary_df_cluster, filename='clustering_summary.pdf')
    
    # Umap
    sc.pl.umap(adata, color=['leiden_R','Sample'], show=False)
    plt.savefig(os.path.join(save_dir, 'clustering_leidens.png'), bbox_inches='tight', dpi=300)
    plt.close()
    
    # Ranking
    sc.tl.rank_genes_groups(adata, 'leiden_R')
    sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=False)
    plt.savefig(os.path.join(save_dir, 'clustering_ranking.png'), bbox_inches='tight', dpi=300)
    plt.close()
    
    # Heatmap
    all_markers = adata.var_names.tolist()
    sc.tl.dendrogram(adata, groupby='leiden_R', n_pcs=npcs)
    sc.pl.matrixplot(adata, all_markers, groupby='leiden_R', dendrogram=True, show=False, 
                     use_raw=False, cmap="vlag", standard_scale=None, title='leiden')
    plt.savefig(os.path.join(save_dir, 'clustering_heatmap.png'), bbox_inches='tight', dpi=300)
    plt.close()
    
    return adata