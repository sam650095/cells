from django.conf import settings

import pandas as pd
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
    
    save_dir = os.path.join(settings.MEDIA_ROOT, 'cluster_result')
    os.makedirs(save_dir, exist_ok=True)
    
    # Summary表格
    summary_df_cluster = summary_cluster(adata)
    dataframe_to_image(summary_df_cluster, os.path.join(save_dir, 'clustering_summary.png'))
    
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
    available_cmaps = plt.colormaps()
    cmap_choice = 'vlag' if 'vlag' in available_cmaps else 'coolwarm'
    # 判斷是否有vlag

    all_markers = adata.var_names.tolist()
    sc.tl.dendrogram(adata, groupby='leiden_R', n_pcs=npcs)
    sc.pl.matrixplot(adata, all_markers, groupby='leiden_R', dendrogram=True, show=False, 
                     use_raw=False, cmap=cmap_choice, standard_scale=None, title='leiden')
    plt.savefig(os.path.join(save_dir, 'clustering_heatmap.png'), bbox_inches='tight', dpi=300)
    plt.close()
    
    return adata
def dataframe_to_image(df, filename):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.axis('off')
    table = pd.plotting.table(ax, df, loc='center', cellLoc='center') 
    table.set_fontsize(60)
    table.scale(3, 3) 
    plt.savefig(filename, bbox_inches='tight')
    plt.close()
def summary_cluster(adata):
    if adata.uns.get('is_merged', False): 
        n_obs, n_vars = adata.shape
        print(f"Merged Data: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
        
        for sample in adata.obs['Sample'].cat.categories:
            sample_adata = adata[adata.obs['Sample'] == sample]
            n_obs, n_vars = sample_adata.shape
            print(f"{sample}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
        
        if 'leiden_R' not in adata.obs.columns:
            adata.obs['leiden_R'] = adata.obs['leiden']
        
        summary_df = pd.DataFrame(columns=['Leidens', 'Count', 'Merged_Proportion'])
        counts = adata.obs['leiden_R'].value_counts()
        proportions = (counts / counts.sum()).round(2)
        summary_df['Leidens'] = counts.index
        summary_df['Count'] = counts.values
        summary_df['Merged_Proportion'] = proportions.values

        for sample in adata.obs['Sample'].cat.categories:
            sample_adata = adata[adata.obs['Sample'] == sample]
            sample_counts = sample_adata.obs['leiden_R'].value_counts()
            sample_proportions = (sample_counts / sample_counts.sum()).round(2) 
            summary_df[f'{sample}_Proportion'] = sample_proportions.reindex(summary_df['Leidens']).fillna(0.0).values
        print(summary_df)
        
    else: 
        n_obs, n_vars = adata.shape 
        sample_name = adata.obs['Sample'].unique()[0]
        print(f"{sample_name}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
        
        if 'leiden_R' not in adata.obs.columns:
            adata.obs['leiden_R'] = adata.obs['leiden']
        
        summary_df = pd.DataFrame(columns=['Leidens', 'Count', f'{sample_name}_Proportion'])
        counts = adata.obs['leiden_R'].value_counts()
        proportions = (counts / counts.sum()).round(2)
        summary_df['Leidens'] = counts.index
        summary_df['Count'] = counts.values
        summary_df[f'{sample_name}_Proportion'] = proportions.values
        print(summary_df)
    
    return summary_df