from django.conf import settings

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import scanpy as sc
import scimap as sm 
import time
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import MiniBatchKMeans
import seaborn as sns
from .saved import *
import matplotlib
matplotlib.use('Agg')

# preprocessing
def produce_and_save_img(adata, saveplace, name):
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)
    save_dir = os.path.join(settings.MEDIA_ROOT, 'qualitycontrol', saveplace)
    os.makedirs(save_dir, exist_ok=True)
    save_name = f'{adata.uns["prefix"]}_{name}.png'
    save_path = os.path.join(save_dir, save_name)
    plt.savefig(save_path)
    plt.close()
    return adata, save_name

def replaceimage():
    origin_dir = os.path.join(settings.MEDIA_ROOT, 'qualitycontrol', 'origin')
    preview_dir = os.path.join(settings.MEDIA_ROOT, 'qualitycontrol', 'preview')
    replaced = []
    for preview_file in os.listdir(preview_dir):
        prefix = preview_file.split('_')[0]
        preview_path = os.path.join(preview_dir, preview_file)
        
        for origin_file in os.listdir(origin_dir):
            if origin_file.startswith(prefix + '_'):
                origin_path = os.path.join(origin_dir, origin_file)
                shutil.copy(preview_path, origin_path)
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
    group_sizes = adata.obs['leiden_R'].value_counts() 
    valid_groups = group_sizes[group_sizes >= 2].index.tolist() 

    if len(valid_groups) > 1:
        sc.tl.rank_genes_groups(adata, 'leiden_R')
        n_cols = 3  
        n_rows = -(-len(valid_groups) // n_cols) 
        fig_width = 6 * n_cols 
        fig_height = 5 * n_rows 
        
        plt.figure(figsize=(fig_width, fig_height))
        
        sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=False,
                                ncols=n_cols,
                                fontsize=12,
                            )
        plt.tight_layout() 
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
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.axis('off')
    table = pd.plotting.table(ax, df, loc='center', cellLoc='center') 
    
    for key, cell in table.get_celld().items():
        cell.set_fontsize(36)
    
    table.scale(3, 3) 
    plt.savefig(filename, bbox_inches='tight')
    plt.close()
def summary_cluster(adata):
    if adata.uns.get('is_merged', False): 
        n_obs, n_vars = adata.shape
        
        for sample in adata.obs['Sample'].cat.categories:
            sample_adata = adata[adata.obs['Sample'] == sample]
            n_obs, n_vars = sample_adata.shape
        
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
        # print(summary_df)
        
    else: 
        n_obs, n_vars = adata.shape 
        sample_name = adata.obs['Sample'].unique()[0]
        # print(f"{sample_name}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
        
        if 'leiden_R' not in adata.obs.columns:
            adata.obs['leiden_R'] = adata.obs['leiden']
        
        summary_df = pd.DataFrame(columns=['Leidens', 'Count', f'{sample_name}_Proportion'])
        counts = adata.obs['leiden_R'].value_counts()
        proportions = (counts / counts.sum()).round(2)
        summary_df['Leidens'] = counts.index
        summary_df['Count'] = counts.values
        summary_df[f'{sample_name}_Proportion'] = proportions.values
        # print(summary_df)
    
    return summary_df
def adding_clusteringumap(adata, chosen_markers=None):
    num_cols = 3
    
    if chosen_markers:
        plot_data = chosen_markers
        color_key = lambda x: x
        title_key = lambda x: x
        filename = 'clustering_markers.png'
    else:
        plot_data = adata.obs['Sample'].cat.categories
        color_key = lambda x: 'leiden_R'
        title_key = lambda x: x
        filename = 'clustering_leidens_bysample.png'
    
    num_plots = len(plot_data)
    num_rows = int(np.ceil(num_plots / num_cols))
    
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 4.5 * num_rows))
    axs = axs.flatten() if num_plots > 1 else [axs]
    
    for i, item in enumerate(plot_data):
        if chosen_markers:
            sc.pl.umap(adata, color=color_key(item), title=title_key(item), ax=axs[i], show=False)
        else:
            sample_adata = adata[adata.obs['Sample'] == item]
            sc.pl.umap(sample_adata, color=color_key(item), title=title_key(item), ax=axs[i], show=False)
    
    for j in range(num_plots, len(axs)):
        fig.delaxes(axs[j])
    
    plt.tight_layout()
    
    if settings and hasattr(settings, 'MEDIA_ROOT'):
        save_dir = os.path.join(settings.MEDIA_ROOT, 'umap_cluster')
    else:
        save_dir = 'umap_cluster'
    os.makedirs(save_dir, exist_ok=True)
    
    save_path = os.path.join(save_dir, filename)
    plt.savefig(save_path)
    plt.close(fig)
    
    return save_path  

# phenotyping
def summary_phenotype(adata, chosen_adata):
    n_obs, n_vars = adata.shape
        
    if len(adata.obs['Sample'].cat.categories) > 1: 
        for sample in adata.obs['Sample'].cat.categories:
            sample_adata = adata[adata.obs['Sample'] == sample]
            n_obs, n_vars = sample_adata.shape
            # print(f"{sample}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")

        summary_df = pd.DataFrame(columns=['Phenotypes', 'Count', 'Merged_Proportion'])
        counts = adata.obs['phenotype'].value_counts()
        proportions = (counts / counts.sum()).round(2)
        summary_df['Phenotypes'] = counts.index
        summary_df['Count'] = counts.values
        summary_df['Merged_Proportion'] = proportions.values

        for sample in adata.obs['Sample'].cat.categories:
            sample_adata = adata[adata.obs['Sample'] == sample]
            sample_counts = sample_adata.obs['phenotype'].value_counts()
            sample_proportions = (sample_counts / sample_counts.sum()).round(2)
            summary_df[f'{sample}_Proportion'] = sample_proportions.reindex(summary_df['Phenotypes']).fillna(0.0).values
        # print(summary_df)
        
    else: 
        sample_name = adata.obs['Sample'].unique()[0]
        summary_df = pd.DataFrame(columns=['Phenotypes', 'Count', f'{sample_name}_Proportion'])
        counts = adata.obs['phenotype'].value_counts()
        proportions = (counts / counts.sum()).round(2)
        summary_df['Phenotypes'] = counts.index
        summary_df['Count'] = counts.values
        summary_df[f'{sample_name}_Proportion'] = proportions.values
        # print(summary_df)
    
    return summary_df, f"{chosen_adata}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"
def phenotype_result(adata, chosen_adata, n_pcs):
    save_dir = os.path.join(settings.MEDIA_ROOT, 'phenotype_result')
    os.makedirs(save_dir, exist_ok=True)
    # Summary表格
    summary_df_phenotype, result_text = summary_phenotype(adata,chosen_adata)
    dataframe_to_image(summary_df_phenotype, os.path.join(save_dir, 'phenotyping_summary.png'))

    # Umap
    fig, axs = plt.subplots(2, 1, figsize=(10, 20)) 
    sc.pl.umap(adata, color='phenotype', ax=axs[0], show=False)
    sc.pl.umap(adata, color='Sample', ax=axs[1], show=False)
    plt.savefig(os.path.join(save_dir, 'phenotyping_leidens.png'), bbox_inches='tight') 
    
    # Ranking
    group_sizes = adata.obs['phenotype'].value_counts()
    valid_groups = group_sizes[group_sizes >= 2].index.tolist()
    
    if len(valid_groups) > 1:
        sc.tl.rank_genes_groups(adata, 'phenotype', groups=valid_groups)
        
        n_cols = 3  
        n_rows = -(-len(valid_groups) // n_cols) 
        fig_width = 6 * n_cols 
        fig_height = 5 * n_rows 
        
        plt.figure(figsize=(fig_width, fig_height))
        
        sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=False,
                                ncols=n_cols,
                                fontsize=12,
                            )
        plt.tight_layout() 
        plt.savefig(os.path.join(save_dir, 'phenotyping_ranking.png'), bbox_inches='tight', dpi=300)
        plt.close()
    
    # Heatmap
    all_markers = adata.var_names.tolist()
    sc.tl.dendrogram(adata, groupby='phenotype', n_pcs=n_pcs)
    sc.pl.matrixplot(adata, all_markers, groupby='phenotype', dendrogram=True, show=False, 
                     use_raw=False, cmap="vlag", standard_scale=None, title='phenotype')
    plt.savefig(os.path.join(save_dir, 'phenotyping_heatmap.png'), bbox_inches='tight')
    
    return adata, result_text
def add_phenotypes_bysample(adata):
    num_samples = len(adata.obs['Sample'].cat.categories)
    fig, axs = plt.subplots(num_samples, 1, figsize=(10, 5 * num_samples))
    if num_samples == 1:
        axs = [axs]
    for i, sample in enumerate(adata.obs['Sample'].cat.categories):
        sample_adata = adata[adata.obs['Sample'] == sample]
        sc.pl.umap(sample_adata, color='phenotype', title=f'{sample}', ax=axs[i], show=False)
    plt.tight_layout()

    if settings and hasattr(settings, 'MEDIA_ROOT'):
        save_dir = os.path.join(settings.MEDIA_ROOT, 'umap_phenotyping')
    else:
        save_dir = 'umap_phenotyping'
    os.makedirs(save_dir, exist_ok=True)
    
    save_path = os.path.join(save_dir, 'phenotypes_bysample.png')
    plt.savefig(save_path)
    plt.close(fig)
    
    return save_path  
def add_phenotypes_markers(adata, chosen_markers):
    num_cols = 3
    num_markers = len(chosen_markers)
    num_rows = np.ceil(num_markers / num_cols).astype(int)
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 4 * num_rows))
    axs = axs.flatten()
    for i, marker in enumerate(chosen_markers):
        sc.pl.umap(adata, color=marker, title=f'{marker}', ax=axs[i], show=False)  
    for j in range(i + 1, len(axs)):
        fig.delaxes(axs[j])
    plt.tight_layout()
    if settings and hasattr(settings, 'MEDIA_ROOT'):
        save_dir = os.path.join(settings.MEDIA_ROOT, 'umap_phenotyping')
    else:
        save_dir = 'umap_phenotyping'
    os.makedirs(save_dir, exist_ok=True)
    
    save_path = os.path.join(save_dir, 'phenotypes_markers.png')
    plt.savefig(save_path)
    plt.close(fig)
    
    return save_path

# spatial analysis
def distances_heatmap(chosen_column):
    adata = read_h5ad_file('adata_spatial_analysis.h5ad')
    save_dir = os.path.join(settings.MEDIA_ROOT, 'spatial_result')
    os.makedirs(save_dir, exist_ok=True)

    if(chosen_column == "leiden_R"):
        sm.pl.spatial_distance(adata, spatial_distance='spatial_distance_leiden_R', phenotype='leiden_R', imageid='Sample', heatmap_summarize=False)
        plt.savefig(os.path.join(save_dir, 'distances_heatmap_leiden_R.png'))
        plt.close()
        return 'distances_heatmap_leiden_R.png'
    else:
        sm.pl.spatial_distance(adata, spatial_distance='spatial_distance_phenotype', phenotype='phenotype', imageid='Sample', heatmap_summarize=False)
        plt.savefig(os.path.join(save_dir,'distances_heatmap_phenotype.png'))
        plt.close()
        return 'distances_heatmap_phenotype.png'
            
def distances_numeric_plot(chosen_column, chosen_cluster): 
    adata = read_h5ad_file('adata_spatial_analysis.h5ad')
    save_dir = os.path.join(settings.MEDIA_ROOT, 'spatial_result')
    os.makedirs(save_dir, exist_ok=True)

    sm.pl.spatial_distance(adata, spatial_distance=f'spatial_distance_{chosen_column}', method='numeric', distance_from=chosen_cluster,
                           imageid='Sample', phenotype=chosen_column, log=True)
    plt.savefig(os.path.join(save_dir,f'distances_numeric_plot_{chosen_column}_{chosen_cluster}.png'))
    plt.close()
    return f'distances_numeric_plot_{chosen_column}_{chosen_cluster}.png'

def interactions_heatmap(chosen_column, chosen_method):
    adata = read_h5ad_file('adata_spatial_analysis.h5ad')
    save_dir = os.path.join(settings.MEDIA_ROOT, 'spatial_result')
    os.makedirs(save_dir, exist_ok=True)
    sm.pl.spatial_interaction(adata, p_val=0.05, summarize_plot=True, row_cluster=True,
                              spatial_interaction=f'spatial_interaction_{chosen_column}_{chosen_method}')
    plt.savefig(os.path.join(save_dir,f'interactions_heatmap_{chosen_column}_{chosen_method}.png'))
    plt.close()
    return f'interactions_heatmap_{chosen_column}_{chosen_method}.png'
    
def interactions_voronoi(chosen_column):
    adata = read_h5ad_file('adata_spatial_analysis.h5ad')
    save_dir = os.path.join(settings.MEDIA_ROOT, 'spatial_result')
    os.makedirs(save_dir, exist_ok=True)
    plt.rcParams['figure.figsize'] = [15, 10]
    sm.pl.voronoi(adata, color_by=chosen_column, voronoi_edge_color = 'black', voronoi_line_width = 0.3, 
                  voronoi_alpha = 0.8, size_max=5000, overlay_points=None, plot_legend=True, legend_size=6)
    plt.savefig(os.path.join(save_dir, f'interactions_voronoi_{chosen_column}.png'))
    plt.close()
    return f'interactions_voronoi_{chosen_column}.png'


# neighbor
def get_windows(job, n_neighbors, tissue_group, exps, X, Y):
    start_time,idx,tissue_name,indices = job
    job_start = time.time()
    
    print ("Starting:", str(idx+1)+'/'+str(len(exps)),': ' + exps[idx])

    tissue = tissue_group.get_group(tissue_name)
    to_fit = tissue.loc[indices][[X,Y]].values

    fit = NearestNeighbors(n_neighbors=n_neighbors).fit(tissue[[X,Y]].values)
    m = fit.kneighbors(to_fit)
    m = m[0], m[1]
    
    # sort_neighbors
    args = m[0].argsort(axis = 1)
    add = np.arange(m[1].shape[0])*m[1].shape[1]
    sorted_indices = m[1].flatten()[args+add[:,None]]

    neighbors = tissue.index.values[sorted_indices]
   
    end_time = time.time()
   
    print ("Finishing:", str(idx+1)+"/"+str(len(exps)),": "+ exps[idx],end_time-job_start,end_time-start_time)
    return neighbors.astype(np.int32)

def preprocessing(adata, chosen_column, k):
    clustering_columns = [col for col in adata.obs.columns if col.startswith(("leiden_R", "phenotype"))]
    neigberhood_df = adata.obs[['X_centroid', 'Y_centroid', 'Sample'] + clustering_columns] 
    neigberhood_df.to_csv('neigberhood_data.csv')
    k = int(k)
    n_neighbors = k
    path = os.path.join(settings.MEDIA_ROOT, 'tempfile/neighbor')
    os.makedirs(path, exist_ok=True)
    path_to_data = os.path.join('data','neigberhood_data.csv')

    X = 'X_centroid'
    Y = 'Y_centroid'
    reg = 'Sample' 
    keep_cols = [X,Y,reg,chosen_column] 

    cells = pd.read_csv(path_to_data)
    cells = pd.concat([cells, pd.get_dummies(cells[chosen_column])], axis=1)
    sum_cols = cells[chosen_column].unique()
    values = cells[sum_cols].values 

    tissue_group = cells[[X, Y, reg]].groupby(reg) 
    exps = list(cells[reg].unique()) 
    tissue_chunks = [(time.time(), exps.index(t), t, a) for t, indices in tissue_group.groups.items() for a in
                     np.array_split(indices, 1)]

    tissues = [get_windows(job, n_neighbors, tissue_group, exps, X, Y) for job in tissue_chunks]

    out_dict = {}
    for neighbors, job in zip(tissues, tissue_chunks):
        chunk = np.arange(len(neighbors))  
        tissue_name = job[2] 
        indices = job[3] 
        window = values[neighbors[chunk, :k].flatten()].reshape(len(chunk), k, len(sum_cols)).sum(axis=1)
        out_dict[(tissue_name, k)] = (window.astype(np.float16), indices)

    windows = {}
    window = pd.concat([pd.DataFrame(out_dict[(exp,k)][0],index = out_dict[(exp,k)][1].astype(int),columns = sum_cols) for exp in exps],axis=0)
    window = window.loc[cells.index.values]
    window = pd.concat([cells[keep_cols], window], axis=1)
    windows[k] = window
    return windows, sum_cols, cells, values, reg

def perform_neighborhoods(adata, chosen_column, k, n_neighborhoods):
    windows, sum_cols, cells, values, reg = preprocessing(adata, chosen_column, k)
    neighborhood_name = "neighborhood" + str(k)
    k_centroids = {}  
    k = int(k)
    n_neighborhoods = int(n_neighborhoods)
    windows2 = windows[k] 
    km = MiniBatchKMeans(n_clusters=n_neighborhoods, random_state=0) 
    labelskm = km.fit_predict(windows2[sum_cols].values) 
    k_centroids[k] = km.cluster_centers_ 
    cells[neighborhood_name] = labelskm 
    cells[neighborhood_name] = cells[neighborhood_name].astype('category')

    # Heatmap
    cell_order = sum_cols
    niche_clusters = k_centroids[k]
    tissue_avgs = values.mean(axis=0)
    fc = np.log2(((niche_clusters + tissue_avgs) / (niche_clusters + tissue_avgs).sum(axis=1, keepdims=True)) / tissue_avgs)
    fc = pd.DataFrame(fc, columns=sum_cols)
    s = sns.clustermap(fc.loc[range(n_neighborhoods), cell_order], vmin=-3, vmax=3, cmap='bwr', row_cluster=False)
    
    save_dir = os.path.join(settings.MEDIA_ROOT, 'neighbor_result')
    os.makedirs(save_dir, exist_ok=True)
    
    plt.savefig(os.path.join(save_dir,"neighborhood_heatmap.png"))

    cells[neighborhood_name] = cells[neighborhood_name].astype('category')
    plt.rcParams['figure.figsize'] = [15, 10]
    sns.lmplot(data = cells, x = 'X_centroid', y='Y_centroid', hue = neighborhood_name, palette = 'bright', 
               height = 15, col = reg, col_wrap = 1, fit_reg = False)
    plt.savefig(os.path.join(save_dir,"neighborhood_lmplot.png"))