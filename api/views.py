import os
import pandas as pd
import scanpy as sc
import numpy as np
import shutil
import bbknn
import logging

from django.conf import settings
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from .serializers import FileUploadSerializer
from .imaging import *
from .saved import *
logger = logging.getLogger(__name__)
def sort_key(filename):
    return int(filename.split('_')[0][1:]) 
def clearmediafiles(temp):
    media_root = settings.MEDIA_ROOT
    temp_dir = os.path.join(media_root, temp)
    messages = []

    if os.path.exists(temp_dir):
        for item in os.listdir(temp_dir):
            item_path = os.path.join(temp_dir, item)
            try:
                if os.path.isfile(item_path) or os.path.islink(item_path):
                    os.unlink(item_path)
                elif os.path.isdir(item_path):
                    shutil.rmtree(item_path)
                messages.append(f'Successfully deleted {item_path}')
            except Exception as e:
                messages.append(f'Failed to delete {item_path}. Reason: {e}')
    else:
        messages.append(f'Directory {temp_dir} does not exist')

    return '\n'.join(messages)    
# Preprocessing     
class FileUploadView(APIView):
    def post(self, request):
        clearmediafiles('tempfile')
        files = request.FILES.getlist('files')
        file_serializers = []

        for file in files:
            serializer = FileUploadSerializer(data={'file': file})
            if serializer.is_valid():
                saved_file = serializer.save()
                file_serializers.append(saved_file)
            else:
                return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
        signal_value_files = sorted(os.listdir(os.path.join(settings.MEDIA_ROOT, 'tempfile', 'signal_value')), key=sort_key)
        metadata_files = sorted(os.listdir(os.path.join(settings.MEDIA_ROOT, 'tempfile', 'metadata')), key=sort_key)
        adata_objects, original_adata_objects, adata_results = self.create_adata(signal_value_files, metadata_files)
        save_data(adata_objects, 'adata_objects')
        save_data(original_adata_objects, 'original_adata_objects')
        return Response({"adata_results": adata_results}, status=status.HTTP_201_CREATED)
    
    def create_adata(self, signal_value_files, metadata_files):
        adata_objects = []
        original_adata_objects = []
        adata_results = []
        signal_value_dir = os.path.join(settings.MEDIA_ROOT, 'tempfile', 'signal_value')
        metadata_dir = os.path.join(settings.MEDIA_ROOT, 'tempfile', 'metadata')

        for signal_value_file in signal_value_files:
            prefix = signal_value_file.split('_')[0]
            matching_metadata_file = next((metadata_file for metadata_file in metadata_files if metadata_file.startswith(prefix)), None)

            if matching_metadata_file:
                data = pd.read_csv(os.path.join(signal_value_dir, signal_value_file), header=0, index_col=0)
                metadata = pd.read_csv(os.path.join(metadata_dir, matching_metadata_file), header=0, index_col=0)
                adata = sc.AnnData(X=data, obs=metadata)
                adata.raw = adata
                adata.uns['prefix'] = prefix
                n_obs, n_vars = adata.shape
                adata_results.append(f"{prefix}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
                adata_objects.append(adata)
                original_adata_objects.append(adata.copy())
        return adata_objects, original_adata_objects, adata_results
class QualityControlView(APIView):
    def post(self, request):
        clear_media = clearmediafiles('tempimage')
        adata_objects, qc_adata_objects, adata_results, save_image_names = self.post_reset_adata()
        save_data(adata_objects, 'adata_objects')
        save_data(qc_adata_objects, 'qc_adata_objects')
        save_data(qc_adata_objects, 'preview_adata_objects')
        return Response({'adata_results': adata_results,'save_image_names': save_image_names}, status=status.HTTP_201_CREATED)
    
    def post_reset_adata(self):
        adata_objects = load_data('original_adata_objects')
        qc_adata_objects = []
        adata_results = []
        save_image_names = []
        for adata in adata_objects:
            n_obs, n_vars = adata.shape
            adata_results.append(f"{adata.uns['prefix']}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
            sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
            adata, imgname = produce_and_save_img(adata, "origin", 'violin_plot')
            qc_adata_objects.append(adata)
            save_image_names.append(imgname)
        
        return adata_objects, qc_adata_objects, adata_results, save_image_names  
class PreviewView(APIView):
    def post(self, request):
        adata_objects = load_data('qc_adata_objects')
        # get data from frontend
        f_sampleSelect = request.data.get('sample')
        minGenes = int(request.data.get('minGenes'))
        filter_method = request.data.get('fmethod')
        lowerlimit = float(request.data.get('lowerlimit')) if request.data.get('lowerlimit') else None
        upperlimit = float(request.data.get('upperlimit')) if request.data.get('upperlimit') else None
        user_inputs = {}
        user_inputs['chosen_adata'] = self.choose_sample(adata_objects,f_sampleSelect)
        user_inputs['min_genes'] = minGenes
        user_inputs['outlier_lims'] = self.input_outliers(user_inputs['chosen_adata'], filter_method, lowerlimit, upperlimit)

        preview_adata, preview_adata_result, preview_image_name = self.preview_filter(user_inputs)
        # chose adata to preview
        save_data(preview_adata, 'preview_adata')
        return Response({'adata_result': preview_adata_result,'save_image_names': preview_image_name}, status=status.HTTP_201_CREATED)
    
    def choose_sample(self,adata_objects, sample_select):
        for adata in adata_objects:
            if adata.uns['prefix'] == sample_select:
                chosen_adata = adata
                return chosen_adata
    def input_outliers(self, chosen_adata, filtermethod, lowerLimit, upperLimit):   
        choice = filtermethod
        if choice == 'Manual':
            lower_lim = lowerLimit
            upper_lim = upperLimit
            return lower_lim, upper_lim
        
        elif choice == 'Quantile':
            lower_lim_quantile = lowerLimit
            upper_lim_quantile = upperLimit
            lower_lim = np.quantile(chosen_adata.X.sum(axis=1), lower_lim_quantile)
            upper_lim = np.quantile(chosen_adata.X.sum(axis=1), upper_lim_quantile)
            return lower_lim, upper_lim
        
        elif choice == '1.5IQR':
            total_counts = chosen_adata.X.sum(axis=1)
            Q1 = np.percentile(total_counts, 25)
            Q3 = np.percentile(total_counts, 75)
            IQR = Q3 - Q1
            lower_lim = Q1 - 1.5 * IQR
            upper_lim = Q3 + 1.5 * IQR
            return lower_lim, upper_lim
    def preview_filter(self, user_inputs):
        chosen_adata = user_inputs['chosen_adata']
        min_genes = user_inputs['min_genes'] 
        lower_lim, upper_lim = user_inputs['outlier_lims']
        sc.pp.filter_cells(chosen_adata, min_genes=min_genes)
        chosen_adata = chosen_adata[(chosen_adata.obs.total_counts > lower_lim) &
                                (chosen_adata.obs.total_counts < upper_lim)]

        chosen_adata, imgname = produce_and_save_img(chosen_adata, 'preview', 'previewimage')
                
        n_obs, n_vars = chosen_adata.shape
        chosen_adata_result = f"{chosen_adata.uns['prefix']}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"
        return chosen_adata, chosen_adata_result, imgname
class ReplaceView(APIView):
    def post(self,request):
        image_names = replaceimage()
        update_adata_objects, filtered_adata_objects, adata_results = self.perform_filter()

        save_data(update_adata_objects, 'adata_objects')
        save_data(filtered_adata_objects, 'preview_adata_objects')

        return Response({'adata_results':adata_results,'save_image_names':image_names}, status=status.HTTP_201_CREATED)
    
    def perform_filter(self):
        preview_adata = load_data('preview_adata')
        preview_adata_objects = load_data('preview_adata_objects')
        updated_adata_objects = []
        filtered_adata_objects = []
        adata_results = []
        for adata in preview_adata_objects:
            if adata.uns['prefix'] == preview_adata.uns['prefix']:
                updated_adata_objects.append(preview_adata)
                filtered_adata_objects.append(preview_adata.copy())
            else:
                updated_adata_objects.append(adata)
                filtered_adata_objects.append(adata.copy())
        for adata in updated_adata_objects:
            n_obs, n_vars = adata.shape
            adata_results.append(f"{adata.uns['prefix']}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
            

        return updated_adata_objects, filtered_adata_objects, adata_results
class ConfirmView(APIView):
    def post(self, request):
        adata_objects = load_data('preview_adata_objects')
        filtered_adata_objects = adata_objects.copy()
        
        save_data(adata_objects, 'adata_objects')
        save_data(filtered_adata_objects, 'filtered_adata_objects')
        return Response({}, status=status.HTTP_201_CREATED)
class NormalizationView(APIView):
    def post(self, request):
        chosen_method = request.data.get('normal_method')
        adata_objects, norm_adata_objects, norm_adata_results = self.perform_normalization(chosen_method)
        save_data(adata_objects, 'adata_objects')
        save_data(norm_adata_objects, 'norm_adata_objects')
        return Response({"adata_results": norm_adata_results}, status=status.HTTP_201_CREATED)
    
    def perform_normalization(self, chosen_method):
        adata_objects = load_data('filtered_adata_objects')
        norm_adata_objects = []
        norm_adata_result = []
        if chosen_method == 'cpm':
            for adata in adata_objects:
                sc.pp.normalize_total(adata, target_sum=1e6)
                sc.pp.log1p(adata)
                norm_adata_result.append(f"Normalization completed for {adata.uns['prefix']}")
                norm_adata_objects.append(adata.copy())
        elif chosen_method == 'clr':
            for adata in adata_objects:
                x_pos = adata.X[adata.X > 0] 
                geo_mean = np.exp(np.sum(np.log1p(x_pos)) / len(x_pos))  
                clr_x = np.log1p(x_pos / geo_mean) 
                adata.X[adata.X > 0] = clr_x 
                norm_adata_result.append(f"Normalization completed for {adata.uns['prefix']}")
                norm_adata_objects.append(adata.copy())
        return adata_objects, norm_adata_objects, norm_adata_result
class MergeView(APIView):
    def post(self, request):
        adata_merged_results = self.perform_merged()
        return Response({'adata_results':adata_merged_results}, status=status.HTTP_201_CREATED)
    def perform_merged(self):
        adata_objects = load_data('norm_adata_objects')
        adata_merged_results = []
        
        if len(adata_objects) == 1:
            adata = adata_objects[0]
            adata.uns['is_merged'] = False
        else:
            adata = sc.concat(adata_objects, axis=0)
            adata.obs_names_make_unique()
            adata.uns['is_merged'] = True
        
        for data in [*adata_objects, adata]:
            prefix = data.uns.get('prefix', 'Merged Data')
            adata_merged_results.append(f"{prefix}: AnnData object after filter and normalization = {data.shape[0]} × {data.shape[1]}")
        
        save_h5ad_file(adata, 'adata_preprocessing.h5ad')
        return adata_merged_results

# Clustering
def load_info(input_file): 
    adata = read_h5ad_file(input_file)
    n_obs, n_vars = adata.shape
    if adata.uns.get('is_merged', False):
        return adata, f"Merged Data: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"
    else:
        return adata, f"{adata.uns['prefix']}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"
    
class PreloadPCAView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_preprocessing.h5ad')
        marker_list = adata.var_names.tolist()
        marker_list.insert(0,"Select All") 
        return Response({'marker_list': marker_list}, status=status.HTTP_201_CREATED)
class PCAView(APIView):
    def post(self, request):
        chosen_markers = request.POST.getlist('markers')
        merged_results, n_pcs_results, save_img_name = self.perform_pca(chosen_markers)
        return Response({"merged_results":merged_results,"n_pcs_results": n_pcs_results,"save_img_name":save_img_name}, status=status.HTTP_201_CREATED)

    def perform_pca(self, chosen_markers):
        adata, merged_results = load_info('adata_preprocessing.h5ad')
        adata = adata[:, chosen_markers]
        
        # /imaging/plot_clustering_pca
        clearmediafiles('pca_img')
        save_img_name = plot_clustering_pca(adata)

        variance_ratio = adata.uns['pca']['variance_ratio']
        cumulative_variance_ratio = np.cumsum(variance_ratio)
        n_pcs = np.where(cumulative_variance_ratio >= 0.95)[0][0] + 1
        n_pcs_results = f'Selected number of PCs: {n_pcs}(cumulative_variance_ratio >= 0.95)'

        save_h5ad_file(adata, 'adata_pca.h5ad')
        save_data(n_pcs, 'n_pcs')
        return merged_results, n_pcs_results, save_img_name
class PreloadCLusteringView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_preprocessing.h5ad')
        if adata.uns.get('is_merged'):
            methods = ['none','harmony', 'combat', 'bbknn']
            return Response({'methods': methods}, status=status.HTTP_201_CREATED)
        else:
            methods = ['none']
            return Response({'methods': methods}, status=status.HTTP_201_CREATED)
class CLusteringView(APIView):
    def post(self, request):
        clearmediafiles('cluster_result')
        chosen_method = request.data.get('method')
        n_neighbors = int(request.data.get('n_neighbors'))
        resolution = float(request.data.get('resolution'))
        n_pcs = load_data('n_pcs')
        self.perform_clustering(chosen_method, n_neighbors, resolution, n_pcs) 
        return Response({}, status=status.HTTP_201_CREATED)
        
    def perform_clustering(self, chosen_method, n_neighbors, resolution, n_pcs):
        adata = read_h5ad_file('adata_pca.h5ad') 
        if chosen_method == 'none':
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=42)
        elif chosen_method == 'bbknn':
            sc.external.pp.bbknn(adata, batch_key='Sample', computation='fast')
        elif chosen_method == 'harmony':
            sc.external.pp.harmony_integrate(adata, key='Sample', random_state=42)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep='X_pca_harmony', random_state=42)
        elif chosen_method == 'combat':
            sc.pp.combat(adata, key='Sample')
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=42)
            
        sc.tl.umap(adata, random_state=42)
        sc.tl.leiden(adata, resolution=resolution, random_state=42)        
        adata = clustering_result(adata, n_pcs) # cluster總結果: Summary, Umap*1, Ranking, Heatmap
        
        save_h5ad_file(adata, 'adata_clustering.h5ad')
        return adata
class PreloadMarkersView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_clustering.h5ad')
        marker_list = adata.var_names.tolist()
        marker_list.insert(0,"Select All") 
        return Response({'marker_list': marker_list}, status=status.HTTP_201_CREATED)
class AddUmapClusterView(APIView):
    def post(self, request, met):
        if met == 'leidens':
            return self.post_leidens(request)
        elif met == 'markers':
            return self.post_markers(request)
        else:
            return Response({"error": "Invalid method"}, status=status.HTTP_400_BAD_REQUEST)

    def post_leidens(self, request):
        adata = read_h5ad_file('adata_clustering.h5ad')
        adding_umap(adata)
        # 處理 leidens 邏輯
        return Response({"message": "Leidens processed"}, status=status.HTTP_201_CREATED)

    def post_markers(self, request):
        adata = read_h5ad_file('adata_clustering.h5ad')
        # 處理 markers 邏輯
        return Response({"message": "Markers processed"}, status=status.HTTP_201_CREATED)