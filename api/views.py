import os
import pandas as pd
import scanpy as sc
import numpy as np
import shutil

from django.conf import settings
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from .serializers import FileUploadSerializer
from .imaging import *
from .redis import *

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
            
class FileUploadView(APIView):
    def post(self, request):
        clear_media = clearmediafiles('tempfile')
        # print(clear_media)
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
        adata_objects, original_adata_objects, adata_result = self.create_adata(signal_value_files, metadata_files)
        save_adata_objects(adata_objects, 'adata_objects')
        save_adata_objects(original_adata_objects, 'original_adata_objects')
        return Response({"adata_result": adata_result}, status=status.HTTP_201_CREATED)
    
    def create_adata(self, signal_value_files, metadata_files):
        adata_objects = []
        original_adata_objects = []
        adata_result = []
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
                adata_result.append(f"{prefix}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
                adata_objects.append(adata)
                original_adata_objects.append(adata.copy())
        return adata_objects, original_adata_objects, adata_result
    
    
class QualityControlView(APIView):
    def post(self, request):
        clear_media = clearmediafiles('tempimage')
        # print(clear_media)
        adata_objects, qc_adata_objects, adata_result, save_image_names = self.post_reset_adata()
        save_adata_objects(adata_objects, 'adata_objects')
        save_adata_objects(qc_adata_objects, 'qc_adata_objects')
        save_adata_objects(qc_adata_objects, 'preview_adata_objects')
        return Response({'adata_result': adata_result,'save_image_names': save_image_names}, status=status.HTTP_201_CREATED)
    
    def post_reset_adata(self):
        adata_objects = load_adata_objects('original_adata_objects')
        qc_adata_objects = []
        adata_result = []
        save_image_names = []
        for adata in adata_objects:
            n_obs, n_vars = adata.shape
            adata_result.append(f"{adata.uns['prefix']}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
            sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
            adata, imgname = produce_and_save_img(adata, "origin", 'violin_plot')
            qc_adata_objects.append(adata)
            save_image_names.append(imgname)
        
        return adata_objects, qc_adata_objects, adata_result, save_image_names
    
class PreviewView(APIView):
    def post(self, request):
        adata_objects = load_adata_objects('qc_adata_objects')
        # get data from frontend
        f_sampleSelect = request.data.get('f_sampleSelect')
        minGenes = int(request.data.get('minGenes'))
        filter_method = request.data.get('filter_method')
        lowerlimit = float(request.data.get('lowerlimit')) if request.data.get('lowerlimit') else None
        upperlimit = float(request.data.get('upperlimit')) if request.data.get('upperlimit') else None
        user_inputs = {}
        user_inputs['chosen_adata'] = self.choose_sample(adata_objects,f_sampleSelect)
        user_inputs['min_genes'] = minGenes
        user_inputs['outlier_lims'] = self.input_outliers(user_inputs['chosen_adata'], filter_method, lowerlimit, upperlimit)

        preview_adata, preview_adata_result, preview_image_name = self.preview_filter(user_inputs)
        # chose adata to preview
        save_adata_objects(preview_adata, 'preview_adata')
        return Response({'adata_result': preview_adata_result,'save_image_names': preview_image_name}, status=status.HTTP_201_CREATED)
    
    def choose_sample(self,adata_objects, sample_select):
        for adata in adata_objects:
            if adata.uns['prefix'] == sample_select:
                chosen_adata = adata
                return chosen_adata

    def input_outliers(self, chosen_adata, filtermethod, lowerLimit, upperLimit):   
        choice = filtermethod
        if choice == 'manual':
            lower_lim = lowerLimit
            upper_lim = upperLimit
            return lower_lim, upper_lim
        
        elif choice == 'quantile':
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
        preview_adata = load_adata_objects('preview_adata')
        preview_adata_objects = load_adata_objects('preview_adata_objects')
        update_adata_objects, filtered_adata_objects, updated_adata_result = self.perform_filter(preview_adata_objects, preview_adata)

        save_adata_objects(update_adata_objects, 'adata_objects')
        save_adata_objects(filtered_adata_objects, 'preview_adata_objects')

        n_obs, n_vars = preview_adata.shape
        preview_adata_result = f"{preview_adata.uns['prefix']}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"
        return Response({'adata':preview_adata_result,'save_image_names':image_names}, status=status.HTTP_201_CREATED)
    
    def perform_filter(self, preview_adata_objects, preview_adata):
        updated_adata_objects = []
        filtered_adata_objects = []
        updated_adata_result = []
        for adata in preview_adata_objects:
            if adata.uns['prefix'] == preview_adata.uns['prefix']:
                updated_adata_objects.append(preview_adata)
                filtered_adata_objects.append(preview_adata.copy())
            else:
                updated_adata_objects.append(adata)
                filtered_adata_objects.append(adata.copy())
        for adata in updated_adata_objects:
            n_obs, n_vars = adata.shape
            updated_adata_result.append(f"{adata.uns['prefix']}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
            sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True)

        return updated_adata_objects, filtered_adata_objects, updated_adata_result
class ConfirmView(APIView):
    def post(self, request):
        return Response({}, status=status.HTTP_201_CREATED)
    
class NormalizationView(APIView):
    def post(self, request):
        adata_objects = load_adata_objects('filtered_adata_objects')
        chosen_method = request.data.get('normal_method')
        adata_objects, norm_adata_objects, norm_adata_result = self.perform_normalization(adata_objects, chosen_method)
        save_adata_objects(norm_adata_objects, 'norm_adata_objects')
        return Response({'message': 'Succeed', 'normal_result':norm_adata_result})
    def perform_normalization(self, adata_objects, chosen_method):
        norm_adata_objects = []
        norm_adata_result = []
        if chosen_method == 'CPM':
            for adata in adata_objects:
                sc.pp.normalize_total(adata, target_sum=1e6)
                sc.pp.log1p(adata)
                norm_adata_result.append(f"Normalization completed for {adata.uns['prefix']}")
                norm_adata_objects.append(adata.copy())
        elif chosen_method == 'CLR':
            for adata in adata_objects:
                x_pos = adata.X[adata.X > 0]  # 提取非零元素 
                geo_mean = np.exp(np.sum(np.log1p(x_pos)) / len(x_pos))  # 計算幾何平均值: exp(非零元素的對數和/元素個數)
                clr_x = np.log1p(x_pos / geo_mean)  # log(x / geo_mean),得到標準化後的值
                adata.X[adata.X > 0] = clr_x # 將clr_x放回原矩陣相應的位置
                norm_adata_result.append(f"Normalization completed for {adata.uns['prefix']}")
                norm_adata_objects.append(adata.copy())
        return adata_objects, norm_adata_objects, norm_adata_result