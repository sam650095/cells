import os
import pandas as pd
import scanpy as sc
import scimap as sm 
import numpy as np
import bbknn
import logging
import json
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
  
# cleardata
class ClearAllDataView(APIView):
    def post(self, request):
        clear_all_data()
        clearmediafiles('')
        clear_all_h5ad_files()
        save_data(0, 'steps')
        return Response({'message':'data are all killed'}, status=status.HTTP_200_OK)
# Preprocessing     
class FileUploadView(APIView):
    def post(self, request):
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
        available_files = {}
        if adata.uns.get('is_merged'):
            available_files['Merged Data'] = adata.copy()
            methods = ['none','harmony', 'combat', 'bbknn']
            save_data(available_files, 'available_files') 
            return Response({'methods': methods}, status=status.HTTP_201_CREATED)
        else:
            methods = ['none']
            available_files[adata.uns['prefix']] = adata.copy()
            save_data(available_files, 'available_files') 
            return Response({'methods': methods}, status=status.HTTP_201_CREATED)
class CLusteringView(APIView):
    def post(self, request):
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
            print("bbknn start")
            sc.external.pp.bbknn(adata, batch_key='Sample', computation='fast')
            print("bbknn finish")
        elif chosen_method == 'harmony':
            sc.external.pp.harmony_integrate(adata, key='Sample', random_state=42)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep='X_pca_harmony', random_state=42)
        elif chosen_method == 'combat':
            sc.pp.combat(adata, key='Sample')
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=42)
            
        sc.tl.umap(adata, random_state=42)
        sc.tl.leiden(adata, resolution=resolution, random_state=42)        
        adata = clustering_result(adata, n_pcs)
        
        save_h5ad_file(adata, 'adata_clustering.h5ad')
        return adata
class PreloadClusterMarkersView(APIView):
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
        adding_clusteringumap(adata)
        return Response({"message": "Leidens processed"}, status=status.HTTP_201_CREATED)

    def post_markers(self, request):
        adata = read_h5ad_file('adata_clustering.h5ad')
        chosen_method = request.data.getlist('markers')
        adding_clusteringumap(adata, chosen_method)
        return Response({"message": "Markers processed"}, status=status.HTTP_201_CREATED)
class GrabClusterNameView(APIView):
    def post(self, request):
        rename_df = self.grabnames()
        return Response({'rename_df': rename_df}, status=status.HTTP_201_CREATED)
    def grabnames(self):
        adata = read_h5ad_file('adata_clustering.h5ad')
        columns = ['CurrentName', 'NewName']
        rename_df = pd.DataFrame(columns=columns)
        rename_df['CurrentName'] = adata.obs['leiden_R'].cat.categories.tolist()
        rename_df['NewName'] = adata.obs['leiden_R'].cat.categories.tolist()
        return rename_df
class RenameClusterView(APIView):
    def post(self, request):
        data = json.loads(request.body)
        self.perform_rename_cluster(pd.DataFrame(data))
        return Response({}, status=status.HTTP_201_CREATED)
    def perform_rename_cluster(self, rename_df):
        adata = read_h5ad_file('adata_clustering.h5ad')
        n_pcs = load_data('n_pcs')

        rename_dict = {}
        for _, row in rename_df.iterrows():
            current_name = str(row['CurrentName'])
            new_name = str(row['NewName'])
            rename_dict[current_name] = new_name

        adata.obs['leiden_R'] = adata.obs['leiden_R'].astype(str).map(rename_dict).astype('category')
        adata.obs['leiden'] = adata.obs['leiden_R']

        clustering_result(adata, n_pcs)
        save_h5ad_file(adata, 'adata_clustering.h5ad')
        return adata
class GrabClustersView(APIView):
    def post(self, request):
        clusters_list = self.grabclusters()
        return Response({'clusters_list': clusters_list}, status=status.HTTP_201_CREATED)
    def grabclusters(self):
        adata = read_h5ad_file('adata_clustering.h5ad')
        cluster_list = adata.obs['leiden'].cat.categories
        return cluster_list
class SubclusterView(APIView):
    def post(self, request):
        chosen_clusters = request.data.getlist('clusters') 
        resolution = float(request.data.get('resolution'))
        self.perform_subclustering(resolution, chosen_clusters)
        return Response({}, status=status.HTTP_201_CREATED)
    def perform_subclustering(self, resolution, chosen_clusters):
        adata = read_h5ad_file('adata_clustering.h5ad')
        n_pcs = load_data('n_pcs')
        sc.tl.leiden(adata, resolution=resolution, restrict_to=('leiden', chosen_clusters))
        clustering_result(adata, n_pcs) 
        save_h5ad_file(adata, 'adata_clustering.h5ad') 
        return adata
class PreloadSubsetView(APIView):
    def post(self, request):
        available_files_result, clustering_columns = self.available_adata()
        return Response({'available_files_result':available_files_result,'clustering_columns':clustering_columns}, status=status.HTTP_201_CREATED)
    def available_adata(self):
        adata = read_h5ad_file('adata_clustering.h5ad')
        clustering_columns = [col for col in adata.obs.columns if col.startswith(("leiden_R", "phenotype", "Sample"))]
        available_files = load_data('available_files')
        available_files_result = []
        if adata.uns.get('is_merged', False):
            available_files['Merged Data'] = adata.copy()
        else:
            available_files[adata.uns['prefix']] = adata.copy()

        for name, adata_item in available_files.items():
            n_obs, n_vars = adata_item.shape
            available_files_result.append(f"{name}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
        return available_files_result, clustering_columns
class GrabClusterSubsetView(APIView):
    def post(self, request):
        chosen_cluster = request.data.get('sample')
        adata = read_h5ad_file('adata_clustering.h5ad')
        cluster = adata.obs[chosen_cluster].cat.categories
        return Response({'cluster':cluster}, status=status.HTTP_201_CREATED)
class SubsetView(APIView):
    def post(self, request, met):
        if(met == "new"):
            available_files_result = self.post_new(request)
        if(met == "merged"):
            available_files_result = self.post_merged()
        return Response({'available_files_result':available_files_result}, status=status.HTTP_201_CREATED)
    def post_merged(self):
        # for no subset
        adata = read_h5ad_file('adata_clustering.h5ad')
        chosen_cluster = "Sample"
        chosen_cluster_names = list(adata.obs[chosen_cluster].cat.categories)
        if(adata.uns['is_merged'] == True):
            subset_name = "Merged Data"
        else:
            subset_name = adata.uns['prefix']
        available_files_result = self.add_subset_cluster(adata, chosen_cluster, chosen_cluster_names, subset_name)
        return available_files_result
    def post_new(self, request):
        # for subset
        adata = read_h5ad_file('adata_clustering.h5ad')
        chosen_cluster = request.data.get('sample')
        chosen_cluster_names = request.data.getlist('clusters')
        subset_name = request.data.get('naming_cluster')
        available_files_result = self.add_subset_cluster(adata, chosen_cluster, chosen_cluster_names, subset_name)
        return available_files_result
    def add_subset_cluster(self, adata, chosen_cluster, chosen_cluster_names, subset_name):
        available_files = load_data('available_files')
        available_files_result = []
        subset_adata = adata[adata.obs[chosen_cluster].isin(chosen_cluster_names)].copy()
        available_files[subset_name] = subset_adata
        for name, adata in available_files.items(): 
            n_obs, n_vars = adata.shape
            available_files_result.append(f"{name}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
        save_data(available_files, 'available_files') 
        return available_files_result    

# Identify the gates
class PreloadIdentifytheGatesView(APIView):
    def post(self, request):
        available_files = load_data('available_files')
        adata_list = list(available_files.keys())
        return Response({'adata_list': adata_list}, status=status.HTTP_201_CREATED)
    
class ChosenAdataResultView(APIView):
    def post(self, request):
        chosen_adata = request.data.get('chosen_adata')
        available_files = load_data('available_files')
        adata = available_files[chosen_adata]
        n_obs, n_vars = adata.shape
        result = f"{chosen_adata}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"
        save_h5ad_file(adata, 'adata_identify_gate.h5ad')
        save_data(chosen_adata, 'chosen_adata') 
        return Response({'result': result}, status=status.HTTP_201_CREATED)
    

class IdentifytheGatesView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_identify_gate.h5ad')
        columns = ['Marker'] + list(adata.obs['Sample'].cat.categories)
        gate_df = pd.DataFrame(columns=columns)
        gate_df['Marker'] = adata.var_names
        save_data(gate_df, 'gate_df') 
        gate_dict = gate_df.applymap(self.clean_value).to_dict(orient='records')
        
        return Response({'gate_df': gate_dict}, status=status.HTTP_201_CREATED)
    def clean_value(self, x):
        if pd.isna(x):
            return None
        if isinstance(x, (int, float)):
            if np.isinf(x):
                return None
            return x
        return str(x) 
class AddValueView(APIView):
    def post(self, request):
        editdata = json.loads(request.POST.get('editdata'))
        gate_df = load_data('gate_df')
        edit_df = pd.DataFrame(editdata)
        edit_df = edit_df.apply(pd.to_numeric, errors='coerce')
        gate_df.set_index('Marker', inplace=True)
        edit_df.set_index('Marker', inplace=True)
        for marker in edit_df.index:
            for column in edit_df.columns:
                if marker in gate_df.index and column in gate_df.columns:
                    new_value = edit_df.loc[marker, column]
                    if pd.notna(new_value):  
                        gate_df.loc[marker, column] = new_value
        
        gate_df.reset_index(inplace=True)
        gate_df = gate_df.replace([np.inf, -np.inf], None)
        gate_df = gate_df.where(pd.notnull(gate_df), None)
        
        result = gate_df.to_dict(orient='records')
        return Response({'result':result}, status=status.HTTP_201_CREATED)
    
# Phenotyping
class PhenotypingView(APIView):
    def post(self, request):
        uploaded_file = request.FILES['file']
        phenotype_dir = os.path.join(settings.MEDIA_ROOT, 'tempfile', 'phenotypefile')
        os.makedirs(phenotype_dir, exist_ok=True)
        file_name = uploaded_file.name
        file_path = os.path.join(phenotype_dir, file_name)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)

        with open(file_path, 'wb+') as destination:
            for chunk in uploaded_file.chunks():
                destination.write(chunk)
        phenotyping_result = self.perform_phenotyping(file_path)
        return Response({'phenotyping_result': phenotyping_result}, status=status.HTTP_201_CREATED)

    def perform_phenotyping(self, file_path):
        adata = read_h5ad_file('adata_identify_gate.h5ad')
        gate_df = load_data('gate_df')
        chosen_adata = load_data('chosen_adata')
        n_pcs = load_data('n_pcs')
        adata = sm.pp.rescale(adata, gate = gate_df, imageid = 'Sample', method = 'by_image')
        phenotype = pd.read_csv(file_path, encoding='utf-8')
        adata = sm.tl.phenotype_cells(adata, phenotype=phenotype,
                                    label="phenotype", imageid = 'Sample')
        adata, phenotyping_result = phenotype_result(adata, chosen_adata, n_pcs) 
        
        save_h5ad_file(adata, 'adata_phenotyping.h5ad')
        return phenotyping_result
    
class AddUmapPhenotypeView(APIView):
    def post(self, request, met):
        if met == 'phenotypes':
            return self.post_phenotypes(request)
        elif met == 'markers':
            return self.post_markers(request)
        else:
            return Response({"error": "Invalid method"}, status=status.HTTP_400_BAD_REQUEST)

    def post_phenotypes(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        add_phenotypes_bysample(adata)
        return Response({"message": "Phenotypes processed"}, status=status.HTTP_201_CREATED)

    def post_markers(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        chosen_method = request.data.getlist('markers')
        add_phenotypes_markers(adata, chosen_method)
        return Response({"message": "Markers processed"}, status=status.HTTP_201_CREATED)
class PreloadPhenotypesMarkersView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        marker_list = adata.var_names.tolist()
        marker_list.insert(0,"Select All") 
        return Response({'marker_list': marker_list}, status=status.HTTP_201_CREATED)
class GrabPhenotypesNameView(APIView):
    def post(self, request):
        rename_df = self.grabnames()
        return Response({'rename_df': rename_df}, status=status.HTTP_201_CREATED)
    def grabnames(self):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        columns = ['CurrentName', 'NewName']
        rename_df = pd.DataFrame(columns=columns)
        rename_df['CurrentName'] = adata.obs['phenotype'].cat.categories.tolist()
        rename_df['NewName'] = adata.obs['phenotype'].cat.categories.tolist()
        return rename_df
class RenamePhenotypeView(APIView):
    def post(self, request):
        data = json.loads(request.body)
        self.perform_rename_cluster(pd.DataFrame(data))
        return Response({}, status=status.HTTP_201_CREATED)
    def perform_rename_cluster(self, rename_df):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        n_pcs = load_data('n_pcs')

        rename_dict = {}
        for _, row in rename_df.iterrows():
            current_name = str(row['CurrentName'])
            new_name = str(row['NewName'])
            rename_dict[current_name] = new_name

        adata.obs['phenotype'] = adata.obs['phenotype'].astype(str).map(rename_dict).astype('category')
        chosen_adata = load_data('chosen_adata')
        n_pcs = load_data('n_pcs')
        phenotype_result(adata, chosen_adata, n_pcs)

        save_h5ad_file(adata, 'adata_phenotyping.h5ad')
        return adata
class GrabDropPhenotypeView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        drop_list = adata.obs['phenotype'].cat.categories.tolist()
        return Response({'drop_list': drop_list}, status=status.HTTP_201_CREATED)
    
class DropPhenotypeView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        drop_groups = request.data.getlist('drops')
        adata = sm.hl.dropFeatures(adata, drop_groups=drop_groups, groups_column='phenotype', subset_raw=False)
        chosen_adata = load_data('chosen_adata')
        n_pcs = load_data('n_pcs')
        phenotype_result(adata, chosen_adata, n_pcs)

        save_h5ad_file(adata, 'adata_phenotyping.h5ad')
        return Response({'message': 'drop phenotype successfully!'}, status=status.HTTP_201_CREATED)
  
# Spatial Analysis
class PreloadSpatialAnalysisView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        chosen_adata = load_data('chosen_adata')

        n_obs, n_vars = adata.shape
        return Response({'preload_result': f"{chosen_adata}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"}, status=status.HTTP_201_CREATED)
class SpatialAnalysisView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        columns_list = [col for col in adata.obs.columns if col.startswith(("leiden_R", "phenotype"))]
        default_chosen_column = 'phenotype'
        cluster_list = adata.obs[default_chosen_column].unique()
        cluster_list_L = adata.obs['leiden_R'].unique()
        default_chosen_cluster = cluster_list[0]
        method_list = ['radius', 'knn']
        default_chosen_method = 'radius'

        for col in columns_list: 
            sm.tl.spatial_distance(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid', 
                                phenotype=col, imageid='Sample', label=f'spatial_distance_{col}')
        for col in columns_list:
            sm.tl.spatial_interaction(adata, phenotype=col, method='knn', radius=30, 
                                    imageid='Sample', label=f'spatial_interaction_{col}_knn', pval_method='zscore')
            sm.tl.spatial_interaction(adata, phenotype=col, method='radius', radius=30, 
                                    imageid='Sample', label=f'spatial_interaction_{col}_radius', pval_method='zscore')
        filename = []
        filename.append(interactions_voronoi(adata, default_chosen_column))
        filename.append(distances_heatmap(adata, default_chosen_column))
        filename.append(distances_numeric_plot(adata, default_chosen_column, default_chosen_cluster))
        filename.append(interactions_heatmap(adata, default_chosen_column, default_chosen_method))
        
        save_h5ad_file(adata, 'adata_spatial_analysis.h5ad')

        return Response({'columns_list': columns_list,'cluster_list':cluster_list, 'cluster_list_L':cluster_list_L, 'method_list':method_list, "filename":filename}, status=status.HTTP_201_CREATED)
class AddSpatial(APIView):
    def post(self, request, met):
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        filename = ''
        if(met == "DH"):
            filename = self.post_dh(request)
        elif(met == "NP"):
            filename = self.post_np(request)
        elif(met == "IH"):
            filename = self.post_ih(request)
        else:
            filename = self.post_vp(request)
        save_h5ad_file(adata, 'adata_spatial_analysis.h5ad')

        return Response({"filename":filename}, status=status.HTTP_201_CREATED)
    def post_dh(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        filename = distances_heatmap(adata, json.loads(request.body).get('dh_ul_input'))
        return filename
    def post_np(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        filename = distances_numeric_plot(adata, json.loads(request.body).get('np_ul_input'), json.loads(request.body).get('npd_ul_input'))
        return filename
    def post_ih(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        filename = interactions_heatmap(adata, json.loads(request.body).get('ih_ul_input'), json.loads(request.body).get('ihm_ul_input'))
        return filename
    def post_vp(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        filename = interactions_voronoi(adata, json.loads(request.body).get('vp_ul_input'))
        return filename

# Neighborhood
class PreloadNeighborView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        n_obs, n_vars = adata.shape
        chosen_adata = load_data('chosen_adata')
        return Response({"preload_result":f"{chosen_adata}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"}, status=status.HTTP_201_CREATED)
        
class NeighborView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        chosen_column = json.loads(request.body).get('n_input')
        k = json.loads(request.body).get('k_neighbor')
        n_neighborhoods =json.loads(request.body).get('n_neighbor')
        perform_neighborhoods(adata, chosen_column, k, n_neighborhoods)
        return Response({"test":"test"}, status=status.HTTP_201_CREATED)
