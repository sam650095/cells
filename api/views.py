import os
import pandas as pd
import scanpy as sc
import scimap as sm 
import numpy as np
import json
from django.conf import settings
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from .imaging import *
from .saved import *
from .models import OperationStep


def sort_key(filename):
    return int(filename.split('_')[0][1:]) 
def SaveSteps(session_id, step, operation_type, input_values, output_values):
    OperationStep.objects.update_or_create(
        session_id=session_id,
        defaults={
            'step':step,
            'operation_type':operation_type,
            'input_values':input_values,
            'output_values':output_values
        }
    )
def GetSteps(session_id):
    steps = OperationStep.objects.filter(session_id=session_id)
    if steps.exists():
        return steps
    else:
        return None

# cleardata
class ClearAllDataView(APIView):
    def post(self, request):
        OperationStep.objects.all().delete()
        
        # clear_all_data()
        clearmediafiles('')
        clear_all_h5ad_files()
        return Response({'message':'data are all killed'}, status=status.HTTP_200_OK)
# Preprocessing     
class FileUploadView(APIView):
    def post(self, request):
        clearmediafiles('tempfile')
        files = request.FILES.getlist('files')
        if(len(files) < 2):
            return Response({'message': "Please uplaod unless 2 files."}, status=status.HTTP_400_BAD_REQUEST)
        file_names = [file.name for file in files]
        file_sizes = [file.size for file in files]
        # file types
        file_type_mapping = {
            'metadata.csv': 'metadata',
            # signal_value
            'value.csv': 'signal_value'
        }
        for file in files:
            file_extension = file.name.split('_')[-1]
            target_dir = file_type_mapping.get(file_extension)
            
            if not target_dir:
                return Response({'message': "Invalid File Type, Please read the rule above and upload again."}, status=status.HTTP_400_BAD_REQUEST)
            
            target_path = os.path.join(settings.MEDIA_ROOT, 'tempfile', target_dir)
            os.makedirs(target_path, exist_ok=True) 
            new_file_path = os.path.join(target_path, file.name)

            with open(new_file_path, 'wb') as f:
                for chunk in file.chunks():
                    f.write(chunk)
            
        signal_value_files = sorted(os.listdir(os.path.join(settings.MEDIA_ROOT, 'tempfile', 'signal_value')), key=sort_key)
        metadata_files = sorted(os.listdir(os.path.join(settings.MEDIA_ROOT, 'tempfile', 'metadata')), key=sort_key)
        result = self.create_adata(signal_value_files, metadata_files)
        if isinstance(result, str):
            return Response({'message': result}, status=status.HTTP_400_BAD_REQUEST)
        adata_objects, original_adata_objects, adata_results = result
        for i, adata in enumerate(original_adata_objects):
            save_h5ad_file(adata, f'original_adata_objects_{i}')
        # saving steps
        SaveSteps(1, 'create_adata', 'file_upload', {'file_names': file_names, 'file_sizes': file_sizes}, {'adata_results': adata_results})
        
        return Response({"adata_results": adata_results}, status=status.HTTP_201_CREATED)
    
    def create_adata(self, signal_value_files, metadata_files):
        adata_objects = []
        original_adata_objects = []
        adata_results = []
        signal_value_dir = os.path.join(settings.MEDIA_ROOT, 'tempfile', 'signal_value')
        metadata_dir = os.path.join(settings.MEDIA_ROOT, 'tempfile', 'metadata')
        if(len(signal_value_files)!=len(metadata_files)):   
            return "File's prefix does not match, Please read the rule above and upload again."
        for signal_value_file in signal_value_files:
            prefix = signal_value_file.split('_')[0]
            matching_metadata_file = next((metadata_file for metadata_file in metadata_files if metadata_file.startswith(prefix)), None)
            if matching_metadata_file!=None:
                data = pd.read_csv(os.path.join(signal_value_dir, signal_value_file), header=0, index_col=0)
                metadata = pd.read_csv(os.path.join(metadata_dir, matching_metadata_file), header=0, index_col=0)
                adata = sc.AnnData(X=data, obs=metadata)
                adata.raw = adata
                adata.uns['prefix'] = prefix
                n_obs, n_vars = adata.shape
                adata_results.append(f"{prefix}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
                adata_objects.append(adata)
                original_adata_objects.append(adata.copy())
            else:
                return "File's prefix does not match, Please read the rule above and upload again."
        return adata_objects, original_adata_objects, adata_results
class QualityControlView(APIView):
    def post(self, request):
        adata_objects, qc_adata_objects, adata_results, save_image_names, marker_list = self.post_reset_adata()
        # save anndata
        for i, adata in enumerate(qc_adata_objects):
            save_h5ad_file(adata, f'qc_adata_objects_{i}')
            save_h5ad_file(adata, f'preview_adata_objects_{i}')
        # saving steps
        SaveSteps(2, 'qualitycontrol', 'process', {}, {'adata_results': adata_results,'save_image_names': save_image_names, 'marker_list': marker_list})
        # default filter
        SaveSteps(3, 'qualitycontrol', 'filter', {}, {})
        # default qui
        SaveSteps(4, 'qualitycontrol', 'qui', {}, {})
        
        return Response({'adata_results': adata_results,'save_image_names': save_image_names, 'marker_list': marker_list}, status=status.HTTP_201_CREATED)
    
    def post_reset_adata(self):
        # read anndata
        files = os.listdir(settings.H5AD_STORAGE_PATH)
        origin_adata_files = [file for file in files if file.startswith('original_adata_objects')]
        adata_objects = [read_h5ad_file(file) for file in origin_adata_files]

        marker_list = []
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
            marker_list.append(adata.var_names.tolist())
        
        return adata_objects, qc_adata_objects, adata_results, save_image_names, marker_list
class PreviewView(APIView):
    def post(self, request):
        # read anndata
        files = os.listdir(settings.H5AD_STORAGE_PATH)
        origin_adata_files = [file for file in files if file.startswith('qc_adata_objects')]
        adata_objects = [read_h5ad_file(file) for file in origin_adata_files]

        # get data from frontend
        f_sampleSelect = request.data.get('sample')
        minGenes = int(request.data.get('minGenes'))
        filter_method = request.data.get('filter_sampleul_input')
        lowerlimit = float(request.data.get('lowerlimit')) if request.data.get('lowerlimit') else None
        upperlimit = float(request.data.get('upperlimit')) if request.data.get('upperlimit') else None
        user_inputs = {}
        user_inputs['chosen_adata'] = self.choose_sample(adata_objects,f_sampleSelect)
        user_inputs['min_genes'] = minGenes
        user_inputs['outlier_lims'] = self.input_outliers(user_inputs['chosen_adata'], filter_method, lowerlimit, upperlimit)

        preview_adata, preview_adata_result, preview_image_name = self.preview_filter(user_inputs)
        # chose adata to preview
        save_h5ad_file(preview_adata, 'preview_adata')
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
        f_sampleSelect = request.data.get('sample')
        minGenes = int(request.data.get('minGenes'))
        filter_method = request.data.get('filter_sampleul_input')
        lowerlimit = float(request.data.get('lowerlimit')) if request.data.get('lowerlimit') else None
        upperlimit = float(request.data.get('upperlimit')) if request.data.get('upperlimit') else None
        filtered_adata_objects, adata_results, changed_adata_result = self.perform_filter()
        
        # save anndata
        for i, adata in enumerate(filtered_adata_objects):
            save_h5ad_file(adata, f'filtered_adata_objects_{i}')
        # update steps
        operation_step = OperationStep.objects.get(session_id=2)
        print(operation_step.output_values.get('adata_results'))
        SaveSteps(2, 'qualitycontrol', 'process', {}, {'adata_results': adata_results,'save_image_names': operation_step.output_values.get('save_image_names'), 'marker_list':operation_step.output_values.get('marker_list')})
       
        SaveSteps(3, 'qualitycontrol', 'filter', {"f_sampleSelect":f_sampleSelect, "minGenes":minGenes, "filter_method":filter_method, "lowerlimit":lowerlimit, "upperlimit":upperlimit}, {'adata_result':changed_adata_result, 'save_image_names':f_sampleSelect+"_previewimage.png"})
        
        return Response({'adata_results':adata_results,'save_image_names':image_names}, status=status.HTTP_201_CREATED)
    
    def perform_filter(self):
        preview_adata = read_h5ad_file('preview_adata')
        # read anndata
        files = os.listdir(settings.H5AD_STORAGE_PATH)
        adata_files = [file for file in files if file.startswith('preview_adata_objects')]
        preview_adata_objects = [read_h5ad_file(file) for file in adata_files]

        updated_adata_objects = []
        filtered_adata_objects = []
        adata_results = []
        changed_adata_result = ""
        for adata in preview_adata_objects:
            if adata.uns['prefix'] == preview_adata.uns['prefix']:
                updated_adata_objects.append(preview_adata)
                filtered_adata_objects.append(preview_adata.copy())
            else:
                updated_adata_objects.append(adata)
                filtered_adata_objects.append(adata.copy())
        for adata in updated_adata_objects:
            n_obs, n_vars = adata.shape
            if adata.uns['prefix'] == preview_adata.uns['prefix']:
                changed_adata_result = f"{adata.uns['prefix']}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"
            adata_results.append(f"{adata.uns['prefix']}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
        print(adata_results)

        return filtered_adata_objects, adata_results, changed_adata_result
class ConfirmView(APIView):
    def post(self, request):
        # read anndata
        files = os.listdir(settings.H5AD_STORAGE_PATH)
        adata_files = [file for file in files if file.startswith('preview_adata_objects')]
        adata_objects = [read_h5ad_file(file) for file in adata_files]
        f_adata_files = [file for file in files if file.startswith('filtered_adata_objects')]
        # print(f_adata_files)
        if(len(f_adata_files) == 0):
            filtered_adata_objects = adata_objects.copy()
            for i, adata in enumerate(filtered_adata_objects):
                save_h5ad_file(adata, f'filtered_adata_objects_{i}')
        return Response({}, status=status.HTTP_201_CREATED)
class NormalizationView(APIView):
    def post(self, request):
        chosen_method = request.data.get('method')
        adata_objects, norm_adata_objects, norm_adata_results = self.perform_normalization(chosen_method)
        # save anndata
        for i, adata in enumerate(norm_adata_objects):
            save_h5ad_file(adata, f'norm_adata_objects_{i}')
            print(f'{adata.shape[0]} × {adata.shape[1]}')
        SaveSteps(5, 'normalization', 'process', {}, {"adata_results": norm_adata_results})
        
        return Response({"adata_results": norm_adata_results}, status=status.HTTP_201_CREATED)
    
    def perform_normalization(self, chosen_method):
        files = os.listdir(settings.H5AD_STORAGE_PATH)
        adata_files = [file for file in files if file.startswith('filtered_adata_objects')]
        adata_objects = [read_h5ad_file(file) for file in adata_files]
        norm_adata_objects = []
        norm_adata_result = []
        if chosen_method == 'CPM':
            for adata in adata_objects:
                sc.pp.normalize_total(adata, target_sum=1e6)
                sc.pp.log1p(adata)
                norm_adata_result.append(f"Normalization completed for {adata.uns['prefix']}")
                print(f"{adata.shape[0]} × {adata.shape[1]}")
                norm_adata_objects.append(adata.copy())
        elif chosen_method == 'CLR':
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
        SaveSteps(6, 'merge', 'process', {}, {"adata_results": adata_merged_results})
        
        return Response({'adata_results':adata_merged_results}, status=status.HTTP_201_CREATED)
    def perform_merged(self):
        files = os.listdir(settings.H5AD_STORAGE_PATH)
        adata_files = [file for file in files if file.startswith('norm_adata_objects')]
        adata_objects = [read_h5ad_file(file) for file in adata_files]
        adata_merged_results = []
        
        if len(adata_objects) == 1:
            adata = adata_objects[0]
            adata.uns['is_merged'] = False
            prefix = adata.uns.get('prefix', 'Merged Data')
            adata_merged_results.append(f"{prefix}: AnnData object after filter and normalization = {adata.shape[0]} × {adata.shape[1]}")
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
        merged_results, n_pcs, n_pcs_results, save_img_name = self.perform_pca(chosen_markers)
        SaveSteps(7, 'pca', 'process', {"chosen_markers": chosen_markers}, {"merged_results":merged_results, "n_pcs": n_pcs, "n_pcs_results": n_pcs_results, "save_img_name": save_img_name})
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
        return merged_results, n_pcs, n_pcs_results, save_img_name
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
        adata = read_h5ad_file('adata_preprocessing.h5ad')
        available_files = {}
        chosen_method = request.data.get('method')
        n_neighbors = int(request.data.get('n_neighbors'))
        resolution = float(request.data.get('resolution'))
        
        n_pcs = GetSteps(7)[0].output_values['n_pcs']
        # n_pcs = load_data('n_pcs')
        if adata.uns.get('is_merged'):
            available_files['Merged Data'] = adata.copy()
            save_h5ad_file(available_files['Merged Data'], 'available_files.h5ad')
        else:
            available_files[adata.uns['prefix']] = adata.copy()
            save_h5ad_file(available_files[adata.uns['prefix']], 'available_files.h5ad')

        self.perform_clustering(chosen_method, n_neighbors, resolution, n_pcs) 
        SaveSteps(8, 'clustering', 'process', {"chosen_method": chosen_method, "resolution":resolution, "n_neighbors": n_neighbors, "n_pcs": n_pcs}, {})
        SaveSteps(9, 'clustering', 'rename', {},{})
        SaveSteps(10, 'clustering', 'subcluster', {},{})
        SaveSteps(11, 'clustering', 'subset', {},{})
        SaveSteps(12, 'clustering', 'umap_add_image', {},{})
        SaveSteps(13, 'clustering', 'umap_add_l_image', {},{})
        return Response({}, status=status.HTTP_201_CREATED)
        
    def perform_clustering(self, chosen_method, n_neighbors, resolution, n_pcs):
        adata = read_h5ad_file('adata_pca.h5ad') 
        if chosen_method == 'none':
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=42)
        elif chosen_method == 'bbknn':
            sc.external.pp.bbknn(adata, batch_key='Sample')
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
        SaveSteps(12, 'clustering', 'umap_add_image', {'chosen_method':chosen_method},{})
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
        SaveSteps(9, 'clustering', 'rename', {},{})
        return Response({}, status=status.HTTP_201_CREATED)
    def perform_rename_cluster(self, rename_df):
        adata = read_h5ad_file('adata_clustering.h5ad')
        n_pcs = GetSteps(7)[0].output_values['n_pcs']

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
        SaveSteps(10, 'clustering', 'subcluster', {"chosen_clusters":chosen_clusters,"resolution":resolution},{})
        return Response({}, status=status.HTTP_201_CREATED)
    def perform_subclustering(self, resolution, chosen_clusters):
        adata = read_h5ad_file('adata_clustering.h5ad')
        n_pcs = GetSteps(7)[0].output_values['n_pcs']
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
        available_files_saved = read_h5ad_file('available_files.h5ad')
        available_files= {}
        available_files_result = []
        if adata.uns.get('is_merged', False):
            available_files['Merged Data'] = available_files_saved
        else:
            available_files[adata.uns['prefix']] = available_files_saved

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
        SaveSteps(11, 'clustering', 'subset', {"chosen_cluster":chosen_cluster,"chosen_cluster_names":chosen_cluster_names,"subset_name":subset_name},{"available_files_result":available_files_result})
        return available_files_result
    def add_subset_cluster(self, adata, chosen_cluster, chosen_cluster_names, subset_name):
        available_files_saved = read_h5ad_file('available_files.h5ad')
        available_files= {}
        if adata.uns.get('is_merged', False):
            available_files['Merged Data'] = available_files_saved
        else:
            available_files[adata.uns['prefix']] = available_files_saved

        available_files_result = []
        subset_adata = adata[adata.obs[chosen_cluster].isin(chosen_cluster_names)].copy()
        available_files[subset_name] = subset_adata
        for name, adata in available_files.items(): 
            n_obs, n_vars = adata.shape
            available_files_result.append(f"{name}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}")
        # save_data(available_files, 'available_files') 
        save_h5ad_file(available_files[subset_name], 'available_files_'+subset_name+'.h5ad')
        return available_files_result    

# Identify the gates
class PreloadIdentifytheGatesView(APIView):
    def post(self, request):
        h5ad_files = [f for f in os.listdir(settings.H5AD_STORAGE_PATH) 
                  if f.startswith('available_files_') and f.endswith('.h5ad')]
    
        adata_list = [file.split('_')[-1].replace('.h5ad', '') for file in h5ad_files]
        
        return Response({'adata_list': adata_list}, status=status.HTTP_201_CREATED)
    
class ChosenAdataResultView(APIView):
    def post(self, request):
        chosen_adata = request.data.get('chosen_adata')
        available_files = read_h5ad_file('available_files_'+chosen_adata+'.h5ad')
        adata = available_files
        n_obs, n_vars = adata.shape
        result = f"{chosen_adata}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"
        save_h5ad_file(adata, 'adata_identify_gate.h5ad')
        SaveSteps(14, 'identifythegates', 'chosen', {'chosen_adata':chosen_adata}, {})
        # save_data(chosen_adata, 'chosen_adata') 
        return Response({'result': result}, status=status.HTTP_201_CREATED)
    

class IdentifytheGatesView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_identify_gate.h5ad')
        columns = ['Marker'] + list(adata.obs['Sample'].cat.categories)
        gate_df = pd.DataFrame(columns=columns)
        gate_df['Marker'] = adata.var_names
        print(gate_df)
        # save_data(gate_df, 'gate_df') 
        SaveSteps(15, 'identifythegates', 'process', {'gate_df': json.dumps(gate_df.to_dict(orient='records'))},{})
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
        
        edit_df = pd.DataFrame(editdata)
        
        edit_df = edit_df.apply(pd.to_numeric, errors='coerce')
        edit_df = edit_df.replace([np.inf, -np.inf, np.nan], None)
        
        for col in edit_df.columns:
            if col != 'Marker': 
                edit_df[col] = edit_df[col].where(
                    (edit_df[col] >= -1e308) & (edit_df[col] <= 1e308), 
                    None
                )
        
        result = edit_df.to_dict(orient='records')
        SaveSteps(15, 'identifythegates', 'process', {'gate_df': json.dumps(result)}, {})
        
        return Response({'result': result}, status=status.HTTP_201_CREATED)
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
        SaveSteps(16, 'phenotyping', 'process', {'file_name': file_name}, {'phenotyping_result': phenotyping_result})
        return Response({'phenotyping_result': phenotyping_result}, status=status.HTTP_201_CREATED)

    def perform_phenotyping(self, file_path):
        adata = read_h5ad_file('adata_identify_gate.h5ad')
        gate_df_data = GetSteps(15)[0].input_values['gate_df']
        if isinstance(gate_df_data, str):
            gate_df_data = json.loads(gate_df_data)
        gate_df = pd.DataFrame(gate_df_data)
        chosen_adata = GetSteps(14)[0].input_values['chosen_adata']
        n_pcs = GetSteps(7)[0].output_values['n_pcs']
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
        SaveSteps(20, 'phenotyping', 'addumapphenotypes', {}, {})
        return Response({"message": "Phenotypes processed"}, status=status.HTTP_201_CREATED)

    def post_markers(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        chosen_method = request.data.getlist('markers')
        add_phenotypes_markers(adata, chosen_method)
        SaveSteps(21, 'phenotyping', 'addumapmarkers', {'chosen_method':chosen_method}, {})
        return Response({"message": "Markers processed"}, status=status.HTTP_201_CREATED)
class PreloadPhenotypesMarkersView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        marker_list = adata.var_names.tolist()
        marker_list.insert(0,"Select All") 
        print(marker_list)
        SaveSteps(17, 'phenotyping', 'preloadmarkers', {}, {'marker_list':marker_list})
        
        print(GetSteps(17)[0])
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
        n_pcs = GetSteps(7)[0].output_values['n_pcs']

        rename_dict = {}
        for _, row in rename_df.iterrows():
            current_name = str(row['CurrentName'])
            new_name = str(row['NewName'])
            rename_dict[current_name] = new_name

        adata.obs['phenotype'] = adata.obs['phenotype'].astype(str).map(rename_dict).astype('category')
        chosen_adata = GetSteps(14)[0].input_values['chosen_adata']
        phenotype_result(adata, chosen_adata, n_pcs)
        SaveSteps(18, 'phenotyping', 'rename', {}, {})
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
        chosen_adata = GetSteps(14)[0].input_values['chosen_adata']
        n_pcs = GetSteps(7)[0].output_values['n_pcs']

        phenotype_result(adata, chosen_adata, n_pcs)
        SaveSteps(19, 'phenotyping', 'drop', {}, {})
        save_h5ad_file(adata, 'adata_phenotyping.h5ad')
        return Response({'message': 'drop phenotype successfully!'}, status=status.HTTP_201_CREATED)
  
# Spatial Analysis
class PreloadSpatialAnalysisView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        chosen_adata = GetSteps(14)[0].input_values['chosen_adata']
        n_obs, n_vars = adata.shape
        return Response({'preload_result': f"{chosen_adata}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"}, status=status.HTTP_201_CREATED)
class SpatialAnalysisView(APIView):
    def post(self, request):
        self.calculate_spatial()
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        
        columns_list = [col for col in adata.obs.columns if col.startswith(("leiden_R", "phenotype"))]
        default_chosen_column = 'phenotype'
        cluster_list = adata.obs[default_chosen_column].unique()
        cluster_list_L = adata.obs['leiden_R'].unique()
        default_chosen_cluster = cluster_list[0]
        method_list = ['radius', 'knn']
        default_chosen_method = 'radius'
        filename = []
        
        DH = distances_heatmap(default_chosen_column)
        DN = distances_numeric_plot(default_chosen_column, default_chosen_cluster)
        IH = interactions_heatmap(default_chosen_column, default_chosen_method)
        IV = interactions_voronoi(default_chosen_column)
        filename.append(DH)
        filename.append(DN)
        filename.append(IH)
        filename.append(IV)
        print(adata.shape)
        SaveSteps(22, 'spatialanalysis', 'process', {}, {"filename":filename})
        return Response({'columns_list': columns_list,'cluster_list':cluster_list, 'cluster_list_L':cluster_list_L, 'method_list':method_list, "filename":filename}, status=status.HTTP_201_CREATED)
    def calculate_spatial(self):
        adata = read_h5ad_file('adata_phenotyping.h5ad')
        columns_list = [col for col in adata.obs.columns if col.startswith(("leiden_R", "phenotype"))]
        for col in columns_list: 
            sm.tl.spatial_distance(adata, x_coordinate='X_centroid', y_coordinate='Y_centroid', 
                                phenotype=col, imageid='Sample', label=f'spatial_distance_{col}')
        for col in columns_list:
            sm.tl.spatial_interaction(adata, phenotype=col, method='knn', radius=30, 
                                    imageid='Sample', label=f'spatial_interaction_{col}_knn', pval_method='zscore')
            sm.tl.spatial_interaction(adata, phenotype=col, method='radius', radius=30, 
                                    imageid='Sample', label=f'spatial_interaction_{col}_radius', pval_method='zscore')
        save_h5ad_file(adata, 'adata_spatial_analysis.h5ad')

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
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        filename = distances_heatmap(json.loads(request.body).get('dh_ul_input'))
        return filename
    def post_np(self, request):
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        filename = distances_numeric_plot(json.loads(request.body).get('np_ul_input'), json.loads(request.body).get('npd_ul_input'))
        return filename
    def post_ih(self, request):
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        filename = interactions_heatmap(json.loads(request.body).get('ih_ul_input'), json.loads(request.body).get('ihm_ul_input'))
        return filename
    def post_vp(self, request):
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        filename = interactions_voronoi(json.loads(request.body).get('vp_ul_input'))
        return filename

# Neighborhood
class PreloadNeighborView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        n_obs, n_vars = adata.shape
        chosen_adata = GetSteps(14)[0].input_values['chosen_adata']
        return Response({"preload_result":f"{chosen_adata}: AnnData object with n_obs × n_vars = {n_obs} × {n_vars}"}, status=status.HTTP_201_CREATED)
        
class NeighborView(APIView):
    def post(self, request):
        adata = read_h5ad_file('adata_spatial_analysis.h5ad')
        chosen_column = json.loads(request.body).get('n_input')
        k = json.loads(request.body).get('k_neighbor')
        n_neighborhoods =json.loads(request.body).get('n_neighbor')
        perform_neighborhoods(adata, chosen_column, k, n_neighborhoods)
        SaveSteps(23, 'neighborhood', 'process', {'chosen_column':chosen_column,'k':k,'n_neighborhoods':n_neighborhoods}, {})
        print(True)
        return Response({"test":"test"}, status=status.HTTP_201_CREATED)
