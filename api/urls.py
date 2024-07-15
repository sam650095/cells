from django.urls import path
from . import views

urlpatterns = [
    path('clear', views.ClearAllDataView.as_view(),name='clear data'),

    path('upload', views.FileUploadView.as_view(), name='file-upload'),
    path('qualitycontrol', views.QualityControlView.as_view(), name='qualicty-control'),
    path('preview', views.PreviewView.as_view(), name="preview-filter"),
    path('replace', views.ReplaceView.as_view(), name="replace-adata"),
    path('confirm', views.ConfirmView.as_view(), name="confirm-filter"),
    path('normal', views.NormalizationView.as_view(), name="normalize"),
    path('merge', views.MergeView.as_view(), name="merged-data"),

    path('preloadpca', views.PreloadPCAView.as_view(), name="PreloadPCA"),
    path('pca', views.PCAView.as_view(), name="PCA"),
    path('preloadclustering', views.PreloadCLusteringView.as_view(), name="PreloadClustering"),
    path('clustering', views.CLusteringView.as_view(), name="Clustering"),
    path('preloadmarkers', views.PreloadMarkersView.as_view(), name="PreloadMarkers"),
    path('addumapcluster/<str:met>', views.AddUmapClusterView.as_view(), name="addingmet"),
    path('grabnames', views.GrabClusterNameView.as_view(), name="grabname"),
    path('rename', views.RenameClusterView.as_view(), name="rename"),
    path('grabclusters', views.GrabClustersView.as_view(), name="grabclusters"),
    path('subclusters', views.SubclusterView.as_view(), name="subclusters"),
    path('preloadsubset', views.PreloadSubsetView.as_view(), name="preloadsubset"),
    path('grabclustersubset', views.GrabClusterSubsetView.as_view(), name="grabclustersubset"),
    path('subset', views.SubsetView.as_view(), name="subset"),

    path('preloadidentifythegates', views.PreloadIdentifytheGatesView.as_view(), name="preloadidentifythegates"),
    path('choseadata', views.ChosenAdataResultView.as_view(),name="choseadata"),
    path('identifythegates', views.IdentifytheGatesView.as_view(), name="identifythegates"),
    path('phenotyping', views.PhenotypingView.as_view(), name="phenotyping")
]