from django.urls import path
from . import views

urlpatterns = [
    path('upload', views.FileUploadView.as_view(), name='file-upload'),
    path('qualitycontrol', views.QualityControlView.as_view(), name='qualicty-control'),
    path('preview', views.PreviewView.as_view(), name="preview-filter"),
    path('replace', views.ReplaceView.as_view(), name="replace-adata"),
    path('confirm', views.ConfirmView.as_view(), name="confirm-filter"),
    path('normal', views.NormalizationView.as_view(), name="normalize"),
    path('merge', views.MergeView.as_view(), name="merged-data"),

    path('preload', views.PreloadView.as_view(), name="Preload"),
    path('pca', views.PCAView.as_view(), name="PCA")
]