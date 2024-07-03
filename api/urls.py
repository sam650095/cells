from django.urls import path
from . import views

urlpatterns = [
    path('upload', views.FileUploadView.as_view(), name='file-upload'),
    path('qualitycontrol',views.QualityControlView.as_view(), name='qualicty-control'),
    path('preview', views.PreviewView.as_view(), name="preview-filter"),
    path('replaceimg', views.ReplaceView.as_view(), name="replace-image"),
    path('confirm', views.ConfirmView.as_view(), name="confirm-filter")
]