from django.urls import path
from . import views
urlpatterns = [
    path('', views.index, name='index'),
    path('<str:process>/<str:method>',views.page, name='page'),
    path('steps/<str:met>', views.StepView.as_view(), name='steps'),
    path('get_image/', views.get_image, name='get_image'),
    path('download_image/', views.download_images, name='download_image'),
]