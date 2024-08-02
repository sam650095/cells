from django.urls import path
from . import views
urlpatterns = [
    path('', views.index, name='index'),
    path('<str:process>/<str:method>',views.page, name='page'),
    path('vkt/gui', views.gui, name='gui'),
    path('get_image/', views.get_image, name='get_image'),
]