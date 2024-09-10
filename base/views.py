from django.shortcuts import render
from django.http import JsonResponse
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from django.conf import settings
import os
from .models import *

from api.saved import *

def index(request):
    return render(request, 'base/index.html')
def page(request, process, method):
    context = {'process': process, 'method': method}
    return render(request, f"{process}/{method}.html", context)
def createadata(request):
    if request.method == 'POST':
        file_name = request.POST.get('file_names')
        UploadData.objects.create(file_name=file_name)
    return True

def get_image(request):
    if request.method == 'POST':
        folder = request.POST.get('folder', '')
        filename = request.POST.get('filename', '')
        
        image_path = os.path.normpath(os.path.join(folder, filename))
        if image_path.startswith('..'):
            return JsonResponse({'error': 'Invalid path'}, status=400)
        
        full_path = os.path.join(settings.MEDIA_ROOT, image_path)
        print(full_path)
        if os.path.exists(full_path) and os.path.isfile(full_path):
            return JsonResponse({'image_path': os.path.join(settings.MEDIA_URL, image_path)})
        else:
            return JsonResponse({'error': 'File not found'}, status=404)
    
    return JsonResponse({'error': 'Invalid request method'}, status=400)
def download_images(request):
    folder = request.POST.get('folder', '')
    print(folder)
    folder_path = os.path.normpath(os.path.join(settings.MEDIA_ROOT, folder))
    if not folder_path.startswith(settings.MEDIA_ROOT):
        return JsonResponse({'error': 'Invalid path'}, status=400)
    
    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        image_files = [f for f in os.listdir(folder_path) if f.lower().endswith(('.png', '.jpg', '.jpeg', '.gif'))]
        image_paths = [os.path.join(settings.MEDIA_URL, folder, f) for f in image_files]
        return JsonResponse({'image_paths': image_paths})
    else:
        return JsonResponse({'error': 'Folder not found'}, status=404)
