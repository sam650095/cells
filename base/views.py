from django.shortcuts import render
from django.shortcuts import render
from django.http import JsonResponse
from django.conf import settings
import os

def index(request):
    return render(request, 'base/index.html')
def page(request, process, method):
    context = {'process': process, 'method': method}
    return render(request, f"{process}/{method}.html", context)
def gui(request):
    context = {}
    return render(request, '/vkt/gui.html', context)
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