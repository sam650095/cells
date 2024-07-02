from django.shortcuts import render

def index(request):
    return render(request, 'base/index.html')
def page(request, process, method):
    context = {'process': process, 'method': method}
    return render(request, f"{process}/{method}.html", context)