import redis
import pickle
from django.conf import settings
import os
import anndata as ad

r = None
def get_redis_connection():
    global r
    if r is None:
        try:
            r = redis.StrictRedis()
        except redis.exceptions.ConnectionError as e:
            print("連線 Redis 發生錯誤:", e)
    return r

def save_adata_objects(adata_objects, key):
    print(key, ':\n', adata_objects)
    r = get_redis_connection()
    if r:
        adata_bytes = pickle.dumps(adata_objects)
        r.set(key, adata_bytes)

def load_adata_objects(key):
    r = get_redis_connection()
    if r:
        adata_bytes = r.get(key)
        if adata_bytes:
            return pickle.loads(adata_bytes)
    return None



def save_h5ad_file(adata, filename):
    os.makedirs(settings.H5AD_STORAGE_PATH, exist_ok=True)
    file_path = os.path.join(settings.H5AD_STORAGE_PATH, filename)
    adata.write_h5ad(file_path)

def read_h5ad_file(filename):
    file_path = os.path.join(settings.H5AD_STORAGE_PATH, filename)
    return ad.read_h5ad(file_path)