import redis
import pickle
from django.conf import settings
import os
import anndata as ad
import shutil

r = None
def get_redis_connection():
    global r
    if r is None:
        try:
            r = redis.StrictRedis()
        except redis.exceptions.ConnectionError as e:
            print("連線 Redis 發生錯誤:", e)
    return r

def save_data(data, key):
    print(key, ':\n', data)
    r = get_redis_connection()
    if r:
        data_bytes = pickle.dumps(data)
        r.set(key, data_bytes)

def load_data(key):
    r = get_redis_connection()
    if r:
        data_bytes = r.get(key)
        if data_bytes:
            return pickle.loads(data_bytes)
    return None

def clear_data(key):
    r = get_redis_connection()
    if r:
        r.delete(key)
    else:
        print("無法連接到 Redis")


def get_all_keys():
    r = get_redis_connection()
    if r:
        keys = r.keys('*')
        result = {}
        for key in keys:
            decoded_key = key.decode('utf-8')
            data = load_data(decoded_key)
            result[decoded_key] = data
            print(f"Key: {decoded_key}")
            print(f"Value: {data}")
            print("-" * 50)  

def clear_all_data():
    r = get_redis_connection()
    if r:
        keys = r.keys('*')
        for key in keys:
            r.delete(key)
        print("所有資料已清除")
    else:
        print("無法連接到 Redis")


def save_h5ad_file(adata, filename):
    os.makedirs(settings.H5AD_STORAGE_PATH, exist_ok=True)
    file_path = os.path.join(settings.H5AD_STORAGE_PATH, filename)
    adata.write_h5ad(file_path)

def read_h5ad_file(filename):
    file_path = os.path.join(settings.H5AD_STORAGE_PATH, filename)
    return ad.read_h5ad(file_path)
def clear_all_h5ad_files():
    # h5ad_files = [f for f in os.listdir(settings.H5AD_STORAGE_PATH) if f.endswith('.h5ad')]
    # [os.remove(os.path.join(settings.H5AD_STORAGE_PATH, file)) for file in h5ad_files]
    shutil.rmtree(settings.H5AD_STORAGE_PATH)
def clearmediafiles(temp):
    media_root = settings.MEDIA_ROOT
    temp_dir = os.path.join(media_root, temp)
    messages = []

    if os.path.exists(temp_dir):
        for item in os.listdir(temp_dir):
            item_path = os.path.join(temp_dir, item)
            try:
                if os.path.isfile(item_path) or os.path.islink(item_path):
                    os.unlink(item_path)
                elif os.path.isdir(item_path):
                    shutil.rmtree(item_path)
                messages.append(f'Successfully deleted {item_path}')
            except Exception as e:
                messages.append(f'Failed to delete {item_path}. Reason: {e}')
    else:
        messages.append(f'Directory {temp_dir} does not exist')

    return '\n'.join(messages)  