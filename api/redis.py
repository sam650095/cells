import redis
import pickle
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
    print(key, adata_objects)
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