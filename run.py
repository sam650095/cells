import subprocess
import sys
import os
import signal
import platform

def kill_port(port):
    if platform.system() == "Windows":
        try:
            output = subprocess.check_output(f"netstat -ano | findstr :{port}", shell=True)
            if output:
                pid = output.split()[-1]
                os.system(f"taskkill /PID {pid} /F")
                print(f"Killed process on port {port}")
        except subprocess.CalledProcessError:
            print(f"No process found on port {port}")
    else:
        try:
            output = subprocess.check_output(["lsof", "-ti", f":{port}"])
            pids = output.decode().split()
            for pid in pids:
                os.kill(int(pid), signal.SIGKILL)
            print(f"Killed process on port {port}")
        except subprocess.CalledProcessError:
            print(f"No process found on port {port}")

def run_servers():
    kill_port(6379)

    if platform.system() == "Windows":
        redis = subprocess.Popen(["wsl", "redis-server"])
    else:
        redis = subprocess.Popen(['redis-server'])

    django = subprocess.Popen([sys.executable, 'manage.py', 'runserver'])
    
    try:
        django.wait()
    finally:
        if platform.system() == "Windows":
            subprocess.run(["wsl", "pkill", "redis-server"])
        else:
            redis.terminate()
        redis.wait()

if __name__ == '__main__':
    run_servers()