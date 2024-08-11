# Django Project

### Url

1. Frontend Slides: https://docs.google.com/presentation/d/1e6EwQYKyvkKRE1Z1Q--5is8NDSv-8fH3/edit#slide=id.p23
2. Backend Codes: https://github.com/sherry-1230/Single-Cell-Analysis/blob/main/%24final.html

## Versions

- Python 3.8

## Installations

1. Create Virtual Environment
   ##### miniconda3
   https://docs.anaconda.com/miniconda/

   ```sh
   conda create --name CellsProjectenv python=3.8
   conda activate CellsProjectenv
   ```

2. Install Requirements:

   ### Install Requirements

   ```sh
   pip install -r requirements.txt
   ```

   ### You mighht need before install requirements.txt

   - In MacOS

   ```sh
   export PYTHONUTF8=1
   ```

   - In Windows

   ```sh
   $env:PYTHONUTF8 = 1
   ```

3. Redis:

   - MacOS

     ```sh
     redis-server
     ```

     ##### # 如果看到這個 port 有被使用

     ```sh
     sudo lsof -i :6379 # 查看port 6379的使用情形
     sudo kill -9 <PID> # 將上面指令看到的PID一欄的數字填上
     sudo systemctl status redis
     sudo systemctl stop redis
     ```

   - Windows

     ##### redis 不支援 Windows，因此只能在 Linux 執行

     1. 進入 Linux

        快速指南 - https://learn.microsoft.com/en-us/windows/wsl/install

        ```sh
        WSL # 進入Linux環境
        WSL --install # 下載Linux
        ```

     2. 下載 Rendis

        快速指南 - https://redis.io/docs/latest/operate/oss_and_stack/install/install-redis/install-redis-on-windows/

        ```sh
        curl -fsSL https://packages.redis.io/gpg | sudo gpg --dearmor -o /usr/share/keyrings/redis-archive-keyring.gpg
        echo "deb [signed-by=/usr/share/keyrings/redis-archive-keyring.gpg] https://packages.redis.io/deb $(lsb_release -cs) main" | sudo tee /etc/apt/sources.list.d/redis.list
        sudo apt-get update
        sudo apt-get install redis
        ```

4. 執行專案

   - 這兩個必須同時啟動
     啟動 redis & 啟動 Django

   ```
   sudo service redis-server start
   python manage.py runserver
   ```