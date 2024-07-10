# Django Project

## 需求

- Python 3.8

## 安裝

1. 虛擬環境：

   ```sh
   conda create --name CellsProjectenv python=3.8
   conda activate CellsProjectenv
   ```

2. 安裝必要套件:

   ### 只能執行到 00 Preprocessing

   ```sh
   pip install -r requirements.txt
   ```

   ### 如果下載有誤，請先打上

   - MacOS

   ```sh
   export PYTHONUTF8=1
   ```

   - Windows

   ```sh
   $env:PYTHONUTF8 = 1
   ```

3. Redis:

   - MacOS

     ##### 如果先前 install requiremnents.txt 是正常，那只要執行以下就不會有意外。

     ```sh
     redis-server
     ```

     ##### # 如果看到這個 port 有被使用

     ```sh
     lsof -i :6379 # 查看port 6379的使用情形
     kill -9 <PID> # 將上面指令看到的PID一欄的數字填上
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

     3. 啟動 redis

        ```sh
        sudo service redis-server start # 沒有報錯應該就沒問題
        ```

4. 執行專案

   ```sh
   python manage.py runserver # 開啟 http://localhost:8000
   ```

- 這兩個必須同時啟動

  ```
  sudo service redis-server start
  python manage.py runserver
  ```
