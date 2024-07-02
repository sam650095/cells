# Django Project

## 需求

- Python 3.8

## 安裝

1. 虛擬環境：

   ```sh
   conda create --name CellsProjectenv python=3.8
   conda activate CellsProjectenv
   ```

2. 安裝:

   ```sh
   pip install -r requirements.txt
   # 只能執行到00 Preprocessing
   ```

3. 建立虛擬環境：

   ```sh
   python -m venv env
   ```

4. 啟動虛擬環境：

   - Windows：

     ```sh
     .\env\Scripts\activate
     ```

   - MacOS/Linux：

     ```sh
     source env/bin/activate
     ```

5. 安裝所需的套件：

   ```sh
   pip install -r requirements.txt
   ```

6. 進行數據庫遷移：

   ```sh
   python manage.py migrate
   ```

7. 啟動開發伺服器：

   ```sh
   python manage.py runserver
   ```

8. 打開瀏覽器，訪問 [http://127.0.0.1:8000/](http://127.0.0.1:8000/) 查看應用程式。
