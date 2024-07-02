# Django Project

## 簡介

這是一個使用 Django 框架開發的網頁應用程式範例。此項目展示了 Django 的基本功能，包括用戶驗證、資料庫操作和模板渲染。

## 特性

- 使用 Django 開發
- 用戶註冊和登入功能
- 資料庫操作（CRUD）
- 模板渲染
- 管理介面

## 需求

- Python 3.x
- Django 3.x 或更新版本
- SQLite（默認資料庫）

## 安裝

1. 克隆此存儲庫：

   ```sh
   git clone https://github.com/yourusername/yourproject.git
   ```

2. 進入項目目錄：

   ```sh
   cd yourproject
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

7. 創建超級用戶（管理員帳號）：

   ```sh
   python manage.py createsuperuser
   ```

8. 啟動開發伺服器：

   ```sh
   python manage.py runserver
   ```

9. 打開瀏覽器，訪問 [http://127.0.0.1:8000/](http://127.0.0.1:8000/) 查看應用程式。

## 使用說明

- 訪問 [http://127.0.0.1:8000/admin/](http://127.0.0.1:8000/admin/) 進入 Django 管理介面。
- 在管理介面中，您可以管理用戶和其他模型數據。

## 貢獻

歡迎任何形式的貢獻！您可以透過提交問題（issues）或拉取請求（pull requests）來貢獻。

## 許可

此項目採用 MIT 許可。詳情請參閱 `LICENSE` 文件。
