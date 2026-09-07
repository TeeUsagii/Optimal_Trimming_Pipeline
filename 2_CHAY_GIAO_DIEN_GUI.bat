@echo off
chcp 65001 > nul
echo ============================================================
echo   OPTIMAL TRIMMING PIPELINE - KHỞI ĐỘNG GIAO DIỆN (GUI)
echo ============================================================
echo.

python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [LOI] Khong tim thay Python tren he thong!
    echo Vui long chay file '1_CAI_DAT_THU_VIEN.bat' truoc de kiem tra.
    echo.
    pause
    exit /b 1
)

echo Dang mo giao dien phan mem...
start "" pythonw app_gui.py
if %errorlevel% neq 0 (
    python app_gui.py
)
