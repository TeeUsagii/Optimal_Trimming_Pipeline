@echo off
chcp 65001 > nul
echo ============================================================
echo   mtDNA SANGER ANALYSIS SUITE - CAI DAT THU VIEN
echo ============================================================
echo.
echo Dang kiem tra Python va cai dat cac thu vien can thiet...
echo.

python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [LOI] Khong tim thay Python tren may nay!
    echo Vui long cai dat Python (khuyen dung Python 3.11) tu:
    echo   https://www.python.org/downloads/
    echo.
    echo * LUU Y QUAN TRONG: Khi cai dat, nho tich vao o:
    echo   [x] "Add Python to PATH" (hoac Add python.exe to PATH)
    echo.
    pause
    exit /b 1
)

echo [+] Da tim thay Python! Dang cai dat biopython, pandas, numpy...
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

if %errorlevel% equ 0 (
    echo.
    echo ============================================================
    echo [THANH CONG] Da cai dat xong tat ca thu vien!
    echo Bay gio ban chi can click dup file '2_CHAY_GIAO_DIEN_GUI.bat'.
    echo ============================================================
) else (
    echo.
    echo [LOI] Co loi trong qua trinh cai dat. Vui long kiem tra ket noi Internet.
)

echo.
pause
