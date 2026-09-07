@echo off
chcp 65001 > nul
echo ============================================================
echo   mtDNA SANGER ANALYSIS SUITE - ĐÓNG GÓI THÀNH FILE .EXE
echo ============================================================
echo.

python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [LỖI] Không tìm thấy Python trên máy tính.
    echo Vui lòng cài đặt Python và tích chọn "Add Python to PATH".
    pause
    exit /b 1
)

echo [1/3] Kiểm tra công cụ PyInstaller...
python -c "import PyInstaller" >nul 2>&1
if %errorlevel% neq 0 (
    echo [+] Đang cài đặt PyInstaller...
    python -m pip install pyinstaller
)

echo.
echo [2/4] Đang đóng gói Giao diện đồ họa (app_gui.py -> mtDNA_Sanger_Suite.exe)...
python -m PyInstaller --noconfirm --onedir --windowed --name "mtDNA_Sanger_Suite" --add-data "data;data" app_gui.py

echo.
echo [3/4] Đang đóng gói Dòng lệnh CLI (run_cli.py -> mtDNA_CLI.exe)...
python -m PyInstaller --noconfirm --onedir --console --name "mtDNA_CLI" run_cli.py
copy /Y dist\mtDNA_CLI\mtDNA_CLI.exe dist\mtDNA_Sanger_Suite\ > nul
rd /S /Q dist\mtDNA_CLI > nul 2>&1

echo.
echo [4/4] Đang sao chép thư mục dữ liệu mẫu và tài liệu vào bản phân phối...
xcopy /E /I /Y data dist\mtDNA_Sanger_Suite\data > nul
copy /Y HUONG_DAN_SU_DUNG.txt dist\mtDNA_Sanger_Suite\ > nul
echo start "" mtDNA_Sanger_Suite.exe > dist\mtDNA_Sanger_Suite\CHAY_GIAO_DIEN.bat
echo mtDNA_CLI.exe --input "data\sample_ab1\real data" --rcrs "data\rCRS.fasta" --output "output_real" > dist\mtDNA_Sanger_Suite\CHAY_CLI.bat
echo pause >> dist\mtDNA_Sanger_Suite\CHAY_CLI.bat

echo.
echo ============================================================
echo [THÀNH CÔNG] ĐÃ ĐÓNG GÓI HOÀN CHỈNH BẢN PORTABLE STANDALONE!
echo Thư mục ứng dụng độc lập:
echo   --> dist\mtDNA_Sanger_Suite\mtDNA_Sanger_Suite.exe
echo.
echo Bạn có thể copy toàn bộ thư mục 'dist\mtDNA_Sanger_Suite' sang
echo bất kỳ máy tính Windows nào khác để chạy ngay mà KHÔNG CẦN
echo cài đặt Python hay bất kỳ thư viện nào!
echo ============================================================
echo.
pause
