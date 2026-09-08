@echo off
chcp 65001 > nul
echo ============================================================
echo   OPTIMAL TRIMMING PIPELINE - ĐÓNG GÓI THÀNH FILE .EXE
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
echo [2/4] Đang đóng gói Giao diện đồ họa (app_gui.py -> Optimal_Trimming_Pipeline.exe)...
python -m PyInstaller --noconfirm --onedir --windowed --name "Optimal_Trimming_Pipeline" --icon "assets\icon.ico" --add-data "data;data" --add-data "assets;assets" --collect-all customtkinter app_gui.py

echo.
echo [3/4] Đang đóng gói Dòng lệnh CLI (run_cli.py -> Optimal_Trimming_CLI.exe)...
python -m PyInstaller --noconfirm --onedir --console --name "Optimal_Trimming_CLI" --icon "assets\icon.ico" --add-data "data;data" run_cli.py
copy /Y dist\Optimal_Trimming_CLI\Optimal_Trimming_CLI.exe dist\Optimal_Trimming_Pipeline\ > nul
rd /S /Q dist\Optimal_Trimming_CLI > nul 2>&1

echo.
echo [4/4] Đang sao chép thư mục dữ liệu mẫu, biểu tượng và tài liệu vào bản phân phối...
xcopy /E /I /Y data dist\Optimal_Trimming_Pipeline\data > nul
xcopy /E /I /Y assets dist\Optimal_Trimming_Pipeline\assets > nul
copy /Y HUONG_DAN_SU_DUNG.txt dist\Optimal_Trimming_Pipeline\ > nul
copy /Y README.md dist\Optimal_Trimming_Pipeline\ > nul
echo @echo off > dist\Optimal_Trimming_Pipeline\CHAY_GIAO_DIEN.bat
echo start "" "%%~dp0Optimal_Trimming_Pipeline.exe" >> dist\Optimal_Trimming_Pipeline\CHAY_GIAO_DIEN.bat
echo @echo off > dist\Optimal_Trimming_Pipeline\CHAY_CLI.bat
echo "%%~dp0Optimal_Trimming_CLI.exe" --input "data\sample_ab1\real data" --rcrs "data\rCRS.fasta" --output "output_real" >> dist\Optimal_Trimming_Pipeline\CHAY_CLI.bat
echo pause >> dist\Optimal_Trimming_Pipeline\CHAY_CLI.bat

echo.
echo ============================================================
echo [THÀNH CÔNG] ĐÃ ĐÓNG GÓI HOÀN CHỈNH BẢN PORTABLE STANDALONE!
echo Thư mục ứng dụng độc lập:
echo   --> dist\Optimal_Trimming_Pipeline\Optimal_Trimming_Pipeline.exe
echo.
echo Bạn có thể copy toàn bộ thư mục 'dist\Optimal_Trimming_Pipeline' sang
echo bất kỳ máy tính Windows nào khác để chạy ngay mà KHÔNG CẦN
echo cài đặt Python hay bất kỳ thư viện nào!
echo ============================================================
echo.
pause
