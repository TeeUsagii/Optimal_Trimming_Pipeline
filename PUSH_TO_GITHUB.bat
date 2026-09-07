@echo off
chcp 65001 > nul
echo ============================================================
echo   OPTIMAL TRIMMING PIPELINE - ĐẨY DỰ ÁN LÊN GITHUB
echo ============================================================
echo.
echo Đang đồng bộ và đẩy mã nguồn lên nhánh main:
echo Remote: https://github.com/Quannguyen513/Optimal_Trimming_Pipeline.git
echo.
git push -u origin main
echo.
if %errorlevel% equ 0 (
    echo [THÀNH CÔNG] Toàn bộ dự án đã được đẩy lên GitHub!
    echo Xem tại: https://github.com/Quannguyen513/Optimal_Trimming_Pipeline
) else (
    echo [CHÚ Ý] Nếu gặp lỗi "Repository not found":
    echo Hãy đảm bảo bạn đã tạo repository tên "Optimal_Trimming_Pipeline"
    echo trên trang: https://github.com/new (tài khoản Quannguyen513).
)
echo.
pause
