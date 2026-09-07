@echo off
chcp 65001 > nul
echo ============================================================
echo   OPTIMAL TRIMMING PIPELINE - CHẠY DÒNG LỆNH (CLI)
echo ============================================================
echo.

set PYTHONIOENCODING=utf-8
python run_cli.py --input "data/sample_ab1/real data" --rcrs "data/rCRS.fasta" --output "output_real"

echo.
echo ============================================================
echo Da hoan thanh! Ban co the kiem tra thu muc output_real.
echo ============================================================
pause
