# Optimal Trimming Pipeline (Human mtDNA Sanger Sequencing)

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/Platform-Windows-lightgrey.svg)]()
[![Sequencher](https://img.shields.io/badge/Integration-Sequencher%20v5.4.6-orange.svg)]()

Hệ thống phân tích và dự đoán vị trí cắt lọc nhiễu tự động (**Non-destructive Noise Boundary Prediction**) cho dữ liệu giải trình tự Sanger mao quản hệ gen ty thể người (**Human mtDNA D-loop: HV1 & HV2/HV3**).

Được tối ưu hóa chuyên sâu cho thiết bị giải trình tự **Applied Biosystems SeqStudio Flex 24** (hóa chất **BigDye Terminator Plus**, mao quản **50 cm**, polymer **POP-7**) và tương thích hoàn toàn với quy trình làm việc trên **Sequencher v5.4.6**.

---

## 🔬 Điểm Nổi Bật & Nguyên Lý Sinh Học Phân Tử

1. **Bảo toàn dữ liệu thô (Non-destructive Approach):**
   - Không can thiệp, không cắt xén cấu trúc file `.ab1` gốc để đảm bảo tính toàn vẹn của điện di đồ huỳnh quang (Electropherogram trace).
   - Dự đoán chính xác ranh giới nhiễu đầu 5' (loại bỏ primer, nhiễu mao quản sớm 10–20 bp) và đầu 3' (sụt giảm huỳnh quang).

2. **Bảo toàn dải sau trượt Poly-C (Homopolymer C-Tract):**
   - Tại vùng trượt Poly-C ở HV1 (16184 - 16193) và HV2 (303 - 315), hiện tượng trượt pha (Stutter) thường gây phân đôi đỉnh sóng.
   - Hệ thống **bảo toàn nguyên vẹn dải sau Poly-C** để làm dữ liệu đối sánh chéo giữa chiều đọc mồi xuôi (Forward) và mồi ngược (Reverse), ưu tiên chiều đọc đối diện sạch để tạo consensus.

3. * Quyết định QC nghiêm ngặt & Nhãn `DENY`:**
   - Tự động phát hiện mẫu suy thoái, chập peak quang học (`?`), hoặc đứt đoạn lớn.
   - Gán nhãn **`DENY` (Khuyến cáo phòng lab làm lại mẫu)** đối với các mẫu không phủ kín >= 85% dải mục tiêu hoặc nghi ngờ nhiễm tạp chéo (>= 20 SNPs).

---

## 📊 Hai Dạng Báo Cáo Phối Hợp (A & C)

- **Dạng A — Bảng Khuyến Nghị Vị Trí Cắt (`trimming_recommendations.csv`):**
  Cung cấp bảng tra cứu cho kỹ thuật viên nhập trực tiếp vào Sequencher:
  `Suggested_5p_Cut | Suggested_3p_Cut | Retained_Length_bp | QC_Status | Lab_Recommendation | Sequencher_Guidance`

- **Dạng C — Điện Di Đồ Tương Tác HTML (`html_reports/*_chromatogram.html`):**
  - Trực quan hóa 4 kênh huỳnh quang thực tế (A: Xanh lá, C: Xanh dương, G: Vàng, T: Đỏ).
  - Vùng cắt 5' và 3' được phủ bóng mờ; dải Poly-C được highlight viền vàng.
  - Chuỗi chuẩn tham chiếu **rCRS (NC_012920.1)** chạy song song dưới từng nucleotide với ký hiệu trùng khớp (`.`) và biến dị SNPs.
  - Banner cảnh báo màu đỏ nổi bật đối với các mẫu bị gán nhãn `DENY`.

---

## 📁 Cấu Trúc Dự Án

```text
Optimal_Trimming_Pipeline/
├── core/
│   ├── ab1_parser.py       # Trích xuất DATA9-12, PLOC, PCON từ .ab1
│   ├── tracy_engine.py     # Nhận diện Clean, Heteroplasmy, Conflict (?), Poly-C, Stutter
│   ├── aligner.py          # Căn gióng rCRS (16024-16365 và 73-576)
│   ├── trimmer.py          # Thuật toán dự đoán ranh giới 5'/3' không phá hủy
│   ├── assembler.py        # Ráp contig 2 chiều & cây quyết định QC (PASS / REVIEW / DENY)
│   └── reporter.py         # Xuất Form A (CSV) và Form C (HTML Chromatogram Viewer)
├── data/
│   ├── rCRS.fasta          # Chuỗi tham chiếu chuẩn quốc tế rCRS (16,569 bp)
│   └── sample_ab1/         # Dữ liệu thực tế 24 file .ab1 từ SeqStudio Flex 24
├── tests/
│   └── test_pipeline.py    # Bộ kiểm thử tự động (Unit Test Suite 6/6 PASS)
├── app_gui.py              # Giao diện đồ họa Desktop (Tkinter)
├── run_cli.py              # Công cụ dòng lệnh xử lý hàng loạt tốc độ cao
├── 1_CAI_DAT_THU_VIEN.bat  # Cài đặt thư viện tự động
├── 2_CHAY_GIAO_DIEN_GUI.bat# Chạy ứng dụng giao diện
├── 3_CHAY_DONG_LENH_CLI.bat# Chạy xử lý dòng lệnh
├── 4_DONG_GOI_EXE.bat      # Script đóng gói thành .EXE độc lập (PyInstaller)
├── requirements.txt        # Danh sách thư viện phụ thuộc
└── HUONG_DAN_SU_DUNG.txt   # Hướng dẫn chi tiết cho phòng lab
```

---

## 🚀 Hướng Dẫn Sử Dụng

### Lựa chọn 1: Sử dụng Giao diện đồ họa (GUI)
1. Chạy file `2_CHAY_GIAO_DIEN_GUI.bat` (hoặc `python app_gui.py`).
2. Chọn thư mục chứa file `.ab1` cần phân tích.
3. Bấm **▶ BẮT ĐẦU PHÂN TÍCH HÀNG LOẠT (1-CLICK RUN)**.
4. Bấm **Mở Bảng Khuyến Nghị Cắt (CSV)** để mở ngay kết quả trong Excel.

### Lựa chọn 2: Sử dụng Dòng lệnh (CLI)
```bash
python run_cli.py --input "data/sample_ab1/real data" --rcrs "data/rCRS.fasta" --output "output"
```

### Lựa chọn 3: Chạy Kiểm thử tự động
```bash
python tests/test_pipeline.py
```

---

## 📋 Tiêu Chuẩn Kiểm Soát Chất Lượng (QC Matrix)

| Trạng Thái QC | Tỷ Lệ Bao Phủ | Mâu Thuẫn (+) | Chập Peak (?) | Số Lượng SNPs | Khuyến Cáo Phòng Lab |
| :---: | :---: | :---: | :---: | :---: | :--- |
| **PASS** | >= 98% | 0 | 0 | < 20 | Đạt chuẩn chất lượng cao, xuất kết quả ngay. |
| **REVIEW** | 88% - 98% | <= 2 | <= 2 | < 20 | Cần chuyên viên soi lại biểu đồ điện di đồ (HTML). |
| **DENY** | < 85% | >= 3 | >= 5 | >= 20 | **Không phân tích — Khuyến cáo làm lại mẫu trong phòng lab!** |


