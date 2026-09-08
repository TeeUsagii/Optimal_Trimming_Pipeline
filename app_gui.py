import os
import sys
import glob
import re
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Dict, Any, List, Optional

# UTF-8 stdout/stderr
if sys.stdout and hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
if sys.stderr and hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8')

# Xác định đường dẫn tương thích cả chế độ script và file .exe đóng gói bởi PyInstaller
if getattr(sys, 'frozen', False):
    EXE_DIR = os.path.dirname(sys.executable)
    INTERNAL_DIR = getattr(sys, '_MEIPASS', EXE_DIR)
else:
    EXE_DIR = os.path.dirname(os.path.abspath(__file__))
    INTERNAL_DIR = EXE_DIR

def find_app_path(rel_path: str) -> str:
    """Tìm đường dẫn tài nguyên theo thứ tự ưu tiên: bên cạnh exe -> trong _internal -> trong _MEIPASS -> relative"""
    candidates = [
        os.path.join(EXE_DIR, rel_path),
        os.path.join(EXE_DIR, "_internal", rel_path),
        os.path.join(INTERNAL_DIR, rel_path),
        os.path.abspath(rel_path)
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return candidates[0]

if INTERNAL_DIR not in sys.path:
    sys.path.insert(0, INTERNAL_DIR)
if EXE_DIR not in sys.path:
    sys.path.insert(0, EXE_DIR)

import customtkinter as ctk
from run_cli import run_pipeline, detect_sample_info

# Thiết lập theme mặc định: Dark Slate Minimalist
ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")


class ConsoleWindow(ctk.CTkToplevel):
    """Cửa sổ Console Log chuyên dụng (ẩn mặc định, mở khi bấm F12 hoặc click nút Xem Log)"""
    def __init__(self, parent):
        super().__init__(parent)
        self.title("Terminal Console - Optimal Trimming Pipeline (Nhấn F12 để đóng/mở)")
        self.geometry("880x520")
        self.minsize(680, 360)
        self.parent = parent

        # Set icon nếu có
        icon_path = find_app_path(os.path.join("assets", "icon.ico"))
        if os.path.exists(icon_path):
            try:
                self.iconbitmap(icon_path)
            except Exception:
                pass

        self._setup_ui()
        self.protocol("WM_DELETE_WINDOW", self.withdraw)
        self.bind("<F12>", lambda e: self.withdraw())

    def _setup_ui(self):
        # Header Toolbar
        head = ctk.CTkFrame(self, fg_color=("white", "#0f172a"), height=48, corner_radius=0)
        head.pack(fill="x", side="top")
        head.pack_propagate(False)

        ctk.CTkLabel(
            head,
            text="💻 NHẬT KÝ CHI TIẾT DÒNG LỆNH (DEVELOPER CONSOLE)",
            font=ctk.CTkFont(family="Segoe UI", size=12, weight="bold"),
            text_color=("#0284c7", "#38bdf8")
        ).pack(side="left", padx=16)

        ctk.CTkButton(
            head,
            text="Đóng (F12)",
            width=80,
            height=28,
            font=ctk.CTkFont(size=11),
            fg_color=("#e2e8f0", "#334155"),
            text_color=("#0f172a", "#f8fafc"),
            hover_color=("#cbd5e1", "#475569"),
            command=self.withdraw
        ).pack(side="right", padx=12)

        ctk.CTkButton(
            head,
            text="Sao Chép",
            width=70,
            height=28,
            font=ctk.CTkFont(size=11),
            fg_color=("#e2e8f0", "#334155"),
            text_color=("#0f172a", "#f8fafc"),
            hover_color=("#cbd5e1", "#475569"),
            command=self._copy_log
        ).pack(side="right", padx=4)

        ctk.CTkButton(
            head,
            text="Xóa Log",
            width=70,
            height=28,
            font=ctk.CTkFont(size=11),
            fg_color=("#e2e8f0", "#334155"),
            text_color=("#0f172a", "#f8fafc"),
            hover_color=("#cbd5e1", "#475569"),
            command=self._clear_log
        ).pack(side="right", padx=4)

        # Log Box
        self.log_text = ctk.CTkTextbox(
            self,
            corner_radius=0,
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color=("#f8fafc", "#0b0f19"),
            text_color=("#0f172a", "#e2e8f0"),
            wrap="char"
        )
        self.log_text.pack(fill="both", expand=True)

    def append_log(self, text: str):
        self.log_text.insert("end", text + "\n")
        self.log_text.see("end")

    def _clear_log(self):
        self.log_text.delete("1.0", "end")

    def _copy_log(self):
        content = self.log_text.get("1.0", "end-1c")
        self.clipboard_clear()
        self.clipboard_append(content)
        messagebox.showinfo("Thông báo", "Đã sao chép toàn bộ nhật ký vào Clipboard.")


class MinimalistTrimmingApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Optimal Trimming Pipeline - Human mtDNA Sanger Sequencing")
        self.geometry("1060x780")
        self.minsize(960, 680)

        # Set App Icon nếu có
        icon_path = find_app_path(os.path.join("assets", "icon.ico"))
        if os.path.exists(icon_path):
            try:
                self.iconbitmap(icon_path)
            except Exception:
                pass

        self.pipeline_running = False
        self.last_results: Optional[Dict[str, Any]] = None
        self.log_lines_count = 0
        self.all_summary_records: List[Dict[str, Any]] = []

        # Khởi tạo Cửa sổ Console Toplevel (ẩn sẵn)
        self.console_win = ConsoleWindow(self)
        self.console_win.withdraw()

        # Phím tắt toàn cục F12 để mở/đóng Console
        self.bind_all("<F12>", self._toggle_console)

        self._setup_ui()
        self._check_initial_data()

    def _setup_ui(self):
        # 1. TOP HEADER BANNER (Minimalist Slate Surface)
        header_frame = ctk.CTkFrame(self, fg_color=("#ffffff", "#0f172a"), corner_radius=0, height=72)
        header_frame.pack(fill="x", side="top")
        header_frame.pack_propagate(False)

        title_container = ctk.CTkFrame(header_frame, fg_color="transparent")
        title_container.pack(side="left", padx=20, pady=12, fill="y")

        app_title = ctk.CTkLabel(
            title_container,
            text="🧬 OPTIMAL TRIMMING PIPELINE",
            font=ctk.CTkFont(family="Segoe UI", size=18, weight="bold"),
            text_color=("#0284c7", "#38bdf8")
        )
        app_title.pack(anchor="w")

        app_subtitle = ctk.CTkLabel(
            title_container,
            text="Hệ thống dự đoán ranh giới cắt lọc nhiễu tự động & Trực quan hóa điện di đồ đối chiếu rCRS",
            font=ctk.CTkFont(family="Segoe UI", size=11),
            text_color=("#475569", "#94a3b8")
        )
        app_subtitle.pack(anchor="w")

        # Nút đổi Dark/Light mode và Console toggle
        right_header = ctk.CTkFrame(header_frame, fg_color="transparent")
        right_header.pack(side="right", padx=20, pady=12)

        self.btn_header_console = ctk.CTkButton(
            right_header,
            text="💻 Console (F12)",
            font=ctk.CTkFont(size=11, weight="bold"),
            width=110,
            height=28,
            fg_color=("#e2e8f0", "#1e293b"),
            text_color=("#0f172a", "#38bdf8"),
            hover_color=("#cbd5e1", "#334155"),
            command=self._toggle_console
        )
        self.btn_header_console.pack(side="right", padx=(10, 0))

        self.theme_switch = ctk.CTkSwitch(
            right_header,
            text="Dark Mode",
            font=ctk.CTkFont(family="Segoe UI", size=11, weight="bold"),
            text_color=("#0f172a", "#f8fafc"),
            command=self._toggle_theme,
            onvalue="dark",
            offvalue="light"
        )
        self.theme_switch.select()
        self.theme_switch.pack(side="right")

        # 2. MAIN TABVIEW
        self.tabview = ctk.CTkTabview(
            self,
            corner_radius=10,
            fg_color=("#ffffff", "#111827"),
            segmented_button_selected_color=("#0284c7", "#0284c7"),
            segmented_button_selected_hover_color=("#0369a1", "#0369a1"),
            segmented_button_unselected_color=("#e2e8f0", "#1e293b"),
            segmented_button_unselected_hover_color=("#cbd5e1", "#334155"),
            segmented_button_fg_color=("#f1f5f9", "#0f172a")
        )
        self.tabview.pack(fill="both", expand=True, padx=16, pady=(10, 14))

        # CẤU HÌNH ĐỘ TƯƠNG PHẢN CHỮ TABVIEW TRONG CẢ LIGHT & DARK MODE
        self.tabview._segmented_button.configure(
            font=ctk.CTkFont(family="Segoe UI", size=12, weight="bold"),
            text_color=("#0f172a", "#f8fafc"),
            selected_color=("#0284c7", "#0284c7"),
            selected_hover_color=("#0369a1", "#0369a1"),
            unselected_color=("#e2e8f0", "#1e293b"),
            unselected_hover_color=("#cbd5e1", "#334155")
        )

        self.tab_analysis = self.tabview.add("⚡ Phân Tích & Giám Sát")
        self.tab_results = self.tabview.add("📊 Kết Quả QC & Sequencher")
        self.tab_settings = self.tabview.add("⚙️ Cài Đặt Thuật Toán Tracy")

        self._build_analysis_tab()
        self._build_results_tab()
        self._build_settings_tab()

    def _toggle_console(self, event=None):
        """Ẩn/hiện cửa sổ Developer Console (F12)"""
        if self.console_win.winfo_viewable():
            self.console_win.withdraw()
        else:
            self.console_win.deiconify()
            self.console_win.lift()

    # -------------------------------------------------------------
    # TAB 1: PHÂN TÍCH & GIÁM SÁT (MINIMALIST & CLEAN, KHÔNG BỊ TRÀN LOG)
    # -------------------------------------------------------------
    def _build_analysis_tab(self):
        tab = self.tab_analysis

        # Section 1: Card Cấu hình dữ liệu
        io_card = ctk.CTkFrame(tab, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        io_card.pack(fill="x", padx=12, pady=(10, 8), ipady=6)

        ctk.CTkLabel(
            io_card,
            text="📁 Cấu Hình Dữ Liệu Vào / Ra",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold"),
            text_color=("#0f172a", "#f8fafc")
        ).pack(anchor="w", padx=14, pady=(6, 4))

        # Input Row
        in_row = ctk.CTkFrame(io_card, fg_color="transparent")
        in_row.pack(fill="x", padx=14, pady=4)
        ctk.CTkLabel(in_row, text="Thư mục file .ab1:", width=160, anchor="w", font=ctk.CTkFont(size=12, weight="bold")).pack(side="left")
        self.input_var = tk.StringVar()
        self.input_entry = ctk.CTkEntry(in_row, textvariable=self.input_var, placeholder_text="Chọn thư mục chứa file .ab1 cần phân tích...")
        self.input_entry.pack(side="left", fill="x", expand=True, padx=8)
        self.input_entry.bind("<KeyRelease>", lambda e: self._on_input_path_changed())
        ctk.CTkButton(in_row, text="Chọn Thư Mục...", width=110, command=self._browse_input).pack(side="right")

        # Scan status badge
        self.scan_status_lbl = ctk.CTkLabel(
            io_card,
            text="",
            font=ctk.CTkFont(family="Segoe UI", size=11, weight="bold"),
            text_color="#10b981",
            anchor="w"
        )
        self.scan_status_lbl.pack(fill="x", padx=180, pady=(0, 4))

        # rCRS Row
        rcrs_row = ctk.CTkFrame(io_card, fg_color="transparent")
        rcrs_row.pack(fill="x", padx=14, pady=4)
        ctk.CTkLabel(rcrs_row, text="Chuỗi chuẩn rCRS (.fasta):", width=160, anchor="w", font=ctk.CTkFont(size=12, weight="bold")).pack(side="left")
        self.rcrs_var = tk.StringVar()
        self.rcrs_entry = ctk.CTkEntry(rcrs_row, textvariable=self.rcrs_var, placeholder_text="Đường dẫn file rCRS.fasta...")
        self.rcrs_entry.pack(side="left", fill="x", expand=True, padx=8)
        ctk.CTkButton(rcrs_row, text="Chọn File...", width=110, command=self._browse_rcrs).pack(side="right")

        # Output Row
        out_row = ctk.CTkFrame(io_card, fg_color="transparent")
        out_row.pack(fill="x", padx=14, pady=4)
        ctk.CTkLabel(out_row, text="Thư mục lưu kết quả:", width=160, anchor="w", font=ctk.CTkFont(size=12, weight="bold")).pack(side="left")
        self.output_var = tk.StringVar(value=os.path.join(EXE_DIR, "output"))
        self.output_entry = ctk.CTkEntry(out_row, textvariable=self.output_var, placeholder_text="Thư mục xuất báo cáo...")
        self.output_entry.pack(side="left", fill="x", expand=True, padx=8)
        ctk.CTkButton(out_row, text="Chọn Thư Mục...", width=110, command=self._browse_output).pack(side="right")

        # Section 2: Action & Live Process Monitor Card (Thay thế khung đen to đùng bằng giao diện tối giản, hiện đại)
        action_card = ctk.CTkFrame(tab, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        action_card.pack(fill="x", padx=12, pady=8, ipady=6)

        ctk.CTkLabel(
            action_card,
            text="⚡ Tiến Trình & Điều Khiển Phân Tích",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold"),
            text_color=("#0f172a", "#f8fafc")
        ).pack(anchor="w", padx=14, pady=(6, 8))

        # Nút bắt đầu phân tích nổi bật
        self.run_btn = ctk.CTkButton(
            action_card,
            text="▶ BẮT ĐẦU PHÂN TÍCH HÀNG LOẠT (1-CLICK RUN)",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold"),
            fg_color="#0284c7",
            hover_color="#0369a1",
            height=44,
            corner_radius=8,
            command=self._start_pipeline
        )
        self.run_btn.pack(fill="x", padx=14, pady=(2, 8))

        # Progress bar
        self.progress_bar = ctk.CTkProgressBar(action_card, height=12, corner_radius=6, progress_color="#38bdf8")
        self.progress_bar.set(0)
        self.progress_bar.pack(fill="x", padx=14, pady=(6, 4))

        # Live Activity Feed Pill (thẻ hiển thị bước phân tích hiện tại mượt mà)
        ticker_frame = ctk.CTkFrame(action_card, fg_color=("#e2e8f0", "#0f172a"), corner_radius=8, height=36)
        ticker_frame.pack(fill="x", padx=14, pady=(4, 8))
        ticker_frame.pack_propagate(False)

        self.status_dot = ctk.CTkLabel(
            ticker_frame,
            text="●",
            font=ctk.CTkFont(size=14),
            text_color="#10b981",
            width=24
        )
        self.status_dot.pack(side="left", padx=(10, 0))

        self.progress_lbl = ctk.CTkLabel(
            ticker_frame,
            text="Sẵn sàng phân tích. Bấm nút phía trên để bắt đầu.",
            font=ctk.CTkFont(family="Segoe UI", size=11, weight="bold"),
            text_color=("#0f172a", "#f8fafc"),
            anchor="w"
        )
        self.progress_lbl.pack(side="left", fill="x", expand=True, padx=4)

        # Hàng nút thao tác nhanh
        action_btn_row = ctk.CTkFrame(action_card, fg_color="transparent")
        action_btn_row.pack(fill="x", padx=14, pady=(4, 6))

        self.btn_open_console = ctk.CTkButton(
            action_btn_row,
            text="💻 Xem Nhật Ký / Console (F12)",
            font=ctk.CTkFont(size=12, weight="bold"),
            fg_color=("#cbd5e1", "#334155"),
            text_color=("#0f172a", "#f8fafc"),
            hover_color=("#94a3b8", "#475569"),
            height=34,
            command=self._toggle_console
        )
        self.btn_open_console.pack(side="left", padx=(0, 6))

        self.open_res_btn = ctk.CTkButton(
            action_btn_row,
            text="📊 Mở Bảng Khuyến Nghị (CSV)",
            font=ctk.CTkFont(size=12),
            fg_color=("#0f766e", "#0f766e"),
            hover_color=("#115e59", "#115e59"),
            height=34,
            command=self._open_trim_csv
        )
        self.open_res_btn.pack(side="left", padx=6)

        self.open_dir_btn = ctk.CTkButton(
            action_btn_row,
            text="📁 Mở Thư Mục Kết Quả",
            font=ctk.CTkFont(size=12),
            fg_color=("#475569", "#475569"),
            hover_color=("#64748b", "#64748b"),
            height=34,
            command=self._open_output_folder
        )
        self.open_dir_btn.pack(side="right", padx=(6, 0))

        # Section 3: Thẻ Xem Nhanh Kết Quả Ngay Trên Tab 1 (Sau khi chạy xong sẽ tự hiện)
        self.quick_summary_card = ctk.CTkFrame(tab, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        self.quick_summary_card.pack(fill="x", padx=12, pady=6, ipady=8)

        self.quick_summary_head = ctk.CTkLabel(
            self.quick_summary_card,
            text="📊 Tóm Tắt Nhanh Lần Chạy Gần Nhất",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold"),
            text_color=("#0f172a", "#f8fafc")
        )
        self.quick_summary_head.pack(anchor="w", padx=14, pady=(6, 6))

        self.quick_metrics_lbl = ctk.CTkLabel(
            self.quick_summary_card,
            text="Chưa có dữ liệu phân tích nào. Hãy bấm 'Bắt đầu phân tích hàng loạt' để tạo báo cáo.",
            font=ctk.CTkFont(family="Segoe UI", size=12),
            text_color=("#64748b", "#94a3b8")
        )
        self.quick_metrics_lbl.pack(anchor="w", padx=14, pady=(0, 6))

        self.quick_nav_btn = ctk.CTkButton(
            self.quick_summary_card,
            text="👉 Xem Chi Tiết Bảng QC & Khuyến Nghị Sequencher (Tab 2)",
            font=ctk.CTkFont(size=12, weight="bold"),
            fg_color="#0284c7",
            hover_color="#0369a1",
            height=32,
            command=lambda: self.tabview.set("📊 Kết Quả QC & Sequencher")
        )
        self.quick_nav_btn.pack(anchor="w", padx=14, pady=(2, 6))

    # -------------------------------------------------------------
    # TAB 2: KẾT QUẢ QC & SEQUENCHER
    # -------------------------------------------------------------
    def _build_results_tab(self):
        tab = self.tab_results

        # 1. KPI Cards Row
        kpi_frame = ctk.CTkFrame(tab, fg_color="transparent")
        kpi_frame.pack(fill="x", padx=10, pady=(10, 6))

        # 4 cards
        self.card_total = self._create_kpi_card(kpi_frame, "TỔNG SỐ PHÂN TÍCH", "0", "#38bdf8")
        self.card_pass = self._create_kpi_card(kpi_frame, "PASS (ĐẠT CHUẨN)", "0", "#10b981")
        self.card_review = self._create_kpi_card(kpi_frame, "REVIEW (SOI LẠI)", "0", "#f59e0b")
        self.card_deny = self._create_kpi_card(kpi_frame, "DENY (LÀM LẠI MẪU)", "0", "#ef4444")

        # 2. Results Table Container
        table_container = ctk.CTkFrame(tab, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        table_container.pack(fill="both", expand=True, padx=10, pady=6)

        table_header = ctk.CTkFrame(table_container, fg_color="transparent")
        table_header.pack(fill="x", padx=14, pady=(8, 2))
        ctk.CTkLabel(
            table_header,
            text="Bảng Tổng Hợp Kiểm Soát Chất Lượng (QC Matrix) & Tọa Độ Cắt Sequencher",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold"),
            text_color=("#0f172a", "#f8fafc")
        ).pack(side="left")

        ctk.CTkLabel(
            table_header,
            text="(Double-click hoặc chọn dòng rồi bấm 'Soi Điện Di Đồ')",
            font=ctk.CTkFont(size=11),
            text_color=("#64748b", "#94a3b8")
        ).pack(side="right")

        # QC Filter & Quick Action Row
        filter_bar = ctk.CTkFrame(table_container, fg_color="transparent")
        filter_bar.pack(fill="x", padx=14, pady=(2, 6))

        ctk.CTkLabel(
            filter_bar,
            text="Lọc theo QC:",
            font=ctk.CTkFont(family="Segoe UI", size=11, weight="bold"),
            text_color=("#475569", "#94a3b8")
        ).pack(side="left", padx=(0, 6))

        self.filter_seg = ctk.CTkSegmentedButton(
            filter_bar,
            values=["Tất Cả", "PASS ✓", "REVIEW 👁️", "DENY ⚠️"],
            command=self._on_filter_changed,
            font=ctk.CTkFont(family="Segoe UI", size=11, weight="bold"),
            selected_color="#0284c7",
            selected_hover_color="#0369a1",
            unselected_color=("#e2e8f0", "#1e293b"),
            unselected_hover_color=("#cbd5e1", "#334155")
        )
        self.filter_seg.set("Tất Cả")
        self.filter_seg.pack(side="left", padx=(0, 10))

        self.btn_open_selected_chroma = ctk.CTkButton(
            filter_bar,
            text="📈 Soi Điện Di Đồ Mẫu Đang Chọn",
            font=ctk.CTkFont(size=11, weight="bold"),
            fg_color="#0284c7",
            hover_color="#0369a1",
            height=28,
            command=self._open_selected_chromatogram
        )
        self.btn_open_selected_chroma.pack(side="right", padx=(6, 0))

        self.btn_copy_coords = ctk.CTkButton(
            filter_bar,
            text="📋 Sao Chép Tọa Độ Sequencher",
            font=ctk.CTkFont(size=11, weight="bold"),
            fg_color=("#e2e8f0", "#334155"),
            text_color=("#0f172a", "#f8fafc"),
            hover_color=("#cbd5e1", "#475569"),
            height=28,
            command=self._copy_selected_trim
        )
        self.btn_copy_coords.pack(side="right", padx=(6, 0))

        # Treeview bọc trong CTkFrame
        tree_frame = ctk.CTkFrame(table_container, fg_color="transparent")
        tree_frame.pack(fill="both", expand=True, padx=12, pady=(0, 8))

        columns = ("Sample", "Region", "Coverage", "SNPs", "Conflict_+", "Peak_?", "QC_Status", "Action_Lab")
        self.tree = ttk.Treeview(tree_frame, columns=columns, show="headings", height=10)
        self.tree.heading("Sample", text="Mẫu (Sample)")
        self.tree.heading("Region", text="Vùng (Region)")
        self.tree.heading("Coverage", text="Độ Phủ (%)")
        self.tree.heading("SNPs", text="SNPs")
        self.tree.heading("Conflict_+", text="Dị Thể (+)")
        self.tree.heading("Peak_?", text="Chập Peak (?)")
        self.tree.heading("QC_Status", text="Trạng Thái QC")
        self.tree.heading("Action_Lab", text="Khuyến Cáo / Hướng Dẫn Sequencher")

        self.tree.column("Sample", width=90, anchor="center")
        self.tree.column("Region", width=90, anchor="center")
        self.tree.column("Coverage", width=90, anchor="center")
        self.tree.column("SNPs", width=70, anchor="center")
        self.tree.column("Conflict_+", width=80, anchor="center")
        self.tree.column("Peak_?", width=90, anchor="center")
        self.tree.column("QC_Status", width=120, anchor="center")
        self.tree.column("Action_Lab", width=360, anchor="w")

        # Treeview scrollbar
        tree_scroll_y = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=tree_scroll_y.set)
        self.tree.pack(side="left", fill="both", expand=True)
        tree_scroll_y.pack(side="right", fill="y")

        self.tree.bind("<Double-1>", self._on_tree_double_click)

        # Style Treeview matching dark/light mode
        self._style_treeview()

        # 3. Action Buttons Toolbar
        action_bar = ctk.CTkFrame(tab, fg_color="transparent")
        action_bar.pack(fill="x", padx=10, pady=(4, 10))

        ctk.CTkButton(
            action_bar,
            text="📊 Mở Bảng Khuyến Nghị (CSV Form A)",
            font=ctk.CTkFont(size=12, weight="bold"),
            fg_color="#0284c7",
            hover_color="#0369a1",
            height=36,
            command=self._open_trim_csv
        ).pack(side="left", padx=(0, 6))

        ctk.CTkButton(
            action_bar,
            text="📑 Mở Tổng Hợp Kiểu Gen (CSV)",
            font=ctk.CTkFont(size=12),
            fg_color=("#334155", "#334155"),
            text_color="#ffffff",
            hover_color=("#475569", "#475569"),
            height=36,
            command=self._open_summary_csv
        ).pack(side="left", padx=6)

        ctk.CTkButton(
            action_bar,
            text="📈 Mở Thư Mục Biểu Đồ HTML (Form C)",
            font=ctk.CTkFont(size=12),
            fg_color=("#0f766e", "#0f766e"),
            hover_color=("#115e59", "#115e59"),
            height=36,
            command=self._open_html_reports
        ).pack(side="left", padx=6)

        ctk.CTkButton(
            action_bar,
            text="📁 Mở Thư Mục Kết Quả",
            font=ctk.CTkFont(size=12),
            fg_color=("#475569", "#475569"),
            hover_color=("#64748b", "#64748b"),
            height=36,
            command=self._open_output_folder
        ).pack(side="right", padx=(6, 0))

    def _create_kpi_card(self, parent, title: str, init_val: str, accent_color: str):
        card = ctk.CTkFrame(parent, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        card.pack(side="left", fill="both", expand=True, padx=4)

        ctk.CTkLabel(
            card,
            text=title,
            font=ctk.CTkFont(family="Segoe UI", size=10, weight="bold"),
            text_color=("#64748b", "#94a3b8")
        ).pack(anchor="center", pady=(8, 0))

        val_lbl = ctk.CTkLabel(
            card,
            text=init_val,
            font=ctk.CTkFont(family="Segoe UI", size=22, weight="bold"),
            text_color=accent_color
        )
        val_lbl.pack(anchor="center", pady=(0, 8))
        return val_lbl

    def _style_treeview(self):
        mode = ctk.get_appearance_mode()
        style = ttk.Style()
        style.theme_use('clam')

        if mode == "Light":
            bg_col = "#ffffff"
            fg_col = "#0f172a"
            head_bg = "#f1f5f9"
            head_fg = "#0284c7"
            sel_bg = "#bae6fd"
            sel_fg = "#0f172a"
        else:
            bg_col = "#1e293b"
            fg_col = "#f8fafc"
            head_bg = "#0f172a"
            head_fg = "#38bdf8"
            sel_bg = "#0369a1"
            sel_fg = "#ffffff"

        style.configure(
            "Treeview",
            background=bg_col,
            foreground=fg_col,
            fieldbackground=bg_col,
            rowheight=28,
            font=("Segoe UI", 10)
        )
        style.configure(
            "Treeview.Heading",
            background=head_bg,
            foreground=head_fg,
            relief="flat",
            font=("Segoe UI", 10, "bold")
        )
        style.map("Treeview", background=[('selected', sel_bg)], foreground=[('selected', sel_fg)])

        if mode == "Light":
            self.tree.tag_configure("PASS", foreground="#059669")
            self.tree.tag_configure("REVIEW", foreground="#d97706")
            self.tree.tag_configure("DENY", foreground="#dc2626")
        else:
            self.tree.tag_configure("PASS", foreground="#34d399")
            self.tree.tag_configure("REVIEW", foreground="#fbbf24")
            self.tree.tag_configure("DENY", foreground="#f87171")

    # -------------------------------------------------------------
    # TAB 3: CÀI ĐẶT THUẬT TOÁN TRACY
    # -------------------------------------------------------------
    def _build_settings_tab(self):
        tab = self.tab_settings

        # Parameters Card
        param_card = ctk.CTkFrame(tab, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        param_card.pack(fill="x", padx=10, pady=(10, 8), ipady=6)

        ctk.CTkLabel(
            param_card,
            text="Ngưỡng Nhận Diện Dị Thể & Lỗi Sóng Quang Học (Tracy Engine)",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold"),
            text_color=("#0f172a", "#f8fafc")
        ).pack(anchor="w", padx=14, pady=(6, 10))

        # Slider 1: Secondary Ratio
        s1_box = ctk.CTkFrame(param_card, fg_color="transparent")
        s1_box.pack(fill="x", padx=14, pady=4)
        ctk.CTkLabel(s1_box, text="Ngưỡng Đỉnh Phụ Tracy (Secondary Peak Ratio):", width=300, anchor="w", font=ctk.CTkFont(size=12, weight="bold")).pack(side="left")
        self.thresh_var = tk.DoubleVar(value=0.25)
        self.thresh_lbl = ctk.CTkLabel(s1_box, text="0.25", width=50, font=ctk.CTkFont(size=12, weight="bold"), text_color="#38bdf8")
        self.thresh_lbl.pack(side="right")
        self.thresh_slider = ctk.CTkSlider(
            s1_box,
            from_=0.10,
            to=0.50,
            number_of_steps=40,
            variable=self.thresh_var,
            command=lambda v: self.thresh_lbl.configure(text=f"{float(v):.2f}")
        )
        self.thresh_slider.pack(side="right", fill="x", expand=True, padx=10)

        ctk.CTkLabel(
            param_card,
            text="  ↳ Nhận diện dị thể (Heteroplasmy) hoặc đỉnh phụ. Mặc định 0.25 (khuyến nghị cho SeqStudio Flex 24).",
            font=ctk.CTkFont(size=11),
            text_color=("#64748b", "#94a3b8")
        ).pack(anchor="w", padx=20, pady=(0, 6))

        # Slider 2: Conflict Threshold
        s2_box = ctk.CTkFrame(param_card, fg_color="transparent")
        s2_box.pack(fill="x", padx=14, pady=4)
        ctk.CTkLabel(s2_box, text="Ngưỡng Chập Peak / Lỗi Cảm Biến (?):", width=300, anchor="w", font=ctk.CTkFont(size=12, weight="bold")).pack(side="left")
        self.conflict_var = tk.DoubleVar(value=0.70)
        self.conflict_lbl = ctk.CTkLabel(s2_box, text="0.70", width=50, font=ctk.CTkFont(size=12, weight="bold"), text_color="#38bdf8")
        self.conflict_lbl.pack(side="right")
        self.conflict_slider = ctk.CTkSlider(
            s2_box,
            from_=0.50,
            to=0.90,
            number_of_steps=40,
            variable=self.conflict_var,
            command=lambda v: self.conflict_lbl.configure(text=f"{float(v):.2f}")
        )
        self.conflict_slider.pack(side="right", fill="x", expand=True, padx=10)

        ctk.CTkLabel(
            param_card,
            text="  ↳ Tỷ lệ đỉnh phụ/đỉnh chính tối thiểu để xác định chập peak quang học ('?'). Mặc định 0.70.",
            font=ctk.CTkFont(size=11),
            text_color=("#64748b", "#94a3b8")
        ).pack(anchor="w", padx=20, pady=(0, 8))

        # Bio Info Card
        bio_card = ctk.CTkFrame(tab, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        bio_card.pack(fill="both", expand=True, padx=10, pady=(4, 10), ipady=6)

        ctk.CTkLabel(
            bio_card,
            text="Quy Chuẩn Phân Tích D-loop Hệ Gen Ty Thể Người (mtDNA Standards)",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold"),
            text_color=("#0f172a", "#f8fafc")
        ).pack(anchor="w", padx=14, pady=(8, 6))

        info_text = (
            "• VÙNG HV1: Tọa độ rCRS 16024 - 16365 (Tổng 342 bp).\n"
            "• VÙNG HV2/HV3: Tọa độ rCRS 73 - 576 (Tổng 504 bp).\n"
            "• BẢO TOÀN DẢI SAU TRƯỢT POLY-C (C-Tract Stutter): Tại HV1 (16184-16193) và HV2 (303-315),\n"
            "  hiện tượng trượt pha do polymerase được bảo toàn để đối sánh chéo giữa mồi xuôi và mồi ngược,\n"
            "  ưu tiên mồi đối diện sạch tín hiệu để xác lập Consensus.\n"
            "• NGUYÊN TẮC NON-DESTRUCTIVE: File .ab1 gốc giữ nguyên 100%. Phần mềm tính toán tọa độ cắt\n"
            "  và xuất bảng Form A (CSV) để kỹ thuật viên nhập trực tiếp vào Sequencher v5.4.6.\n"
            "• PHÂN LOẠI QC & CẢNH BÁO DENY:\n"
            "  - PASS: Độ phủ >= 98%, 0 điểm mâu thuẫn (+), 0 điểm chập peak (?).\n"
            "  - REVIEW: Độ phủ 88% - 98%, hoặc có 1-2 điểm dị thể/nhiễu nhẹ.\n"
            "  - DENY (Cảnh báo đỏ): Độ phủ < 85%, hoặc đứt đoạn >= 20 bp, chập peak >= 5 điểm,\n"
            "    hoặc >= 20 SNPs (nghi ngờ tạp nhiễm). Khuyến cáo kỹ thuật viên KHÔNG import Sequencher\n"
            "    mà làm lại mẫu trong phòng thí nghiệm!"
        )

        ctk.CTkLabel(
            bio_card,
            text=info_text,
            font=ctk.CTkFont(family="Segoe UI", size=11),
            justify="left",
            text_color=("#334155", "#cbd5e1")
        ).pack(anchor="w", padx=14, pady=(0, 8))

    # -------------------------------------------------------------
    # LOGIC & SỰ KIỆN TƯƠNG TÁC
    # -------------------------------------------------------------
    def _check_initial_data(self):
        """Tự động kiểm tra và điền đường dẫn mặc định khi khởi chạy"""
        # 1. Tìm rCRS.fasta
        default_rcrs = find_app_path(os.path.join("data", "rCRS.fasta"))
        if os.path.exists(default_rcrs):
            self.rcrs_var.set(default_rcrs)

        # 2. Tìm thư mục dữ liệu mẫu
        sample_candidates = [
            find_app_path(os.path.join("data", "sample_ab1", "real data")),
            find_app_path(os.path.join("data", "sample_ab1")),
        ]
        for sc in sample_candidates:
            if os.path.exists(sc) and glob.glob(os.path.join(sc, "**", "*.ab1"), recursive=True):
                self.input_var.set(sc)
                break

        self._on_input_path_changed()

        # 3. Tự động nạp dữ liệu kết quả phân tích gần nhất nếu có sẵn
        out_cand = self.output_var.get().strip()
        summary_f = os.path.join(out_cand, "mtDNA_batch_summary.csv")
        if not os.path.exists(summary_f):
            demo_cand = find_app_path(os.path.join("output_demo", "mtDNA_batch_summary.csv"))
            if os.path.exists(demo_cand):
                summary_f = demo_cand
                self.output_var.set(find_app_path("output_demo"))

        if os.path.exists(summary_f):
            try:
                import csv
                loaded_recs = []
                with open(summary_f, "r", encoding="utf-8-sig") as f:
                    reader = csv.DictReader(f)
                    for r in reader:
                        qc = r.get("QC_Status", "REVIEW")
                        if "REVIEW" in qc:
                            qc = "REVIEW"
                        r["QC_Status"] = qc
                        if "Coverage_%" in r:
                            try:
                                r["Coverage_%"] = float(r["Coverage_%"])
                            except Exception:
                                pass
                        loaded_recs.append(r)
                if loaded_recs:
                    self.all_summary_records = loaded_recs
                    total = len(loaded_recs)
                    pass_cnt = sum(1 for r in loaded_recs if r.get("QC_Status") == "PASS")
                    rev_cnt = sum(1 for r in loaded_recs if r.get("QC_Status") == "REVIEW")
                    deny_cnt = sum(1 for r in loaded_recs if r.get("QC_Status") == "DENY")
                    self.card_total.configure(text=str(total))
                    self.card_pass.configure(text=str(pass_cnt))
                    self.card_review.configure(text=str(rev_cnt))
                    self.card_deny.configure(text=str(deny_cnt))
                    self._populate_treeview(loaded_recs)
                    self.quick_metrics_lbl.configure(
                        text=f"✓ Đã nạp kết quả gần nhất: {total} mẫu | PASS: {pass_cnt} | REVIEW: {rev_cnt} | DENY: {deny_cnt}",
                        text_color=("#0f172a", "#38bdf8")
                    )
            except Exception:
                pass

    def _on_input_path_changed(self):
        """Quét và đếm nhanh số lượng file .ab1 và mẫu trong thư mục đã chọn"""
        in_path = self.input_var.get().strip()
        if not in_path or not os.path.isdir(in_path):
            self.scan_status_lbl.configure(text="")
            return

        ab1s = glob.glob(os.path.join(in_path, "**", "*.ab1"), recursive=True)
        if not ab1s:
            self.scan_status_lbl.configure(
                text="⚠️ Chưa tìm thấy file .ab1 nào trong thư mục này.",
                text_color="#f59e0b"
            )
            return

        # Thống kê mẫu
        samples = set()
        for f in ab1s:
            sid, _, _, _ = detect_sample_info(f)
            if sid:
                samples.add(sid)

        sample_preview = ", ".join(list(samples)[:6])
        if len(samples) > 6:
            sample_preview += f", ... (+{len(samples)-6} mẫu)"

        self.scan_status_lbl.configure(
            text=f"✓ Đã phát hiện {len(ab1s)} file .ab1 ({len(samples)} mẫu: {sample_preview})",
            text_color="#10b981"
        )

    def _browse_input(self):
        d = filedialog.askdirectory(title="Chọn thư mục chứa các file .ab1")
        if d:
            self.input_var.set(d)
            self._on_input_path_changed()

    def _browse_rcrs(self):
        f = filedialog.askopenfilename(
            title="Chọn file chuẩn rCRS FASTA",
            filetypes=[("FASTA files", "*.fasta;*.fa;*.txt"), ("All files", "*.*")]
        )
        if f:
            self.rcrs_var.set(f)

    def _browse_output(self):
        d = filedialog.askdirectory(title="Chọn thư mục lưu kết quả phân tích")
        if d:
            self.output_var.set(d)

    def _toggle_theme(self):
        mode = self.theme_switch.get()
        ctk.set_appearance_mode(mode)
        self.theme_switch.configure(text="Dark Mode" if mode == "dark" else "Light Mode")
        # Đồng bộ màu cho segmented buttons
        self.tabview._segmented_button.configure(
            text_color=("#0f172a", "#f8fafc")
        )
        self._style_treeview()

    def _log(self, text: str):
        self.log_lines_count += 1
        self.console_win.append_log(text)
        # Cập nhật nhãn nút console hiển thị số dòng
        self.btn_header_console.configure(text=f"💻 Console ({self.log_lines_count})")
        self.btn_open_console.configure(text=f"💻 Xem Nhật Ký / Console ({self.log_lines_count} dòng - F12)")

    def _start_pipeline(self):
        if self.pipeline_running:
            return

        in_dir = self.input_var.get().strip()
        rcrs_f = self.rcrs_var.get().strip()
        out_dir = self.output_var.get().strip()

        if not os.path.exists(in_dir):
            messagebox.showerror("Lỗi", f"Thư mục chứa file .ab1 không tồn tại:\n{in_dir}")
            return
        if not os.path.exists(rcrs_f):
            messagebox.showerror("Lỗi", f"File chuẩn rCRS.fasta không tồn tại:\n{rcrs_f}")
            return

        self.pipeline_running = True
        self.run_btn.configure(state="disabled", text="⏳ ĐANG PHÂN TÍCH DỮ LIỆU...")
        self.progress_bar.set(0)
        self.status_dot.configure(text_color="#f59e0b")
        self.progress_lbl.configure(text="Đang chuẩn bị nạp dữ liệu và kiểm tra...")
        self.log_lines_count = 0
        self.console_win._clear_log()

        def worker():
            class GuiStdoutRedirector:
                def __init__(self, app_inst):
                    self.app = app_inst
                def write(self, s):
                    v = s.strip()
                    if v:
                        self.app.after(0, self.app._log, v)
                def flush(self):
                    pass

            old_stdout = sys.stdout
            sys.stdout = GuiStdoutRedirector(self)

            def progress_cb(current, total, msg):
                frac = current / total if total > 0 else 0
                pct = int(frac * 100)
                self.after(0, lambda: self.progress_bar.set(frac))
                self.after(0, lambda: self.progress_lbl.configure(text=f"[{pct}%] {msg} ({current}/{total})"))

            try:
                res = run_pipeline(
                    input_dir=in_dir,
                    rcrs_path=rcrs_f,
                    output_dir=out_dir,
                    peak_thresh=self.thresh_var.get(),
                    conflict_thresh=self.conflict_var.get(),
                    progress_callback=progress_cb
                )
                self.after(0, lambda: self._on_pipeline_finished(res))
            except Exception as e:
                import traceback
                err_str = traceback.format_exc()
                self.after(0, lambda: self._on_pipeline_error(str(e), err_str))
            finally:
                sys.stdout = old_stdout
                self.after(0, self._reset_run_btn)

        threading.Thread(target=worker, daemon=True).start()

    def _reset_run_btn(self):
        self.pipeline_running = False
        self.run_btn.configure(state="normal", text="▶ BẮT ĐẦU PHÂN TÍCH HÀNG LOẠT (1-CLICK RUN)")

    def _on_pipeline_error(self, err_msg: str, full_trace: str):
        self.status_dot.configure(text_color="#ef4444")
        self.progress_lbl.configure(text=f"❌ Có lỗi xảy ra: {err_msg}")
        self._log(f"\n[LỖI NGHIÊM TRỌNG]: {err_msg}\n{full_trace}")
        messagebox.showerror("Lỗi Phân Tích", f"Quá trình phân tích gặp lỗi:\n{err_msg}\n\n(Bấm F12 để xem chi tiết log lỗi)")

    def _on_pipeline_finished(self, res: Optional[Dict[str, Any]]):
        self.progress_bar.set(1.0)
        self.status_dot.configure(text_color="#10b981")
        self.progress_lbl.configure(text="✓ Phân tích hoàn tất thành công 100%!")
        self.last_results = res

        if not res:
            return

        summary_records = res.get("summary_records", [])
        self.all_summary_records = summary_records

        # Cập nhật KPI
        total = len(summary_records)
        pass_cnt = sum(1 for r in summary_records if r.get("QC_Status") == "PASS")
        rev_cnt = sum(1 for r in summary_records if r.get("QC_Status") == "REVIEW")
        deny_cnt = sum(1 for r in summary_records if r.get("QC_Status") == "DENY")

        self.card_total.configure(text=str(total))
        self.card_pass.configure(text=str(pass_cnt))
        self.card_review.configure(text=str(rev_cnt))
        self.card_deny.configure(text=str(deny_cnt))

        # Cập nhật Thẻ Tóm Tắt Nhanh trên Tab 1
        summary_text = (
            f"✓ Đã hoàn tất: {total} mẫu phân tích  |  "
            f"PASS (Đạt chuẩn): {pass_cnt}  |  "
            f"REVIEW (Soi lại): {rev_cnt}  |  "
            f"DENY (Cảnh báo làm lại): {deny_cnt}"
        )
        self.quick_metrics_lbl.configure(text=summary_text, text_color=("#0f172a", "#38bdf8"))

        # Cập nhật Treeview theo bộ lọc hiện tại
        self.filter_seg.set("Tất Cả")
        self._populate_treeview(summary_records)

        # Chuyển sang Tab Kết Quả để kỹ thuật viên xem ngay
        self.tabview.set("📊 Kết Quả QC & Sequencher")
        messagebox.showinfo(
            "Phân Tích Thành Công",
            f"Đã xử lý xong toàn bộ {total} mẫu phân tích!\n\n"
            f"• PASS (Đạt chuẩn): {pass_cnt}\n"
            f"• REVIEW (Cần soi lại): {rev_cnt}\n"
            f"• DENY (Cảnh báo làm lại): {deny_cnt}\n\n"
            "Bạn có thể mở ngay Bảng Khuyến Nghị Sequencher (CSV) hoặc Điện di đồ tương tác (HTML)."
        )

    def _on_filter_changed(self, value=None):
        """Lọc danh sách mẫu trong bảng theo trạng thái QC"""
        if not self.all_summary_records:
            return
        val = value or self.filter_seg.get()
        if "PASS" in val:
            filtered = [r for r in self.all_summary_records if r.get("QC_Status") == "PASS"]
        elif "REVIEW" in val:
            filtered = [r for r in self.all_summary_records if r.get("QC_Status") == "REVIEW"]
        elif "DENY" in val:
            filtered = [r for r in self.all_summary_records if r.get("QC_Status") == "DENY"]
        else:
            filtered = self.all_summary_records
        self._populate_treeview(filtered)

    def _populate_treeview(self, records: List[Dict[str, Any]]):
        """Hiển thị danh sách mẫu lên bảng Treeview"""
        for item in self.tree.get_children():
            self.tree.delete(item)

        for rec in records:
            qc = rec.get("QC_Status", "REVIEW")
            tag = qc if qc in ["PASS", "REVIEW", "DENY"] else "REVIEW"
            rec_text = rec.get("Lab_Recommendation", "")
            if not rec_text and qc == "DENY":
                rec_text = "DENY: Khuyến cáo phòng lab làm lại mẫu!"
            elif not rec_text and qc == "PASS":
                rec_text = "Đạt chuẩn chất lượng cao, xuất kết quả."

            self.tree.insert(
                "",
                "end",
                values=(
                    rec.get("Sample_ID", ""),
                    rec.get("Region", ""),
                    f"{rec.get('Coverage_%', 0):.1f}%",
                    rec.get("Variants_Count", 0),
                    rec.get("Ambiguity_Count_+", 0),
                    rec.get("Conflict_Count_?", 0),
                    qc,
                    rec_text
                ),
                tags=(tag,)
            )

    def _open_selected_chromatogram(self):
        """Mở biểu đồ điện di đồ HTML của mẫu đang được chọn trong bảng"""
        sel = self.tree.selection()
        if not sel:
            children = self.tree.get_children()
            if children:
                self.tree.selection_set(children[0])
                sel = (children[0],)
            else:
                messagebox.showinfo("Thông báo", "Chưa có mẫu nào trong bảng kết quả để soi.")
                return

        item = self.tree.item(sel[0])
        vals = item.get("values", [])
        if len(vals) < 2:
            return

        sid = vals[0]
        reg = vals[1]
        out_dir = self.output_var.get().strip()
        html_dir = os.path.join(out_dir, "html_reports")

        if os.path.exists(html_dir):
            pat = os.path.join(html_dir, f"{sid}_{reg}_*.html")
            matches = glob.glob(pat)
            if matches:
                os.startfile(matches[0])
                return
            pat_sample = os.path.join(html_dir, f"{sid}_*.html")
            matches_s = glob.glob(pat_sample)
            if matches_s:
                os.startfile(matches_s[0])
                return
            os.startfile(html_dir)
        else:
            messagebox.showinfo("Thông báo", f"Thư mục html_reports chưa tồn tại trong:\n{out_dir}")

    def _copy_selected_trim(self):
        """Sao chép thông số cắt Sequencher của mẫu đang chọn vào clipboard"""
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("Thông báo", "Vui lòng click chọn một mẫu trong bảng trước khi sao chép.")
            return

        item = self.tree.item(sel[0])
        vals = item.get("values", [])
        if len(vals) < 2:
            return

        sid = vals[0]
        reg = vals[1]
        qc = vals[6] if len(vals) > 6 else ""

        out_dir = self.output_var.get().strip()
        trim_csv = os.path.join(out_dir, "trimming_recommendations.csv")
        found_rows = []
        if os.path.exists(trim_csv):
            try:
                import csv
                with open(trim_csv, "r", encoding="utf-8-sig") as f:
                    reader = csv.DictReader(f)
                    for r in reader:
                        if r.get("Sample_ID") == sid and r.get("Region") == reg:
                            found_rows.append(r)
            except Exception:
                pass

        if found_rows:
            lines = [f"Mẫu: {sid} | Vùng: {reg} | Trạng thái QC: {qc}"]
            for r in found_rows:
                lines.append(f"• File {r.get('Direction', '')}: 5' Cut = {r.get('Suggested_5p_Cut_1based', '')}, 3' Cut = {r.get('Suggested_3p_Cut_1based', '')} (Giữ lại {r.get('Retained_Length_bp', '')} bp)")
            text_to_copy = "\n".join(lines)
        else:
            text_to_copy = f"Mẫu: {sid} | Vùng: {reg} | QC: {qc}"

        self.clipboard_clear()
        self.clipboard_append(text_to_copy)
        messagebox.showinfo("Đã Sao Chép Tọa Độ", f"Đã lưu tọa độ cắt Sequencher của {sid} ({reg}) vào bộ nhớ tạm!\n\n{text_to_copy}")

    def _on_tree_double_click(self, event):
        """Double click vào một dòng trong Treeview để mở Điện Di Đồ HTML tương ứng"""
        self._open_selected_chromatogram()

    def _open_trim_csv(self):
        out_dir = self.output_var.get()
        p = os.path.join(out_dir, "trimming_recommendations.csv")
        if os.path.exists(p):
            os.startfile(p)
        else:
            messagebox.showinfo("Thông báo", f"Chưa tìm thấy file trimming_recommendations.csv trong {out_dir}.\nVui lòng bấm bắt đầu phân tích trước.")

    def _open_summary_csv(self):
        out_dir = self.output_var.get()
        p = os.path.join(out_dir, "mtDNA_batch_summary.csv")
        if os.path.exists(p):
            os.startfile(p)
        else:
            messagebox.showinfo("Thông báo", "Chưa tìm thấy file mtDNA_batch_summary.csv.")

    def _open_html_reports(self):
        out_dir = self.output_var.get()
        p = os.path.join(out_dir, "html_reports")
        if os.path.exists(p):
            os.startfile(p)
        else:
            messagebox.showinfo("Thông báo", "Thư mục html_reports chưa được tạo.")

    def _open_output_folder(self):
        out_dir = self.output_var.get()
        if os.path.exists(out_dir):
            os.startfile(out_dir)
        else:
            messagebox.showinfo("Thông báo", "Thư mục kết quả chưa tồn tại.")


if __name__ == "__main__":
    app = MinimalistTrimmingApp()
    app.mainloop()
