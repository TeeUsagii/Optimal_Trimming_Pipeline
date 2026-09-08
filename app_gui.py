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


class MinimalistTrimmingApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Optimal Trimming Pipeline - Human mtDNA Sanger Sequencing")
        self.geometry("1060x780")
        self.minsize(940, 680)

        # Set App Icon nếu có
        icon_path = find_app_path(os.path.join("assets", "icon.ico"))
        if os.path.exists(icon_path):
            try:
                self.iconbitmap(icon_path)
            except Exception:
                pass

        self.pipeline_running = False
        self.last_results: Optional[Dict[str, Any]] = None

        self._setup_ui()
        self._check_initial_data()

    def _setup_ui(self):
        # 1. TOP HEADER BANNER (Minimalist Slate Surface)
        header_frame = ctk.CTkFrame(self, fg_color=("#f1f5f9", "#0f172a"), corner_radius=0, height=72)
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
            text="Dự đoán ranh giới cắt lọc nhiễu tự động & Trực quan hóa điện di đồ đối chiếu rCRS (Non-destructive)",
            font=ctk.CTkFont(family="Segoe UI", size=11),
            text_color=("#64748b", "#94a3b8")
        )
        app_subtitle.pack(anchor="w")

        # Nút đổi Dark/Light mode và badge phiên bản
        right_header = ctk.CTkFrame(header_frame, fg_color="transparent")
        right_header.pack(side="right", padx=20, pady=12)

        ver_badge = ctk.CTkLabel(
            right_header,
            text="v2.5 Standalone",
            font=ctk.CTkFont(family="Segoe UI", size=11, weight="bold"),
            fg_color=("#e2e8f0", "#1e293b"),
            corner_radius=6,
            padx=8,
            pady=2,
            text_color=("#0369a1", "#7dd3fc")
        )
        ver_badge.pack(side="right", padx=(10, 0))

        self.theme_switch = ctk.CTkSwitch(
            right_header,
            text="Dark Mode",
            font=ctk.CTkFont(family="Segoe UI", size=11),
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
            segmented_button_fg_color=("#cbd5e1", "#0f172a")
        )
        self.tabview.pack(fill="both", expand=True, padx=16, pady=(10, 14))

        self.tab_analysis = self.tabview.add("⚡ Phân Tích & Giám Sát")
        self.tab_results = self.tabview.add("📊 Kết Quả QC & Sequencher")
        self.tab_settings = self.tabview.add("⚙️ Cài Đặt Thuật Toán Tracy")

        self._build_analysis_tab()
        self._build_results_tab()
        self._build_settings_tab()

    # -------------------------------------------------------------
    # TAB 1: PHÂN TÍCH & GIÁM SÁT
    # -------------------------------------------------------------
    def _build_analysis_tab(self):
        tab = self.tab_analysis

        # Section 1: Card Cấu hình dữ liệu
        io_card = ctk.CTkFrame(tab, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        io_card.pack(fill="x", padx=10, pady=(10, 8), ipady=6)

        ctk.CTkLabel(
            io_card,
            text="Cấu Hình Dữ Liệu Vào / Ra",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold"),
            text_color=("#0f172a", "#f8fafc")
        ).pack(anchor="w", padx=14, pady=(6, 4))

        # Input Row
        in_row = ctk.CTkFrame(io_card, fg_color="transparent")
        in_row.pack(fill="x", padx=14, pady=3)
        ctk.CTkLabel(in_row, text="Thư mục file .ab1:", width=160, anchor="w", font=ctk.CTkFont(size=12)).pack(side="left")
        self.input_var = tk.StringVar()
        self.input_entry = ctk.CTkEntry(in_row, textvariable=self.input_var, placeholder_text="Chọn thư mục chứa file .ab1 cần phân tích...")
        self.input_entry.pack(side="left", fill="x", expand=True, padx=8)
        self.input_entry.bind("<KeyRelease>", lambda e: self._on_input_path_changed())
        ctk.CTkButton(in_row, text="Chọn Thư Mục...", width=110, command=self._browse_input).pack(side="right")

        # Scan status badge
        self.scan_status_lbl = ctk.CTkLabel(
            io_card,
            text="",
            font=ctk.CTkFont(family="Segoe UI", size=11),
            text_color="#10b981",
            anchor="w"
        )
        self.scan_status_lbl.pack(fill="x", padx=180, pady=(0, 4))

        # rCRS Row
        rcrs_row = ctk.CTkFrame(io_card, fg_color="transparent")
        rcrs_row.pack(fill="x", padx=14, pady=3)
        ctk.CTkLabel(rcrs_row, text="Chuỗi chuẩn rCRS (.fasta):", width=160, anchor="w", font=ctk.CTkFont(size=12)).pack(side="left")
        self.rcrs_var = tk.StringVar()
        self.rcrs_entry = ctk.CTkEntry(rcrs_row, textvariable=self.rcrs_var, placeholder_text="Đường dẫn file rCRS.fasta...")
        self.rcrs_entry.pack(side="left", fill="x", expand=True, padx=8)
        ctk.CTkButton(rcrs_row, text="Chọn File...", width=110, command=self._browse_rcrs).pack(side="right")

        # Output Row
        out_row = ctk.CTkFrame(io_card, fg_color="transparent")
        out_row.pack(fill="x", padx=14, pady=3)
        ctk.CTkLabel(out_row, text="Thư mục lưu kết quả:", width=160, anchor="w", font=ctk.CTkFont(size=12)).pack(side="left")
        self.output_var = tk.StringVar(value=os.path.join(EXE_DIR, "output"))
        self.output_entry = ctk.CTkEntry(out_row, textvariable=self.output_var, placeholder_text="Thư mục xuất báo cáo...")
        self.output_entry.pack(side="left", fill="x", expand=True, padx=8)
        ctk.CTkButton(out_row, text="Chọn Thư Mục...", width=110, command=self._browse_output).pack(side="right")

        # Section 2: Action & Progress Card
        action_card = ctk.CTkFrame(tab, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        action_card.pack(fill="x", padx=10, pady=6, ipady=4)

        btn_row = ctk.CTkFrame(action_card, fg_color="transparent")
        btn_row.pack(fill="x", padx=14, pady=(8, 4))

        self.run_btn = ctk.CTkButton(
            btn_row,
            text="▶ BẮT ĐẦU PHÂN TÍCH HÀNG LOẠT (1-CLICK RUN)",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold"),
            fg_color="#0284c7",
            hover_color="#0369a1",
            height=40,
            command=self._start_pipeline
        )
        self.run_btn.pack(side="left", fill="x", expand=True, padx=(0, 8))

        self.open_res_btn = ctk.CTkButton(
            btn_row,
            text="📊 Mở Bảng Khuyến Nghị (CSV)",
            font=ctk.CTkFont(size=12),
            fg_color=("#0f766e", "#0f766e"),
            hover_color=("#115e59", "#115e59"),
            height=40,
            command=self._open_trim_csv
        )
        self.open_res_btn.pack(side="right", padx=(0, 0))

        # Progress bar
        self.progress_bar = ctk.CTkProgressBar(action_card, height=10, corner_radius=5, progress_color="#38bdf8")
        self.progress_bar.set(0)
        self.progress_bar.pack(fill="x", padx=14, pady=(6, 2))

        self.progress_lbl = ctk.CTkLabel(
            action_card,
            text="Sẵn sàng thực hiện phân tích",
            font=ctk.CTkFont(size=11),
            text_color=("#64748b", "#94a3b8")
        )
        self.progress_lbl.pack(anchor="w", padx=14, pady=(0, 6))

        # Section 3: Log Console Card
        log_card = ctk.CTkFrame(tab, corner_radius=10, fg_color=("#f8fafc", "#1e293b"), border_width=1, border_color=("#e2e8f0", "#334155"))
        log_card.pack(fill="both", expand=True, padx=10, pady=(6, 8))

        log_head = ctk.CTkFrame(log_card, fg_color="transparent")
        log_head.pack(fill="x", padx=14, pady=(8, 4))

        ctk.CTkLabel(
            log_head,
            text="Nhật Ký Xử Lý Trực Tiếp (Real-time Console)",
            font=ctk.CTkFont(family="Segoe UI", size=12, weight="bold")
        ).pack(side="left")

        ctk.CTkButton(
            log_head,
            text="Xóa Log",
            width=70,
            height=26,
            font=ctk.CTkFont(size=11),
            fg_color=("#cbd5e1", "#334155"),
            text_color=("#0f172a", "#f8fafc"),
            hover_color=("#94a3b8", "#475569"),
            command=self._clear_log
        ).pack(side="right")

        self.log_box = ctk.CTkTextbox(
            log_card,
            corner_radius=8,
            font=ctk.CTkFont(family="Consolas", size=11),
            fg_color=("#0f172a", "#0b0f19"),
            text_color="#e2e8f0",
            wrap="char"
        )
        self.log_box.pack(fill="both", expand=True, padx=12, pady=(0, 10))

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
        table_header.pack(fill="x", padx=14, pady=(8, 4))
        ctk.CTkLabel(
            table_header,
            text="Bảng Tổng Hợp Kiểm Soát Chất Lượng (QC Matrix) & Tọa Độ Cắt Sequencher",
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold")
        ).pack(side="left")

        ctk.CTkLabel(
            table_header,
            text="(Double click vào dòng bất kỳ để mở Biểu Đồ Điện Di Đồ HTML)",
            font=ctk.CTkFont(size=11),
            text_color=("#64748b", "#94a3b8")
        ).pack(side="right")

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
        style = ttk.Style()
        style.theme_use('clam')
        style.configure(
            "Treeview",
            background="#1e293b",
            foreground="#f8fafc",
            fieldbackground="#1e293b",
            rowheight=28,
            font=("Segoe UI", 10)
        )
        style.configure(
            "Treeview.Heading",
            background="#0f172a",
            foreground="#38bdf8",
            relief="flat",
            font=("Segoe UI", 10, "bold")
        )
        style.map("Treeview", background=[('selected', '#0369a1')])
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
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold")
        ).pack(anchor="w", padx=14, pady=(6, 10))

        # Slider 1: Secondary Ratio
        s1_box = ctk.CTkFrame(param_card, fg_color="transparent")
        s1_box.pack(fill="x", padx=14, pady=4)
        ctk.CTkLabel(s1_box, text="Ngưỡng Đỉnh Phụ Tracy (Secondary Peak Ratio):", width=300, anchor="w", font=ctk.CTkFont(size=12)).pack(side="left")
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
        ctk.CTkLabel(s2_box, text="Ngưỡng Chập Peak / Lỗi Cảm Biến (?):", width=300, anchor="w", font=ctk.CTkFont(size=12)).pack(side="left")
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
            font=ctk.CTkFont(family="Segoe UI", size=13, weight="bold")
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

    def _log(self, text: str):
        self.log_box.insert("end", text + "\n")
        self.log_box.see("end")

    def _clear_log(self):
        self.log_box.delete("1.0", "end")

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
        self.progress_lbl.configure(text="Đang chuẩn bị nạp dữ liệu và kiểm tra...")
        self._clear_log()

        # Switch to analysis tab
        self.tabview.set("⚡ Phân Tích & Giám Sát")

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
        self.progress_lbl.configure(text=f"❌ Có lỗi xảy ra: {err_msg}", text_color="#ef4444")
        self._log(f"\n[LỖI NGHIÊM TRỌNG]: {err_msg}\n{full_trace}")
        messagebox.showerror("Lỗi Phân Tích", f"Quá trình phân tích gặp lỗi:\n{err_msg}")

    def _on_pipeline_finished(self, res: Optional[Dict[str, Any]]):
        self.progress_bar.set(1.0)
        self.progress_lbl.configure(text="✓ Phân tích hoàn tất thành công 100%!", text_color="#10b981")
        self.last_results = res

        if not res:
            return

        summary_records = res.get("summary_records", [])

        # Cập nhật KPI
        total = len(summary_records)
        pass_cnt = sum(1 for r in summary_records if r.get("QC_Status") == "PASS")
        rev_cnt = sum(1 for r in summary_records if r.get("QC_Status") == "REVIEW")
        deny_cnt = sum(1 for r in summary_records if r.get("QC_Status") == "DENY")

        self.card_total.configure(text=str(total))
        self.card_pass.configure(text=str(pass_cnt))
        self.card_review.configure(text=str(rev_cnt))
        self.card_deny.configure(text=str(deny_cnt))

        # Cập nhật Treeview
        for item in self.tree.get_children():
            self.tree.delete(item)

        for rec in summary_records:
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

    def _on_tree_double_click(self, event):
        """Double click vào một dòng trong Treeview để mở Điện Di Đồ HTML tương ứng"""
        sel = self.tree.selection()
        if not sel or not self.last_results:
            return
        item = self.tree.item(sel[0])
        vals = item.get("values", [])
        if len(vals) < 2:
            return

        sid = vals[0]
        reg = vals[1]
        html_dir = self.last_results.get("html_dir")
        if html_dir and os.path.exists(html_dir):
            # Tìm file chromatogram của mẫu này
            pat = os.path.join(html_dir, f"{sid}_{reg}_*.html")
            matches = glob.glob(pat)
            if matches:
                os.startfile(matches[0])
                return

            # Thử mở thư mục html
            os.startfile(html_dir)

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
