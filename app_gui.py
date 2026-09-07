import os
import sys
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

if sys.stdout and hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
if sys.stderr and hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8')

# Xác định thư mục ứng dụng (hỗ trợ cả chạy script trực tiếp và chạy file .exe đóng gói)
if getattr(sys, 'frozen', False):
    EXE_DIR = os.path.dirname(sys.executable)
    INTERNAL_DIR = getattr(sys, '_MEIPASS', EXE_DIR)
else:
    EXE_DIR = os.path.dirname(os.path.abspath(__file__))
    INTERNAL_DIR = EXE_DIR

def find_app_path(rel_path: str) -> str:
    p1 = os.path.join(EXE_DIR, rel_path)
    if os.path.exists(p1):
        return p1
    p2 = os.path.join(INTERNAL_DIR, rel_path)
    if os.path.exists(p2):
        return p2
    return p1

CURRENT_DIR = EXE_DIR
if INTERNAL_DIR not in sys.path:
    sys.path.insert(0, INTERNAL_DIR)
if EXE_DIR not in sys.path:
    sys.path.insert(0, EXE_DIR)

from run_cli import run_pipeline

class OptimalTrimmingGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Optimal Trimming Pipeline - Human mtDNA Sanger Sequencing")
        self.root.geometry("820x680")
        self.root.minsize(760, 600)

        # Style cấu hình
        self.style = ttk.Style()
        self.style.theme_use('clam')

        self.setup_ui()

    def setup_ui(self):
        # 1. Header Banner
        header_frame = tk.Frame(self.root, bg="#1e293b", padx=16, pady=12)
        header_frame.pack(fill=tk.X)

        title_label = tk.Label(header_frame, text="OPTIMAL TRIMMING PIPELINE", font=("Segoe UI", 15, "bold"), fg="#38bdf8", bg="#1e293b")
        title_label.pack(anchor=tk.W)

        subtitle_label = tk.Label(
            header_frame,
            text="Dự đoán vị trí loại bỏ vùng nhiễu (Non-destructive Trim Guidance) & Trực quan hóa điện di đồ đối chiếu rCRS",
            font=("Segoe UI", 9),
            fg="#94a3b8",
            bg="#1e293b"
        )
        subtitle_label.pack(anchor=tk.W, pady=(2, 0))

        # Main Container
        main_frame = ttk.Frame(self.root, padding=16)
        main_frame.pack(fill=tk.BOTH, expand=True)

        # 2. File & Directory Selection Frame
        io_frame = ttk.LabelFrame(main_frame, text="Cấu Hình Dữ Liệu Vào / Ra", padding=12)
        io_frame.pack(fill=tk.X, pady=(0, 10))

        # Input Dir
        default_input = find_app_path(os.path.join("data", "sample_ab1"))
        ttk.Label(io_frame, text="Thư mục chứa file .ab1:").grid(row=0, column=0, sticky=tk.W, pady=4)
        self.input_var = tk.StringVar(value=default_input if os.path.exists(default_input) else "")
        ttk.Entry(io_frame, textvariable=self.input_var, width=54).grid(row=0, column=1, padx=6, pady=4, sticky=tk.EW)
        ttk.Button(io_frame, text="Chọn thư mục...", command=self.browse_input).grid(row=0, column=2, pady=4)

        # rCRS File
        ttk.Label(io_frame, text="File tham chiếu rCRS (.fasta):").grid(row=1, column=0, sticky=tk.W, pady=4)
        default_rcrs = find_app_path(os.path.join("data", "rCRS.fasta"))
        self.rcrs_var = tk.StringVar(value=default_rcrs if os.path.exists(default_rcrs) else "")
        ttk.Entry(io_frame, textvariable=self.rcrs_var, width=54).grid(row=1, column=1, padx=6, pady=4, sticky=tk.EW)
        ttk.Button(io_frame, text="Chọn file...", command=self.browse_rcrs).grid(row=1, column=2, pady=4)

        # Output Dir
        ttk.Label(io_frame, text="Thư mục lưu kết quả:").grid(row=2, column=0, sticky=tk.W, pady=4)
        default_output = os.path.join(EXE_DIR, "output")
        self.output_var = tk.StringVar(value=default_output)
        ttk.Entry(io_frame, textvariable=self.output_var, width=54).grid(row=2, column=1, padx=6, pady=4, sticky=tk.EW)
        ttk.Button(io_frame, text="Chọn thư mục...", command=self.browse_output).grid(row=2, column=2, pady=4)

        io_frame.columnconfigure(1, weight=1)

        # 3. Parameters Frame
        param_frame = ttk.LabelFrame(main_frame, text="Thông Số Thuật Toán Tracy & Cắt Gọt", padding=12)
        param_frame.pack(fill=tk.X, pady=(0, 10))

        # Sliders
        ttk.Label(param_frame, text="Ngưỡng Đỉnh Phụ Tracy (Secondary Ratio):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.thresh_var = tk.DoubleVar(value=0.25)
        self.thresh_slider = ttk.Scale(param_frame, from_=0.10, to=0.50, variable=self.thresh_var, orient=tk.HORIZONTAL)
        self.thresh_slider.grid(row=0, column=1, padx=8, pady=2, sticky=tk.EW)
        self.thresh_lbl = ttk.Label(param_frame, text="0.25", width=6)
        self.thresh_lbl.grid(row=0, column=2, sticky=tk.W)
        self.thresh_slider.configure(command=lambda v: self.thresh_lbl.configure(text=f"{float(v):.2f}"))

        ttk.Label(param_frame, text="Ngưỡng Chồng Peak / Lỗi Sensor (?):").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.conflict_var = tk.DoubleVar(value=0.70)
        self.conflict_slider = ttk.Scale(param_frame, from_=0.50, to=0.90, variable=self.conflict_var, orient=tk.HORIZONTAL)
        self.conflict_slider.grid(row=1, column=1, padx=8, pady=2, sticky=tk.EW)
        self.conflict_lbl = ttk.Label(param_frame, text="0.70", width=6)
        self.conflict_lbl.grid(row=1, column=2, sticky=tk.W)
        self.conflict_slider.configure(command=lambda v: self.conflict_lbl.configure(text=f"{float(v):.2f}"))

        param_frame.columnconfigure(1, weight=1)

        # 4. Action Buttons & Progress Bar
        action_frame = ttk.Frame(main_frame)
        action_frame.pack(fill=tk.X, pady=(0, 8))

        self.run_btn = tk.Button(
            action_frame,
            text="▶ BẮT ĐẦU PHÂN TÍCH HÀNG LOẠT (1-CLICK RUN)",
            font=("Segoe UI", 10, "bold"),
            bg="#0284c7",
            fg="white",
            activebackground="#0369a1",
            activeforeground="white",
            padx=16,
            pady=8,
            cursor="hand2",
            command=self.start_processing
        )
        self.run_btn.pack(side=tk.LEFT)

        self.open_out_btn = ttk.Button(action_frame, text="Mở Thư Mục Kết Quả", command=self.open_output_folder)
        self.open_out_btn.pack(side=tk.RIGHT, padx=4)

        self.open_trim_btn = ttk.Button(action_frame, text="Mở Bảng Khuyến Nghị Cắt (CSV)", command=self.open_trim_table)
        self.open_trim_btn.pack(side=tk.RIGHT, padx=4)

        # 5. Log Console Window
        log_frame = ttk.LabelFrame(main_frame, text="Nhật Ký Xử Lý Trực Tiếp (Real-time Log)", padding=6)
        log_frame.pack(fill=tk.BOTH, expand=True)

        self.log_text = tk.Text(log_frame, bg="#0f172a", fg="#e2e8f0", font=("Consolas", 9), wrap=tk.WORD, height=12)
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        scrollbar = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.log_text.config(yscrollcommand=scrollbar.set)

    def log(self, message: str):
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)

    def browse_input(self):
        d = filedialog.askdirectory(title="Chọn thư mục chứa file .ab1")
        if d:
            self.input_var.set(d)

    def browse_rcrs(self):
        f = filedialog.askopenfilename(title="Chọn file rCRS FASTA", filetypes=[("FASTA files", "*.fasta;*.fa;*.txt"), ("All files", "*.*")])
        if f:
            self.rcrs_var.set(f)

    def browse_output(self):
        d = filedialog.askdirectory(title="Chọn thư mục xuất kết quả")
        if d:
            self.output_var.set(d)

    def open_output_folder(self):
        out_dir = self.output_var.get()
        if os.path.exists(out_dir):
            os.startfile(out_dir)
        else:
            messagebox.showinfo("Thông báo", "Thư mục kết quả chưa được tạo.")

    def open_trim_table(self):
        out_dir = self.output_var.get()
        csv_path = os.path.join(out_dir, "trimming_recommendations.csv")
        if os.path.exists(csv_path):
            os.startfile(csv_path)
        else:
            messagebox.showinfo("Thông báo", "Chưa tìm thấy file trimming_recommendations.csv. Vui lòng bấm bắt đầu phân tích trước.")

    def start_processing(self):
        input_dir = self.input_var.get()
        rcrs_file = self.rcrs_var.get()
        output_dir = self.output_var.get()

        if not os.path.exists(input_dir):
            messagebox.showerror("Lỗi", "Thư mục input không tồn tại!")
            return
        if not os.path.exists(rcrs_file):
            messagebox.showerror("Lỗi", "File rCRS.fasta không tồn tại! Vui lòng kiểm tra lại.")
            return

        self.run_btn.config(state=tk.DISABLED, text="Đang xử lý dữ liệu...")
        self.log_text.delete(1.0, tk.END)

        def worker():
            # Chuyển hướng stdout để hiển thị trên text console
            class TextRedirector:
                def __init__(self, gui_instance):
                    self.gui = gui_instance
                def write(self, str_val):
                    val = str_val.strip()
                    if val:
                        self.gui.root.after(0, self.gui.log, val)
                def flush(self):
                    pass

            old_stdout = sys.stdout
            sys.stdout = TextRedirector(self)
            try:
                run_pipeline(
                    input_dir=input_dir,
                    rcrs_path=rcrs_file,
                    output_dir=output_dir,
                    peak_thresh=self.thresh_var.get(),
                    conflict_thresh=self.conflict_var.get()
                )
                self.root.after(0, lambda: messagebox.showinfo("Thành công", "Đã xử lý xong toàn bộ các mẫu!"))
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Lỗi", f"Có lỗi xảy ra: {e}"))
            finally:
                sys.stdout = old_stdout
                self.root.after(0, lambda: self.run_btn.config(state=tk.NORMAL, text="▶ BẮT ĐẦU PHÂN TÍCH HÀNG LOẠT (1-CLICK RUN)"))

        t = threading.Thread(target=worker, daemon=True)
        t.start()

if __name__ == "__main__":
    root = tk.Tk()
    app = OptimalTrimmingGUI(root)
    root.mainloop()
