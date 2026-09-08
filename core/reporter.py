import os
import json
import csv
from typing import Dict, List, Any, Optional

try:
    import pandas as pd
except ImportError:
    pd = None

class SuiteReporter:
    """
    Module xuất báo cáo cho Optimal Trimming Pipeline:
    1. Xuất file FASTA bảo toàn chuỗi gốc kèm tọa độ cắt đề xuất trong Header.
    2. Xuất Bảng khuyến nghị vị trí cắt chi tiết (trimming_recommendations.csv) cho Sequencher 5.4.6 (Dạng A).
    3. Xuất bảng tổng hợp kết quả đột biến và chất lượng (mtDNA_batch_summary.csv).
    4. Xuất giao diện HTML trực quan hóa đồ thị điện di đồ chromatogram 4 màu thực tế,
       vạch chỉ thị cắt 5'/3', cụm Poly-C, và dòng đối chiếu chuỗi chuẩn rCRS (Dạng C).
    """

    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self.fasta_dir = os.path.join(output_dir, "trimmed_fasta")
        self.html_dir = os.path.join(output_dir, "html_reports")
        os.makedirs(self.fasta_dir, exist_ok=True)
        os.makedirs(self.html_dir, exist_ok=True)

    def export_trimmed_fasta(self, sample_id: str, region: str, direction: str, trim_result: Dict[str, Any]) -> str:
        """
        Xuất file FASTA: Giữ nguyên chuỗi gốc, chú thích dải đề xuất giữ lại trong Header.
        """
        file_name = f"{sample_id}_{region}_{direction}_trimmed.fasta"
        file_path = os.path.join(self.fasta_dir, file_name)

        raw_seq = trim_result.get("raw_sequence", trim_result.get("trimmed_seq", ""))
        keep_5p = trim_result.get("suggested_5p_cut_1based", trim_result.get("trim_start_idx", 0) + 1)
        keep_3p = trim_result.get("suggested_3p_cut_1based", trim_result.get("trim_end_idx", len(raw_seq)))
        r_span = trim_result.get("rcrs_span", (0, 0))
        poly_c = trim_result.get("poly_c_status", "NONE")

        header = (f">{sample_id}_{region}_{direction}|RawLen:{len(raw_seq)}|"
                  f"SuggestedKeepRange:{keep_5p}-{keep_3p}|rCRS:{r_span[0]}-{r_span[1]}|"
                  f"PolyC:{poly_c}|{trim_result.get('trim_reason', '')}")

        with open(file_path, "w", encoding="utf-8") as f:
            f.write(f"{header}\n{raw_seq}\n")
        return file_path

    def export_consensus_fasta(self, sample_id: str, region: str, contig_res: Dict[str, Any]) -> str:
        """Xuất file FASTA chuỗi consensus và so khớp rCRS"""
        file_name = f"{sample_id}_{region}_Consensus.fasta"
        file_path = os.path.join(self.fasta_dir, file_name)

        with open(file_path, "w", encoding="utf-8") as f:
            f.write(f">{sample_id}_{region}_Consensus|Range:{contig_res['target_range']}|QC:{contig_res['qc_status']}\n")
            f.write(f"{contig_res['consensus_sequence']}\n")
            f.write(f">{sample_id}_{region}_rCRS_Match ('.' = Match)\n")
            f.write(f"{contig_res['rcrs_comparison']}\n")
        return file_path

    def export_trim_recommendations_csv(self, recommendations: List[Dict[str, Any]]) -> str:
        """
        DẠNG A: Xuất bảng tổng hợp khuyến nghị vị trí cắt chi tiết cho từng file .ab1
        để kỹ thuật viên mở song song khi thao tác trên Sequencher 5.4.6.
        """
        csv_path = os.path.join(self.output_dir, "trimming_recommendations.csv")
        if pd is not None:
            try:
                df = pd.DataFrame(recommendations)
                df.to_csv(csv_path, index=False, encoding="utf-8-sig")
                return csv_path
            except Exception:
                pass

        if recommendations:
            fieldnames = list(recommendations[0].keys())
            with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(recommendations)
        else:
            with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
                f.write("")
        return csv_path

    def export_summary_csv(self, all_results: List[Dict[str, Any]]) -> str:
        """Xuất bảng tổng kết toàn bộ các mẫu ra file CSV"""
        csv_path = os.path.join(self.output_dir, "mtDNA_batch_summary.csv")
        if pd is not None:
            try:
                df = pd.DataFrame(all_results)
                df.to_csv(csv_path, index=False, encoding="utf-8-sig")
                return csv_path
            except Exception:
                pass

        if all_results:
            fieldnames = list(all_results[0].keys())
            with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_results)
        else:
            with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
                f.write("")
        return csv_path

    def export_chromatogram_html(self,
                                 sample_id: str,
                                 region: str,
                                 direction: str,
                                 parsed_ab1: Dict[str, Any],
                                 analyzed_data: Dict[str, Any],
                                 trim_res: Dict[str, Any],
                                 alignment_info: Optional[Dict[str, Any]],
                                 rcrs_seq: str,
                                 contig_res: Optional[Dict[str, Any]] = None) -> str:
        """
        DẠNG C: Trực quan hóa điện di đồ Chromatogram 4 màu thực tế,
        vạch ranh giới cắt 5' và 3', cụm Poly-C, và chuỗi tham chiếu rCRS chạy song song.
        """
        html_file = os.path.join(self.html_dir, f"{sample_id}_{region}_{direction}_chromatogram.html")

        # Chuẩn bị dữ liệu tích hợp rCRS cho từng base
        query_to_rcrs = alignment_info.get("query_to_rcrs", {}) if alignment_info else {}
        bases_info = analyzed_data.get("bases", [])

        viewer_bases = []
        for i, b in enumerate(bases_info):
            called = b.get("called_base", b.get("primary_base", "N"))
            ploc = b.get("trace_pos", 0)
            q = b.get("phred_score", 10)
            r_coord = query_to_rcrs.get(i)

            ref_char = ""
            is_match = False
            is_variant = False

            if r_coord is not None and 1 <= r_coord <= len(rcrs_seq):
                ref_char = rcrs_seq[r_coord - 1]
                # Nếu là Reverse read, so với chuỗi đối ứng đã đảo bổ sung
                check_char = called
                if direction.upper() in ["R", "REV", "REVERSE"]:
                    comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                    check_char = comp_map.get(called, called)

                if check_char == ref_char:
                    is_match = True
                elif check_char in ["A", "C", "G", "T"] and ref_char in ["A", "C", "G", "T"]:
                    is_variant = True

            viewer_bases.append({
                "idx": i + 1,
                "ploc": ploc,
                "called": called,
                "q": q,
                "ratio": b.get("secondary_ratio", 0.0),
                "status": b.get("status", "CLEAN"),
                "probs": b.get("probabilities", {}),
                "is_conflict": b.get("is_conflict", False),
                "is_post_polyc": b.get("is_post_polyc", False),
                "rcrs_coord": r_coord,
                "rcrs_base": ref_char,
                "is_match": is_match,
                "is_variant": is_variant
            })

        # Dữ liệu đóng gói gửi sang client
        payload = {
            "sample_id": sample_id,
            "region": region,
            "direction": direction,
            "file_name": parsed_ab1.get("file_name", ""),
            "traces": parsed_ab1.get("channel_traces", {}),
            "bases": viewer_bases,
            "trim": {
                "suggested_5p_cut_1based": trim_res.get("suggested_5p_cut_1based", 1),
                "suggested_3p_cut_1based": trim_res.get("suggested_3p_cut_1based", len(viewer_bases)),
                "retained_length": trim_res.get("retained_length", len(viewer_bases)),
                "rcrs_span": trim_res.get("rcrs_span", (0, 0)),
                "poly_c_status": trim_res.get("poly_c_status", "NONE"),
                "noise_5p_reason": trim_res.get("noise_5p_reason", ""),
                "noise_3p_reason": trim_res.get("noise_3p_reason", "")
            },
            "qc_status": contig_res.get("qc_status", "N/A") if contig_res else "N/A",
            "variants_count": contig_res.get("variant_count", 0) if contig_res else 0
        }

        def _safe_json(o):
            if hasattr(o, 'item'):
                return o.item()
            if hasattr(o, '__iter__') and not isinstance(o, (str, bytes, dict)):
                return list(o)
            return str(o)

        json_data = json.dumps(payload, default=_safe_json)

        qc_status = contig_res.get("qc_status", "N/A") if contig_res else "N/A"
        qc_recommendation = contig_res.get("qc_recommendation", "") if contig_res else ""
        qc_details = contig_res.get("qc_details", "") if contig_res else ""

        deny_alert_html = ""
        qc_badge_html = ""
        if qc_status == "DENY":
            deny_alert_html = f"""
    <!-- DENY ALERT BANNER -->
    <div class="bg-rose-950/90 border-2 border-rose-500 rounded-xl p-4 shadow-xl flex items-start gap-3.5">
      <span class="text-2xl">🚫</span>
      <div class="flex-1">
        <div class="flex items-center gap-2">
          <span class="px-2.5 py-0.5 rounded bg-rose-600 text-white font-bold text-xs uppercase tracking-wider">CẢNH BÁO QC: DENY</span>
          <h2 class="text-sm font-bold text-rose-200">KHUYẾN CÁO PHÒNG LAB: KHÔNG THỂ PHÂN TÍCH — CẦN LÀM LẠI MẪU</h2>
        </div>
        <p class="text-xs text-rose-300 mt-1.5"><b>Lý do:</b> {qc_details}</p>
        <p class="text-[11px] text-rose-400 mt-0.5 italic">Các đỉnh sóng chồng chéo không rõ nguyên nhân hoặc dải giữ lại không phủ kín dải mục tiêu. Không được đưa vào báo cáo!</p>
      </div>
    </div>
"""
            qc_badge_html = f'<span class="px-2.5 py-1 rounded bg-rose-600/30 text-rose-300 border border-rose-500 font-bold">QC: DENY ⚠️</span>'
        elif qc_status == "PASS":
            qc_badge_html = f'<span class="px-2.5 py-1 rounded bg-emerald-500/20 text-emerald-400 border border-emerald-500/40 font-bold">QC: PASS ✓</span>'
        else:
            qc_badge_html = f'<span class="px-2.5 py-1 rounded bg-amber-500/20 text-amber-400 border border-amber-500/40 font-bold">QC: REVIEW 👁️</span>'

        html_content = f"""<!DOCTYPE html>
<html lang="vi">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{sample_id} - {region} ({direction}) Chromatogram Inspector</title>
  <style>
    :root {{
      --bg-main: #090d16;
      --bg-card: #111827;
      --bg-card-sub: #1e293b;
      --border-color: #334155;
      --text-main: #f8fafc;
      --text-sub: #94a3b8;
      --accent: #38bdf8;
      --accent-hover: #0284c7;
      --canvas-bg: #090d16;
      --hud-bg: rgba(17, 24, 39, 0.94);
      --hud-border: #38bdf8;
      --font-mono: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
    }}
    body.light {{
      --bg-main: #f1f5f9;
      --bg-card: #ffffff;
      --bg-card-sub: #f8fafc;
      --border-color: #cbd5e1;
      --text-main: #0f172a;
      --text-sub: #64748b;
      --accent: #0284c7;
      --accent-hover: #0369a1;
      --canvas-bg: #ffffff;
      --hud-bg: rgba(255, 255, 255, 0.96);
      --hud-border: #0284c7;
    }}
    * {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background-color: var(--bg-main);
      color: var(--text-main);
      padding: 10px 14px;
      line-height: 1.4;
      transition: background-color 0.2s ease, color 0.2s ease;
      overflow-x: hidden;
    }}
    .container {{ max-width: 1440px; margin: 0 auto; display: flex; flex-direction: column; gap: 8px; }}
    .card {{
      background: var(--bg-card);
      border: 1px solid var(--border-color);
      border-radius: 8px;
      padding: 10px 14px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
    }}
    .btn {{
      display: inline-flex; align-items: center; justify-content: center; gap: 4px;
      padding: 5px 10px; font-size: 11px; font-weight: 600; border-radius: 6px;
      border: 1px solid var(--border-color); background: var(--bg-card-sub);
      color: var(--text-main); cursor: pointer; user-select: none; transition: all 0.15s ease;
      white-space: nowrap;
    }}
    .btn:hover {{ background: var(--border-color); }}
    .btn.active {{ border-color: var(--accent); background: var(--accent); color: #ffffff !important; }}
    .badge {{
      display: inline-flex; align-items: center; padding: 3px 8px; border-radius: 6px;
      font-size: 11px; font-weight: bold;
    }}
    .badge-pass {{ background: rgba(16, 185, 129, 0.18); color: #10b981; border: 1px solid rgba(16, 185, 129, 0.4); }}
    .badge-review {{ background: rgba(245, 158, 11, 0.18); color: #f59e0b; border: 1px solid rgba(245, 158, 11, 0.4); }}
    .badge-deny {{ background: rgba(239, 68, 68, 0.18); color: #ef4444; border: 1px solid rgba(239, 68, 68, 0.4); }}
    
    /* Mini-Map Styles */
    .minimap-container {{
      position: relative; width: 100%; height: 38px; background: var(--canvas-bg);
      border: 1px solid var(--border-color); border-radius: 6px; overflow: hidden; cursor: pointer;
    }}
    #minimap-canvas {{ width: 100%; height: 100%; display: block; }}
    .minimap-lens {{
      position: absolute; top: 0; height: 100%;
      background: rgba(56, 189, 248, 0.22); border: 2px solid var(--accent);
      border-radius: 4px; cursor: grab; pointer-events: auto;
    }}
    .minimap-lens:active {{ cursor: grabbing; }}

    /* Main Viewport Single-Canvas */
    .viewport-box {{
      position: relative; width: 100%; height: 380px; background: var(--canvas-bg);
      border: 1px solid var(--border-color); border-radius: 8px; overflow-x: auto; overflow-y: hidden;
      scroll-behavior: auto;
    }}
    #chroma-canvas {{ display: block; background: var(--canvas-bg); cursor: crosshair; }}
    
    /* Control Toolbar */
    .toolbar {{ display: flex; flex-wrap: wrap; align-items: center; justify-content: space-between; gap: 8px; font-size: 11px; }}
    .tool-group {{ display: flex; align-items: center; gap: 5px; }}
    .slider-wrap {{ display: flex; align-items: center; gap: 5px; }}
    .slider-wrap input[type=range] {{ width: 75px; cursor: pointer; accent-color: var(--accent); }}
    
    /* Floating Crosshair HUD */
    #floating-hud {{
      position: absolute; display: none; z-index: 50; pointer-events: none;
      background: var(--hud-bg); border: 1.5px solid var(--hud-border);
      border-radius: 8px; padding: 8px 12px; font-size: 11px;
      box-shadow: 0 8px 24px rgba(0,0,0,0.4);
      backdrop-filter: blur(8px);
      min-width: 260px;
    }}

    /* Toast Alert */
    #toast {{
      position: fixed; bottom: 20px; right: 20px; z-index: 100;
      background: #0284c7; color: #ffffff; padding: 9px 16px;
      border-radius: 6px; font-size: 12px; font-weight: bold;
      box-shadow: 0 4px 16px rgba(0,0,0,0.3);
      opacity: 0; transform: translateY(10px); transition: all 0.25s ease;
      pointer-events: none;
    }}
    #toast.show {{ opacity: 1; transform: translateY(0); }}

    /* Custom Scrollbar */
    ::-webkit-scrollbar {{ height: 9px; width: 8px; }}
    ::-webkit-scrollbar-track {{ background: var(--bg-card); }}
    ::-webkit-scrollbar-thumb {{ background: var(--border-color); border-radius: 5px; }}
    ::-webkit-scrollbar-thumb:hover {{ background: #64748b; }}
  </style>
</head>
<body class="dark">
  <div class="container">

    {deny_alert_html}

    <!-- Header Card -->
    <div class="card" style="display: flex; flex-wrap: wrap; align-items: center; justify-content: space-between; gap: 10px;">
      <div style="display: flex; align-items: center; gap: 10px;">
        <span style="font-size: 24px;">🧬</span>
        <div>
          <div style="display: flex; align-items: center; gap: 8px;">
            <h1 style="font-size: 16px; font-weight: bold; letter-spacing: 0.5px;">{sample_id}</h1>
            <span style="font-family: var(--font-mono); font-size: 12px; color: var(--accent); font-weight: bold;">[{region} - Mồi {direction}]</span>
            <span style="font-size: 10px; padding: 2px 6px; border-radius: 4px; background: var(--bg-card-sub); color: var(--text-sub); border: 1px solid var(--border-color);">SeqStudio Flex 24</span>
            {qc_badge_html}
          </div>
          <p style="font-size: 11px; color: var(--text-sub); font-family: var(--font-mono); margin-top: 1px;">File: {parsed_ab1.get("file_name", "")}</p>
        </div>
      </div>

      <!-- Quick Metrics & Sequencher Actions -->
      <div style="display: flex; flex-wrap: wrap; align-items: center; gap: 6px;">
        <!-- 1-Click Copy Sequencher Coordinates -->
        <button id="btn-copy-sequencher" class="btn" style="background: #0284c7; color: #ffffff; border-color: #0369a1; font-weight: bold;" title="Sao chép tọa độ cắt 5' và 3' để nhập thẳng vào Sequencher 5.4.6">
          📋 Sao Chép Vị Trí Cắt Sequencher
        </button>

        <div class="badge" style="background: var(--bg-card-sub); border: 1px solid var(--border-color);">
          <span style="color: var(--text-sub);">Giữ lại:</span>
          <b style="color: #10b981; margin-left: 4px;">Base {trim_res.get('suggested_5p_cut_1based', 1)} → {trim_res.get('suggested_3p_cut_1based', len(viewer_bases))}</b>
          <span style="color: var(--text-sub); margin-left: 3px; font-size: 10px;">({trim_res.get('retained_length', len(viewer_bases))} bp)</span>
        </div>

        <div class="badge" style="background: var(--bg-card-sub); border: 1px solid var(--border-color);">
          <span style="color: var(--text-sub);">rCRS:</span>
          <b style="color: #38bdf8; margin-left: 4px;">{trim_res.get('rcrs_span', (0, 0))[0]} - {trim_res.get('rcrs_span', (0, 0))[1]}</b>
        </div>

        <div class="badge" style="background: rgba(168, 85, 247, 0.15); color: #c084fc; border: 1px solid rgba(168, 85, 247, 0.4);">
          <span>Poly-C:</span>
          <b style="margin-left: 4px;">{trim_res.get('poly_c_status', 'NONE')}</b>
        </div>

        <button id="btn-theme-toggle" class="btn" title="Chuyển đổi Dark Mode / Light Mode">☀️ Light</button>
        <button id="btn-snapshot" class="btn" style="border-color: var(--border-color);">📸 Chụp PNG</button>
      </div>
    </div>

    <!-- Rationales & Recommendations Strip -->
    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 6px;">
      <div class="card" style="padding: 6px 12px; font-size: 11px; display: flex; align-items: center;">
        <b style="color: #f43f5e; margin-right: 6px;">✂ Ranh giới 5':</b>
        <span style="color: var(--text-sub);">{trim_res.get('noise_5p_reason', 'Loại bỏ primer ngoài và nhiễu mao quản sớm')}</span>
      </div>
      <div class="card" style="padding: 6px 12px; font-size: 11px; display: flex; align-items: center;">
        <b style="color: #f59e0b; margin-right: 6px;">✂ Ranh giới 3':</b>
        <span style="color: var(--text-sub);">{trim_res.get('noise_3p_reason', 'Giới hạn theo khung rCRS và sụt giảm huỳnh quang')}</span>
      </div>
    </div>

    <!-- Mini-Map Card (Tổng quan toàn cảnh 0-100%) -->
    <div class="card" style="padding: 6px 12px;">
      <div style="display: flex; align-items: center; justify-content: space-between; margin-bottom: 3px; font-size: 10px; color: var(--text-sub); font-weight: bold;">
        <span>🗺️ ĐỒ THỊ TOÀN CẢNH (MINI-MAP OVERVIEW)</span>
        <span>(Click hoặc kéo khung xanh để định vị nhanh tức thì)</span>
      </div>
      <div class="minimap-container" id="minimap-wrapper">
        <canvas id="minimap-canvas" height="38"></canvas>
        <div class="minimap-lens" id="minimap-lens"></div>
      </div>
    </div>

    <!-- Main Trace Viewer Card (Unified Single-Canvas) -->
    <div class="card" style="display: flex; flex-direction: column; gap: 6px; padding: 8px 12px;">
      
      <!-- Toolbar Controls -->
      <div class="toolbar">
        <!-- Channel Filters -->
        <div class="tool-group">
          <span style="font-weight: bold; margin-right: 2px; color: var(--text-sub);">Kênh:</span>
          <button id="ch-A" class="btn active" style="border-color: #10b981; color: #10b981;" title="Phím tắt: A">● A</button>
          <button id="ch-C" class="btn active" style="border-color: #2563eb; color: #2563eb;" title="Phím tắt: C">● C</button>
          <button id="ch-G" class="btn active" style="border-color: #f59e0b; color: #f59e0b;" title="Phím tắt: G">● G</button>
          <button id="ch-T" class="btn active" style="border-color: #ef4444; color: #ef4444;" title="Phím tắt: T">● T</button>
        </div>

        <!-- Zoom & Amplitude Gain -->
        <div class="tool-group">
          <!-- Horizontal Zoom -->
          <div class="slider-wrap">
            <span style="color: var(--text-sub);">Thu/Phóng X:</span>
            <button id="btn-zoom-out" class="btn" style="padding: 2px 6px;">-</button>
            <input type="range" id="zoom-range" min="0.6" max="3.5" step="0.05" value="1.0">
            <button id="btn-zoom-in" class="btn" style="padding: 2px 6px;">+</button>
          </div>

          <div style="width: 1px; height: 16px; background: var(--border-color); margin: 0 4px;"></div>

          <!-- Vertical Amplitude Gain -->
          <div class="slider-wrap">
            <span style="color: var(--text-sub);" title="Tăng chiều cao sóng để soi rõ dị thể hoặc đỉnh thấp">Biên độ Y:</span>
            <input type="range" id="gain-range" min="0.5" max="6.0" step="0.25" value="1.0">
            <button id="btn-auto-gain" class="btn" style="padding: 2px 6px; color: var(--accent); font-weight: bold;" title="Tự động cân bằng độ cao sóng dựa trên vùng đang xem">⚡ Auto</button>
            <button id="btn-gain-reset" class="btn" style="padding: 2px 6px;" title="Đặt lại độ cao mặc định 1x">1x</button>
          </div>
        </div>

        <!-- Smart Clinical Jumps & Search -->
        <div class="tool-group">
          <button id="btn-jump-5p" class="btn" style="color: #f43f5e; border-color: #f43f5e;">✂ 5' Cắt</button>
          <button id="btn-jump-polyc" class="btn" style="color: #c084fc; border-color: #c084fc;">🧬 Poly-C</button>
          <button id="btn-jump-3p" class="btn" style="color: #f59e0b; border-color: #f59e0b;">✂ 3' Cắt</button>
          <button id="btn-jump-snp" class="btn" style="color: #38bdf8; border-color: #38bdf8;" title="Nhảy qua các vị trí biến dị SNPs (Phím Space)">⚡ SNP</button>
          <button id="btn-jump-hetero" class="btn" style="color: #eab308; border-color: #eab308;" title="Nhảy qua các vị trí có đỉnh phụ >= 15%">⚠️ Đỉnh Kép</button>
          
          <input type="text" id="input-search" placeholder="Nhảy tới #150 hoặc 16189..." 
                 style="padding: 4px 8px; font-size: 11px; border-radius: 6px; border: 1px solid var(--border-color); background: var(--bg-card-sub); color: var(--text-main); width: 145px;">
        </div>
      </div>

      <!-- Viewport Canvas Container (Single Unified Canvas) -->
      <div class="viewport-box" id="scroll-wrapper">
        <canvas id="chroma-canvas" height="380"></canvas>
        
        <!-- Floating HUD Tooltip -->
        <div id="floating-hud">
          <div id="hud-title" style="font-weight: bold; font-family: var(--font-mono); font-size: 12px; margin-bottom: 4px; display: flex; justify-content: space-between;"></div>
          <div id="hud-rcrs" style="margin-bottom: 4px; font-size: 11px;"></div>
          <div id="hud-quality" style="margin-bottom: 4px; font-size: 11px;"></div>
          <div id="hud-rfu" style="font-family: var(--font-mono); font-size: 10px; color: var(--text-sub); border-top: 1px solid var(--border-color); padding-top: 4px;"></div>
        </div>
      </div>

      <!-- Footer Help Hints -->
      <div style="display: flex; align-items: center; justify-content: space-between; font-size: 10px; color: var(--text-sub); padding: 2px 4px;">
        <div>
          <span>💡 Phím tắt: <b>← / →</b> cuộn dải | <b>Ctrl + Cuộn chuột</b> thu phóng | <b>A, C, G, T</b> bật/tắt kênh | <b>Space</b> nhảy SNP tiếp theo</span>
        </div>
        <div>
          <span>Thước đo: <span style="color: #10b981;">● Base Index</span> | <span style="color: #38bdf8;">● rCRS Coordinates</span> | <span style="color: #ef4444;">✂ 5'/3' Trim Boundaries</span></span>
        </div>
      </div>

    </div>

    <!-- Toast Notification Element -->
    <div id="toast"></div>

  </div>

  <script>
    const DATA = {json_data};

    const traces = DATA.traces || {{}};
    const bases = DATA.bases || [];
    const trim = DATA.trim || {{}};

    // State Variables
    let zoomScale = 1.0;
    let gainScale = 1.0;
    const channelsActive = {{ 'A': true, 'C': true, 'G': true, 'T': true }};

    const BASE_COLORS = {{
      'A': '#10b981', // Emerald Green
      'C': '#2563eb', // Royal Blue
      'G': '#f59e0b', // Amber Gold
      'T': '#ef4444', // Crimson Red
      'N': '#64748b'
    }};

    // DOM Elements
    const canvas = document.getElementById('chroma-canvas');
    const ctx = canvas.getContext('2d');
    const minimapCanvas = document.getElementById('minimap-canvas');
    const miniCtx = minimapCanvas.getContext('2d');
    const minimapWrapper = document.getElementById('minimap-wrapper');
    const minimapLens = document.getElementById('minimap-lens');
    const scrollWrapper = document.getElementById('scroll-wrapper');
    const floatingHud = document.getElementById('floating-hud');
    const toast = document.getElementById('toast');

    // Scan points calculation
    const numScanPoints = (traces['A'] && traces['A'].length) || (bases.length * 12);
    const plocs = bases.map(b => b.ploc);
    const minPloc = plocs[0] || 0;
    const maxPloc = plocs[plocs.length - 1] || numScanPoints;

    // Maximum RFU across channels (exclude extreme artifacts)
    let maxRfu = 100;
    ['A', 'C', 'G', 'T'].forEach(c => {{
      if (traces[c]) {{
        const start = Math.floor(traces[c].length * 0.05);
        const end = Math.floor(traces[c].length * 0.95);
        for (let i = start; i < end; i += 4) {{
          if (traces[c][i] > maxRfu) maxRfu = traces[c][i];
        }}
      }}
    }});
    if (maxRfu <= 100) maxRfu = 1000;

    // Mouse & Hover state
    let hoveredBase = null;
    let mouseCanvasX = -1;

    function getBaseSpacing() {{
      return 22 * zoomScale;
    }}

    function getX(scanIdx) {{
      const baseSpacing = getBaseSpacing();
      return 45 + ((scanIdx - minPloc) / (maxPloc - minPloc)) * (bases.length * baseSpacing);
    }}

    function getScanFromX(x) {{
      const baseSpacing = getBaseSpacing();
      const frac = (x - 45) / (bases.length * baseSpacing);
      return minPloc + frac * (maxPloc - minPloc);
    }}

    // Toast Alert Helper
    function showToast(msg) {{
      toast.textContent = msg;
      toast.classList.add('show');
      setTimeout(() => toast.classList.remove('show'), 2800);
    }}

    // Render Mini-Map Overview (Run once or on window resize)
    function renderMinimap() {{
      const w = minimapWrapper.clientWidth;
      const h = minimapCanvas.height = 38;
      minimapCanvas.width = w;
      miniCtx.clearRect(0, 0, w, h);

      if (bases.length === 0) return;

      // 1. Shaded trim 5' & 3'
      const cut5p = Math.max(0, (trim.suggested_5p_cut_1based - 1) / bases.length * w);
      const cut3p = Math.min(w, (trim.suggested_3p_cut_1based - 1) / bases.length * w);

      miniCtx.fillStyle = 'rgba(239, 68, 68, 0.22)';
      miniCtx.fillRect(0, 0, cut5p, h);

      miniCtx.fillStyle = 'rgba(245, 158, 11, 0.22)';
      miniCtx.fillRect(cut3p, 0, w - cut3p, h);

      // 2. Draw compressed waveforms
      ['A', 'C', 'G', 'T'].forEach(ch => {{
        if (!channelsActive[ch] || !traces[ch]) return;
        miniCtx.strokeStyle = BASE_COLORS[ch];
        miniCtx.lineWidth = 1;
        miniCtx.beginPath();
        for (let x = 0; x < w; x++) {{
          const scanIdx = Math.floor(x / w * traces[ch].length);
          const val = traces[ch][scanIdx] || 0;
          const y = h - 2 - (val / (maxRfu * 1.1)) * (h - 4);
          if (x === 0) miniCtx.moveTo(x, y);
          else miniCtx.lineTo(x, y);
        }}
        miniCtx.stroke();
      }});

      // 3. Mark SNPs as vibrant cyan markers
      bases.forEach((b, i) => {{
        if (b.is_variant) {{
          const x = (i / bases.length) * w;
          miniCtx.fillStyle = '#38bdf8';
          miniCtx.beginPath();
          miniCtx.arc(x, 5, 2.5, 0, Math.PI * 2);
          miniCtx.fill();
        }}
      }});

      updateMinimapLens();
    }}

    function updateMinimapLens() {{
      const totalWidth = canvas.width;
      const viewWidth = scrollWrapper.clientWidth;
      const scrollLeft = scrollWrapper.scrollLeft;

      const minimapW = minimapWrapper.clientWidth;
      const lensW = Math.max(16, (viewWidth / totalWidth) * minimapW);
      const lensX = (scrollLeft / totalWidth) * minimapW;

      minimapLens.style.width = lensW + 'px';
      minimapLens.style.left = Math.min(lensX, minimapW - lensW) + 'px';
    }}

    // =========================================================================
    // UNIFIED SINGLE-CANVAS RENDERER (60 FPS, SẮC NÉT, KHÔNG CHỒNG LẤN CHỮ)
    // =========================================================================
    function renderMain() {{
      const baseSpacing = getBaseSpacing();
      const totalWidth = Math.max(bases.length * baseSpacing + 120, scrollWrapper.clientWidth);

      // Handle DevicePixelRatio for Crisp Retina Rendering
      const dpr = window.devicePixelRatio || 1;
      canvas.width = totalWidth * dpr;
      canvas.height = 380 * dpr;
      canvas.style.width = totalWidth + 'px';
      canvas.style.height = '380px';

      ctx.resetTransform ? ctx.resetTransform() : ctx.setTransform(1, 0, 0, 1, 0, 0);
      ctx.scale(dpr, dpr);

      const isLight = document.body.classList.contains('light');

      // Colors configuration
      const colTextMain = isLight ? '#0f172a' : '#f8fafc';
      const colTextSub = isLight ? '#64748b' : '#94a3b8';
      const colBorder = isLight ? '#cbd5e1' : '#334155';
      const colBaseline = isLight ? '#e2e8f0' : '#1e293b';
      const colOverlay = isLight ? 'rgba(226, 232, 240, 0.72)' : 'rgba(15, 23, 42, 0.76)';

      ctx.clearRect(0, 0, totalWidth, 380);

      // Section Geometry:
      // Y: 0 -> 24: Top Base Index Ruler
      // Y: 25 -> 54: Called Base Letters
      // Y: 55 -> 72: Phred Q-score Micro-Bars
      // Y: 75 -> 315: Chromatogram Waveform Area (Baseline at Y: 305)
      // Y: 320 -> 375: rCRS Reference Track & Major Coordinates

      const waveTop = 75;
      const waveHeight = 230;
      const waveBaseY = waveTop + waveHeight;

      // 1. Shaded 5' and 3' Trim Zones
      const cut5pIdx = trim.suggested_5p_cut_1based - 1;
      const cut3pIdx = trim.suggested_3p_cut_1based - 1;
      const x5p = cut5pIdx >= 0 && cut5pIdx < bases.length ? getX(bases[cut5pIdx].ploc) : 0;
      const x3p = cut3pIdx >= 0 && cut3pIdx < bases.length ? getX(bases[cut3pIdx].ploc) : totalWidth;

      ctx.fillStyle = colOverlay;
      ctx.fillRect(0, 0, x5p, 380);
      ctx.fillRect(x3p, 0, totalWidth - x3p, 380);

      // Trim Boundary Lines & Flags
      ctx.setLineDash([5, 4]);
      ctx.lineWidth = 1.5;
      
      // 5' Boundary
      ctx.strokeStyle = '#f43f5e';
      ctx.beginPath(); ctx.moveTo(x5p, 0); ctx.lineTo(x5p, 380); ctx.stroke();

      // 3' Boundary
      ctx.strokeStyle = '#f59e0b';
      ctx.beginPath(); ctx.moveTo(x3p, 0); ctx.lineTo(x3p, 380); ctx.stroke();
      ctx.setLineDash([]);

      // Trim Labels
      ctx.font = 'bold 10px sans-serif';
      ctx.fillStyle = '#f43f5e';
      ctx.fillText(`✂ 5' CẮT (#${{trim.suggested_5p_cut_1based}})`, Math.max(10, x5p - 80), 16);

      ctx.fillStyle = '#f59e0b';
      ctx.fillText(`✂ 3' CẮT (#${{trim.suggested_3p_cut_1based}})`, x3p + 8, 16);

      // 2. Track Dividers & Baseline
      ctx.strokeStyle = colBorder;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(0, 24); ctx.lineTo(totalWidth, 24); // Under Ruler
      ctx.moveTo(0, 72); ctx.lineTo(totalWidth, 72); // Under Phred Bars
      ctx.moveTo(0, 318); ctx.lineTo(totalWidth, 318); // Above rCRS
      ctx.stroke();

      // Soft Baseline
      ctx.strokeStyle = colBaseline;
      ctx.beginPath(); ctx.moveTo(0, waveBaseY); ctx.lineTo(totalWidth, waveBaseY); ctx.stroke();

      // 3. Draw Waveforms (A, C, G, T)
      ['A', 'C', 'G', 'T'].forEach(ch => {{
        if (!channelsActive[ch] || !traces[ch]) return;

        ctx.strokeStyle = BASE_COLORS[ch];
        ctx.lineWidth = 1.7;
        ctx.lineJoin = 'round';
        ctx.beginPath();

        let started = false;
        const step = zoomScale > 1.3 ? 1 : 2;
        for (let i = 0; i < traces[ch].length; i += step) {{
          const x = getX(i);
          if (x < -20 || x > totalWidth + 20) continue;

          const val = traces[ch][i] * gainScale;
          const y = waveBaseY - (val / (maxRfu * 1.1)) * (waveHeight - 15);

          if (!started) {{ ctx.moveTo(x, y); started = true; }}
          else {{ ctx.lineTo(x, y); }}
        }}
        ctx.stroke();
      }});

      // 4. Draw Basecalls, Phred Bars, and Smart rCRS Track
      const showEveryBaseIndex = zoomScale >= 1.6;
      const stepTick = zoomScale >= 1.2 ? 5 : 10;

      bases.forEach((b, idx) => {{
        const x = getX(b.ploc);
        if (x < -30 || x > totalWidth + 30) return;

        // A. Base Index Number (Y: 10 -> 20)
        if (showEveryBaseIndex || b.idx === 1 || b.idx % stepTick === 0) {{
          ctx.font = '9px var(--font-mono)';
          ctx.fillStyle = colTextSub;
          ctx.textAlign = 'center';
          ctx.fillText(b.idx, x, 16);

          // Small tick
          ctx.strokeStyle = colBorder;
          ctx.beginPath(); ctx.moveTo(x, 19); ctx.lineTo(x, 24); ctx.stroke();
        }}

        // B. Primary Base Letter (Y: 46)
        let baseColor = BASE_COLORS[b.called] || colTextSub;
        ctx.font = 'bold 15px var(--font-mono)';
        ctx.textAlign = 'center';

        // Check if conflict / ambiguous
        if (b.is_conflict) {{
          // Ambiguous pill background
          ctx.fillStyle = '#9333ea';
          ctx.fillRect(x - 8, 30, 16, 18);
          ctx.fillStyle = '#ffffff';
          ctx.fillText(b.called, x, 44);
        }} else {{
          ctx.fillStyle = baseColor;
          ctx.fillText(b.called, x, 44);
        }}

        // C. Phred Quality (Q) Micro-Bar (Y: 56 -> 70)
        const qH = Math.min(14, Math.max(3, (b.q / 60) * 14));
        const qY = 70 - qH;
        ctx.fillStyle = b.q >= 30 ? '#10b981' : (b.q >= 20 ? '#38bdf8' : '#ef4444');
        ctx.fillRect(x - 3, qY, 6, qH);

        // D. rCRS Reference Alignment (Y: 335 -> 370)
        if (b.rcrs_coord) {{
          if (b.is_variant) {{
            // Variant SNP Pill
            ctx.fillStyle = 'rgba(56, 189, 248, 0.22)';
            ctx.strokeStyle = '#38bdf8';
            ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.roundRect(x - 9, 324, 18, 18, 4);
            ctx.fill();
            ctx.stroke();

            ctx.font = 'bold 12px var(--font-mono)';
            ctx.fillStyle = '#38bdf8';
            ctx.textAlign = 'center';
            ctx.fillText(b.rcrs_base, x, 337);

            // Print exact coordinate under variant
            ctx.font = 'bold 9px var(--font-mono)';
            ctx.fillStyle = '#38bdf8';
            ctx.fillText(b.rcrs_coord, x, 355);
          }} else {{
            // Match: subtle dot (·)
            ctx.font = 'bold 13px sans-serif';
            ctx.fillStyle = colTextSub;
            ctx.textAlign = 'center';
            ctx.fillText('·', x, 336);

            // Print rCRS coordinate only every 10 bp so numbers NEVER overlap!
            if (b.rcrs_coord % 10 === 0 || b.idx === 1) {{
              ctx.font = '9px var(--font-mono)';
              ctx.fillStyle = colTextSub;
              ctx.fillText(b.rcrs_coord, x, 355);

              ctx.strokeStyle = colBorder;
              ctx.beginPath(); ctx.moveTo(x, 342); ctx.lineTo(x, 346); ctx.stroke();
            }}
          }}
        }} else {{
          // Outside rCRS target span
          ctx.font = '9px sans-serif';
          ctx.fillStyle = colTextSub;
          ctx.textAlign = 'center';
          ctx.fillText('-', x, 336);
        }}
      }});

      // 5. Draw Crosshair & Column Highlight on Hover
      if (mouseCanvasX >= 0 && hoveredBase) {{
        const hx = getX(hoveredBase.ploc);

        // Column highlight
        ctx.fillStyle = isLight ? 'rgba(2, 132, 199, 0.08)' : 'rgba(56, 189, 248, 0.12)';
        ctx.fillRect(hx - baseSpacing / 2, 0, baseSpacing, 380);

        // Vertical Hairline Crosshair
        ctx.strokeStyle = 'rgba(56, 189, 248, 0.85)';
        ctx.lineWidth = 1.2;
        ctx.setLineDash([4, 3]);
        ctx.beginPath(); ctx.moveTo(hx, 0); ctx.lineTo(hx, 380); ctx.stroke();
        ctx.setLineDash([]);
      }}

      updateMinimapLens();
    }}

    // =========================================================================
    // INTERACTIVE FLOATING HUD & MOUSE INSPECTOR
    // =========================================================================
    function findNearestBase(canvasX) {{
      let nearest = null;
      let minDiff = Infinity;
      for (let i = 0; i < bases.length; i++) {{
        const bx = getX(bases[i].ploc);
        const diff = Math.abs(bx - canvasX);
        if (diff < minDiff) {{
          minDiff = diff;
          nearest = bases[i];
        }}
      }}
      return minDiff < getBaseSpacing() * 1.5 ? nearest : null;
    }}

    canvas.addEventListener('mousemove', (e) => {{
      const rect = canvas.getBoundingClientRect();
      mouseCanvasX = e.clientX - rect.left;
      hoveredBase = findNearestBase(mouseCanvasX);

      renderMain();

      if (hoveredBase) {{
        showHud(e.clientX, e.clientY, hoveredBase);
      }} else {{
        floatingHud.style.display = 'none';
      }}
    }});

    canvas.addEventListener('mouseleave', () => {{
      mouseCanvasX = -1;
      hoveredBase = null;
      floatingHud.style.display = 'none';
      renderMain();
    }});

    function showHud(clientX, clientY, b) {{
      const p = b.probs || {{}};
      const qColor = b.q >= 30 ? '#10b981' : (b.q >= 20 ? '#38bdf8' : '#ef4444');
      const qText = b.q >= 30 ? 'Rất cao' : (b.q >= 20 ? 'Trung bình' : 'Kém');

      // Title line
      const titleEl = document.getElementById('hud-title');
      titleEl.innerHTML = `
        <span style="color: ${{BASE_COLORS[b.called] || '#94a3b8'}};">Base #${{b.idx}} [${{b.called}}]</span>
        <span style="color: var(--text-sub); font-size: 10px;">Scan: ${{b.ploc}}</span>
      `;

      // rCRS status line
      const rcrsEl = document.getElementById('hud-rcrs');
      if (b.rcrs_coord) {{
        if (b.is_variant) {{
          rcrsEl.innerHTML = `<b style="color: #38bdf8;">⚡ rCRS ${{b.rcrs_coord}}: ${{b.rcrs_base}} → ${{b.called}} (ĐỘT BIẾN)</b>`;
        }} else {{
          rcrsEl.innerHTML = `<span style="color: var(--text-sub);">rCRS: ${{b.rcrs_coord}} (${{b.rcrs_base}} Match ✓)</span>`;
        }}
      }} else {{
        rcrsEl.innerHTML = `<span style="color: var(--text-sub);">Ngoài khung tham chiếu rCRS</span>`;
      }}

      // Quality & Secondary Ratio
      const qEl = document.getElementById('hud-quality');
      const heteroPct = (b.ratio * 100).toFixed(1);
      const heteroTag = b.ratio >= 0.15 ? `<b style="color: #f59e0b; margin-left: 6px;">[⚠️ Đỉnh kép: ${{heteroPct}}%]</b>` : '';
      qEl.innerHTML = `<span>Phred Q: <b style="color: ${{qColor}};">${{b.q}}</b> (${{qText}})</span>${{heteroTag}}`;

      // Channel RFUs
      const rfuEl = document.getElementById('hud-rfu');
      rfuEl.innerHTML = `
        <span style="color: #10b981;">A: ${{p.A||0}}%</span> &nbsp;|&nbsp; 
        <span style="color: #2563eb;">C: ${{p.C||0}}%</span> &nbsp;|&nbsp; 
        <span style="color: #f59e0b;">G: ${{p.G||0}}%</span> &nbsp;|&nbsp; 
        <span style="color: #ef4444;">T: ${{p.T||0}}%</span>
      `;

      floatingHud.style.display = 'block';

      // Position HUD avoiding edges
      const hudW = 270;
      const hudH = 110;
      let posX = clientX + 15;
      let posY = clientY + 15;

      if (posX + hudW > window.innerWidth) posX = clientX - hudW - 15;
      if (posY + hudH > window.innerHeight) posY = clientY - hudH - 15;

      floatingHud.style.left = posX + 'px';
      floatingHud.style.top = posY + 'px';
      floatingHud.style.position = 'fixed';
    }}

    // =========================================================================
    // ZOOM & AUTO-GAIN CONTROLS
    // =========================================================================
    function setZoom(newZoom, anchorX = null) {{
      const oldZoom = zoomScale;
      zoomScale = Math.max(0.6, Math.min(3.5, newZoom));
      document.getElementById('zoom-range').value = zoomScale;

      const currentScroll = scrollWrapper.scrollLeft;
      const clientW = scrollWrapper.clientWidth;
      const targetCanvasX = anchorX !== null ? anchorX : (currentScroll + clientW / 2);

      // Find anchor scan point
      const scanAnchor = getScanFromX(targetCanvasX);

      renderMain();

      // Recenter back to anchor scan point
      const newCanvasX = getX(scanAnchor);
      scrollWrapper.scrollLeft = newCanvasX - (anchorX !== null ? (anchorX - currentScroll) : (clientW / 2));
    }}

    document.getElementById('zoom-range').addEventListener('input', (e) => setZoom(parseFloat(e.target.value)));
    document.getElementById('btn-zoom-in').addEventListener('click', () => setZoom(zoomScale + 0.2));
    document.getElementById('btn-zoom-out').addEventListener('click', () => setZoom(zoomScale - 0.2));

    // Smooth Ctrl + MouseWheel Zoom
    scrollWrapper.addEventListener('wheel', (e) => {{
      if (e.ctrlKey) {{
        e.preventDefault();
        const delta = e.deltaY < 0 ? 0.15 : -0.15;
        const rect = canvas.getBoundingClientRect();
        const mouseX = e.clientX - rect.left;
        setZoom(zoomScale + delta, mouseX);
      }}
    }}, {{ passive: false }});

    // Vertical Amplitude Gain
    const gainRange = document.getElementById('gain-range');
    gainRange.addEventListener('input', (e) => {{
      gainScale = parseFloat(e.target.value);
      renderMain();
      renderMinimap();
    }});
    document.getElementById('btn-gain-reset').addEventListener('click', () => {{
      gainScale = 1.0;
      gainRange.value = 1.0;
      renderMain();
      renderMinimap();
    }});

    // Auto-Gain: Tự động cân bằng chiều cao sóng trong vùng xem hiện tại
    document.getElementById('btn-auto-gain').addEventListener('click', () => {{
      const scrollLeft = scrollWrapper.scrollLeft;
      const clientW = scrollWrapper.clientWidth;
      const startScan = Math.max(0, getScanFromX(scrollLeft));
      const endScan = Math.min(numScanPoints, getScanFromX(scrollLeft + clientW));

      let localMax = 10;
      ['A', 'C', 'G', 'T'].forEach(ch => {{
        if (channelsActive[ch] && traces[ch]) {{
          const s = Math.floor(startScan);
          const e = Math.floor(endScan);
          for (let i = s; i <= e; i += 2) {{
            if (traces[ch][i] > localMax) localMax = traces[ch][i];
          }}
        }}
      }});

      if (localMax > 10) {{
        const targetWaveHeight = 180; // pixels
        const idealGain = ((maxRfu * 1.1) / localMax) * (targetWaveHeight / 215);
        gainScale = Math.max(0.5, Math.min(6.0, idealGain));
        gainRange.value = gainScale.toFixed(2);
        renderMain();
        renderMinimap();
        showToast(`⚡ Auto-Gain: ${{gainScale.toFixed(2)}}x`);
      }}
    }});

    // =========================================================================
    // 1-CLICK COPY SEQUENCHER COORDINATES
    // =========================================================================
    document.getElementById('btn-copy-sequencher').addEventListener('click', () => {{
      const cut5 = trim.suggested_5p_cut_1based || 1;
      const cut3 = trim.suggested_3p_cut_1based || bases.length;
      const len = trim.retained_length || (cut3 - cut5 + 1);
      const textToCopy = `Sample: ${{DATA.sample_id}}_${{DATA.region}}_${{DATA.direction}} | 5' Trim: ${{cut5}} | 3' Trim: ${{cut3}} | Retained: ${{len}} bp`;
      
      navigator.clipboard.writeText(textToCopy).then(() => {{
        showToast(`✓ Đã sao chép: 5' Cut = ${{cut5}}, 3' Cut = ${{cut3}}`);
      }}).catch(() => {{
        prompt('Sao chép thông số Sequencher:', textToCopy);
      }});
    }});

    // Channel Toggle Handlers
    function toggleChannel(ch) {{
      channelsActive[ch] = !channelsActive[ch];
      const btn = document.getElementById(`ch-${{ch}}`);
      btn.classList.toggle('active', channelsActive[ch]);
      btn.style.opacity = channelsActive[ch] ? '1' : '0.4';
      renderMain();
      renderMinimap();
    }}

    ['A', 'C', 'G', 'T'].forEach(ch => {{
      document.getElementById(`ch-${{ch}}`).addEventListener('click', () => toggleChannel(ch));
    }});

    // Mini-Map Interaction
    let isDraggingLens = false;
    minimapWrapper.addEventListener('mousedown', (e) => {{
      isDraggingLens = true;
      moveLensToMouse(e);
    }});
    window.addEventListener('mousemove', (e) => {{
      if (isDraggingLens) moveLensToMouse(e);
    }});
    window.addEventListener('mouseup', () => {{ isDraggingLens = false; }});

    function moveLensToMouse(e) {{
      const rect = minimapWrapper.getBoundingClientRect();
      const mouseX = Math.max(0, Math.min(e.clientX - rect.left, rect.width));
      const ratio = mouseX / rect.width;
      scrollWrapper.scrollLeft = ratio * canvas.width - scrollWrapper.clientWidth / 2;
    }}

    scrollWrapper.addEventListener('scroll', updateMinimapLens);

    // =========================================================================
    // SMART CLINICAL JUMPS
    // =========================================================================
    function jumpToBaseIndex(idx) {{
      if (idx >= 0 && idx < bases.length) {{
        const targetX = getX(bases[idx].ploc);
        scrollWrapper.scrollTo({{ left: targetX - scrollWrapper.clientWidth / 2, behavior: 'smooth' }});
        hoveredBase = bases[idx];
        mouseCanvasX = targetX;
        renderMain();
      }}
    }}

    document.getElementById('btn-jump-5p').addEventListener('click', () => {{
      jumpToBaseIndex((trim.suggested_5p_cut_1based || 1) - 1);
    }});

    document.getElementById('btn-jump-3p').addEventListener('click', () => {{
      jumpToBaseIndex((trim.suggested_3p_cut_1based || bases.length) - 1);
    }});

    document.getElementById('btn-jump-polyc').addEventListener('click', () => {{
      const pIdx = bases.findIndex(b => b.status === 'POLYC' || b.status === 'POLYC_END');
      if (pIdx !== -1) {{
        jumpToBaseIndex(pIdx);
        showToast('🧬 Đã định vị vùng Poly-C tract');
      }} else {{
        showToast('Không phát hiện dải Poly-C trên mẫu này');
      }}
    }});

    // Jump to next SNP
    let currentSnpIdx = -1;
    function jumpNextSnp() {{
      const variantIndices = bases.map((b, i) => b.is_variant ? i : -1).filter(i => i !== -1);
      if (variantIndices.length === 0) {{
        showToast('Mẫu hoàn toàn trùng khớp với chuẩn rCRS (0 SNPs)');
        return;
      }}
      currentSnpIdx = (currentSnpIdx + 1) % variantIndices.length;
      const bIdx = variantIndices[currentSnpIdx];
      jumpToBaseIndex(bIdx);
      const b = bases[bIdx];
      showToast(`⚡ SNP [${{currentSnpIdx + 1}}/${{variantIndices.length}}]: rCRS ${{b.rcrs_coord}} (${{b.rcrs_base}} → ${{b.called}})`);
    }}
    document.getElementById('btn-jump-snp').addEventListener('click', jumpNextSnp);

    // Jump to next Heteroplasmy / Double Peak
    let currentHeteroIdx = -1;
    document.getElementById('btn-jump-hetero').addEventListener('click', () => {{
      const heteroIndices = bases.map((b, i) => (b.ratio >= 0.15 || b.is_conflict) ? i : -1).filter(i => i !== -1);
      if (heteroIndices.length === 0) {{
        showToast('Không phát hiện vị trí đỉnh phụ >= 15%');
        return;
      }}
      currentHeteroIdx = (currentHeteroIdx + 1) % heteroIndices.length;
      const bIdx = heteroIndices[currentHeteroIdx];
      jumpToBaseIndex(bIdx);
      const b = bases[bIdx];
      const ratioPct = (b.ratio * 100).toFixed(1);
      showToast(`⚠️ Đỉnh kép: Base #${{b.idx}} [${{b.called}}] - Tỷ lệ: ${{ratioPct}}%`);
    }});

    // Search Input (#150 or 16189)
    document.getElementById('input-search').addEventListener('keydown', (e) => {{
      if (e.key === 'Enter') {{
        const q = e.target.value.trim();
        if (q.startsWith('#')) {{
          const bNum = parseInt(q.slice(1), 10);
          jumpToBaseIndex(bNum - 1);
        }} else if (!isNaN(q)) {{
          const coord = parseInt(q, 10);
          const found = bases.findIndex(b => b.rcrs_coord === coord);
          if (found !== -1) jumpToBaseIndex(found);
          else {{
            const bNum = parseInt(q, 10);
            jumpToBaseIndex(bNum - 1);
          }}
        }}
      }}
    }});

    // Keyboard Shortcuts (Clinical efficiency)
    window.addEventListener('keydown', (e) => {{
      if (e.target.tagName === 'INPUT') return;

      if (e.key === 'a' || e.key === 'A') toggleChannel('A');
      else if (e.key === 'c' || e.key === 'C') toggleChannel('C');
      else if (e.key === 'g' || e.key === 'G') toggleChannel('G');
      else if (e.key === 't' || e.key === 'T') toggleChannel('T');
      else if (e.key === ' ') {{
        e.preventDefault();
        jumpNextSnp();
      }} else if (e.key === 'ArrowRight') {{
        scrollWrapper.scrollLeft += 40;
      }} else if (e.key === 'ArrowLeft') {{
        scrollWrapper.scrollLeft -= 40;
      }} else if (e.key === '+' || e.key === '=') {{
        setZoom(zoomScale + 0.2);
      }} else if (e.key === '-' || e.key === '_') {{
        setZoom(zoomScale - 0.2);
      }}
    }});

    // Dark/Light Theme Toggle
    const themeBtn = document.getElementById('btn-theme-toggle');
    themeBtn.addEventListener('click', () => {{
      const isDark = document.body.classList.contains('dark');
      document.body.classList.toggle('dark', !isDark);
      document.body.classList.toggle('light', isDark);
      themeBtn.textContent = isDark ? '🌙 Dark' : '☀️ Light';
      renderMain();
      renderMinimap();
    }});

    // Export Snapshot as PNG
    document.getElementById('btn-snapshot').addEventListener('click', () => {{
      const link = document.createElement('a');
      link.download = `${{DATA.sample_id}}_${{DATA.region}}_${{DATA.direction}}_chromatogram.png`;
      link.href = canvas.toDataURL('image/png');
      link.click();
      showToast('📸 Đã tải xuống ảnh chụp điện di đồ PNG');
    }});

    // Initial Launch
    window.addEventListener('load', () => {{
      renderMinimap();
      renderMain();
    }});
    window.addEventListener('resize', () => {{
      renderMinimap();
      renderMain();
    }});
  </script>
</body>
</html>
"""
        with open(html_file, "w", encoding="utf-8") as f:
            f.write(html_content)
        return html_file

    def export_interactive_html(self, sample_id: str, region: str, analyzed_bases: List[Dict[str, Any]], contig_res: Optional[Dict[str, Any]] = None) -> str:
        """Fallback tương thích ngược"""
        return ""
