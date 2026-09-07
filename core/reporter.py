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
    Module xuất báo cáo đa luồng cho mtDNA Sanger Analysis Suite:
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
  <script src="https://www.gstatic.com/antigravity/web/dev/tailwindcss.min.js"></script>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; }}
    .trace-canvas {{ background-color: #0f172a; cursor: crosshair; }}
    .base-col:hover {{ background-color: rgba(56, 189, 248, 0.15); }}
    /* Custom scrollbar */
    ::-webkit-scrollbar {{ height: 10px; width: 8px; }}
    ::-webkit-scrollbar-track {{ background: #1e293b; }}
    ::-webkit-scrollbar-thumb {{ background: #475569; border-radius: 5px; }}
    ::-webkit-scrollbar-thumb:hover {{ background: #64748b; }}
  </style>
</head>
<body class="bg-slate-950 text-slate-100 p-3 min-h-screen">
  <div class="max-w-7xl mx-auto flex flex-col gap-3">

    {deny_alert_html}

    <!-- Top Header Banner -->
    <div class="bg-slate-900 border border-slate-800 rounded-xl p-4 shadow-xl flex flex-wrap items-center justify-between gap-4">
      <div>
        <div class="flex items-center gap-3">
          <span class="inline-block w-3.5 h-3.5 rounded-full bg-cyan-400 shadow-[0_0_10px_#22d3ee]"></span>
          <h1 class="text-xl font-bold text-white tracking-wide">
            {sample_id} <span class="text-cyan-400 font-mono">[{region} - Chiều {direction}]</span>
          </h1>
          <span class="text-xs px-2.5 py-0.5 rounded-full bg-slate-800 text-slate-400 border border-slate-700">SeqStudio Flex 24</span>
          {qc_badge_html}
        </div>
        <p class="text-xs text-slate-400 mt-1 font-mono">File: {parsed_ab1.get("file_name", "")}</p>
      </div>

      <!-- Recommendation Summary Badges -->
      <div class="flex flex-wrap items-center gap-2 text-xs">
        <div class="px-3 py-1.5 rounded-lg bg-slate-800 border border-slate-700">
          <span class="text-slate-400">Ranh giới giữ lại:</span>
          <span class="font-bold text-emerald-400 ml-1">Base {trim_res.get('suggested_5p_cut_1based', 1)} → {trim_res.get('suggested_3p_cut_1based', len(viewer_bases))}</span>
          <span class="text-slate-500 text-[11px]">({trim_res.get('retained_length', len(viewer_bases))} bp)</span>
        </div>
        <div class="px-3 py-1.5 rounded-lg bg-slate-800 border border-slate-700">
          <span class="text-slate-400">rCRS Span:</span>
          <span class="font-bold text-sky-400 ml-1">{trim_res.get('rcrs_span', (0, 0))[0]} - {trim_res.get('rcrs_span', (0, 0))[1]}</span>
        </div>
        <div class="px-3 py-1.5 rounded-lg bg-purple-950/60 border border-purple-800/60 text-purple-300">
          <span>Poly-C:</span>
          <span class="font-bold ml-1">{trim_res.get('poly_c_status', 'NONE')}</span>
        </div>
      </div>
    </div>

    <!-- Trimming Rationale Card -->
    <div class="bg-slate-900/90 border border-slate-800/80 rounded-xl p-3.5 text-xs grid grid-cols-1 md:grid-cols-2 gap-3">
      <div class="flex items-start gap-2 bg-slate-950/60 p-2.5 rounded-lg border border-slate-800">
        <span class="text-rose-400 font-bold">✂ Cắt 5' đề xuất:</span>
        <span class="text-slate-300">{trim_res.get('noise_5p_reason', 'Loại bỏ mồi ngoài và nhiễu 10-20 bp đầu')}</span>
      </div>
      <div class="flex items-start gap-2 bg-slate-950/60 p-2.5 rounded-lg border border-slate-800">
        <span class="text-amber-400 font-bold">✂ Cắt 3' đề xuất:</span>
        <span class="text-slate-300">{trim_res.get('noise_3p_reason', 'Giới hạn theo khung rCRS và sụt giảm tín hiệu cuối read')}</span>
      </div>
    </div>

    <!-- Main Chromatogram & Track Viewer -->
    <div class="bg-slate-900 border border-slate-800 rounded-xl p-4 shadow-xl flex flex-col gap-3">
      
      <!-- Toolbar & Legend -->
      <div class="flex flex-wrap items-center justify-between gap-3 text-xs border-b border-slate-800 pb-3">
        <!-- 4-Channel Legend -->
        <div class="flex items-center gap-4 font-bold">
          <span class="flex items-center gap-1.5"><span class="w-3 h-3 rounded-sm bg-emerald-500"></span> A (Adenine)</span>
          <span class="flex items-center gap-1.5"><span class="w-3 h-3 rounded-sm bg-blue-500"></span> C (Cytosine)</span>
          <span class="flex items-center gap-1.5"><span class="w-3 h-3 rounded-sm bg-amber-400"></span> G (Guanine)</span>
          <span class="flex items-center gap-1.5"><span class="w-3 h-3 rounded-sm bg-red-500"></span> T (Thymine)</span>
        </div>

        <!-- Jump to markers & Controls -->
        <div class="flex items-center gap-2">
          <button id="btn-jump-5p" class="px-2.5 py-1 rounded bg-slate-800 hover:bg-slate-700 text-slate-200 border border-slate-700">Nhảy tới 5' Cut</button>
          <button id="btn-jump-polyc" class="px-2.5 py-1 rounded bg-purple-900/60 hover:bg-purple-800/80 text-purple-200 border border-purple-700">Nhảy tới Poly-C</button>
          <button id="btn-jump-3p" class="px-2.5 py-1 rounded bg-slate-800 hover:bg-slate-700 text-slate-200 border border-slate-700">Nhảy tới 3' Cut</button>
          <div class="h-4 w-px bg-slate-700 mx-1"></div>
          <button id="btn-zoom-in" class="px-2 py-1 rounded bg-slate-800 hover:bg-slate-700 text-slate-200 font-bold border border-slate-700" title="Phóng to">+</button>
          <button id="btn-zoom-out" class="px-2 py-1 rounded bg-slate-800 hover:bg-slate-700 text-slate-200 font-bold border border-slate-700" title="Thu nhỏ">-</button>
        </div>
      </div>

      <!-- Horizontal Scrollable Container -->
      <div id="scroll-wrapper" class="overflow-x-auto overflow-y-hidden border border-slate-800 rounded-lg relative bg-slate-950">
        
        <!-- Basecalls Track (Above waves) -->
        <div id="basecall-track" class="h-10 relative select-none border-b border-slate-800/80 bg-slate-900/40"></div>

        <!-- Chromatogram Wave Canvas -->
        <canvas id="chroma-canvas" class="trace-canvas block" height="240"></canvas>

        <!-- rCRS Alignment Track (Directly below waves) -->
        <div id="rcrs-track" class="h-12 relative select-none border-t border-slate-800/80 bg-slate-900/60"></div>
      </div>

      <!-- Interactive Hover Tooltip Box -->
      <div id="hover-box" class="bg-slate-950 p-3 rounded-lg border border-slate-800 text-xs flex flex-wrap items-center justify-between gap-3 min-h-[42px]">
        <div class="text-slate-400">Rà chuột trên đồ thị điện di để soi chi tiết từng nucleotide, điểm Phred và so sánh chuỗi chuẩn rCRS...</div>
      </div>

    </div>

  </div>

  <script>
    const DATA = {json_data};

    // Color definitions
    const COLORS = {{
      'A': '#10b981', // Emerald
      'C': '#2563eb', // Blue
      'G': '#f59e0b', // Amber/Gold (thay màu đen để tương phản sắc nét trên nền tối)
      'T': '#ef4444'  // Red
    }};

    let zoomScale = 1.0;
    const canvas = document.getElementById('chroma-canvas');
    const ctx = canvas.getContext('2d');
    const baseTrack = document.getElementById('basecall-track');
    const rcrsTrack = document.getElementById('rcrs-track');
    const scrollWrapper = document.getElementById('scroll-wrapper');
    const hoverBox = document.getElementById('hover-box');

    const traces = DATA.traces;
    const bases = DATA.bases;
    const trim = DATA.trim;

    const numScanPoints = (traces['A'] && traces['A'].length) || (bases.length * 12);
    
    // Tìm RFU cực đại để chuẩn hóa chiều cao sóng
    let maxRfu = 100;
    ['A', 'C', 'G', 'T'].forEach(c => {{
      if (traces[c]) {{
        for (let i = 0; i < traces[c].length; i += 5) {{
          if (traces[c][i] > maxRfu) maxRfu = traces[c][i];
        }}
      }}
    }});

    function renderAll() {{
      const baseSpacing = 16 * zoomScale;
      const totalWidth = Math.max(bases.length * baseSpacing + 100, scrollWrapper.clientWidth);

      canvas.width = totalWidth;
      canvas.height = 240;
      baseTrack.style.width = totalWidth + 'px';
      rcrsTrack.style.width = totalWidth + 'px';

      // 1. Vẽ đồ thị sóng 4 màu
      ctx.clearRect(0, 0, canvas.width, canvas.height);

      // Tính tọa độ X cho từng điểm quét
      const plocs = bases.map(b => b.ploc);
      const minPloc = plocs[0] || 0;
      const maxPloc = plocs[plocs.length - 1] || numScanPoints;

      function getX(scanIdx) {{
        return 40 + ((scanIdx - minPloc) / (maxPloc - minPloc)) * (bases.length * baseSpacing);
      }}

      // Vẽ vùng bóng mờ Trim 5' và 3'
      const cut5pIdx = trim.suggested_5p_cut_1based - 1;
      const cut3pIdx = trim.suggested_3p_cut_1based - 1;

      const x5p = cut5pIdx >= 0 && cut5pIdx < bases.length ? getX(bases[cut5pIdx].ploc) : 0;
      const x3p = cut3pIdx >= 0 && cut3pIdx < bases.length ? getX(bases[cut3pIdx].ploc) : totalWidth;

      // Vùng 5' mờ
      ctx.fillStyle = 'rgba(15, 23, 42, 0.75)';
      ctx.fillRect(0, 0, x5p, canvas.height);

      // Vùng 3' mờ
      ctx.fillRect(x3p, 0, totalWidth - x3p, canvas.height);

      // Vạch ranh giới nét đứt
      ctx.setLineDash([5, 4]);
      ctx.lineWidth = 1.5;
      
      // Vạch 5'
      ctx.strokeStyle = '#f43f5e';
      ctx.beginPath();
      ctx.moveTo(x5p, 0); ctx.lineTo(x5p, canvas.height);
      ctx.stroke();

      // Vạch 3'
      ctx.strokeStyle = '#f59e0b';
      ctx.beginPath();
      ctx.moveTo(x3p, 0); ctx.lineTo(x3p, canvas.height);
      ctx.stroke();

      ctx.setLineDash([]); // Reset nét liền

      // Vẽ 4 kênh sóng
      const channels = ['A', 'C', 'G', 'T'];
      channels.forEach(ch => {{
        const trace = traces[ch];
        if (!trace) return;

        ctx.strokeStyle = COLORS[ch];
        ctx.lineWidth = 1.5;
        ctx.beginPath();

        let started = false;
        // Quét bước nhảy theo zoom để tối ưu hiệu năng render
        const step = zoomScale > 1.2 ? 1 : 2;

        for (let i = 0; i < trace.length; i += step) {{
          const x = getX(i);
          if (x < -50 || x > totalWidth + 50) continue;

          // Chuẩn hóa chiều cao (scale 85% chiều cao canvas)
          const val = trace[i];
          const y = canvas.height - 10 - (val / (maxRfu * 1.1)) * (canvas.height - 20);

          if (!started) {{
            ctx.moveTo(x, y);
            started = true;
          }} else {{
            ctx.lineTo(x, y);
          }}
        }}
        ctx.stroke();
      }});

      // 2. Render Basecalls Track & rCRS Track
      baseTrack.innerHTML = '';
      rcrsTrack.innerHTML = '';

      bases.forEach((b, idx) => {{
        const x = getX(b.ploc);

        // A. Basecall Element
        const bEl = document.createElement('div');
        bEl.className = 'absolute flex flex-col items-center cursor-pointer transition-transform hover:scale-125';
        bEl.style.left = (x - 9) + 'px';
        bEl.style.top = '2px';
        bEl.style.width = '18px';

        let color = COLORS[b.called] || '#e2e8f0';
        let bgStyle = '';
        if (b.is_conflict) {{
          color = '#ffffff';
          bgStyle = 'background: #9333ea; border-radius: 4px; padding: 0 3px;';
        }}

        bEl.innerHTML = `
          <span class="text-[10px] text-slate-500 font-mono">${{b.idx}}</span>
          <span class="font-bold text-sm leading-none font-mono" style="color: ${{color}}; ${{bgStyle}}">${{b.called}}</span>
        `;

        bEl.addEventListener('mouseenter', () => showDetail(b));
        baseTrack.appendChild(bEl);

        // B. rCRS Alignment Track Element
        const rEl = document.createElement('div');
        rEl.className = 'absolute flex flex-col items-center cursor-pointer text-xs font-mono';
        rEl.style.left = (x - 12) + 'px';
        rEl.style.top = '4px';
        rEl.style.width = '24px';

        if (b.rcrs_coord) {{
          let matchHtml = '';
          if (b.is_match) {{
            matchHtml = '<span class="text-slate-400 font-bold">.</span>';
          }} else if (b.is_variant) {{
            matchHtml = `<span class="px-1 rounded bg-amber-500/20 text-amber-400 font-bold border border-amber-500/40">${{b.rcrs_base}}</span>`;
          }} else {{
            matchHtml = `<span class="text-slate-400 font-bold">${{b.rcrs_base || '?'}}</span>`;
          }}

          rEl.innerHTML = `
            ${{matchHtml}}
            <span class="text-[9px] text-slate-500 leading-tight mt-0.5">${{b.rcrs_coord}}</span>
          `;
        }} else {{
          rEl.innerHTML = `<span class="text-slate-600">-</span>`;
        }}

        rEl.addEventListener('mouseenter', () => showDetail(b));
        rcrsTrack.appendChild(rEl);
      }});
    }}

    function showDetail(b) {{
      const p = b.probs || {{}};
      const probStr = `A:${{p.A||0}}% | C:${{p.C||0}}% | G:${{p.G||0}}% | T:${{p.T||0}}%`;
      const rcrsText = b.rcrs_coord 
        ? `<span class="text-sky-400 font-bold">rCRS: ${{b.rcrs_coord}} (${{b.rcrs_base}} → ${{b.called}})</span>` 
        : `<span class="text-slate-500">Ngoài dải rCRS</span>`;
      
      const polyNote = b.is_post_polyc 
        ? `<span class="px-1.5 py-0.5 rounded bg-purple-900/60 text-purple-300 border border-purple-700">Dữ liệu đối sánh sau Poly-C</span>` 
        : '';

      hoverBox.innerHTML = `
        <div class="flex items-center gap-3">
          <span class="font-bold text-sm text-white font-mono">Base #${{b.idx}} [${{b.called}}]</span>
          <span class="text-slate-400">Scan: <b class="text-slate-200">${{b.ploc}}</b></span>
          <span class="text-slate-400">Phred Q: <b class="${{b.q >= 30 ? 'text-emerald-400' : (b.q >= 20 ? 'text-cyan-400' : 'text-rose-400')}}">${{b.q}}</b></span>
          <span class="text-slate-400">Đỉnh phụ: <b class="text-slate-200">${{(b.ratio * 100).toFixed(1)}}%</b></span>
        </div>
        <div class="flex items-center gap-3">
          ${{rcrsText}}
          <span class="text-slate-400 font-mono text-[11px]">${{probStr}}</span>
          ${{polyNote}}
        </div>
      `;
    }}

    // Event handlers
    document.getElementById('btn-zoom-in').addEventListener('click', () => {{
      zoomScale = Math.min(zoomScale * 1.25, 3.5);
      renderAll();
    }});

    document.getElementById('btn-zoom-out').addEventListener('click', () => {{
      zoomScale = Math.max(zoomScale / 1.25, 0.6);
      renderAll();
    }});

    document.getElementById('btn-jump-5p').addEventListener('click', () => {{
      const cut5pIdx = trim.suggested_5p_cut_1based - 1;
      if (cut5pIdx >= 0 && cut5pIdx < bases.length) {{
        const baseSpacing = 16 * zoomScale;
        scrollWrapper.scrollTo({{ left: cut5pIdx * baseSpacing - 100, behavior: 'smooth' }});
      }}
    }});

    document.getElementById('btn-jump-3p').addEventListener('click', () => {{
      const cut3pIdx = trim.suggested_3p_cut_1based - 1;
      if (cut3pIdx >= 0 && cut3pIdx < bases.length) {{
        const baseSpacing = 16 * zoomScale;
        scrollWrapper.scrollTo({{ left: cut3pIdx * baseSpacing - 200, behavior: 'smooth' }});
      }}
    }});

    document.getElementById('btn-jump-polyc').addEventListener('click', () => {{
      // Tìm vị trí Poly-C đầu tiên
      const pIdx = bases.findIndex(b => b.status === 'POLYC' || b.status === 'POLYC_END');
      if (pIdx !== -1) {{
        const baseSpacing = 16 * zoomScale;
        scrollWrapper.scrollTo({{ left: pIdx * baseSpacing - 150, behavior: 'smooth' }});
      }} else {{
        alert('Không phát hiện dải Poly-C trên mẫu này.');
      }}
    }});

    // Khởi chạy
    window.addEventListener('load', renderAll);
    window.addEventListener('resize', renderAll);
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
