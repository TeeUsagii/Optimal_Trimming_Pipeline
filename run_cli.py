import os
import sys
import glob
import re
import argparse
from typing import Dict, List, Tuple, Optional, Callable, Any

if sys.stdout and hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
if sys.stderr and hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8')
from core.ab1_parser import AB1Parser
from core.tracy_engine import TracyEngine
from core.aligner import RCRSAligner
from core.trimmer import SmartTrimmer
from core.assembler import ContigAssembler
from core.reporter import SuiteReporter

PRIMER_CONFIG = {
    'F15971': {'region': 'HV1', 'direction': 'F'},
    'R16410': {'region': 'HV1', 'direction': 'R'},
    'F15':    {'region': 'HV2_HV3', 'direction': 'F'},
    'R639':   {'region': 'HV2_HV3', 'direction': 'R'},
}

def detect_sample_info(file_path: str) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    """
    Tự động phân tích tên file .ab1 theo chuẩn phòng xét nghiệm:
    Định dạng: [Well]_[SampleID]-[Primer]_[Date]_[Time].ab1
    Ví dụ: H05_L2-R639_20260718_113234.ab1 -> SampleID='L2', Primer='R639', Region='HV2_HV3', Direction='R'

    Từ dữ liệu mồi, suy biến chính xác:
    - F15971 -> Vùng HV1, Chiều F
    - R16410 -> Vùng HV1, Chiều R
    - F15    -> Vùng HV2_HV3, Chiều F
    - R639   -> Vùng HV2_HV3, Chiều R
    """
    fname = os.path.splitext(os.path.basename(file_path))[0]

    # 1. Tìm mồi (ưu tiên F15971 trước F15)
    matched_primer = None
    for p in ['F15971', 'R16410', 'R639', 'F15']:
        pattern = rf'[-_]{p}(?:_|$)'
        if re.search(pattern, fname, re.IGNORECASE):
            matched_primer = p
            break

    # Nếu không thấy mồi chuẩn trong tên, thử fallback theo quy ước cũ
    if not matched_primer:
        # Fallback theo từ khóa HV1 / HV2
        reg = "HV1" if "HV1" in fname.upper() else ("HV2_HV3" if any(k in fname.upper() for k in ["HV2", "HV3"]) else None)
        direct = "R" if re.search(r'(_|-)(R|REV)($|_)', fname, re.IGNORECASE) else "F"
        parts = re.split(r'[_.-]', fname)
        sid = parts[0] if parts else fname
        return (sid, reg, direct, None) if reg else (None, None, None, None)

    cfg = PRIMER_CONFIG[matched_primer]
    parts = fname.split('_')
    sample_id = None

    if len(parts) >= 2:
        # Nếu parts[0] là tọa độ giếng trên plate (ví dụ: A04, H05: 1 chữ cái + 2 số)
        if re.match(r'^[A-H]\d{2}$', parts[0], re.IGNORECASE):
            sub = parts[1]
            sample_id = re.sub(rf'[-_]{matched_primer}$', '', sub, flags=re.IGNORECASE)
        else:
            sub = parts[0]
            sample_id = re.sub(rf'[-_]{matched_primer}$', '', sub, flags=re.IGNORECASE)
    else:
        sample_id = re.sub(rf'[-_]{matched_primer}.*', '', fname, flags=re.IGNORECASE)

    return sample_id, cfg['region'], cfg['direction'], matched_primer

def run_pipeline(input_dir: str,
                 rcrs_path: str,
                 output_dir: str,
                 peak_thresh: float = 0.25,
                 conflict_thresh: float = 0.70,
                 progress_callback: Optional[Callable[[int, int, str], None]] = None) -> Dict[str, Any]:
    print("=" * 70)
    print("      OPTIMAL TRIMMING PIPELINE (HIGH-THROUGHPUT ENGINE)")
    print("=" * 70)
    print(f"[+] Thư mục dữ liệu thô: {input_dir}")
    print(f"[+] Chuỗi tham chiếu rCRS: {rcrs_path}")
    print(f"[+] Thư mục xuất kết quả: {output_dir}")
    print(f"[+] Ngưỡng Tracy Peak Threshold: {peak_thresh}")
    print(f"[+] Ngưỡng Chồng Peak (?): {conflict_thresh}")

    # Khởi tạo các Core Engines
    aligner = RCRSAligner(rcrs_path)
    tracy = TracyEngine(peak_threshold=peak_thresh, conflict_threshold=conflict_thresh)
    trimmer = SmartTrimmer()
    assembler = ContigAssembler(aligner.rcrs_seq)
    reporter = SuiteReporter(output_dir)

    # 1. Thu thập và gom nhóm các file .ab1
    ab1_files = glob.glob(os.path.join(input_dir, "**", "*.ab1"), recursive=True)
    if not ab1_files:
        print(f"[-] Không tìm thấy file .ab1 nào trong: {input_dir}")
        return

    print(f"[+] Tìm thấy tổng cộng {len(ab1_files)} file .ab1. Đang gom nhóm theo Mẫu & Vùng...")

    samples_map: Dict[str, Dict[str, Dict[str, str]]] = {}
    for f in ab1_files:
        sid, reg, direct, primer = detect_sample_info(f)
        if not sid or not reg or not direct:
            continue
        if sid not in samples_map:
            samples_map[sid] = {}
        if reg not in samples_map[sid]:
            samples_map[sid][reg] = {}
        samples_map[sid][reg][direct] = f

    print(f"[+] Tổng số mẫu phân tích: {len(samples_map)} mẫu.")
    print("-" * 70)

    summary_records = []
    trim_recommendations = []
    processed_count = 0
    total_tasks = sum(len(regions) for regions in samples_map.values())

    # 2. Xử lý từng mẫu
    for sid, regions in samples_map.items():
        for reg, files in regions.items():
            processed_count += 1
            if progress_callback and total_tasks > 0:
                progress_callback(processed_count, total_tasks, f"Đang phân tích {sid} ({reg})...")
            fwd_file = files.get("F")
            rev_file = files.get("R")

            fwd_data = None
            fwd_aln = None
            fwd_trim = None

            rev_data = None
            rev_aln = None
            rev_trim = None

            # 1. Phân tích Forward (nếu có)
            if fwd_file:
                parsed_f = AB1Parser.parse(fwd_file)
                if parsed_f:
                    fwd_data = tracy.analyze(parsed_f)
                    fwd_aln = aligner.align(fwd_data["primary_sequence"], reg, is_reverse=False)
                    if fwd_aln:
                        fwd_trim = trimmer.evaluate_boundaries(fwd_data, fwd_aln, reg, is_reverse=False)
                        reporter.export_trimmed_fasta(sid, reg, "F", fwd_trim)

            # 2. Phân tích Reverse (nếu có)
            if rev_file:
                parsed_r = AB1Parser.parse(rev_file)
                if parsed_r:
                    rev_data = tracy.analyze(parsed_r)
                    rev_aln = aligner.align(rev_data["primary_sequence"], reg, is_reverse=True)
                    if rev_aln:
                        rev_trim = trimmer.evaluate_boundaries(rev_data, rev_aln, reg, is_reverse=True)
                        reporter.export_trimmed_fasta(sid, reg, "R", rev_trim)

            # 3. Ráp Contig 2 chiều & Đánh giá Tiêu chuẩn QC (PASS / REVIEW / DENY)
            contig_res = assembler.assemble_region(
                region=reg,
                fwd_data=fwd_data,
                fwd_trim=fwd_trim,
                fwd_aln=fwd_aln,
                rev_data=rev_data,
                rev_trim=rev_trim,
                rev_aln=rev_aln
            )

            # Xuất FASTA Consensus
            reporter.export_consensus_fasta(sid, reg, contig_res)

            # 4. Xuất Báo cáo HTML Chromatogram (kèm đối chiếu rCRS và cảnh báo DENY nếu có)
            if fwd_file and parsed_f and fwd_data and fwd_trim:
                reporter.export_chromatogram_html(
                    sample_id=sid,
                    region=reg,
                    direction="F",
                    parsed_ab1=parsed_f,
                    analyzed_data=fwd_data,
                    trim_res=fwd_trim,
                    alignment_info=fwd_aln,
                    rcrs_seq=aligner.rcrs_seq,
                    contig_res=contig_res
                )
            if rev_file and parsed_r and rev_data and rev_trim:
                reporter.export_chromatogram_html(
                    sample_id=sid,
                    region=reg,
                    direction="R",
                    parsed_ab1=parsed_r,
                    analyzed_data=rev_data,
                    trim_res=rev_trim,
                    alignment_info=rev_aln,
                    rcrs_seq=aligner.rcrs_seq,
                    contig_res=contig_res
                )

            # 5. Ghi nhận Bảng khuyến nghị vị trí cắt chi tiết (Dạng A cho Sequencher v5.4.6)
            is_deny = (contig_res["qc_status"] == "DENY")
            if fwd_file:
                _, _, _, primer_f = detect_sample_info(fwd_file)
                fname_f = parsed_f["file_name"] if parsed_f else os.path.basename(fwd_file)
                raw_len_f = len(parsed_f["raw_sequence"]) if parsed_f else 0
                if fwd_trim:
                    seq_guidance_f = (
                        f"⚠️ DENY: Khuyến cáo phòng lab làm lại mẫu ({contig_res.get('qc_details', '')}). Không import vào Sequencher."
                        if is_deny else
                        f"Trim 5' bases 1-{fwd_trim['suggested_5p_cut_1based']-1}, Trim 3' from {fwd_trim['suggested_3p_cut_1based']+1}"
                    )
                    trim_recommendations.append({
                        "File_Name": fname_f,
                        "Sample_ID": sid,
                        "Region": reg,
                        "Primer": primer_f or "F",
                        "Direction": "Forward",
                        "Raw_Length_bp": raw_len_f,
                        "Suggested_5p_Cut": fwd_trim["suggested_5p_cut_1based"],
                        "Suggested_3p_Cut": fwd_trim["suggested_3p_cut_1based"],
                        "Retained_Range": f"{fwd_trim['suggested_5p_cut_1based']} - {fwd_trim['suggested_3p_cut_1based']}",
                        "Retained_Length_bp": fwd_trim["retained_length"],
                        "rCRS_Span_Covered": f"{fwd_trim['rcrs_span'][0]} - {fwd_trim['rcrs_span'][1]}",
                        "PolyC_Status": fwd_trim["poly_c_status"],
                        "Noise_5p_Reason": fwd_trim["noise_5p_reason"],
                        "Noise_3p_Reason": fwd_trim["noise_3p_reason"],
                        "QC_Status": contig_res["qc_status"],
                        "Lab_Recommendation": contig_res.get("qc_recommendation", ""),
                        "Sequencher_Guidance": seq_guidance_f
                    })
                else:
                    trim_recommendations.append({
                        "File_Name": fname_f,
                        "Sample_ID": sid,
                        "Region": reg,
                        "Primer": primer_f or "F",
                        "Direction": "Forward",
                        "Raw_Length_bp": raw_len_f,
                        "Suggested_5p_Cut": "N/A",
                        "Suggested_3p_Cut": "N/A",
                        "Retained_Range": "N/A",
                        "Retained_Length_bp": 0,
                        "rCRS_Span_Covered": "N/A",
                        "PolyC_Status": "NONE",
                        "Noise_5p_Reason": "Không căn gióng được với rCRS",
                        "Noise_3p_Reason": "Không căn gióng được với rCRS",
                        "QC_Status": "DENY",
                        "Lab_Recommendation": "DENY (Khuyến cáo làm lại mẫu trong phòng thí nghiệm)",
                        "Sequencher_Guidance": "⚠️ DENY: Mẫu quá nhiễu hoặc sai primer, không căn gióng được rCRS. Khuyến cáo làm lại mẫu!"
                    })

            if rev_file:
                _, _, _, primer_r = detect_sample_info(rev_file)
                fname_r = parsed_r["file_name"] if parsed_r else os.path.basename(rev_file)
                raw_len_r = len(parsed_r["raw_sequence"]) if parsed_r else 0
                if rev_trim:
                    seq_guidance_r = (
                        f"⚠️ DENY: Khuyến cáo phòng lab làm lại mẫu ({contig_res.get('qc_details', '')}). Không import vào Sequencher."
                        if is_deny else
                        f"Trim 5' bases 1-{rev_trim['suggested_5p_cut_1based']-1}, Trim 3' from {rev_trim['suggested_3p_cut_1based']+1}"
                    )
                    trim_recommendations.append({
                        "File_Name": fname_r,
                        "Sample_ID": sid,
                        "Region": reg,
                        "Primer": primer_r or "R",
                        "Direction": "Reverse",
                        "Raw_Length_bp": raw_len_r,
                        "Suggested_5p_Cut": rev_trim["suggested_5p_cut_1based"],
                        "Suggested_3p_Cut": rev_trim["suggested_3p_cut_1based"],
                        "Retained_Range": f"{rev_trim['suggested_5p_cut_1based']} - {rev_trim['suggested_3p_cut_1based']}",
                        "Retained_Length_bp": rev_trim["retained_length"],
                        "rCRS_Span_Covered": f"{rev_trim['rcrs_span'][0]} - {rev_trim['rcrs_span'][1]}",
                        "PolyC_Status": rev_trim["poly_c_status"],
                        "Noise_5p_Reason": rev_trim["noise_5p_reason"],
                        "Noise_3p_Reason": rev_trim["noise_3p_reason"],
                        "QC_Status": contig_res["qc_status"],
                        "Lab_Recommendation": contig_res.get("qc_recommendation", ""),
                        "Sequencher_Guidance": seq_guidance_r
                    })
                else:
                    trim_recommendations.append({
                        "File_Name": fname_r,
                        "Sample_ID": sid,
                        "Region": reg,
                        "Primer": primer_r or "R",
                        "Direction": "Reverse",
                        "Raw_Length_bp": raw_len_r,
                        "Suggested_5p_Cut": "N/A",
                        "Suggested_3p_Cut": "N/A",
                        "Retained_Range": "N/A",
                        "Retained_Length_bp": 0,
                        "rCRS_Span_Covered": "N/A",
                        "PolyC_Status": "NONE",
                        "Noise_5p_Reason": "Không căn gióng được với rCRS",
                        "Noise_3p_Reason": "Không căn gióng được với rCRS",
                        "QC_Status": "DENY",
                        "Lab_Recommendation": "DENY (Khuyến cáo làm lại mẫu trong phòng thí nghiệm)",
                        "Sequencher_Guidance": "⚠️ DENY: Mẫu quá nhiễu hoặc sai primer, không căn gióng được rCRS. Khuyến cáo làm lại mẫu!"
                    })

            # Ghi nhận kết quả tóm tắt
            summary_records.append({
                "Sample_ID": sid,
                "Region": reg,
                "Target_Range": contig_res["target_range"],
                "Coverage_%": contig_res["coverage_pct"],
                "QC_Status": contig_res["qc_status"],
                "Lab_Recommendation": contig_res.get("qc_recommendation", ""),
                "QC_Details": contig_res.get("qc_details", ""),
                "Variants_Count": contig_res["variant_count"],
                "Variants_List": "; ".join(contig_res["variants"]),
                "Ambiguity_Count_+": contig_res["ambiguity_count"],
                "Conflict_Count_?": contig_res["conflict_question_count"],
                "Suggested_Trim_Forward": f"Cut 5': {fwd_trim['suggested_5p_cut_1based']}, Cut 3': {fwd_trim['suggested_3p_cut_1based']} ({fwd_trim['retained_length']} bp, {fwd_trim['poly_c_status']})" if fwd_trim else "Missing",
                "Suggested_Trim_Reverse": f"Cut 5': {rev_trim['suggested_5p_cut_1based']}, Cut 3': {rev_trim['suggested_3p_cut_1based']} ({rev_trim['retained_length']} bp, {fwd_trim['poly_c_status'] if fwd_trim else 'NONE'})" if rev_trim else "Missing"
            })

            qc_str = contig_res['qc_status']
            if qc_str == "DENY":
                qc_display = "DENY ⚠️ (CẦN LÀM LẠI MẪU)"
            elif qc_str == "PASS":
                qc_display = "PASS ✓"
            else:
                qc_display = "REVIEW 👁️"

            print(f"[{processed_count:03d}] Mẫu: {sid:<8} | Vùng: {reg:<7} | QC: {qc_display:<26} | SNPs: {contig_res['variant_count']:<3} | (+): {contig_res['ambiguity_count']}")

    # 4. Xuất Bảng khuyến nghị vị trí cắt chi tiết (Dạng A cho Sequencher)
    trim_csv_out = reporter.export_trim_recommendations_csv(trim_recommendations)

    # 5. Xuất file tổng hợp CSV
    csv_out = reporter.export_summary_csv(summary_records)
    print("-" * 70)
    print(f"[✓] HOÀN TẤT PHÂN TÍCH!")
    print(f"[✓] Bảng khuyến nghị vị trí cắt (Dạng A cho Sequencher): {trim_csv_out}")
    print(f"[✓] Bảng tổng kết kiểu gen và SNPs: {csv_out}")
    print(f"[✓] Điện di đồ Chromatogram tương tác đối chiếu rCRS (Dạng C): {reporter.html_dir}")
    print(f"[✓] File FASTA chuỗi bảo toàn kèm tọa độ đề xuất: {reporter.fasta_dir}")
    print("=" * 70)

    return {
        "summary_records": summary_records,
        "trim_recommendations": trim_recommendations,
        "trim_csv_path": trim_csv_out,
        "summary_csv_path": csv_out,
        "html_dir": reporter.html_dir,
        "fasta_dir": reporter.fasta_dir,
        "total_samples": len(samples_map),
        "processed_count": processed_count
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Optimal Trimming Pipeline for Human mtDNA Sanger Sequencing")
    parser.add_argument("--input", "-i", default="./data/sample_ab1", help="Thư mục chứa file .ab1")
    parser.add_argument("--rcrs", "-r", default="./data/rCRS.fasta", help="Đường dẫn file rCRS.fasta")
    parser.add_argument("--output", "-o", default="./output", help="Thư mục lưu kết quả")
    parser.add_argument("--threshold", "-t", type=float, default=0.25, help="Ngưỡng tỷ lệ đỉnh phụ Tracy")
    parser.add_argument("--conflict", "-c", type=float, default=0.70, help="Ngưỡng nhận diện chồng peak (?)")

    args = parser.parse_args()
    run_pipeline(args.input, args.rcrs, args.output, args.threshold, args.conflict)
