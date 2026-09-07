from typing import Dict, List, Any, Optional, Tuple
from Bio.Seq import Seq

IUPAC_MAP = {
    tuple(sorted(['A', 'G'])): 'R',
    tuple(sorted(['C', 'T'])): 'Y',
    tuple(sorted(['G', 'C'])): 'S',
    tuple(sorted(['A', 'T'])): 'W',
    tuple(sorted(['G', 'T'])): 'K',
    tuple(sorted(['A', 'C'])): 'M',
}

class ContigAssembler:
    """
    Ráp contig 2 chiều (Forward & Reverse) và đối chiếu với rCRS:
    - Giải quyết xung đột giữa 2 chiều dựa trên điểm chất lượng Phred & tỷ lệ đỉnh phụ.
    - Tự động triệt tiêu các mâu thuẫn rác (do trượt Poly-C, dye blob).
    - So sánh với rCRS: Trùng khớp ký hiệu '.', biến dị ghi nhận SNP, xung đột không giải quyết được ghi dấu '+'.
    """

    def __init__(self, rcrs_seq: str):
        self.rcrs_seq = rcrs_seq.upper()

    def assemble_region(self,
                        region: str,
                        fwd_data: Optional[Dict[str, Any]],
                        fwd_trim: Optional[Dict[str, Any]],
                        fwd_aln: Optional[Dict[str, Any]],
                        rev_data: Optional[Dict[str, Any]],
                        rev_trim: Optional[Dict[str, Any]],
                        rev_aln: Optional[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Ráp contig cho 1 vùng (HV1 hoặc HV2_HV3).
        """
        is_hv1 = "HV1" in region.upper()
        target_start = 16024 if is_hv1 else 73
        target_end = 16365 if is_hv1 else 576
        target_len = target_end - target_start + 1

        # Lập chỉ mục rCRS -> base_info cho Forward
        fwd_map = {}
        if fwd_trim and fwd_aln:
            f_bases = fwd_data["bases"]
            for q_idx in range(fwd_trim["trim_start_idx"], fwd_trim["trim_end_idx"]):
                if q_idx in fwd_aln["query_to_rcrs"]:
                    r_pos = fwd_aln["query_to_rcrs"][q_idx]
                    fwd_map[r_pos] = f_bases[q_idx]

        # Lập chỉ mục rCRS -> base_info cho Reverse
        rev_map = {}
        if rev_trim and rev_aln:
            r_bases = rev_data["bases"]
            for q_idx in range(rev_trim["trim_start_idx"], rev_trim["trim_end_idx"]):
                if q_idx in rev_aln["query_to_rcrs"]:
                    r_pos = rev_aln["query_to_rcrs"][q_idx]
                    # Chú ý: Reverse read trên rCRS cần lấy base bổ sung (Reverse Complement)
                    b_obj = dict(r_bases[q_idx])
                    # Đảo bổ sung base
                    orig_base = b_obj["primary_base"]
                    b_obj["rc_base"] = str(Seq(orig_base).reverse_complement()) if orig_base in "ACGT" else orig_base
                    rev_map[r_pos] = b_obj

        consensus_chars = []
        rcrs_comparison = []
        variants_list = []
        ambiguity_count = 0
        conflict_question_count = 0
        covered_positions = 0

        # Rà soát từng nucleotide trong dải mục tiêu
        for r_pos in range(target_start, target_end + 1):
            ref_base = self.rcrs_seq[r_pos - 1] if r_pos <= len(self.rcrs_seq) else "N"
            f_info = fwd_map.get(r_pos)
            r_info = rev_map.get(r_pos)

            call_base = "N"
            status_note = ""

            if f_info and r_info:
                covered_positions += 1
                f_base = f_info["primary_base"]
                r_base = r_info["rc_base"]

                if f_base == r_base:
                    call_base = f_base
                else:
                    # Xung đột giữa Forward và Reverse -> Phân giải theo Phred & Tracy Purity
                    f_q = f_info["phred_score"]
                    r_q = r_info["phred_score"]
                    f_ratio = f_info["secondary_ratio"]
                    r_ratio = r_info["secondary_ratio"]

                    # Ưu tiên chiều sạch đối diện bù vào dải sau Poly-C bị trượt
                    if f_info.get("is_post_polyc") and not r_info.get("is_post_polyc") and r_q >= 20:
                        call_base = r_base
                        status_note = "Resolved (Rev clean overrides Fwd post-PolyC stutter; Fwd kept for comparison)"
                    elif r_info.get("is_post_polyc") and not f_info.get("is_post_polyc") and f_q >= 20:
                        call_base = f_base
                        status_note = "Resolved (Fwd clean overrides Rev post-PolyC stutter; Rev kept for comparison)"
                    # Trường hợp 1 bên sạch, 1 bên dính nhiễu/trượt
                    elif f_q >= 25 and f_ratio < 0.15 and (r_q < 18 or r_ratio > 0.35):
                        call_base = f_base
                        status_note = "Resolved (Fwd clean overrides Rev noise)"
                    elif r_q >= 25 and r_ratio < 0.15 and (f_q < 18 or f_ratio > 0.35):
                        call_base = r_base
                        status_note = "Resolved (Rev clean overrides Fwd noise)"
                    elif f_ratio >= 0.20 and r_ratio >= 0.20:
                        # Dị thể điểm thực sự xuất hiện trên cả 2 chiều
                        pair = tuple(sorted([f_base, r_base]))
                        call_base = IUPAC_MAP.get(pair, "+")
                        status_note = "True Heteroplasmy (Dual Peak Confirmed)"
                    else:
                        # Không thể phân giải an toàn
                        call_base = "+"
                        ambiguity_count += 1
                        status_note = "Ambiguity (+)"

            elif f_info:
                covered_positions += 1
                if f_info.get("called_base") == "?":
                    conflict_question_count += 1
                    call_base = "?"
                else:
                    call_base = f_info["primary_base"]

            elif r_info:
                covered_positions += 1
                if r_info.get("called_base") == "?":
                    conflict_question_count += 1
                    call_base = "?"
                else:
                    call_base = r_info["rc_base"]

            else:
                call_base = "N"

            consensus_chars.append(call_base)

            # So sánh với rCRS
            if call_base == ref_base:
                rcrs_comparison.append(".")
            elif call_base in ["N", "?"]:
                rcrs_comparison.append(call_base)
            elif call_base == "+":
                rcrs_comparison.append("+")
            else:
                rcrs_comparison.append(call_base)
                # Ghi nhận biến dị SNP theo chuẩn quốc tế
                variants_list.append(f"{r_pos}{ref_base}>{call_base}")

        consensus_seq = "".join(consensus_chars)
        match_dots_seq = "".join(rcrs_comparison)
        coverage_pct = round((covered_positions / target_len) * 100, 1)

        # Đánh giá QC Status theo tiêu chuẩn nghiêm ngặt của phòng xét nghiệm
        deny_reasons = []

        # 1. Kiểm tra độ bao phủ dải mục tiêu (HV1: 16024-16365, HV2: 73-576)
        if coverage_pct < 85.0:
            deny_reasons.append(f"Độ bao phủ dải mục tiêu không đạt ({coverage_pct}% < 85%), mồi xuôi và mồi ngược không phủ kín")

        # 2. Kiểm tra khoảng đứt đoạn lớn (consecutive Ns)
        max_consecutive_n = 0
        current_n = 0
        for ch in consensus_chars:
            if ch == "N":
                current_n += 1
                if current_n > max_consecutive_n:
                    max_consecutive_n = current_n
            else:
                current_n = 0
        if max_consecutive_n >= 20:
            deny_reasons.append(f"Bị đứt đoạn {max_consecutive_n} nucleotide liên tiếp trong dải mục tiêu (mồi xuôi và ngược không giao nhau)")

        # 3. Kiểm tra chập peak, nhiễu cảm biến/quang học (?) và mâu thuẫn 2 chiều (+)
        if conflict_question_count >= 5:
            deny_reasons.append(f"Peak chồng chéo/lỗi cảm biến quá nặng ({conflict_question_count} điểm '?')")
        if ambiguity_count >= 3:
            deny_reasons.append(f"Mâu thuẫn 2 chiều không thể phân giải ({ambiguity_count} điểm '+')")

        # 4. Kiểm tra số lượng biến dị dị thường (dấu hiệu mẫu nhiễm tạp hoặc chạy hỏng)
        if len(variants_list) >= 20:
            deny_reasons.append(f"Số lượng đột biến dị thường ({len(variants_list)} SNPs), nghi ngờ mẫu nhiễm tạp hoặc chạy hỏng")

        if deny_reasons:
            qc_status = "DENY"
            qc_recommendation = "DENY (Khuyến cáo làm lại mẫu trong phòng thí nghiệm)"
            qc_details = " | ".join(deny_reasons)
        elif ambiguity_count == 0 and conflict_question_count == 0 and coverage_pct >= 98.0 and len(variants_list) < 20:
            qc_status = "PASS"
            qc_recommendation = "PASS (Đạt tiêu chuẩn chất lượng cao)"
            qc_details = "Chất lượng giải trình tự tuyệt đối"
        elif ambiguity_count <= 2 and conflict_question_count <= 2 and coverage_pct >= 88.0 and len(variants_list) < 20:
            qc_status = "REVIEW_NEEDED"
            qc_recommendation = "REVIEW_NEEDED (Cần chuyên viên soi lại biểu đồ điện di đồ)"
            qc_details = "Tồn tại điểm dị thể hoặc độ phủ chấp nhận được"
        else:
            qc_status = "DENY"
            qc_recommendation = "DENY (Khuyến cáo làm lại mẫu trong phòng thí nghiệm)"
            qc_details = "Không đạt tiêu chuẩn kiểm soát chất lượng tối thiểu"

        return {
            "region": region,
            "target_range": f"{target_start}-{target_end}",
            "target_length": target_len,
            "covered_positions": covered_positions,
            "coverage_pct": coverage_pct,
            "consensus_sequence": consensus_seq,
            "rcrs_comparison": match_dots_seq,
            "variants": variants_list,
            "variant_count": len(variants_list),
            "ambiguity_count": ambiguity_count,
            "conflict_question_count": conflict_question_count,
            "qc_status": qc_status,
            "qc_recommendation": qc_recommendation,
            "qc_details": qc_details
        }
