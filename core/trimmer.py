from typing import Dict, Any, Tuple, Optional, List
import numpy as np

# Khung tọa độ chuẩn theo hướng dẫn quốc tế (ISFG / SWGDAM mtDNA guidelines)
TARGET_BOUNDS = {
    "HV1": {
        "start": 16024,
        "end": 16365,
        "poly_c_start": 16184,
        "poly_c_end": 16193
    },
    "HV2_HV3": {
        "start": 73,
        "end": 576,
        "poly_c_start": 303,
        "poly_c_end": 315
    }
}

class SmartTrimmer:
    """
    Bộ thuật toán dự đoán vị trí loại bỏ vùng nhiễu không phá hủy (Non-destructive Noise Boundary Prediction):
    1. Dự đoán điểm cắt đầu 5': Dò tìm điểm thoát khỏi vùng nhiễu mồi/injection noise (10-20 bp đầu)
       dựa trên cửa sổ trượt điểm Phred (Q >= 20) và độ ổn định khoảng cách đỉnh (PLOC).
    2. Dự đoán điểm cắt đuôi 3': Dò tìm điểm rơi chất lượng cuối mao quản hoặc giới hạn rCRS.
    3. Xử lý chuyên biệt dải Poly-C: Bảo toàn toàn bộ chuỗi sau Poly-C để làm dữ liệu đối sánh
       giữa mồi xuôi và mồi ngược (Dual-strand comparative data), KHÔNG cắt cụt.
    """

    def __init__(self):
        pass

    @staticmethod
    def _get_phred(b: Dict[str, Any]) -> float:
        if "phred_score" in b:
            return float(b["phred_score"])
        return 35.0 if b.get("status") == "CLEAN" else 10.0

    @classmethod
    def _find_5prime_quality_boundary(cls, bases: List[Dict[str, Any]], window_size: int = 5, min_q: float = 20.0) -> int:
        """
        Tìm chỉ số base đầu tiên thoát khỏi vùng nhiễu đầu 5' bằng cửa sổ trượt Phred Q >= 20.
        """
        n = len(bases)
        if n < window_size:
            return 0

        quals = [cls._get_phred(b) for b in bases]
        for i in range(n - window_size + 1):
            window_avg = sum(quals[i:i + window_size]) / window_size
            # Đảm bảo cửa sổ đạt chuẩn và base hiện tại không phải là N/? với Q < 15
            if window_avg >= min_q and quals[i] >= 15:
                return i
        return 0

    @classmethod
    def _find_3prime_quality_boundary(cls, bases: List[Dict[str, Any]], start_check: int, min_q: float = 15.0, consecutive_low: int = 15) -> int:
        """
        Dò tìm điểm sụt giảm chất lượng dốc đứng ở đuôi 3' (ddNTP depletion & POP-7 band broadening).
        """
        n = len(bases)
        quals = [cls._get_phred(b) for b in bases]

        low_count = 0
        cutoff_idx = n

        for i in range(start_check, n):
            if quals[i] < min_q:
                low_count += 1
                if low_count >= consecutive_low:
                    # Chốt điểm cắt tại vị trí bắt đầu chuỗi sụt giảm
                    cutoff_idx = i - consecutive_low + 1
                    break
            else:
                low_count = 0

        return max(start_check, cutoff_idx)

    @classmethod
    def evaluate_boundaries(cls,
                            analyzed_data: Dict[str, Any],
                            alignment_info: Optional[Dict[str, Any]],
                            region: str,
                            is_reverse: bool) -> Dict[str, Any]:
        """
        Dự đoán vị trí loại bỏ vùng nhiễu và vùng mồi, bảo toàn 100% dữ liệu gốc.
        """
        reg_key = "HV1" if "HV1" in region.upper() else "HV2_HV3"
        bounds = TARGET_BOUNDS[reg_key]
        target_start = bounds["start"]
        target_end = bounds["end"]
        poly_c_start = bounds["poly_c_start"]
        poly_c_end = bounds["poly_c_end"]

        raw_seq = analyzed_data["primary_sequence"]
        bases = analyzed_data["bases"]
        num_bases = len(raw_seq)

        if not alignment_info or "rcrs_to_query" not in alignment_info:
            # Fallback nếu chưa align được: dùng thuần chất lượng sóng
            q_5p = cls._find_5prime_quality_boundary(bases, window_size=5, min_q=20.0)
            q_3p = cls._find_3prime_quality_boundary(bases, start_check=max(q_5p + 50, num_bases // 2))
            return {
                "raw_sequence": raw_seq,
                "raw_bases": bases,
                "trimmed_seq": raw_seq[q_5p:q_3p],
                "trimmed_bases": bases[q_5p:q_3p],
                "trim_start_idx": q_5p,
                "trim_end_idx": q_3p,
                "suggested_5p_cut_1based": q_5p + 1,
                "suggested_3p_cut_1based": q_3p,
                "retained_length": q_3p - q_5p,
                "rcrs_span": (0, 0),
                "poly_c_detected": False,
                "poly_c_status": "NONE",
                "noise_5p_reason": f"Quality-based 5' trim (Q >= 20 window at base {q_5p + 1})",
                "noise_3p_reason": f"Quality-based 3' trim (Signal decay at base {q_3p})",
                "trim_reason": "UNALIGNED_QUALITY_TRIM",
                "is_reverse": is_reverse,
                "stutter_detected": False
            }

        rcrs_to_query = alignment_info["rcrs_to_query"]
        query_to_rcrs = alignment_info.get("query_to_rcrs", {})

        # 1. Phát hiện tình trạng Poly-C và trượt Stutter
        stutter_detected = any(b.get("status") in ["STUTTER", "POLYC_END"] for b in bases)
        poly_c_present = any(b.get("status") in ["POLYC", "POLYC_END"] for b in bases)

        if stutter_detected:
            poly_c_status = "STUTTER_RETAINED_FOR_COMPARISON"
        elif poly_c_present:
            poly_c_status = "CLEAN"
        else:
            poly_c_status = "NONE"

        # 2. Xác định ranh giới chất lượng sóng
        q_5p = cls._find_5prime_quality_boundary(bases, window_size=5, min_q=20.0)
        q_3p = cls._find_3prime_quality_boundary(bases, start_check=max(q_5p + 100, num_bases // 2))

        # 3. Kết hợp ranh giới chất lượng với Tọa độ mục tiêu rCRS
        if not is_reverse:  # CHUỖI XUÔI (FORWARD)
            # Đầu 5': Cắt bỏ mồi và nhiễu đầu (chọn điểm bắt đầu vào dải target_start hoặc điểm thoát nhiễu Q >= 20)
            if target_start in rcrs_to_query:
                rcrs_5p = rcrs_to_query[target_start]
                # Chọn điểm bảo đảm cả 2: vừa vào dải rCRS vừa thoát nhiễu mồi
                cut_5p = max(rcrs_5p, q_5p)
                noise_5p_reason = f"Primer & 5' early noise trimmed (Bases 1-{cut_5p}, Target rCRS {target_start})"
            else:
                cut_5p = q_5p
                noise_5p_reason = f"5' early noise trimmed (Bases 1-{cut_5p}, Q < 20)"

            # Đuôi 3': GIỮ LẠI DẢI SAU POLYC ĐỂ ĐỐI SÁNH!
            # Chỉ cắt tại giới hạn rCRS target_end hoặc khi tín hiệu sụt giảm thực sự (q_3p)
            if target_end in rcrs_to_query:
                rcrs_3p = rcrs_to_query[target_end] + 1
                cut_3p = min(rcrs_3p, q_3p)
                noise_3p_reason = f"Target rCRS {target_end} boundary reached"
            else:
                cut_3p = q_3p
                noise_3p_reason = f"Signal decay 3' boundary (Base {cut_3p})"

            if stutter_detected:
                noise_3p_reason += " | Poly-C stutter retained for dual-strand comparison"

        else:  # CHUỖI NGƯỢC (REVERSE)
            # Trên rCRS: Chiều Reverse đọc ngược từ target_end về target_start
            # Đầu 5' của Reverse read tương ứng với vùng target_end trên rCRS
            if target_end in rcrs_to_query:
                rcrs_5p = rcrs_to_query[target_end]
                cut_5p = max(rcrs_5p, q_5p)
                noise_5p_reason = f"Reverse primer & 5' noise trimmed (Target rCRS {target_end})"
            else:
                cut_5p = q_5p
                noise_5p_reason = f"Reverse 5' noise trimmed (Bases 1-{cut_5p}, Q < 20)"

            # Đuôi 3' của Reverse read đọc về phía target_start
            if target_start in rcrs_to_query:
                rcrs_3p = rcrs_to_query[target_start] + 1
                cut_3p = min(rcrs_3p, q_3p)
                noise_3p_reason = f"Target rCRS {target_start} boundary reached"
            else:
                cut_3p = q_3p
                noise_3p_reason = f"Signal decay 3' boundary (Base {cut_3p})"

            if stutter_detected:
                noise_3p_reason += " | Reverse strand provides high-confidence consensus over Poly-C"

        # Đảm bảo thứ tự chỉ số hợp lệ
        trim_start_idx = min(cut_5p, cut_3p)
        trim_end_idx = max(cut_5p, cut_3p)

        # Tránh trường hợp cắt hết sạch chuỗi
        if trim_end_idx <= trim_start_idx or (trim_end_idx - trim_start_idx) < 30:
            trim_start_idx = min(15, num_bases // 4)
            trim_end_idx = max(num_bases - 15, num_bases * 3 // 4)

        # Tính khoảng tọa độ rCRS thực tế được bao phủ
        mapped_rcrs = [query_to_rcrs[i] for i in range(trim_start_idx, trim_end_idx) if i in query_to_rcrs]
        rcrs_span = (int(min(mapped_rcrs)), int(max(mapped_rcrs))) if mapped_rcrs else (0, 0)

        # Đánh dấu thuộc tính is_comparative cho các base sau Poly-C
        for idx in range(num_bases):
            base_info = bases[idx]
            if idx in query_to_rcrs:
                r_pos = query_to_rcrs[idx]
                if not is_reverse:
                    base_info["is_post_polyc"] = bool((r_pos > poly_c_end) and stutter_detected)
                else:
                    base_info["is_post_polyc"] = bool((r_pos < poly_c_start) and stutter_detected)
            else:
                base_info["is_post_polyc"] = False

        return {
            "raw_sequence": raw_seq,
            "raw_bases": bases,
            "trimmed_seq": raw_seq[trim_start_idx:trim_end_idx],
            "trimmed_bases": bases[trim_start_idx:trim_end_idx],
            "trim_start_idx": trim_start_idx,
            "trim_end_idx": trim_end_idx,
            "suggested_5p_cut_1based": trim_start_idx + 1,
            "suggested_3p_cut_1based": trim_end_idx,
            "retained_length": trim_end_idx - trim_start_idx,
            "trimmed_length": trim_end_idx - trim_start_idx,
            "rcrs_span": rcrs_span,
            "poly_c_detected": poly_c_present or stutter_detected,
            "poly_c_status": poly_c_status,
            "noise_5p_reason": noise_5p_reason,
            "noise_3p_reason": noise_3p_reason,
            "trim_reason": f"5' Trim: {noise_5p_reason} | 3' Trim: {noise_3p_reason}",
            "is_reverse": is_reverse,
            "stutter_detected": stutter_detected
        }

    # Phương thức tương thích ngược cho code cũ
    @classmethod
    def trim(cls, analyzed_data: Dict[str, Any], alignment_info: Dict[str, Any], region: str, is_reverse: bool) -> Dict[str, Any]:
        return cls.evaluate_boundaries(analyzed_data, alignment_info, region, is_reverse)
