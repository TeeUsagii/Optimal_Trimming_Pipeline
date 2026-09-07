from typing import Dict, List, Any, Tuple
import math

IUPAC_MAP = {
    tuple(sorted(['A', 'G'])): 'R',
    tuple(sorted(['C', 'T'])): 'Y',
    tuple(sorted(['G', 'C'])): 'S',
    tuple(sorted(['A', 'T'])): 'W',
    tuple(sorted(['G', 'T'])): 'K',
    tuple(sorted(['A', 'C'])): 'M',
}

class TracyEngine:
    """
    Thuật toán Basecaller mô phỏng Tracy (Gear Genomics):
    - Quét từng basecall interval quanh vị trí đỉnh (PLOC).
    - Tính tỷ lệ tín hiệu huỳnh quang cực đại của 4 kênh (A, C, G, T).
    - Nhận diện trạng thái '?' (Peak Collision / Sensor Crosstalk) kèm % xác suất khả thi.
    - Phát hiện hiện tượng trượt Poly-C Stutter và đánh dấu điểm cắt lý tưởng.
    """

    def __init__(self, peak_threshold: float = 0.25, conflict_threshold: float = 0.70):
        self.peak_threshold = peak_threshold
        self.conflict_threshold = conflict_threshold

    def analyze(self, ab1_data: Dict[str, Any]) -> Dict[str, Any]:
        """Phân tích chi tiết từng nucleotide trong file .ab1"""
        channel_traces = ab1_data["channel_traces"]
        ploc = ab1_data["peak_locations"]
        qualities = ab1_data["phred_qualities"]
        num_bases = len(ploc)

        analyzed_bases = []
        consecutive_c = 0
        in_poly_c = False
        stutter_mode = False

        for i in range(num_bases):
            peak_pos = ploc[i]
            
            # 1. Tính toán Interval cửa sổ xung quanh đỉnh
            left = (ploc[i-1] + peak_pos) // 2 if i > 0 else max(0, peak_pos - 6)
            right = (peak_pos + ploc[i+1]) // 2 if i < num_bases - 1 else peak_pos + 6

            # 2. Tìm RFU cực đại của 4 kênh trong interval
            signals = {}
            for b in ['A', 'C', 'G', 'T']:
                trace_slice = channel_traces[b][left:right+1]
                signals[b] = max(trace_slice) if len(trace_slice) > 0 else 0

            total_rfu = sum(signals.values())
            
            # Sắp xếp theo cường độ tín hiệu giảm dần
            sorted_sig = sorted(signals.items(), key=lambda x: x[1], reverse=True)
            top1_base, top1_rfu = sorted_sig[0]
            top2_base, top2_rfu = sorted_sig[1]

            # 3. Tính tỷ lệ đỉnh phụ & độ tinh khiết
            ratio = round(top2_rfu / top1_rfu, 3) if top1_rfu > 0 else 0.0
            purity = max(0, min(100, int((1.0 - ratio) * 100)))
            q_score = qualities[i] if i < len(qualities) else 10

            # 4. Tính % khả thi của từng base
            probs = {}
            for b in ['A', 'C', 'G', 'T']:
                p = round((signals[b] / total_rfu) * 100, 1) if total_rfu > 0 else 25.0
                probs[b] = p

            # 5. Phát hiện trạng thái '?' (Chồng peak / lỗi sensor)
            # Điều kiện: 2 đỉnh ngang ngửa nhau (ratio >= 0.70 và Phred < 18) HOẶC không có đỉnh nào vượt trội (p_top < 45%)
            is_conflict = (ratio >= self.conflict_threshold and q_score < 18) or (probs[top1_base] < 45.0)

            # 6. Phát hiện Poly-C và trượt Stutter
            if top1_base == 'C' and not stutter_mode:
                if ratio < 0.30 and q_score >= 20:
                    consecutive_c += 1
                    if consecutive_c >= 5:
                        in_poly_c = True
                elif in_poly_c and (ratio >= 0.30 or q_score < 20):
                    # Tín hiệu C bắt đầu bị vỡ / xuất hiện peak phụ cao -> Bắt đầu trượt Stutter!
                    stutter_mode = True
                    in_poly_c = False
            else:
                if in_poly_c:
                    stutter_mode = True  # Đã đi qua cụm Poly-C sang base khác
                    in_poly_c = False
                consecutive_c = 0

            # Xác định nhãn trạng thái (Status)
            if stutter_mode and (ratio > 0.30 or q_score < 22):
                status = "STUTTER"
                called_char = top1_base
                pair = tuple(sorted([top1_base, top2_base]))
                iupac = IUPAC_MAP.get(pair, top1_base)
            elif is_conflict:
                called_char = "?"
                status = "CONFLICT_SENSOR"
                iupac = "?"
            elif in_poly_c:
                status = "POLYC"
                iupac = "C"
                called_char = "C"
            elif ratio >= self.peak_threshold and q_score >= 18:
                status = "HETEROPLASMY"
                pair = tuple(sorted([top1_base, top2_base]))
                iupac = IUPAC_MAP.get(pair, top1_base)
                called_char = top1_base
            else:
                status = "CLEAN"
                iupac = top1_base
                called_char = top1_base

            # Tạo danh sách gợi ý ứng viên cho vị trí ?
            candidates = sorted(probs.items(), key=lambda x: x[1], reverse=True)

            analyzed_bases.append({
                "index": i + 1,
                "trace_pos": peak_pos,
                "called_base": called_char,
                "primary_base": top1_base,
                "secondary_base": top2_base if ratio >= self.peak_threshold else "-",
                "iupac": iupac,
                "phred_score": q_score,
                "secondary_ratio": ratio,
                "purity_pct": purity,
                "signals": signals,
                "probabilities": probs,
                "candidates": candidates,
                "status": status,
                "is_conflict": is_conflict
            })

        # 7. Đánh dấu điểm cắt tại C cuối cùng của Poly-C trước khi vỡ tín hiệu
        for i in range(len(analyzed_bases) - 1):
            if analyzed_bases[i]["status"] == "POLYC" and analyzed_bases[i+1]["status"] in ["STUTTER", "CONFLICT_SENSOR"]:
                analyzed_bases[i]["status"] = "POLYC_END"

        primary_sequence = "".join([b["primary_base"] for b in analyzed_bases])
        called_sequence = "".join([b["called_base"] for b in analyzed_bases])
        iupac_sequence = "".join([b["iupac"] for b in analyzed_bases])

        return {
            "file_name": ab1_data["file_name"],
            "total_bases": num_bases,
            "primary_sequence": primary_sequence,
            "called_sequence": called_sequence,
            "iupac_sequence": iupac_sequence,
            "bases": analyzed_bases
        }
