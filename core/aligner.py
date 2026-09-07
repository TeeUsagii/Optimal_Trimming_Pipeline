import os
from typing import Dict, Any, Optional, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

class RCRSAligner:
    """
    Module căn chỉnh chuỗi đọc Sanger vào hệ tọa độ chuẩn rCRS (NC_012920.1).
    Sử dụng thuật toán quy hoạch động cục bộ (Smith-Waterman) qua PairwiseAligner.
    Tạo bản đồ ánh xạ song ánh giữa vị trí đọc (Read Index) và tọa độ rCRS chuẩn.
    """

    def __init__(self, rcrs_path: str):
        if not os.path.exists(rcrs_path):
            raise FileNotFoundError(f"Không tìm thấy file rCRS tại: {rcrs_path}")

        record = SeqIO.read(rcrs_path, "fasta")
        self.rcrs_seq = str(record.seq).upper()
        self.rcrs_len = len(self.rcrs_seq)

        # Cấu hình Aligner tối ưu cho Sanger mtDNA
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'local'
        self.aligner.match_score = 2.0
        self.aligner.mismatch_score = -1.0
        self.aligner.open_gap_score = -3.0
        self.aligner.extend_gap_score = -1.0

    def get_target_window(self, region: str) -> Tuple[int, int, str]:
        """Lấy phân đoạn rCRS tương ứng với vùng khuếch đại để tăng tốc độ align gấp 100 lần"""
        reg = region.upper()
        if "HV1" in reg:
            # Mồi F15971 / R16410 -> Lấy cửa sổ 15900 đến 16500
            start = 15900
            end = 16500
            subseq = self.rcrs_seq[start - 1 : end]
            return start, end, subseq
        else:
            # HV2/HV3 (Mồi F15 / R639) -> Lấy cửa sổ 1 đến 700
            start = 1
            end = 700
            subseq = self.rcrs_seq[start - 1 : end]
            return start, end, subseq

    def align(self, query_seq: str, region: str, is_reverse: bool = False) -> Optional[Dict[str, Any]]:
        """
        Căn chỉnh query_seq vào rCRS:
        - Nếu is_reverse=True: đảo ngược bổ sung (reverse complement) để align cùng chiều rCRS.
        - Trả về bản đồ ánh xạ rcrs_to_query và query_to_rcrs.
        """
        ref_start_pos, ref_end_pos, ref_subseq = self.get_target_window(region)

        # Nếu là Reverse, lấy Reverse Complement để so cùng chiều với rCRS
        align_query = str(Seq(query_seq).reverse_complement()) if is_reverse else query_seq

        # 1. Kiểm tra sơ bộ điểm alignment để loại bỏ mẫu trắng/âm tính (Negative Control / Primer Dimer)
        try:
            score = self.aligner.score(ref_subseq, align_query)
            if score < 100.0:  # Không đủ độ dài hoặc tương đồng với rCRS
                return None
            
            # 2. Lấy alignment tối ưu đầu tiên mà không gọi len(alignments) để tránh OverflowError
            alignments = self.aligner.align(ref_subseq, align_query)
            best_aln = next(iter(alignments), None)
            if best_aln is None:
                return None
        except Exception:
            return None

        ref_blocks, query_blocks = best_aln.aligned

        rcrs_to_query = {}
        query_to_rcrs = {}

        for (r_start, r_end), (q_start, q_end) in zip(ref_blocks, query_blocks):
            for i in range(r_end - r_start):
                r_coord = ref_start_pos + r_start + i
                q_idx = q_start + i
                
                # Nếu là Reverse, ánh xạ lại chỉ số gốc trên file read ban đầu
                orig_read_idx = (len(query_seq) - 1 - q_idx) if is_reverse else q_idx

                rcrs_to_query[r_coord] = orig_read_idx
                query_to_rcrs[orig_read_idx] = r_coord

        return {
            "score": best_aln.score,
            "ref_start_mapped": min(rcrs_to_query.keys()) if rcrs_to_query else None,
            "ref_end_mapped": max(rcrs_to_query.keys()) if rcrs_to_query else None,
            "rcrs_to_query": rcrs_to_query,
            "query_to_rcrs": query_to_rcrs,
            "is_reverse": is_reverse,
            "alignment_obj": best_aln
        }
