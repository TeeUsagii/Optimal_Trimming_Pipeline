import os
from typing import Dict, List, Optional, Any
from Bio import SeqIO

class AB1Parser:
    """
    Parser chuyên dụng bóc tách dữ liệu nhị phân ABIF (.ab1) từ máy mao quản (SeqStudio / ABI 3500).
    Giải mã Filter Wheel Order (FWO_1) để đảm bảo 4 kênh huỳnh quang A, C, G, T luôn chính xác.
    """

    @staticmethod
    def parse(file_path: str) -> Optional[Dict[str, Any]]:
        """
        Bóc tách toàn bộ thông tin từ file .ab1:
        - raw_sequence: Chuỗi base KB Basecaller đã gọi sẵn
        - phred_qualities: Mảng điểm chất lượng Phred (0-60)
        - peak_locations: Mảng vị trí đỉnh (PLOC) trên đồ thị điện di
        - channel_traces: Dữ liệu sóng huỳnh quang của 4 kênh {'A': [...], 'C': [...], 'G': [...], 'T': [...]}
        - fwo: Thứ tự màng lọc huỳnh quang (Filter Wheel Order)
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File không tồn tại: {file_path}")

        try:
            record = SeqIO.read(file_path, "abi")
        except Exception as e:
            print(f"[-] Lỗi đọc file ABIF {file_path}: {e}")
            return None

        raw = record.annotations.get("abif_raw", {})
        if not raw:
            print(f"[-] File {file_path} không chứa thẻ abif_raw hợp lệ.")
            return None

        # 1. Giải mã FWO_1 (Filter Wheel Order)
        fwo_tag = raw.get('FWO_1', b'GATC')
        fwo = fwo_tag.decode('utf-8') if isinstance(fwo_tag, bytes) else str(fwo_tag)
        fwo = fwo.upper().strip()
        if len(fwo) < 4:
            fwo = 'GATC'

        # 2. Ánh xạ DATA9-12 vào đúng nucleotide
        # DATA9 tương ứng với fwo[0], DATA10 với fwo[1], DATA11 với fwo[2], DATA12 với fwo[3]
        channel_traces = {
            fwo[0]: list(raw.get('DATA9', [])),
            fwo[1]: list(raw.get('DATA10', [])),
            fwo[2]: list(raw.get('DATA11', [])),
            fwo[3]: list(raw.get('DATA12', []))
        }

        # 3. Lấy vị trí đỉnh (PLOC2 là analyzed peak locations, fallback PLOC1)
        ploc = list(raw.get('PLOC2', raw.get('PLOC1', [])))

        # 4. Lấy chuỗi gọi sẵn (PBAS2 là analyzed basecalls, fallback PBAS1)
        pbas = raw.get('PBAS2', raw.get('PBAS1', b''))
        if isinstance(pbas, bytes):
            pbas = pbas.decode('utf-8')
        raw_sequence = str(pbas).upper().strip()

        # 5. Lấy điểm chất lượng Phred
        qualities = record.letter_annotations.get("phred_quality", [])
        if not qualities:
            pcon = raw.get('PCON2', raw.get('PCON1', b''))
            qualities = list(pcon)

        # 6. Đảm bảo tính đồng bộ về độ dài
        min_len = min(len(raw_sequence), len(ploc))
        if len(qualities) < min_len:
            qualities.extend([10] * (min_len - len(qualities)))

        return {
            "file_name": os.path.basename(file_path),
            "file_path": file_path,
            "raw_sequence": raw_sequence[:min_len],
            "phred_qualities": qualities[:min_len],
            "peak_locations": ploc[:min_len],
            "channel_traces": channel_traces,
            "fwo": fwo,
            "total_bases": min_len
        }
