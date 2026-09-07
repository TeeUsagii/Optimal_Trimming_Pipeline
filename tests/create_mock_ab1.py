import os
import struct
from typing import Dict, List
from Bio.Seq import Seq

def write_abif(file_path: str, fwo: str, pbas: str, pcon: List[int], ploc: List[int], traces: Dict[str, List[int]]):
    """
    Tạo file nhị phân ABIF (.ab1) hợp lệ chuẩn ABI 3130/3500/SeqStudio tương thích 100% Biopython AbiIO.
    """
    entries = []
    data_bytes = bytearray()

    def add_entry(name: str, num: int, elem_type: int, elem_size: int, num_elems: int, raw_val_bytes: bytes):
        nonlocal data_bytes
        data_len = len(raw_val_bytes)
        if data_len <= 4:
            # Dữ liệu nhỏ hơn hoặc bằng 4 bytes được lưu trực tiếp vào trường data_offset (byte 20-23)
            offset_bytes = raw_val_bytes.ljust(4, b'\x00')
            offset_val = struct.unpack(">I", offset_bytes)[0]
        else:
            offset_val = 30 + len(data_bytes) # 30 bytes là header (4 magic + 26 headfmt)
            data_bytes.extend(raw_val_bytes)
        
        entries.append({
            "name": name.encode('ascii')[:4].ljust(4, b' '),
            "num": num,
            "type": elem_type,
            "elem_size": elem_size,
            "num_elems": num_elems,
            "data_size": data_len,
            "offset": offset_val
        })

    # 1. FWO_1
    add_entry("FWO_", 1, 2, 1, 4, fwo.encode('ascii'))
    
    # 2. PBAS1 & PBAS2 (Type 2 = pString/ASCII)
    pbas_bytes = pbas.encode('ascii')
    add_entry("PBAS", 1, 2, 1, len(pbas), pbas_bytes)
    add_entry("PBAS", 2, 2, 1, len(pbas), pbas_bytes)

    # 3. PCON1 & PCON2 (Type 2 = bytes quality scores)
    pcon_bytes = bytes(pcon)
    add_entry("PCON", 1, 2, 1, len(pcon), pcon_bytes)
    add_entry("PCON", 2, 2, 1, len(pcon), pcon_bytes)

    # 4. PLOC1 & PLOC2 (Type 4 = 16-bit short)
    ploc_bytes = struct.pack(f">{len(ploc)}h", *ploc)
    add_entry("PLOC", 1, 4, 2, len(ploc), ploc_bytes)
    add_entry("PLOC", 2, 4, 2, len(ploc), ploc_bytes)

    # 5. DATA9 - DATA12 (Type 4 = 16-bit short)
    trace_len = len(traces[fwo[0]])
    for i in range(4):
        b_char = fwo[i]
        tr_bytes = struct.pack(f">{trace_len}h", *traces[b_char])
        add_entry("DATA", 9 + i, 4, 2, trace_len, tr_bytes)

    dir_offset = 30 + len(data_bytes)

    # Header: 4 bytes magic + 26 bytes _HEADFMT (>H4sI2H3I - 8 items)
    header = b'ABIF' + struct.pack(
        ">H4sI2H3I",
        101,              # version (H)
        b'tdir',          # tag_name (4s)
        1,                # tag_number (I)
        1023,             # elem_type (H)
        28,               # elem_size (H)
        len(entries),     # num_elems (I)
        len(entries) * 28,# data_size (I)
        dir_offset        # data_offset (I)
    )

    dir_data = bytearray()
    for e in entries:
        dir_data.extend(struct.pack(
            ">4sI2H4I",
            e["name"],
            e["num"],
            e["type"],
            e["elem_size"],
            e["num_elems"],
            e["data_size"],
            e["offset"],
            0
        ))

    with open(file_path, "wb") as f:
        f.write(header)
        f.write(data_bytes)
        f.write(dir_data)

def generate_sample_dataset(output_dir: str, rcrs_path: str):
    os.makedirs(output_dir, exist_ok=True)
    with open(rcrs_path, "r") as f:
        lines = [line.strip() for line in f if not line.startswith(">")]
        rcrs_seq = "".join(lines).upper()

    print(f"[+] Đang tạo bộ dữ liệu thử nghiệm .ab1 mô phỏng tại: {output_dir}")

    samples = ["Sample01", "Sample02"]
    fwo = "GATC"

    for sid in samples:
        # 1. HV1 Forward (F15971): 15971 đến 16390 (~420 bp)
        hv1_sub = rcrs_seq[15970:16390]
        # Thêm biến dị 16189T>C cho Sample01
        if sid == "Sample01" and len(hv1_sub) > (16189 - 15971):
            mut_idx = 16189 - 15971
            hv1_sub = hv1_sub[:mut_idx] + "C" + hv1_sub[mut_idx+1:]

        num_b = len(hv1_sub)
        ploc = [30 + i * 11 for i in range(num_b)]
        total_trace_pts = ploc[-1] + 30
        traces = {b: [15] * total_trace_pts for b in ['A', 'C', 'G', 'T']}
        quals = []

        for i, base in enumerate(hv1_sub):
            p = ploc[i]
            # Mồi nhiễu 30 bp đầu
            if i < 30:
                quals.append(12)
                traces[base][p] = 500
                traces['T' if base != 'T' else 'A'][p] = 400
            elif 16184 <= (15971 + i) <= 16193:
                # Poly-C
                quals.append(38)
                traces['C'][p] = 1400
            elif (15971 + i) > 16193 and sid == "Sample01":
                # Stutter do 16189T>C!
                quals.append(14)
                traces[base][p] = 450
                traces['C'][p] = 380
            else:
                quals.append(45)
                traces[base][p] = 1600

        write_abif(os.path.join(output_dir, f"{sid}_HV1_F.ab1"), fwo, hv1_sub, quals, ploc, traces)

        # 2. HV1 Reverse (R16410): từ 16410 lùi về 15990
        hv1_rev_region = rcrs_seq[15990:16410]
        hv1_rev_sub = str(Seq(hv1_rev_region).reverse_complement())
        num_b_r = len(hv1_rev_sub)
        ploc_r = [30 + i * 11 for i in range(num_b_r)]
        total_trace_pts_r = ploc_r[-1] + 30
        traces_r = {b: [15] * total_trace_pts_r for b in ['A', 'C', 'G', 'T']}
        quals_r = [45] * num_b_r
        for i, base in enumerate(hv1_rev_sub):
            p = ploc_r[i]
            traces_r[base][p] = 1500

        write_abif(os.path.join(output_dir, f"{sid}_HV1_R.ab1"), fwo, hv1_rev_sub, quals_r, ploc_r, traces_r)

        # 3. HV2-HV3 Forward (F15): 15 đến 620
        hv2_f_sub = rcrs_seq[14:620]
        num_b2 = len(hv2_f_sub)
        ploc2 = [30 + i * 11 for i in range(num_b2)]
        traces2 = {b: [15] * (ploc2[-1] + 30) for b in ['A', 'C', 'G', 'T']}
        quals2 = []
        for i, base in enumerate(hv2_f_sub):
            p = ploc2[i]
            if i == 80: # Điểm nhiễu ?
                quals2.append(8)
                traces2['A'][p] = 600
                traces2['G'][p] = 580
            else:
                quals2.append(42)
                traces2[base][p] = 1400

        write_abif(os.path.join(output_dir, f"{sid}_HV2_F.ab1"), fwo, hv2_f_sub, quals2, ploc2, traces2)

        # 4. HV2-HV3 Reverse (R639): 639 lùi về 20
        hv2_r_sub = str(Seq(rcrs_seq[19:639]).reverse_complement())
        num_b2_r = len(hv2_r_sub)
        ploc2_r = [30 + i * 11 for i in range(num_b2_r)]
        traces2_r = {b: [15] * (ploc2_r[-1] + 30) for b in ['A', 'C', 'G', 'T']}
        quals2_r = [44] * num_b2_r
        for i, base in enumerate(hv2_r_sub):
            p = ploc2_r[i]
            traces2_r[base][p] = 1500

        write_abif(os.path.join(output_dir, f"{sid}_HV2_R.ab1"), fwo, hv2_r_sub, quals2_r, ploc2_r, traces2_r)

    print(f"[PASS] Đã tạo thành công dữ liệu mẫu cho {len(samples)} mẫu (8 file .ab1 chuẩn ABIF)!")

if __name__ == "__main__":
    out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "sample_ab1")
    rcrs = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "rCRS.fasta")
    generate_sample_dataset(out_dir, rcrs)
