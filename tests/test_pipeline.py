import os
import sys
import unittest

CURRENT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

if sys.stdout and hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
if sys.stderr and hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8')

from core.tracy_engine import TracyEngine
from core.aligner import RCRSAligner
from core.trimmer import SmartTrimmer
from core.assembler import ContigAssembler

class TestOptimalTrimmingPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rcrs_path = os.path.join(CURRENT_DIR, "data", "rCRS.fasta")
        cls.aligner = RCRSAligner(cls.rcrs_path)
        cls.tracy = TracyEngine(peak_threshold=0.25, conflict_threshold=0.70)
        cls.trimmer = SmartTrimmer()
        cls.assembler = ContigAssembler(cls.aligner.rcrs_seq)

    def test_rcrs_loaded(self):
        """Kiểm tra độ dài rCRS chuẩn 16,569 bp"""
        self.assertEqual(self.aligner.rcrs_len, 16569)
        print("[PASS] Test 1: File rCRS.fasta nạp thành công (16,569 bp).")

    def test_tracy_engine_and_conflict_detection(self):
        """Kiểm tra TracyEngine nhận diện đúng Clean, Poly-C, Stutter, và Conflicted (?)"""
        # Giả lập dữ liệu từ AB1Parser
        # Vị trí 1: Clean A (Q=45)
        # Vị trí 2: Conflicted ? (A=500, G=480, Q=10)
        # Vị trí 3-8: 6 base C (Poly-C tract)
        # Vị trí 9+: Stutter (C=400, T=350, Q=12)
        num_bases = 12
        ploc = [100 + i * 12 for i in range(num_bases)]
        quals = [45, 10, 42, 40, 38, 35, 30, 22, 14, 12, 10, 8]
        
        # 4 traces
        traces = {'A': [0]*250, 'C': [0]*250, 'G': [0]*250, 'T': [0]*250}
        
        # Base 0: Clean A
        traces['A'][ploc[0]] = 1500; traces['C'][ploc[0]] = 20; traces['G'][ploc[0]] = 30; traces['T'][ploc[0]] = 10
        # Base 1: Conflict A vs G (?)
        traces['A'][ploc[1]] = 520; traces['G'][ploc[1]] = 500; traces['C'][ploc[1]] = 40; traces['T'][ploc[1]] = 30
        # Base 2-7: 6 C's
        for idx in range(2, 8):
            p = ploc[idx]
            traces['C'][p] = 1200; traces['A'][p] = 20; traces['G'][p] = 30; traces['T'][p] = 40
        # Base 8-11: Stutter (C vs T)
        for idx in range(8, 12):
            p = ploc[idx]
            traces['C'][p] = 450; traces['T'][p] = 380; traces['A'][p] = 180; traces['G'][p] = 150

        mock_ab1 = {
            "file_name": "test_mock.ab1",
            "peak_locations": ploc,
            "phred_qualities": quals,
            "channel_traces": traces
        }

        result = self.tracy.analyze(mock_ab1)
        bases = result["bases"]

        # Kiểm tra Base 0: Clean A
        self.assertEqual(bases[0]["status"], "CLEAN")
        self.assertEqual(bases[0]["called_base"], "A")

        # Kiểm tra Base 1: Conflict (?)
        self.assertTrue(bases[1]["is_conflict"])
        self.assertEqual(bases[1]["called_base"], "?")
        # Kiểm tra % xác suất
        pA = bases[1]["probabilities"]["A"]
        pG = bases[1]["probabilities"]["G"]
        self.assertAlmostEqual(pA + pG, 93.5, delta=3.0)

        # Kiểm tra Base 7: POLYC_END (Điểm cắt đề xuất trước khi vỡ tín hiệu)
        self.assertEqual(bases[7]["status"], "POLYC_END")

        # Kiểm tra Base 8: STUTTER
        self.assertEqual(bases[8]["status"], "STUTTER")
        print("[PASS] Test 2: TracyEngine phát hiện chính xác Clean, '?' (kèm % xác suất), Poly-C và Stutter.")

    def test_hv1_alignment_and_trimming_boundaries(self):
        """Kiểm tra Tiêu chuẩn 1 & 2: HV1 bắt đầu chính xác tại 16024 và bảo toàn rCRS"""
        # Lấy đoạn rCRS HV1 từ 15980 đến 16380 (bao gồm cả mồi ngoài)
        hv1_seq = self.aligner.rcrs_seq[15979:16380] # 15980 đến 16380
        
        # Align Forward
        aln = self.aligner.align(hv1_seq, "HV1", is_reverse=False)
        self.assertIsNotNone(aln)
        self.assertIn(16024, aln["rcrs_to_query"])
        self.assertIn(16365, aln["rcrs_to_query"])

        # Mock analyzed data
        num_bases = len(hv1_seq)
        mock_bases = [{"status": "CLEAN"} for _ in range(num_bases)]
        analyzed_data = {"primary_sequence": hv1_seq, "bases": mock_bases}

        trim_res = self.trimmer.trim(analyzed_data, aln, "HV1", is_reverse=False)

        # Kiểm tra rCRS span sau khi trim: Luôn bắt đầu từ 16024 đến 16365
        self.assertEqual(trim_res["rcrs_span"][0], 16024)
        self.assertEqual(trim_res["rcrs_span"][1], 16365)
        self.assertEqual(trim_res["trimmed_length"], 342)
        print(f"[PASS] Test 3: Tiêu chuẩn 1 & 2 được đảm bảo: HV1 Trimmed Length = {trim_res['trimmed_length']} bp (16024 -> 16365).")

    def test_hv2_hv3_alignment_and_trimming_boundaries(self):
        """Kiểm tra HV2-HV3 dải đích 73 - 576 (504 bp)"""
        hv2_seq = self.aligner.rcrs_seq[20:600] # từ 21 đến 600
        aln = self.aligner.align(hv2_seq, "HV2_HV3", is_reverse=False)
        self.assertIsNotNone(aln)
        self.assertIn(73, aln["rcrs_to_query"])
        self.assertIn(576, aln["rcrs_to_query"])

        mock_bases = [{"status": "CLEAN"} for _ in range(len(hv2_seq))]
        analyzed_data = {"primary_sequence": hv2_seq, "bases": mock_bases}

        trim_res = self.trimmer.trim(analyzed_data, aln, "HV2_HV3", is_reverse=False)
        self.assertEqual(trim_res["rcrs_span"][0], 73)
        self.assertEqual(trim_res["rcrs_span"][1], 576)
        self.assertEqual(trim_res["trimmed_length"], 504)
        print(f"[PASS] Test 4: HV2 & HV3 Trimmed Length = {trim_res['trimmed_length']} bp (73 -> 576).")

    def test_polyc_stutter_preserved_for_comparison(self):
        """Kiểm tra: Khi gặp Poly-C stutter, dữ liệu sau Poly-C vẫn được bảo toàn để làm dữ liệu đối sánh mồi xuôi & ngược"""
        hv1_seq = self.aligner.rcrs_seq[15979:16380] # 15980 đến 16380
        aln = self.aligner.align(hv1_seq, "HV1", is_reverse=False)
        self.assertIsNotNone(aln)

        # Giả lập base có Poly-C và Stutter ở đuôi
        mock_bases = [{"status": "CLEAN", "phred_score": 35} for _ in range(len(hv1_seq))]
        # Gán vị trí 16184-16193 là POLYC và sau đó là STUTTER
        for i, b in enumerate(mock_bases):
            r_coord = aln["query_to_rcrs"].get(i)
            if r_coord and 16184 <= r_coord <= 16193:
                b["status"] = "POLYC"
            elif r_coord and r_coord > 16193:
                b["status"] = "STUTTER"
                b["phred_score"] = 18

        analyzed_data = {"primary_sequence": hv1_seq, "bases": mock_bases}
        res = self.trimmer.evaluate_boundaries(analyzed_data, aln, "HV1", is_reverse=False)

        # Kiểm tra trạng thái Poly-C
        self.assertEqual(res["poly_c_status"], "STUTTER_RETAINED_FOR_COMPARISON")
        # Điểm cắt 3' không bị chặt tại 16193 mà tiếp tục bao phủ tới rCRS 16365
        self.assertIn(16365, aln["rcrs_to_query"])
        self.assertGreater(res["suggested_3p_cut_1based"], aln["rcrs_to_query"][16193])
        print("[PASS] Test 5: Dải sau Poly-C được bảo toàn trọn vẹn làm dữ liệu đối sánh mồi xuôi & mồi ngược.")

    def test_qc_deny_logic(self):
        """Kiểm tra logic phân loại QC: PASS, REVIEW_NEEDED và gán nhãn DENY khi dữ liệu quá xấu hoặc không bao phủ"""
        # 1. Test case DENY do đứt đoạn hoặc không bao phủ (< 85%)
        contig_fail_cov = self.assembler.assemble_region(
            region="HV1",
            fwd_data=None,
            fwd_trim=None,
            fwd_aln=None,
            rev_data=None,
            rev_trim=None,
            rev_aln=None
        )
        self.assertEqual(contig_fail_cov["qc_status"], "DENY")
        self.assertIn("DENY", contig_fail_cov["qc_recommendation"])
        self.assertIn("Độ bao phủ dải mục tiêu không đạt", contig_fail_cov["qc_details"])

        # 2. Test case PASS hoàn hảo
        hv1_seq = self.aligner.rcrs_seq[16023:16365] # 342 bp chuẩn
        aln = self.aligner.align(hv1_seq, "HV1", is_reverse=False)
        mock_bases = [{"called_base": b, "primary_base": b, "phred_score": 40, "is_conflict": False, "status": "CLEAN"} for b in hv1_seq]
        f_data = {"bases": mock_bases}
        f_trim = {"trim_start_idx": 0, "trim_end_idx": len(hv1_seq)}
        
        contig_pass = self.assembler.assemble_region(
            region="HV1",
            fwd_data=f_data,
            fwd_trim=f_trim,
            fwd_aln=aln,
            rev_data=None,
            rev_trim=None,
            rev_aln=None
        )
        self.assertEqual(contig_pass["qc_status"], "PASS")
        self.assertEqual(contig_pass["coverage_pct"], 100.0)

        # 3. Test case DENY do quá nhiều peak chồng chéo / lỗi sensor (?) >= 5
        mock_bases_bad = list(mock_bases)
        for i in range(6):
            mock_bases_bad[i] = {"called_base": "?", "primary_base": "N", "phred_score": 5, "is_conflict": True, "status": "CONFLICT_SENSOR"}
        f_bad_data = {"bases": mock_bases_bad}
        contig_deny_sensor = self.assembler.assemble_region(
            region="HV1",
            fwd_data=f_bad_data,
            fwd_trim=f_trim,
            fwd_aln=aln,
            rev_data=None,
            rev_trim=None,
            rev_aln=None
        )
        self.assertEqual(contig_deny_sensor["qc_status"], "DENY")
        self.assertIn("Peak chồng chéo/lỗi cảm biến quá nặng", contig_deny_sensor["qc_details"])

        print("[PASS] Test 6: Hệ thống phân loại QC chính xác; gán nhãn DENY và khuyến cáo làm lại mẫu đối với dữ liệu hỏng.")

if __name__ == "__main__":
    unittest.main()

