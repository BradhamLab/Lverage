import tempfile
import unittest
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from benchmarks.orf_finder import benchmark


class BenchmarkCoordinateTests(unittest.TestCase):

    def test_locus_coordinates_convert_both_strands_to_half_open(self):
        forward = {"start": 1, "end": 34, "sense": "+", "length": 33, "trailing": False}
        reverse = {"start": 34, "end": 1, "sense": "-", "length": 33, "trailing": False}

        self.assertEqual(benchmark._locus_interval(forward, 33), (0, 33, "+", 33))
        self.assertEqual(benchmark._locus_interval(reverse, 33), (0, 33, "-", 33))

    def test_interval_iou_uses_half_open_bounds(self):
        self.assertEqual(benchmark.interval_iou((0, 10), (5, 15)), 1 / 3)
        self.assertEqual(benchmark.interval_iou((0, 10), (10, 20)), 0.0)

    def test_incomplete_terminal_locus_uses_actual_span(self):
        locus = {"start": 1, "end": 75, "sense": "+", "length": 74, "trailing": True}

        self.assertEqual(benchmark._locus_interval(locus, 75), (0, 75, "+", 75))

    def test_reverse_incomplete_terminal_locus_uses_original_input_bounds(self):
        locus = {"start": 76, "end": 1, "sense": "-", "length": 74, "trailing": True}

        self.assertEqual(benchmark._locus_interval(locus, 75), (0, 75, "-", 75))


class BenchmarkPairingTests(unittest.TestCase):

    @staticmethod
    def candidate(candidate_id, start, end, strand="+"):
        return benchmark.Candidate(
            candidate_id, "test", "TX.1", start, end, strand, end - start,
            "MTEST", False, False, 0.0, False,
        )

    def test_all_positive_overlaps_are_retained_and_exact_matches_separate(self):
        left = [self.candidate("a", 0, 10), self.candidate("b", 5, 15)]
        right = [self.candidate("x", 0, 10), self.candidate("y", 8, 18), self.candidate("z", 1, 4, "-")]

        result = benchmark.compare_candidates(left, right)

        self.assertEqual([(a.candidate_id, b.candidate_id) for a, b, _ in result["exact"]], [("a", "x")])
        self.assertEqual(
            {(a.candidate_id, b.candidate_id) for a, b, _ in result["overlaps"]},
            {("a", "x"), ("a", "y"), ("b", "x"), ("b", "y")},
        )
        self.assertEqual(sum(iou == 1.0 for _, _, iou in result["overlaps"]), 1)
        self.assertEqual(result["unmatched_left"], [])
        self.assertEqual(result["unmatched_right"], ["z"])

    def test_threshold_pairing_is_deterministic_and_one_to_one(self):
        left = [self.candidate("left-b", 0, 10), self.candidate("left-a", 0, 10)]
        right = [self.candidate("right-b", 0, 10), self.candidate("right-a", 0, 10)]

        result = benchmark.compare_candidates(left, right)

        self.assertEqual(
            [(a.candidate_id, b.candidate_id) for a, b, _ in result["exact"]],
            [("left-a", "right-a"), ("left-b", "right-b")],
        )
        self.assertEqual(result["threshold_pairs"], [])

    def test_greedy_iou_tie_breaks_by_stable_candidate_ids(self):
        left = [self.candidate("A", 0, 10), self.candidate("B", 4, 14)]
        right = [self.candidate("X", 2, 12), self.candidate("Y", 2, 12)]

        result = benchmark.compare_candidates(left, right)

        self.assertEqual(
            [(a.candidate_id, b.candidate_id) for a, b, _ in result["threshold_pairs"]],
            [("A", "X"), ("B", "Y")],
        )


class BenchmarkNcbiParserTests(unittest.TestCase):

    def parse(self, record_id, interval_header, orf_tag, cds_sequence, protein, input_sequence):
        record = SeqRecord(Seq(input_sequence), id=record_id)
        record.features = [SeqFeature(FeatureLocation(0, len(input_sequence), strand=1), type="CDS", qualifiers={"gene": ["TEST"], "translation": ["M"]})]
        with tempfile.TemporaryDirectory() as temp_dir:
            cds_path = Path(temp_dir) / "cds.fasta"
            protein_path = Path(temp_dir) / "protein.fasta"
            cds_path.write_text(f">lcl|{record_id}:{interval_header} {orf_tag}\n{cds_sequence}\n")
            protein_path.write_text(f">lcl|{orf_tag}\n{protein}\n")
            return benchmark.parse_ncbi_outputs(cds_path, protein_path, {record_id: record})

    def test_plus_strand_location_stop_inclusion_and_protein_crosscheck(self):
        cds = "ATG" + "GCT" * 9 + "TAA"
        protein = "M" + "A" * 9
        candidates = self.parse("TX.1", "1-33", "ORF1_TX.1:0:32", cds, protein, cds)

        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0].interval, (0, 33, "+"))
        self.assertEqual(candidates[0].nucleotide_length, 33)
        self.assertTrue(candidates[0].terminal_stop_in_nt)
        self.assertEqual(candidates[0].protein, protein)

    def test_minus_strand_location_normalizes_on_original_input(self):
        oriented_cds = "ATG" + "GCT" * 9 + "TAA"
        input_sequence = str(Seq(oriented_cds).reverse_complement())
        protein = "M" + "A" * 9
        candidates = self.parse("TX.1", "c33-1", "ORF1_TX.1:32:0", oriented_cds, protein, input_sequence)

        self.assertEqual(candidates[0].interval, (0, 33, "-"))

    def test_mismatched_ncbi_cds_and_protein_is_rejected(self):
        cds = "ATG" + "GCT" * 9 + "TAA"
        with self.assertRaisesRegex(ValueError, "does not match translated CDS"):
            self.parse("TX.1", "1-33", "ORF1_TX.1:0:32", cds, "M", cds)


class BenchmarkAdapterAssociationTests(unittest.TestCase):

    def setUp(self):
        self.record_id = "SYNTHETIC.1"
        self.sequence = "ATG" + "GCT" * 24
        record = SeqRecord(Seq(self.sequence), id=self.record_id)
        record.features = [SeqFeature(FeatureLocation(0, len(self.sequence), strand=1), type="CDS", qualifiers={"gene": ["SYN"], "translation": ["M" + "A" * 24]})]
        self.old_records = benchmark._RECORDS_BY_ID
        benchmark._RECORDS_BY_ID = {self.record_id: record}

    def tearDown(self):
        benchmark._RECORDS_BY_ID = self.old_records

    def test_adapter_corrects_package_one_nt_terminal_threshold_behavior(self):
        corrected = benchmark.corrected_adapter_candidates(self.record_id, self.sequence)
        raw = benchmark.python_raw_candidates(self.record_id, self.sequence)

        self.assertEqual(len(corrected), 1)
        self.assertEqual(corrected[0].interval, (0, 75, "+"))
        self.assertEqual(corrected[0].protein, "M" + "A" * 24)
        self.assertEqual(raw, [])

    def test_corrected_candidate_proteins_match_production_adapter_with_duplicates(self):
        orf = "ATG" + "GCT" * 23 + "TAA"
        sequence = orf + orf
        record = SeqRecord(Seq(sequence), id=self.record_id)
        record.features = [SeqFeature(FeatureLocation(0, len(sequence), strand=1), type="CDS", qualifiers={"gene": ["SYN"], "translation": ["M" + "A" * 23]})]
        benchmark._RECORDS_BY_ID[self.record_id] = record

        candidates = benchmark.corrected_adapter_candidates(self.record_id, sequence)

        proteins = [candidate.protein for candidate in candidates]
        self.assertGreaterEqual(Counter(proteins).most_common(1)[0][1], 2)


if __name__ == "__main__":
    unittest.main()
