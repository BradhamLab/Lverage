import unittest

from Bio.Seq import Seq

from lverage.domain_scanner import DomainScannerTemplate
from lverage.motif_database import MotifDBTemplate
from lverage.orf_searcher import OrfSearcherTemplate, OrffinderOrfSearcher
from lverage.ortholog_searcher import OrthologSearcherTemplate
from lverage.pipeline import Lverage, LverageCode


class IncompleteOrfSearcher(OrfSearcherTemplate):
    pass


class ConcreteOrfSearcher(OrfSearcherTemplate):

    def get_orfs(self, sequence : str) -> list[str]:
        return [sequence]


class StubMotifDatabase(MotifDBTemplate):

    name = "Stub Motif Database"

    def search(self, request):
        return []

    def check_species_validity(self, species_tax_id: int):
        return True


class RecordingDomainScanner(DomainScannerTemplate):

    def __init__(self):
        self.scanned = []

    def get_domains(self, sequence: str):
        self.scanned.append(sequence)
        return []


class StubOrthologSearcher(OrthologSearcherTemplate):

    def get_orthologs(self, sequence: str):
        return []


class OrfSearcherTemplateTests(unittest.TestCase):

    def test_template_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            OrfSearcherTemplate()

    def test_incomplete_subclass_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            IncompleteOrfSearcher()

    def test_concrete_subclass_returns_orfs(self):
        searcher = ConcreteOrfSearcher()

        self.assertEqual(searcher.get_orfs("ATGAAATAG"), ["ATGAAATAG"])


class OrffinderOrfSearcherTests(unittest.TestCase):

    def test_finds_forward_and_reverse_strand_proteins(self):
        searcher = OrffinderOrfSearcher(minimum_length=9)

        self.assertEqual(searcher.get_orfs("ATGAAATAG"), ["MK"])
        self.assertEqual(searcher.get_orfs("CTATTTCAT"), ["MK"])

    def test_minimum_length_is_inclusive_and_measured_in_nucleotides(self):
        self.assertEqual(OrffinderOrfSearcher(minimum_length=9).get_orfs("ATGAAATAG"), ["MK"])
        self.assertEqual(OrffinderOrfSearcher(minimum_length=10).get_orfs("ATGAAATAG"), [])

    def test_incomplete_terminal_orfs_follow_trimming_and_length_settings(self):
        self.assertEqual(OrffinderOrfSearcher(minimum_length=7).get_orfs("ATGAAATG"), ["MK"])
        self.assertEqual(OrffinderOrfSearcher(minimum_length=8).get_orfs("ATGAAATG"), ["MK"])
        self.assertEqual(OrffinderOrfSearcher(minimum_length=9).get_orfs("ATGAAATG"), [])
        self.assertEqual(
            OrffinderOrfSearcher(minimum_length=7, trim_trailing=True).get_orfs("ATGAAATG"),
            [],
        )

    def test_incomplete_orfs_preserve_full_codons_for_each_suffix_remainder(self):
        cases = (
            ("ATGAAAGGG", 0, "MKG"),
            ("ATGAAAGGGA", 1, "MKG"),
            ("ATGAAAGGGAA", 2, "MKG"),
            ("ATGAAAGGGAAA", 0, "MKGK"),
        )
        for nucleotide_sequence, expected_remainder, expected_protein in cases:
            with self.subTest(nucleotide_sequence=nucleotide_sequence):
                self.assertEqual(len(nucleotide_sequence) % 3, expected_remainder)
                self.assertEqual(
                    OrffinderOrfSearcher(minimum_length=9).get_orfs(nucleotide_sequence),
                    [expected_protein],
                )
                self.assertEqual(
                    OrffinderOrfSearcher(minimum_length=9).get_orfs(
                        str(Seq(nucleotide_sequence).reverse_complement())
                    ),
                    [expected_protein],
                )

    def test_incomplete_orf_minimum_length_uses_actual_nucleotide_span(self):
        self.assertEqual(OrffinderOrfSearcher(minimum_length=9).get_orfs("ATGAAAGGG"), ["MKG"])
        self.assertEqual(OrffinderOrfSearcher(minimum_length=10).get_orfs("ATGAAAGGG"), [])
        self.assertEqual(OrffinderOrfSearcher(minimum_length=10).get_orfs("ATGAAAGGGA"), ["MKG"])
        self.assertEqual(OrffinderOrfSearcher(minimum_length=11).get_orfs("ATGAAAGGGA"), [])

    def test_start_codon_alone_produces_one_residue_not_an_empty_candidate(self):
        self.assertEqual(OrffinderOrfSearcher(minimum_length=1).get_orfs("ATG"), ["M"])

    def test_removes_terminal_stop_but_preserves_internal_translation(self):
        searcher = OrffinderOrfSearcher(minimum_length=9)

        self.assertEqual(searcher.get_orfs("ATGAAATAG"), ["MK"])
        self.assertEqual(searcher.get_orfs("ATGNNNTAA"), ["MX"])

    def test_uses_the_installed_standard_translation_table(self):
        self.assertEqual(OrffinderOrfSearcher(minimum_length=6).get_orfs("ATGTGATAG"), ["M"])

    def test_returns_empty_list_when_no_orf_is_found(self):
        searcher = OrffinderOrfSearcher(minimum_length=9)

        self.assertEqual(searcher.get_orfs("CCCCCCCCCC"), [])

    def test_remove_nested_false_disables_the_nested_candidate_filter(self):
        sequence = "AATATCAATGGGTCAGCGGTACGCCGGGGGTACCGTAAGGCACATAGTATTAGAATGATA"

        unfiltered = OrffinderOrfSearcher(minimum_length=9).get_orfs(sequence)
        filtered = OrffinderOrfSearcher(minimum_length=9, remove_nested=True).get_orfs(sequence)

        self.assertEqual(unfiltered, ["MGQRYAGGTVRHIVLE", "MCLTVPPAYR"])
        self.assertEqual(filtered, ["MGQRYAGGTVRHIVLE"])

    def test_rejects_invalid_sequences(self):
        searcher = OrffinderOrfSearcher(minimum_length=9)

        for invalid in ("", " ", "ATGUAA", "ATG-AA"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    searcher.get_orfs(invalid)
        for invalid in (None, 123):
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    searcher.get_orfs(invalid)

    def test_validates_configuration(self):
        with self.assertRaises(TypeError):
            OrffinderOrfSearcher(minimum_length=True)
        with self.assertRaises(ValueError):
            OrffinderOrfSearcher(minimum_length=0)
        with self.assertRaises(ValueError):
            OrffinderOrfSearcher(start_codons=["GTG"])
        with self.assertRaises(ValueError):
            OrffinderOrfSearcher(start_codons=[])
        with self.assertRaises(TypeError):
            OrffinderOrfSearcher(start_codons="ATG")
        with self.assertRaises(TypeError):
            OrffinderOrfSearcher(remove_nested=1)

    def test_adapter_passes_protein_candidates_through_pipeline(self):
        domain_scanner = RecordingDomainScanner()
        pipeline = Lverage(
            [StubMotifDatabase()],
            OrffinderOrfSearcher(minimum_length=9),
            domain_scanner,
            StubOrthologSearcher(),
        )

        self.assertEqual(pipeline.run("ATGAAATAG"), [])
        self.assertEqual(domain_scanner.scanned, ["MK"])
        self.assertEqual(pipeline.lverage_code, LverageCode.NO_VALID_DOMAIN)


if __name__ == "__main__":
    unittest.main()
