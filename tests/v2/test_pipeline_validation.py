import unittest

from lverage.domain_scanner import DomainRecord, DomainScannerTemplate
from lverage.motif_database import MotifDBTemplate, MotifSearchRequest
from lverage.orf_searcher import OrfSearcherTemplate
from lverage.ortholog_searcher import OrthologSearcherTemplate
from lverage.pipeline import Lverage


class StubMotifDB(MotifDBTemplate):

    name = "Stub Motif Database"

    def search(self, request : MotifSearchRequest):
        return []

    def check_species_validity(self, species_tax_id : int):
        return True


class StubOrfSearcher(OrfSearcherTemplate):

    def get_orfs(self, sequence : str) -> list[str]:
        return [sequence]


class StubDomainScanner(DomainScannerTemplate):

    def get_domains(self, sequence : str):
        return []


class StubOrthologSearcher(OrthologSearcherTemplate):

    def get_orthologs(self, sequence : str):
        return []


class LverageValidationTests(unittest.TestCase):

    def build_lverage(self, valid_pfam_list=None, motif_database_list=None, **kwargs):
        if motif_database_list is None:
            motif_database_list = [StubMotifDB()]

        return Lverage(
            motif_database_list=motif_database_list,
            orf_searcher=kwargs.get("orf_searcher", StubOrfSearcher()),
            domain_scanner=kwargs.get("domain_scanner", StubDomainScanner()),
            ortholog_searcher=kwargs.get("ortholog_searcher", StubOrthologSearcher()),
            valid_pfam_list=valid_pfam_list,
            dbd_identity_thresh=kwargs.get("dbd_identity_thresh", 0.7),
        )

    def test_empty_valid_pfam_list_is_accepted(self):
        self.assertEqual(self.build_lverage([]).valid_pfam_list, [])

    def test_valid_pfam_list_rejects_non_string_values(self):
        with self.assertRaises(TypeError):
            self.build_lverage([46])

    def test_constructor_copies_configuration_lists(self):
        motif_database = StubMotifDB()
        motif_database_list = [motif_database]
        valid_pfam_list = ["PF00046"]
        lverage = self.build_lverage(valid_pfam_list, motif_database_list)
        motif_database_list.clear()
        valid_pfam_list.clear()

        self.assertEqual(lverage.motif_database_list, [motif_database])
        self.assertEqual(lverage.valid_pfam_list, ["PF00046"])

    def test_constructor_requires_an_ortholog_searcher(self):
        with self.assertRaises(TypeError):
            self.build_lverage(ortholog_searcher=object())

    def test_constructor_does_not_query_motif_databases(self):
        class FailingMotifDB(StubMotifDB):

            def check_species_validity(self, species_tax_id : int):
                raise AssertionError("constructor contacted a remote service")

        self.build_lverage(motif_database_list=[FailingMotifDB()])

    def test_identity_threshold_is_a_ratio(self):
        self.assertEqual(self.build_lverage(dbd_identity_thresh=0).dbd_identity_thresh, 0)
        self.assertEqual(self.build_lverage(dbd_identity_thresh=1).dbd_identity_thresh, 1)
        with self.assertRaises(ValueError):
            self.build_lverage(dbd_identity_thresh=1.01)

    def test_orf_search_matches_versionless_pfam_accession(self):
        lverage = self.build_lverage(["PF00046"])
        lverage.domain_scanner.get_domains = lambda sequence: [
            DomainRecord("Homeobox", "PF00046.1", 1, 5)
        ]

        orf, domains = lverage._select_query_orf(["DNA"])

        self.assertEqual(orf, "DNA")
        self.assertEqual(domains[0].accession, "PF00046.1")

    def test_orf_search_requires_exact_versioned_pfam_accession(self):
        lverage = self.build_lverage(["PF00046.2"])
        lverage.domain_scanner.get_domains = lambda sequence: [
            DomainRecord("Homeobox", "PF00046.1", 1, 5)
        ]

        self.assertEqual(lverage._select_query_orf(["DNA"]), (None, []))


if __name__ == "__main__":
    unittest.main()
