import sys
import types
import unittest


ete3_stub = types.ModuleType("ete3")


class StubNCBITaxa:

    def get_rank(self, tax_ids):
        return {tax_id: "species" for tax_id in tax_ids}


ete3_stub.NCBITaxa = StubNCBITaxa
sys.modules.setdefault("ete3", ete3_stub)

validate_email_stub = types.ModuleType("validate_email")
validate_email_stub.validate_email = lambda email: True
sys.modules.setdefault("validate_email", validate_email_stub)

from lverage.domain_scanner import DomainRecord, DomainScannerTemplate
from lverage.motif_database import MotifDBTemplate, MotifSearchRequest
from lverage.orf_searcher import OrfSearcherTemplate
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


class LverageValidationTests(unittest.TestCase):

    def build_lverage(self, valid_pfam_list=None, motif_database_list=None, ortholog_species_list=None):
        if motif_database_list is None:
            motif_database_list = [StubMotifDB()]

        return Lverage(
            motif_database_list=motif_database_list,
            orf_searcher=StubOrfSearcher(),
            domain_scanner=StubDomainScanner(),
            ortholog_species_list=ortholog_species_list,
            valid_pfam_list=valid_pfam_list,
            email="user@example.com"
        )

    def test_empty_valid_pfam_list_is_accepted(self):
        lverage = self.build_lverage([])

        self.assertEqual(lverage.valid_pfam_list, [])

    def test_valid_pfam_list_accepts_strings(self):
        lverage = self.build_lverage(["PF00046"])

        self.assertEqual(lverage.valid_pfam_list, ["PF00046"])

    def test_valid_pfam_list_rejects_non_string_values(self):
        with self.assertRaises(TypeError):
            self.build_lverage([46])

    def test_default_lists_are_unique_to_each_instance(self):
        first_lverage = self.build_lverage()
        second_lverage = self.build_lverage()

        self.assertIsNot(first_lverage.ortholog_species_list, second_lverage.ortholog_species_list)
        self.assertIsNot(first_lverage.valid_pfam_list, second_lverage.valid_pfam_list)

    def test_constructor_copies_configuration_lists(self):
        motif_database = StubMotifDB()
        motif_database_list = [motif_database]
        ortholog_species_list = [9606]
        valid_pfam_list = ["PF00046"]

        lverage = self.build_lverage(
            valid_pfam_list=valid_pfam_list,
            motif_database_list=motif_database_list,
            ortholog_species_list=ortholog_species_list
        )

        motif_database_list.clear()
        ortholog_species_list.clear()
        valid_pfam_list.clear()

        self.assertEqual(lverage.motif_database_list, [motif_database])
        self.assertEqual(lverage.ortholog_species_list, [9606])
        self.assertEqual(lverage.valid_pfam_list, ["PF00046"])

    def test_orf_search_matches_versionless_pfam_accession(self):
        lverage = Lverage.__new__(Lverage)
        lverage.tf_sequences = ["DNA"]
        lverage.orf_searcher = StubOrfSearcher()
        lverage.domain_scanner = StubDomainScanner()
        lverage.domain_scanner.get_domains = lambda sequence: [
            DomainRecord("Homeobox", "PF00046.1", 1, 5)
        ]
        lverage.valid_pfam_list = ["PF00046"]

        lverage._Lverage__search_orfs()

        self.assertEqual(lverage.orf, "DNA")
        self.assertEqual(lverage.valid_domains[0].accession, "PF00046.1")

    def test_orf_search_requires_exact_versioned_pfam_accession(self):
        lverage = Lverage.__new__(Lverage)
        lverage.tf_sequences = ["DNA"]
        lverage.orf_searcher = StubOrfSearcher()
        lverage.domain_scanner = StubDomainScanner()
        lverage.domain_scanner.get_domains = lambda sequence: [
            DomainRecord("Homeobox", "PF00046.1", 1, 5)
        ]
        lverage.valid_pfam_list = ["PF00046.2"]

        lverage._Lverage__search_orfs()

        self.assertIsNone(lverage.orf)


if __name__ == "__main__":
    unittest.main()
