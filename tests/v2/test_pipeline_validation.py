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

from lverage.domain_scanner import DomainScannerTemplate
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

    def build_lverage(self, valid_pfam_list):
        return Lverage(
            motif_database_list=[StubMotifDB()],
            orf_searcher=StubOrfSearcher(),
            domain_scanner=StubDomainScanner(),
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


if __name__ == "__main__":
    unittest.main()
