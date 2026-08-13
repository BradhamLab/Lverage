import unittest

from lverage.domain_scanner import DomainRecord
from lverage.motif_database import MotifSearchRequest


class MotifSearchRequestTests(unittest.TestCase):

    def test_request_stores_motif_search_evidence(self):
        query_domain = DomainRecord("Homeobox", "PF00046.1", 10, 70)
        ortholog_domain = DomainRecord("Homeobox", "PF00046.1", 20, 80)
        request = MotifSearchRequest(
            query_sequence="QUERYPROTEIN",
            query_domain=query_domain,
            ortholog_sequence="ORTHOLOGPROTEIN",
            ortholog_domain=ortholog_domain,
            ortholog_species_tax_id=9606
        )

        self.assertEqual(request.query_sequence, "QUERYPROTEIN")
        self.assertEqual(request.query_domain, query_domain)
        self.assertEqual(request.ortholog_sequence, "ORTHOLOGPROTEIN")
        self.assertEqual(request.ortholog_domain, ortholog_domain)
        self.assertEqual(request.ortholog_species_tax_id, 9606)


if __name__ == "__main__":
    unittest.main()
