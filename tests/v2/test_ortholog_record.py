import unittest

from lverage.ortholog_searcher import OrthologRecord


class OrthologRecordTests(unittest.TestCase):

    def test_record_stores_ortholog_evidence(self):
        record = OrthologRecord(
            "NP_000001.1", "Protein", "Homo sapiens", 9606,
            "MPEPTIDE", 1e-20, 0.75, 0.80
        )

        self.assertEqual(record.accession, "NP_000001.1")
        self.assertEqual(record.species_tax_id, 9606)
        self.assertEqual(record.sequence, "MPEPTIDE")
        self.assertEqual(record.identity, 0.75)
        self.assertEqual(record.query_coverage, 0.80)


if __name__ == "__main__":
    unittest.main()
