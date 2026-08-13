import unittest

from lverage.domain_scanner import DomainRecord


class DomainRecordTests(unittest.TestCase):

    def test_record_stores_domain_values(self):
        domain = DomainRecord("Homeobox", "PF00046.1", 2, 5)

        self.assertEqual(domain.name, "Homeobox")
        self.assertEqual(domain.accession, "PF00046.1")
        self.assertEqual(domain.start, 2)
        self.assertEqual(domain.end, 5)

    def test_coordinates_support_python_sequence_slicing(self):
        domain = DomainRecord("Homeobox", "PF00046.1", 2, 5)
        protein_sequence = "ABCDEFGHIJ"

        self.assertEqual(protein_sequence[domain.start:domain.end], "CDE")


if __name__ == "__main__":
    unittest.main()
