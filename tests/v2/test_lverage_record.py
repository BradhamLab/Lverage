import unittest

from lverage.domain_scanner import DomainRecord
from lverage.motif_database import MotifDBRecordTemplate
from lverage.ortholog_searcher import OrthologRecord
from lverage.records import LverageRecord


class ExampleMotifRecord(MotifDBRecordTemplate):

    headers = ("Matrix ID", "Motif Name")

    def __init__(self, matrix_id, motif_name):
        self.matrix_id = matrix_id
        self.motif_name = motif_name

    def get_values(self):
        return [self.matrix_id, self.motif_name]


class LverageRecordTests(unittest.TestCase):

    def make_record(self):
        return LverageRecord(
            query_domain=DomainRecord("Homeobox", "PF00046.32", 3, 9),
            ortholog=OrthologRecord(
                "NP_000001.1",
                "Homeobox protein [Homo sapiens]",
                "Homo sapiens",
                9606,
                "MPEPTIDE",
                1e-30,
                0.9,
                0.8,
            ),
            ortholog_domain=DomainRecord("Homeobox", "PF00046.30", 1, 7),
            domain_identity=0.75,
            motif_database_name="ExampleDB",
            motif_record=ExampleMotifRecord("MA0001.1", "Example motif"),
        )

    def test_headers_flatten_core_and_motif_schema(self):
        record = self.make_record()

        self.assertEqual(record.get_headers()[-3:], ["Motif Database", "Matrix ID", "Motif Name"])

    def test_values_flatten_nested_evidence(self):
        record = self.make_record()
        values = record.get_values()

        self.assertEqual(values[:4], ["Homeobox", "PF00046.32", 3, 9])
        self.assertEqual(values[4], "NP_000001.1")
        self.assertEqual(values[12:18], ["Homeobox", "PF00046.30", 1, 7, 0.75, "ExampleDB"])
        self.assertEqual(values[-2:], ["MA0001.1", "Example motif"])

    def test_headers_are_returned_as_a_new_list(self):
        record = self.make_record()
        headers = record.get_headers()
        headers.append("Extra")

        self.assertNotIn("Extra", record.get_headers())


if __name__ == "__main__":
    unittest.main()
