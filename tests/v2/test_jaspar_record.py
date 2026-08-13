import unittest

from lverage.jaspar import JasparRecord


class JasparRecordTests(unittest.TestCase):

    def test_record_serializes_in_header_order(self):
        pfm = {"A": [1, 2], "C": [3, 4], "G": [4, 3], "T": [2, 1]}
        record = JasparRecord(
            "MA0001.1",
            "Example motif",
            pfm,
            "https://jaspar.elixir.no/matrix/MA0001.1/",
            "Homeobox",
            1e-8,
        )

        self.assertEqual(record.get_headers(), [
            "Matrix ID", "Motif Name", "PFM", "Motif URL",
            "Motif Class", "Inference E-value",
        ])
        self.assertEqual(record.get_values(), [
            "MA0001.1",
            "Example motif",
            pfm,
            "https://jaspar.elixir.no/matrix/MA0001.1/",
            "Homeobox",
            1e-8,
        ])

    def test_headers_are_immutable_and_returned_as_a_list_copy(self):
        record = JasparRecord("MA0001.1", "Motif", {}, "URL", "Class", 1e-8)
        headers = record.get_headers()
        headers.append("Extra")

        self.assertIsInstance(JasparRecord.headers, tuple)
        self.assertNotIn("Extra", record.get_headers())


if __name__ == "__main__":
    unittest.main()
