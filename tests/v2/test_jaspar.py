import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from lverage.domain_scanner import DomainRecord
from lverage.jaspar import Jaspar2026MotifDB
from lverage.motif_database import MotifSearchRequest


def make_request(sequence="ORTHOLOG", tax_id=9606, domain=None):
    domain = domain or DomainRecord("Homeobox", "PF00046.1", 0, len(sequence))
    return MotifSearchRequest(
        "QUERY", DomainRecord("Homeobox", "PF00046.1", 0, 5),
        sequence, domain, tax_id,
    )


class LocalJasparTests(unittest.TestCase):

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.database_path = self.root / "JASPAR2026.sqlite"
        with sqlite3.connect(self.database_path) as connection:
            connection.executescript("""
                CREATE TABLE MATRIX (ID INTEGER, BASE_ID TEXT, VERSION INTEGER, COLLECTION TEXT, NAME TEXT);
                CREATE TABLE MATRIX_SPECIES (ID INTEGER, TAX_ID TEXT);
                CREATE TABLE MATRIX_ANNOTATION (ID INTEGER, TAG TEXT, VAL TEXT);
                CREATE TABLE MATRIX_DATA (ID INTEGER, ROW TEXT, COL INTEGER, VAL REAL);
            """)
            for matrix_id, internal_id, name, species in (
                ("MA0001.1", 1, "First motif", [9606]),
                ("MA0001.2", 4, "First motif newer", [9606]),
                ("MA0002.3", 2, "Second motif", [10090, 9606]),
                ("MA0003.1", 3, "Mouse-only motif", [10090]),
            ):
                base_id, version = matrix_id.split(".")
                connection.execute(
                    "INSERT INTO MATRIX VALUES (?, ?, ?, 'CORE', ?)",
                    (internal_id, base_id, int(version), name),
                )
                connection.executemany(
                    "INSERT INTO MATRIX_SPECIES VALUES (?, ?)",
                    [(internal_id, str(tax_id)) for tax_id in species],
                )
                connection.execute(
                    "INSERT INTO MATRIX_ANNOTATION VALUES (?, 'class', 'Homeobox')",
                    (internal_id,),
                )
                for row_index, base in enumerate("ACGT"):
                    connection.executemany(
                        "INSERT INTO MATRIX_DATA VALUES (?, ?, ?, ?)",
                        [(internal_id, base, column, row_index + column)
                         for column in range(2)],
                    )
        self.adapter = Jaspar2026MotifDB(
            data_dir=self.root,
            sqlite_db=self.database_path,
            blastp_path="blastp-test",
            hmmscan_path="hmmscan-test",
            hmmalign_path="hmmalign-test",
        )
        self.adapter._validated = True

    def tearDown(self):
        self.temporary.cleanup()

    def test_local_species_filter_threshold_and_limit_preserve_metrics(self):
        self.adapter.n_hits = 1
        self.adapter.escore_threshold = 1e-6
        rows = [
            {"matrix_id": "MA0003.1", "tf_name": "Mouse", "evalue": 1e-12,
             "query_region": "1-30", "tf_region": "4-32", "dbd_identity": 0.99},
            {"matrix_id": "MA0002.3", "tf_name": "Second", "evalue": 1e-10,
             "query_region": "2-30", "tf_region": "3-31", "dbd_identity": 0.81},
            {"matrix_id": "MA0001.1", "tf_name": "First", "evalue": 1e-5,
             "query_region": "1-30", "tf_region": "1-30", "dbd_identity": 0.7},
        ]
        with patch("lverage.jaspar.load_cisbp_models", return_value={None: {}}), \
            patch("lverage.jaspar.load_jaspar_data", return_value={}), \
            patch("lverage.jaspar.infer_sequence", return_value=rows):
            records = self.adapter.search(make_request())

        self.assertEqual([record.matrix_id for record in records], ["MA0002.3"])
        self.assertEqual(records[0].inference_evalue, 1e-10)
        self.assertEqual(records[0].inference_dbd_identity, 0.81)
        self.assertEqual(records[0].motif_name, "Second motif")
        self.assertEqual(records[0].pfm["A"], [0.0, 1.0])

    def test_exact_matrix_version_is_required(self):
        with self.assertRaisesRegex(ValueError, "MA0001.3"):
            self.adapter._resolve_matrix({
                "matrix_id": "MA0001.3", "tf_name": "Wrong version",
                "evalue": 1e-10, "query_region": "1-10",
                "tf_region": "1-10", "dbd_identity": 0.8,
            })

    def test_search_retains_multiple_returned_versions_of_one_matrix(self):
        rows = [
            {"matrix_id": "MA0001.1", "tf_name": "First", "evalue": 1e-10,
             "query_region": "1-30", "tf_region": "1-30", "dbd_identity": 0.8},
            {"matrix_id": "MA0001.2", "tf_name": "First newer", "evalue": 2e-10,
             "query_region": "1-30", "tf_region": "1-30", "dbd_identity": 0.85},
        ]
        with patch("lverage.jaspar.load_cisbp_models", return_value={None: {}}), \
            patch("lverage.jaspar.load_jaspar_data", return_value={}), \
            patch("lverage.jaspar.infer_sequence", return_value=rows):
            records = self.adapter.search(make_request())

        self.assertEqual(
            [record.matrix_id for record in records], ["MA0001.1", "MA0001.2"]
        )

    def test_species_check_uses_exact_local_taxonomy_ids(self):
        self.adapter._validated = False
        with patch.object(self.adapter, "_validate_configuration"):
            self.assertTrue(self.adapter.check_species_validity(9606))
            self.assertTrue(self.adapter.check_species_validity(10090))
            self.assertFalse(self.adapter.check_species_validity(7227))

    def test_search_calls_in_package_inference_with_window_and_configured_executables(self):
        sequence = "A" * 2000
        prediction = {
            "matrix_id": "MA0001.1", "tf_name": "Example TF", "evalue": 1e-12,
            "query_region": "1-40", "tf_region": "2-41", "dbd_identity": 0.875,
        }
        with patch("lverage.jaspar.load_cisbp_models", return_value={None: {}}), \
                patch("lverage.jaspar.load_jaspar_data", return_value={}), \
                patch("lverage.jaspar.infer_sequence", return_value=[prediction]) as infer:
            records = self.adapter.search(make_request(sequence=sequence))

        self.assertEqual([record.matrix_id for record in records], ["MA0001.1"])
        self.assertEqual(infer.call_args.args[0], sequence)
        self.assertEqual(infer.call_args.kwargs["query_id"], "lverage_query")
        self.assertEqual(infer.call_args.kwargs["blastp_path"], "blastp-test")
        self.assertEqual(infer.call_args.kwargs["hmmscan_path"], "hmmscan-test")
        self.assertEqual(infer.call_args.kwargs["hmmalign_path"], "hmmalign-test")
        self.assertEqual(infer.call_args.kwargs["data_dir"], self.root)
        self.assertFalse(Path(infer.call_args.kwargs["work_dir"]).exists())

    def test_existing_online_behavior_window_is_unchanged(self):
        sequence = "A" * 5000
        domain = DomainRecord("Homeobox", "PF00046.1", 2400, 2450)
        self.assertEqual(
            self.adapter._window_ortholog_sequence(sequence, domain),
            sequence[1425:3425],
        )


if __name__ == "__main__":
    unittest.main()