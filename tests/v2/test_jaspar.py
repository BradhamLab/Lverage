import sys
import types
import unittest
from unittest.mock import Mock, patch


requests_stub = types.ModuleType("requests")
requests_stub.get = Mock()
sys.modules.setdefault("requests", requests_stub)

from lverage.jaspar import Jaspar2024MotifDB


class Jaspar2024MotifDBTests(unittest.TestCase):

    def test_import_does_not_request_species(self):
        requests_stub.get.assert_not_called()

    @patch("lverage.jaspar.requests.get")
    def test_species_are_loaded_once_and_cached(self, mock_get):
        response = Mock()
        response.json.return_value = {
            "results": [
                {"tax_id": 9606}
            ]
        }
        mock_get.return_value = response
        motif_database = Jaspar2024MotifDB()

        self.assertTrue(motif_database.check_species_validity(9606))
        self.assertFalse(motif_database.check_species_validity(10090))

        mock_get.assert_called_once_with(
            motif_database.jaspar_rest_species_url,
            params=motif_database.species_params
        )


if __name__ == "__main__":
    unittest.main()
