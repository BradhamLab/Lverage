import unittest
from unittest.mock import Mock, patch

from lverage.jaspar import Jaspar2024MotifDB
from lverage.domain_scanner import DomainRecord
from lverage.motif_database import MotifSearchRequest


class Jaspar2024MotifDBTests(unittest.TestCase):

    def test_search_uses_motif_request_contract(self):
        request = MotifSearchRequest(
            "QUERY", DomainRecord("Homeobox", "PF00046.1", 0, 5),
            "ORTHOLOG", DomainRecord("Homeobox", "PF00046.1", 0, 5), 9606
        )

        with self.assertRaises(NotImplementedError):
            Jaspar2024MotifDB().search(request)

    def test_import_does_not_request_species(self):
        with patch("lverage.jaspar.requests.get") as mock_get:
            Jaspar2024MotifDB()

        mock_get.assert_not_called()

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
