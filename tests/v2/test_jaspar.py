import unittest
from unittest.mock import Mock, call, patch

import requests

from lverage.domain_scanner import DomainRecord
from lverage.jaspar import Jaspar2024MotifDB
from lverage.motif_database import MotifSearchRequest


def make_response(payload):
    response = Mock()
    response.json.return_value = payload
    return response


def make_request(sequence="ORTHOLOG", tax_id=9606):
    return MotifSearchRequest(
        "QUERY",
        DomainRecord("Homeobox", "PF00046.1", 0, 5),
        sequence,
        DomainRecord("Homeobox", "PF00046.2", 0, 5),
        tax_id,
    )


def make_motif(matrix_id, species_ids, motif_class=None):
    if motif_class is None:
        motif_class = ["Homeobox"]
    return {
        "matrix_id": matrix_id,
        "name": f"Motif {matrix_id}",
        "pfm": {"A": [1], "C": [2], "G": [3], "T": [4]},
        "class": motif_class,
        "species": [{"tax_id": species_id} for species_id in species_ids],
    }


class Jaspar2024MotifDBTests(unittest.TestCase):

    def test_import_and_constructor_do_not_make_requests(self):
        with patch("lverage.jaspar.requests.get") as mock_get:
            Jaspar2024MotifDB()

        mock_get.assert_not_called()

    @patch("lverage.jaspar.requests.get")
    def test_search_uses_windowed_ortholog_sequence_and_timeout(self, mock_get):
        mock_get.return_value = make_response({"results": []})
        sequence = "A" * 5000
        request = MotifSearchRequest(
            "QUERY",
            DomainRecord("Homeobox", "PF00046.1", 0, 5),
            sequence,
            DomainRecord("Homeobox", "PF00046.2", 2400, 2450),
            9606,
        )
        database = Jaspar2024MotifDB(request_timeout_seconds=7)
        expected_window = database._window_ortholog_sequence(sequence, request.ortholog_domain)

        self.assertEqual(database.search(request), [])
        mock_get.assert_called_once_with(
            f"https://jaspar.elixir.no/api/v1/infer/{expected_window}/",
            timeout=7,
        )
        mock_get.return_value.raise_for_status.assert_called_once_with()

    @patch("lverage.jaspar.requests.get")
    def test_results_are_sorted_filtered_and_limited_after_acceptance(self, mock_get):
        infer_response = make_response({"results": [
            {"evalue": 1e-4, "url": "https://jaspar.elixir.no/api/v1/matrix/HIGH/"},
            {"evalue": 1e-10, "url": "https://jaspar.elixir.no/api/v1/matrix/WRONG/"},
            {"evalue": 1e-9, "url": "https://jaspar.elixir.no/api/v1/matrix/FIRST/"},
            {"evalue": 1e-8, "url": "https://jaspar.elixir.no/api/v1/matrix/SECOND/"},
            {"evalue": 1e-7, "url": "https://jaspar.elixir.no/api/v1/matrix/THIRD/"},
        ]})
        mock_get.side_effect = [
            infer_response,
            make_response(make_motif("WRONG", [10090])),
            make_response(make_motif("FIRST", [9606])),
            make_response(make_motif("SECOND", [10090, 9606])),
        ]
        database = Jaspar2024MotifDB(n_hits=2, escore_threshold=1e-6)

        records = database.search(make_request())

        self.assertEqual([record.matrix_id for record in records], ["FIRST", "SECOND"])
        self.assertEqual([record.inference_evalue for record in records], [1e-9, 1e-8])
        self.assertEqual(mock_get.call_count, 4)

    @patch("lverage.jaspar.requests.get")
    def test_record_parses_current_motif_fields_and_url(self, mock_get):
        mock_get.side_effect = [
            make_response({"results": [{"evalue": "1e-12", "url": "DETAIL"}]}),
            make_response(make_motif("MA0001.1", [9606], "Homeobox")),
        ]

        record = Jaspar2024MotifDB().search(make_request())[0]

        self.assertEqual(record.matrix_id, "MA0001.1")
        self.assertEqual(record.motif_name, "Motif MA0001.1")
        self.assertEqual(record.motif_class, "Homeobox")
        self.assertEqual(record.motif_url, "https://jaspar.elixir.no/matrix/MA0001.1/")
        self.assertEqual(record.inference_evalue, 1e-12)

    @patch("lverage.jaspar.requests.get")
    def test_species_pagination_is_loaded_once_and_cached(self, mock_get):
        mock_get.side_effect = [
            make_response({
                "results": [{"tax_id": 9606}],
                "next": "https://jaspar.elixir.no/api/v1/species/?page=2",
            }),
            make_response({"results": [{"tax_id": "10090"}], "next": None}),
        ]
        database = Jaspar2024MotifDB(request_timeout_seconds=11)

        self.assertTrue(database.check_species_validity(9606))
        self.assertTrue(database.check_species_validity(10090))
        self.assertFalse(database.check_species_validity(7227))
        self.assertEqual(mock_get.call_args_list, [
            call(
                database.jaspar_rest_species_url,
                params=database.species_params,
                timeout=11,
            ),
            call(
                "https://jaspar.elixir.no/api/v1/species/?page=2",
                params=None,
                timeout=11,
            ),
        ])

    @patch("lverage.jaspar.requests.get")
    def test_http_errors_propagate(self, mock_get):
        response = make_response({"results": []})
        response.raise_for_status.side_effect = requests.HTTPError("service failure")
        mock_get.return_value = response

        with self.assertRaises(requests.HTTPError):
            Jaspar2024MotifDB().search(make_request())

    @patch("lverage.jaspar.requests.get")
    def test_request_timeouts_propagate(self, mock_get):
        mock_get.side_effect = requests.Timeout("service timeout")

        with self.assertRaises(requests.Timeout):
            Jaspar2024MotifDB().search(make_request())

    @patch("lverage.jaspar.requests.get")
    def test_malformed_responses_raise_value_error(self, mock_get):
        mock_get.return_value = make_response({"unexpected": []})

        with self.assertRaises(ValueError):
            Jaspar2024MotifDB().search(make_request())

    def test_search_validates_request_contract(self):
        with self.assertRaises(TypeError):
            Jaspar2024MotifDB().search(object())


if __name__ == "__main__":
    unittest.main()
