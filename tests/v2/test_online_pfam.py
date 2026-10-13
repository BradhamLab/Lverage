import unittest
from types import SimpleNamespace
from unittest.mock import patch

import requests

from lverage.domain_scanner import DomainRecord
from lverage.pfam import OnlinePfamScanner


class OnlinePfamScannerTests(unittest.TestCase):

    def setUp(self):
        self.post_patcher = patch("lverage.pfam.requests.post")
        self.get_patcher = patch("lverage.pfam.requests.get")
        self.sleep_patcher = patch("lverage.pfam.time.sleep")
        self.mock_post = self.post_patcher.start()
        self.mock_get = self.get_patcher.start()
        self.sleep_patcher.start()
        self.addCleanup(self.post_patcher.stop)
        self.addCleanup(self.get_patcher.stop)
        self.addCleanup(self.sleep_patcher.stop)

    def make_response(self, text="", json_data=None):
        response = SimpleNamespace(text=text, json=lambda: json_data)
        response.raise_for_status = unittest.mock.Mock()
        return response

    def make_scanner(self, **kwargs):
        return OnlinePfamScanner("researcher@example.org", **kwargs)

    def test_scan_submits_polls_and_parses_domains(self):
        self.mock_post.return_value = self.make_response("job-123\n")
        self.mock_get.side_effect = [
            self.make_response("RUNNING"),
            self.make_response("FINISHED"),
            self.make_response(json_data=[
                {"identifier": "txt", "mediaType": "text/plain"},
                {"identifier": "json", "mediaType": "application/json"},
            ]),
            self.make_response(json_data=[{
                "name": "Homeobox",
                "acc": "PF00046.1",
                "env": {"from": 2, "to": 7},
            }]),
        ]

        domains = self.make_scanner(time_interval=2).get_domains("MPEPTIDE")

        self.assertEqual(domains, [DomainRecord("Homeobox", "PF00046.1", 1, 7)])
        self.mock_post.assert_called_once_with(
            "https://www.ebi.ac.uk/Tools/services/rest/pfamscan/run/",
            data={
                "email": "researcher@example.org",
                "sequence": "MPEPTIDE",
                "format": "json",
            },
            timeout=30,
        )
        self.assertEqual(self.mock_get.call_count, 4)
        self.assertEqual(self.mock_get.call_args_list[0].args[0].rsplit("/", 1)[-1], "job-123")
        self.assertEqual(self.mock_get.call_args_list[2].args[0].rsplit("/", 1)[-1], "job-123")
        self.assertEqual(self.mock_get.call_args_list[3].args[0].rsplit("/", 1)[-1], "json")

    def test_failed_job_raises(self):
        for status in ("ERROR", "FAILURE", "NOT_FOUND"):
            with self.subTest(status=status):
                self.mock_post.return_value = self.make_response("job-123")
                self.mock_get.return_value = self.make_response(status)

                with self.assertRaises(RuntimeError):
                    self.make_scanner().get_domains("QUERY")

    def test_job_timeout_raises(self):
        self.mock_post.return_value = self.make_response("job-123")
        self.mock_get.return_value = self.make_response("RUNNING")

        with self.assertRaises(TimeoutError):
            self.make_scanner(try_count=2).get_domains("QUERY")

    def test_unknown_job_status_raises(self):
        self.mock_post.return_value = self.make_response("job-123")
        self.mock_get.return_value = self.make_response("PAUSED")

        with self.assertRaises(RuntimeError):
            self.make_scanner().get_domains("QUERY")

    def test_http_error_propagates(self):
        self.mock_post.return_value = self.make_response("job-123")
        self.mock_post.return_value.raise_for_status.side_effect = requests.HTTPError("unavailable")

        with self.assertRaises(requests.HTTPError):
            self.make_scanner().get_domains("QUERY")

    def test_empty_sequence_is_rejected(self):
        with self.assertRaises(ValueError):
            self.make_scanner().get_domains(" ")


if __name__ == "__main__":
    unittest.main()