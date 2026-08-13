import subprocess
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import requests

from lverage.blast import NCBI_EFETCH_URL, RemoteBlastSearcher


def make_alignment(accession="NP_000001.1", species="Homo sapiens"):
    return SimpleNamespace(
        accession=accession,
        hit_def=f"Homeobox protein [{species}]",
        hit_id=f"ref|{accession}|",
        hsps=[SimpleNamespace(expect=1e-20, identities=80, align_length=100)],
    )


class RemoteBlastSearcherTests(unittest.TestCase):

    def setUp(self):
        self.which_patcher = patch("lverage.blast.shutil.which", return_value="/resolved/blastp")
        self.run_patcher = patch("lverage.blast.subprocess.run")
        self.get_patcher = patch("lverage.blast.requests.get")
        self.parse_patcher = patch("lverage.blast.NCBIXML.parse")
        self.seq_parse_patcher = patch("lverage.blast.SeqIO.parse")
        self.which_patcher.start()
        self.mock_run = self.run_patcher.start()
        self.mock_get = self.get_patcher.start()
        self.mock_parse = self.parse_patcher.start()
        self.mock_seq_parse = self.seq_parse_patcher.start()
        self.addCleanup(self.which_patcher.stop)
        self.addCleanup(self.run_patcher.stop)
        self.addCleanup(self.get_patcher.stop)
        self.addCleanup(self.parse_patcher.stop)
        self.addCleanup(self.seq_parse_patcher.stop)

        self.mock_run.return_value = SimpleNamespace(stdout="<xml>")
        self.mock_parse.return_value = iter([SimpleNamespace(alignments=[make_alignment()])])
        self.response = SimpleNamespace(text=">NP_000001.1 protein\nMPEPTIDE\n")
        self.response.raise_for_status = unittest.mock.Mock()
        self.mock_get.return_value = self.response
        self.mock_seq_parse.return_value = [SimpleNamespace(id="NP_000001.1", seq="MPEPTIDE")]

    def make_searcher(self, **kwargs):
        return RemoteBlastSearcher(
            {"Homo sapiens": 9606, "Mus musculus": 10090},
            "researcher@example.org",
            **kwargs,
        )

    def test_search_builds_remote_taxonomic_command_and_timeout(self):
        searcher = self.make_searcher(timeout_seconds=17)

        records = searcher.get_orthologs("A" * 200)

        command = self.mock_run.call_args.args[0]
        self.assertEqual(command, [
            "/resolved/blastp", "-remote", "-db", "nr", "-entrez_query",
            "(txid9606[ORGN] OR txid10090[ORGN])", "-outfmt", "5",
            "-evalue", "1e-06", "-max_target_seqs", "20",
        ])
        self.assertEqual(self.mock_run.call_args.kwargs["timeout"], 17)
        self.assertEqual(records[0].sequence, "MPEPTIDE")

    def test_complete_proteins_are_batch_fetched(self):
        searcher = self.make_searcher(request_timeout_seconds=9)

        searcher.get_orthologs("QUERY")

        self.mock_get.assert_called_once_with(
            NCBI_EFETCH_URL,
            params={
                "db": "protein",
                "id": "NP_000001.1",
                "rettype": "fasta",
                "retmode": "text",
                "email": "researcher@example.org",
                "tool": "Lverage",
            },
            timeout=9,
        )
        self.response.raise_for_status.assert_called_once_with()

    def test_missing_individual_protein_is_skipped(self):
        searcher = self.make_searcher()
        self.mock_seq_parse.return_value = []

        self.assertEqual(searcher.get_orthologs("QUERY"), [])

    def test_blast_timeout_propagates(self):
        searcher = self.make_searcher()
        self.mock_run.side_effect = subprocess.TimeoutExpired("blastp", 3600)

        with self.assertRaises(subprocess.TimeoutExpired):
            searcher.get_orthologs("QUERY")

    def test_http_error_propagates(self):
        searcher = self.make_searcher()
        self.response.raise_for_status.side_effect = requests.HTTPError("service unavailable")

        with self.assertRaises(requests.HTTPError):
            searcher.get_orthologs("QUERY")

    def test_blast_result_parse_failure_propagates(self):
        searcher = self.make_searcher()
        self.mock_parse.side_effect = ValueError("invalid XML")

        with self.assertRaises(ValueError):
            searcher.get_orthologs("QUERY")


if __name__ == "__main__":
    unittest.main()
