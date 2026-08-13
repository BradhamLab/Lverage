import subprocess
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from lverage.blast import DEFAULT_EXCLUDED_TERMS, LocalBlastSearcher


def make_alignment(description="Homeobox protein [Homo sapiens]", hsps=None):
    if hsps is None:
        hsps = [SimpleNamespace(expect=1e-20, identities=80, align_length=100)]
    return SimpleNamespace(
        accession="NP_000001.1",
        hit_def=description,
        hit_id="ref|NP_000001.1|",
        hsps=hsps,
    )


class LocalBlastSearcherTests(unittest.TestCase):

    def setUp(self):
        self.which_patcher = patch("lverage.blast.shutil.which", side_effect=lambda path: f"/resolved/{path}")
        self.run_patcher = patch("lverage.blast.subprocess.run")
        self.mock_which = self.which_patcher.start()
        self.mock_run = self.run_patcher.start()
        self.addCleanup(self.which_patcher.stop)
        self.addCleanup(self.run_patcher.stop)
        self.mock_run.return_value = SimpleNamespace(stdout="database information")

    def make_searcher(self, **kwargs):
        return LocalBlastSearcher("proteins", {"Homo sapiens": 9606}, **kwargs)

    def test_constructor_resolves_tools_and_validates_database_prefix(self):
        searcher = self.make_searcher(blastp_path="custom-blastp", blastdbcmd_path="custom-blastdbcmd")

        self.assertEqual(searcher.blastp_path, "/resolved/custom-blastp")
        self.assertEqual(searcher.blastdbcmd_path, "/resolved/custom-blastdbcmd")
        self.mock_run.assert_called_once_with(
            ["/resolved/custom-blastdbcmd", "-db", "proteins", "-info"],
            capture_output=True,
            check=True,
            text=True,
        )

    def test_constructor_copies_configuration(self):
        species_map = {"Homo sapiens": 9606}
        excluded_terms = []
        searcher = LocalBlastSearcher("proteins", species_map, excluded_terms=excluded_terms)
        species_map["Mus musculus"] = 10090
        excluded_terms.append("predicted")

        self.assertEqual(searcher.species_map, {"Homo sapiens": 9606})
        self.assertEqual(searcher.excluded_terms, [])
        self.assertEqual(LocalBlastSearcher("proteins", species_map).excluded_terms, list(DEFAULT_EXCLUDED_TERMS))

    def test_search_builds_command_and_uses_best_hsp_metrics(self):
        searcher = self.make_searcher(evalue_threshold=1e-8, top_n=7)
        alignment = make_alignment(hsps=[
            SimpleNamespace(expect=1e-5, identities=50, align_length=75),
            SimpleNamespace(expect=1e-30, identities=90, align_length=100),
        ])
        self.mock_run.side_effect = [
            SimpleNamespace(stdout="<xml>"),
            SimpleNamespace(stdout="MPEPTIDE\nSEQUENCE\n"),
        ]

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([SimpleNamespace(alignments=[alignment])])):
            records = searcher.get_orthologs("A" * 200)

        self.assertEqual(records[0].sequence, "MPEPTIDESEQUENCE")
        self.assertEqual(records[0].evalue, 1e-30)
        self.assertEqual(records[0].identity, 0.9)
        self.assertEqual(records[0].query_coverage, 0.5)
        self.assertEqual(records[0].species_tax_id, 9606)
        self.assertEqual(self.mock_run.call_args_list[1].args[0], [
            "/resolved/blastp", "-db", "proteins", "-outfmt", "5",
            "-evalue", "1e-08", "-max_target_seqs", "7",
        ])
        self.assertEqual(self.mock_run.call_args_list[2].args[0], [
            "/resolved/blastdbcmd", "-db", "proteins", "-entry",
            "ref|NP_000001.1|", "-outfmt", "%s",
        ])

    def test_species_resolution_is_case_insensitive(self):
        searcher = self.make_searcher()
        alignment = make_alignment("Homeobox protein [hOmO SaPiEnS]")
        self.mock_run.side_effect = [SimpleNamespace(stdout="<xml>"), SimpleNamespace(stdout="SEQUENCE")]

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([SimpleNamespace(alignments=[alignment])])):
            records = searcher.get_orthologs("QUERY")

        self.assertEqual(records[0].species_name, "Homo sapiens")

    def test_excluded_and_unresolved_hits_are_skipped(self):
        searcher = self.make_searcher()
        alignments = [
            make_alignment("Hypothetical protein [Homo sapiens]"),
            make_alignment("Homeobox protein [Mus musculus]"),
            make_alignment("Homeobox protein without species"),
        ]
        self.mock_run.return_value = SimpleNamespace(stdout="<xml>")

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([SimpleNamespace(alignments=alignments)])):
            records = searcher.get_orthologs("QUERY")

        self.assertEqual(records, [])

    def test_empty_exclusion_list_disables_description_filtering(self):
        searcher = self.make_searcher(excluded_terms=[])
        alignment = make_alignment("Hypothetical protein [Homo sapiens]")
        self.mock_run.side_effect = [SimpleNamespace(stdout="<xml>"), SimpleNamespace(stdout="SEQUENCE")]

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([SimpleNamespace(alignments=[alignment])])):
            records = searcher.get_orthologs("QUERY")

        self.assertEqual(len(records), 1)

    def test_blast_process_failure_propagates(self):
        searcher = self.make_searcher()
        self.mock_run.side_effect = subprocess.CalledProcessError(2, "blastp")

        with self.assertRaises(subprocess.CalledProcessError):
            searcher.get_orthologs("QUERY")

    def test_whole_document_parse_failure_propagates(self):
        searcher = self.make_searcher()
        self.mock_run.return_value = SimpleNamespace(stdout="invalid")

        with patch("lverage.blast.NCBIXML.parse", side_effect=ValueError("invalid XML")):
            with self.assertRaises(ValueError):
                searcher.get_orthologs("QUERY")


if __name__ == "__main__":
    unittest.main()
