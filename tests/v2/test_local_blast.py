import subprocess
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from lverage.blast import DEFAULT_EXCLUDED_TERMS, LocalBlastSearcher, _parse_blast_hits
from lverage.evidence import EvidenceCollector


FIXTURE_DIR = Path(__file__).parent / "fixtures"


def make_alignment(description="Homeobox protein [Homo sapiens]", hsps=None):
    if hsps is None:
        hsps = [SimpleNamespace(
            expect=1e-20,
            identities=80,
            align_length=100,
            query_start=1,
            query_end=100,
        )]
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
            SimpleNamespace(
                expect=1e-5,
                identities=50,
                align_length=75,
                query_start=1,
                query_end=75,
            ),
            SimpleNamespace(
                expect=1e-30,
                identities=90,
                align_length=100,
                query_start=51,
                query_end=150,
            ),
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
            "/resolved/blastp", "-db", "proteins", "-query", "-", "-outfmt", "5",
            "-evalue", "1e-08", "-max_target_seqs", "7",
        ])
        self.assertEqual(
            self.mock_run.call_args_list[1].kwargs["input"],
            ">query\n" + "A" * 200 + "\n",
        )
        self.assertEqual(self.mock_run.call_args_list[2].args[0], [
            "/resolved/blastdbcmd", "-db", "proteins", "-entry",
            "ref|NP_000001.1|", "-outfmt", "%s",
        ])

    def test_coverage_uses_query_coordinates_for_gapped_alignment(self):
        searcher = self.make_searcher()
        alignment = make_alignment(hsps=[SimpleNamespace(
            expect=1e-20,
            identities=8,
            align_length=12,
            query_start=11,
            query_end=18,
        )])
        self.mock_run.side_effect = [
            SimpleNamespace(stdout="<xml>"),
            SimpleNamespace(stdout="SEQUENCE"),
        ]

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([
            SimpleNamespace(alignments=[alignment])
        ])):
            records = searcher.get_orthologs("A" * 20)

        self.assertEqual(records[0].query_coverage, 0.4)

    def test_sequence_retrieval_falls_back_to_accession(self):
        searcher = self.make_searcher()
        alignment = make_alignment()
        self.mock_run.side_effect = [
            SimpleNamespace(stdout="<xml>"),
            subprocess.CalledProcessError(1, "blastdbcmd"),
            SimpleNamespace(stdout="SEQUENCE"),
        ]

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([
            SimpleNamespace(alignments=[alignment])
        ])):
            records = searcher.get_orthologs("QUERY")

        self.assertEqual(records[0].sequence, "SEQUENCE")
        self.assertEqual(self.mock_run.call_args_list[3].args[0][-3], "NP_000001.1")

    def test_sequence_retrieval_failure_for_both_identifiers_skips_hit(self):
        searcher = self.make_searcher()
        alignment = make_alignment()
        self.mock_run.side_effect = [
            SimpleNamespace(stdout="<xml>"),
            subprocess.CalledProcessError(1, "blastdbcmd"),
            subprocess.CalledProcessError(1, "blastdbcmd"),
        ]

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([
            SimpleNamespace(alignments=[alignment])
        ])):
            records = searcher.get_orthologs("QUERY")

        self.assertEqual(records, [])

    def test_incremental_hits_remain_when_later_hit_retrieval_raises(self):
        searcher = self.make_searcher()
        alignments = []
        for index in range(1, 5):
            alignment = make_alignment()
            alignment.accession = f"NP_{index:06d}.1"
            alignment.hit_id = f"ref|NP_{index:06d}.1|"
            alignments.append(alignment)
        original_error = RuntimeError("unexpected retrieval failure")
        self.mock_run.side_effect = [
            SimpleNamespace(stdout="<xml>"),
            SimpleNamespace(stdout="SEQUENCE1"),
            SimpleNamespace(stdout="SEQUENCE2"),
            SimpleNamespace(stdout="SEQUENCE3"),
            original_error,
        ]
        evidence = EvidenceCollector()

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([
            SimpleNamespace(alignments=alignments)
        ])):
            try:
                searcher.get_orthologs("QUERY", evidence_collector=evidence)
            except RuntimeError as error:
                self.assertIs(error, original_error)
                traceback_names = []
                traceback = error.__traceback__
                while traceback is not None:
                    traceback_names.append(traceback.tb_frame.f_code.co_name)
                    traceback = traceback.tb_next
                self.assertIn("_execute_mock_call", traceback_names)
                self.assertIn("__get_complete_sequence", traceback_names)
            else:
                self.fail("expected the injected BLAST retrieval exception")

        retained = [event for event in evidence.events if event["event"] == "ortholog"]
        failures = [event for event in evidence.events if event["event"] == "ortholog_retrieval" and event["status"] == "failed"]
        summaries = [event for event in evidence.events if event["event"] == "blast_search_summary"]
        self.assertEqual([event["ortholog"].accession for event in retained], [
            "NP_000001.1", "NP_000002.1", "NP_000003.1",
        ])
        self.assertEqual(len(failures), 1)
        self.assertEqual(failures[0]["hit_id"], "ref|NP_000004.1|")
        self.assertEqual(len(summaries), 1)
        self.assertEqual(summaries[0]["status"], "retrieval-failed")
        self.assertEqual(summaries[0]["retrieved_count"], 3)

    def test_fourth_hit_skip_keeps_first_three_and_reports_partial_retrieval(self):
        searcher = self.make_searcher()
        alignments = []
        for index in range(1, 5):
            alignment = make_alignment()
            alignment.accession = f"NP_{index:06d}.1"
            alignment.hit_id = f"ref|NP_{index:06d}.1|"
            alignments.append(alignment)
        self.mock_run.side_effect = [
            SimpleNamespace(stdout="<xml>"),
            SimpleNamespace(stdout="SEQUENCE1"),
            SimpleNamespace(stdout="SEQUENCE2"),
            SimpleNamespace(stdout="SEQUENCE3"),
            subprocess.CalledProcessError(1, "blastdbcmd", stderr="missing hit 4"),
            subprocess.CalledProcessError(1, "blastdbcmd", stderr="missing accession 4"),
        ]
        evidence = EvidenceCollector()

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([
            SimpleNamespace(alignments=alignments)
        ])):
            records = searcher.get_orthologs("QUERY", evidence_collector=evidence)

        self.assertEqual([record.accession for record in records], [
            "NP_000001.1", "NP_000002.1", "NP_000003.1",
        ])
        ortholog_events = [event for event in evidence.events if event["event"] == "ortholog"]
        skipped = [event for event in evidence.events if event["event"] == "ortholog_retrieval" and event["status"] == "skipped"]
        summaries = [event for event in evidence.events if event["event"] == "blast_search_summary"]
        self.assertEqual(len(ortholog_events), 3)
        self.assertEqual(len(skipped), 1)
        self.assertEqual(skipped[0]["hit_id"], "ref|NP_000004.1|")
        self.assertEqual(len(summaries), 1)
        self.assertEqual(summaries[0]["status"], "partial-retrieval-failure")

    def test_evidence_callback_failure_is_not_a_retrieval_failure(self):
        searcher = self.make_searcher()
        alignment = make_alignment()
        self.mock_run.side_effect = [
            SimpleNamespace(stdout="<xml>"),
            SimpleNamespace(stdout="SEQUENCE"),
        ]

        def fail_ortholog_row(event):
            if event["event"] == "ortholog":
                raise OSError("bundle write failed")

        evidence = EvidenceCollector(callback=fail_ortholog_row)
        with patch("lverage.blast.NCBIXML.parse", return_value=iter([
            SimpleNamespace(alignments=[alignment])
        ])):
            with self.assertRaisesRegex(OSError, "bundle write failed"):
                searcher.get_orthologs("QUERY", evidence_collector=evidence)

        self.assertFalse(any(
            event["event"] == "ortholog_retrieval"
            and event["status"] in {"failed", "skipped"}
            for event in evidence.events
        ))
        summaries = [event for event in evidence.events if event["event"] == "blast_search_summary"]
        self.assertEqual(len(summaries), 1)
        self.assertEqual(summaries[0]["status"], "evidence-write-failed")

    def test_ortholog_record_construction_failure_has_distinct_summary(self):
        searcher = self.make_searcher()
        original_error = ValueError("invalid ortholog record")
        self.mock_run.side_effect = [
            SimpleNamespace(stdout="<xml>"),
            SimpleNamespace(stdout="SEQUENCE"),
        ]
        evidence = EvidenceCollector()

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([
            SimpleNamespace(alignments=[make_alignment()])
        ])), patch("lverage.blast._make_ortholog", side_effect=original_error):
            with self.assertRaises(ValueError) as raised:
                searcher.get_orthologs("QUERY", evidence_collector=evidence)

        self.assertIs(raised.exception, original_error)
        summaries = [event for event in evidence.events if event["event"] == "blast_search_summary"]
        self.assertEqual(len(summaries), 1)
        self.assertEqual(summaries[0]["status"], "record-construction-failed")
        self.assertFalse(any(
            event["event"] == "blast_search_summary"
            and event["status"] == "evidence-write-failed"
            for event in evidence.events
        ))

    def test_no_alignments_emit_explicit_no_hits_with_reporting(self):
        searcher = self.make_searcher()
        self.mock_run.return_value = SimpleNamespace(stdout="<xml>")
        evidence = EvidenceCollector()

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([
            SimpleNamespace(alignments=[])
        ])):
            self.assertEqual(searcher.get_orthologs("QUERY", evidence_collector=evidence), [])

        searches = [event for event in evidence.events if event["event"] == "blast_search"]
        self.assertEqual(len(searches), 1)
        self.assertEqual(searches[0]["status"], "no-hits")

    def test_callers_without_evidence_reporting_keep_existing_signature_behavior(self):
        searcher = self.make_searcher()
        self.mock_run.side_effect = [
            SimpleNamespace(stdout="<xml>"),
            SimpleNamespace(stdout="SEQUENCE"),
        ]
        with patch("lverage.blast.NCBIXML.parse", return_value=iter([
            SimpleNamespace(alignments=[make_alignment()])
        ])):
            records = searcher.get_orthologs("QUERY")
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].hit_id, "ref|NP_000001.1|")

    def test_species_resolution_is_case_insensitive(self):
        searcher = self.make_searcher()
        alignment = make_alignment("Homeobox protein [hOmO SaPiEnS]")
        self.mock_run.side_effect = [SimpleNamespace(stdout="<xml>"), SimpleNamespace(stdout="SEQUENCE")]

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([SimpleNamespace(alignments=[alignment])])):
            records = searcher.get_orthologs("QUERY")

        self.assertEqual(records[0].species_name, "Homo sapiens")

    def test_species_is_read_from_the_last_bracketed_field(self):
        searcher = self.make_searcher()
        alignment = make_alignment("superoxide dismutase [Cu-Zn] [Homo sapiens]")
        self.mock_run.side_effect = [SimpleNamespace(stdout="<xml>"), SimpleNamespace(stdout="SEQUENCE")]

        with patch("lverage.blast.NCBIXML.parse", return_value=iter([SimpleNamespace(alignments=[alignment])])):
            records = searcher.get_orthologs("QUERY")

        self.assertEqual(len(records), 1)
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

    def test_biopython_parses_blast_xml_fixture(self):
        blast_xml = (FIXTURE_DIR / "blastp.xml").read_text()

        hits = _parse_blast_hits(
            blast_xml,
            4,
            {"homo sapiens": ("Homo sapiens", 9606)},
            [],
        )

        self.assertEqual(hits[0]["accession"], "NP_000001.1")
        self.assertEqual(hits[0]["evalue"], 1e-20)
        self.assertEqual(hits[0]["identity"], 0.75)
        self.assertEqual(hits[0]["query_coverage"], 1.0)


if __name__ == "__main__":
    unittest.main()
