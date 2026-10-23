import csv
import hashlib
import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from lverage.cli import RESULT_HEADERS, _build_pipeline, _validate_inputs, build_parser, main
from lverage.domain_scanner import DomainRecord, DomainScannerTemplate
from lverage.jaspar import JasparRecord
from lverage.motif_database import MotifDBTemplate
from lverage.orf_searcher import OrfSearcherTemplate
from lverage.ortholog_searcher import OrthologRecord, OrthologSearcherTemplate
from lverage.pipeline import Lverage


class StubOrfSearcher(OrfSearcherTemplate):

    def __init__(self, results):
        self.results = results

    def get_orfs(self, sequence):
        return list(self.results.get(sequence, []))


class StubDomainScanner(DomainScannerTemplate):

    def __init__(self, results=None, fail_on=None):
        self.results = results or {}
        self.fail_on = fail_on
        self.scanned = []

    def get_domains(self, sequence):
        self.scanned.append(sequence)
        if sequence == self.fail_on:
            raise RuntimeError("simulated Pfam failure")
        return list(self.results.get(sequence, []))


class StubOrthologSearcher(OrthologSearcherTemplate):

    def __init__(self, orthologs=None, fail_on=None):
        self.orthologs = orthologs or []
        self.fail_on = fail_on

    def get_orthologs(self, sequence):
        if sequence == self.fail_on:
            raise RuntimeError("simulated ortholog search failure")
        return list(self.orthologs)


class IncrementalPartialSearcher(StubOrthologSearcher):
    supports_evidence_collector = True

    def get_orthologs(self, sequence, *, evidence_collector=None):
        evidence_collector.emit(
            "blast_search", status="hits-found", hit_count=4,
            retained_hit_count=4, retrieved_count=0, skipped_count=0,
        )
        for ortholog in self.orthologs:
            evidence_collector.emit("ortholog", ortholog=ortholog, status="retained", error="")
        evidence_collector.emit(
            "ortholog_retrieval", status="skipped", hit_id="native|hit4|",
            accession="NP_4", attempted_identifier="", error="simulated retrieval failure",
        )
        evidence_collector.emit(
            "blast_search_summary", status="partial-retrieval-failure",
            hit_count=4, retrieved_count=3, skipped_count=1,
        )
        return list(self.orthologs)


class StubMotifDatabase(MotifDBTemplate):

    name = "StubJASPAR"

    def __init__(self, records=None, fail_on_call=None):
        self.records = records or []
        self.fail_on_call = fail_on_call
        self.calls = 0

    def check_species_validity(self, species_tax_id):
        return True

    def search(self, request):
        self.calls += 1
        if self.calls == self.fail_on_call:
            raise RuntimeError("simulated JASPAR failure")
        return list(self.records)


def make_motif(matrix_id, name):
    return JasparRecord(
        matrix_id, name, {"A": [1, 2], "C": [3, 4]},
        f"https://example.org/matrix/{matrix_id}/", "Homeobox", 1e-8,
    )


def make_ortholog(accession="NP_1", sequence="AAAA", description=None):
    return OrthologRecord(
        accession, description or "Homeobox protein [Homo sapiens]",
        "Homo sapiens", 9606, sequence, 1e-20, 0.8, 0.9,
        hit_id=f"native|{accession}|",
    )


class LverageCliTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        self.input_dir = self.root / "inputs"
        self.input_dir.mkdir()
        self.pfam_dir = self.root / "pfam"
        self.pfam_dir.mkdir()
        self.species_map = self.root / "species.json"
        self.species_map.write_text('{"Homo sapiens": 9606}', encoding="utf-8")

    def tearDown(self):
        self.temp_dir.cleanup()

    def write_fasta(self, name, sequence="ACGT"):
        path = self.input_dir / name
        path.write_text(f">{Path(name).stem}\n{sequence}\n", encoding="utf-8")
        return path

    def arguments(self, output_name="bundle", *extra):
        return [
            "--input-dir", str(self.input_dir),
            "--blast-db", str(self.root / "blast"),
            "--species-map", str(self.species_map),
            "--pfam-db", str(self.pfam_dir),
            "--output-dir", str(self.root / output_name),
            *extra,
        ]

    def make_pipeline(self, *, orfs=None, domains=None, orthologs=None, database=None,
                      scanner_fail_on=None, ortholog_fail_on=None):
        return Lverage(
            motif_database_list=[database or StubMotifDatabase()],
            orf_searcher=StubOrfSearcher(orfs or {}),
            domain_scanner=StubDomainScanner(domains, fail_on=scanner_fail_on),
            ortholog_searcher=StubOrthologSearcher(orthologs, fail_on=ortholog_fail_on),
        )

    def run_cli(self, pipeline, *arguments):
        stderr = io.StringIO()
        with patch("lverage.cli._build_pipeline", return_value=pipeline):
            with redirect_stdout(io.StringIO()), redirect_stderr(stderr):
                code = main(list(arguments))
        return code, stderr.getvalue()

    def read_rows(self, output_dir):
        with (output_dir / "results.tsv").open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_bundle_contains_all_three_files_and_unique_orfs_with_gene_associations(self):
        self.write_fasta("alpha.fa")
        self.write_fasta("beta.fa", "TGCA")
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["AAAA", "CCCC"], "TGCA": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[],
        )
        destination = self.root / "bundle"

        code, _ = self.run_cli(pipeline, *self.arguments())

        self.assertEqual(code, 0)
        self.assertEqual({path.name for path in destination.iterdir()}, {
            "results.tsv", "orfs.fasta", "metadata.json",
        })
        fasta = (destination / "orfs.fasta").read_text(encoding="utf-8")
        self.assertEqual(fasta.count(">"), 2)
        shared_id = "orf_" + hashlib.sha256(b"AAAA").hexdigest()
        self.assertIn(f">{shared_id}\nAAAA\n", fasta)
        rows = self.read_rows(destination)
        associations = {
            row["Input Gene"] for row in rows
            if row["Row Type"] == "orf" and row["Query ORF ID"] == shared_id
        }
        self.assertEqual(associations, {"alpha", "beta"})
        self.assertEqual(len(RESULT_HEADERS), len(next(csv.reader(
            io.StringIO((destination / "results.tsv").read_text(encoding="utf-8")), delimiter="\t"
        ))))
        metadata = json.loads((destination / "metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["summary"]["genes_completed"], 2)
        self.assertEqual(metadata["completion"], "completed")

    def test_no_orfs_and_no_blast_hits_are_biological_no_hits(self):
        self.write_fasta("empty_orfs.fa")
        no_orf_pipeline = self.make_pipeline(orfs={"ACGT": []})
        no_orf_dir = self.root / "no-orfs"
        no_orf_code, _ = self.run_cli(no_orf_pipeline, *self.arguments("no-orfs"))
        no_orf_rows = self.read_rows(no_orf_dir)
        self.assertEqual(no_orf_code, 0)
        self.assertTrue(any(row["Row Type"] == "input" for row in no_orf_rows))
        self.assertTrue(any(row["Status"] == "biological-no-hit" for row in no_orf_rows))
        self.assertTrue(any(row["Stage"] == "blast" and row["Status"] == "not-reached" for row in no_orf_rows))

        self.write_fasta("no_hits.fa")
        no_hit_pipeline = self.make_pipeline(
            orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[],
        )
        no_hit_dir = self.root / "no-hits"
        no_hit_code, _ = self.run_cli(no_hit_pipeline, *self.arguments("no-hits"))
        no_hit_rows = self.read_rows(no_hit_dir)
        self.assertEqual(no_hit_code, 0)
        self.assertTrue(any(row["Stage"] == "blast" and row["Status"] == "biological-no-hit" for row in no_hit_rows))

    def test_below_threshold_comparison_is_exported_with_distinct_status(self):
        self.write_fasta("below.fa")
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["AAAA"]},
            domains={
                "AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)],
                "TTTT": [DomainRecord("Homeobox", "PF00046.2", 0, 4)],
            },
            orthologs=[make_ortholog("NP_low", "TTTT")],
        )
        destination = self.root / "below"

        code, _ = self.run_cli(pipeline, *self.arguments("below"))

        rows = self.read_rows(destination)
        comparisons = [row for row in rows if row["Row Type"] == "domain_comparison"]
        self.assertEqual(code, 0)
        self.assertEqual(len(comparisons), 1)
        self.assertEqual(comparisons[0]["Status"], "below-threshold")
        self.assertEqual(comparisons[0]["Threshold Passed"], "false")
        self.assertEqual(comparisons[0]["Domain Pair Identity Ratio"], "0.0")
        self.assertTrue(any(row["Stage"] == "motif" and row["Status"] == "not-reached" for row in rows))

    def test_hits_one_to_three_are_exported_once_when_hit_four_is_skipped(self):
        self.write_fasta("blast_partial.fa")
        retained = [make_ortholog(f"NP_{index}", "AAAA") for index in range(1, 4)]
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
        )
        pipeline.ortholog_searcher = IncrementalPartialSearcher(retained)
        destination = self.root / "blast-partial"

        code, _ = self.run_cli(pipeline, *self.arguments("blast-partial"))

        rows = self.read_rows(destination)
        ortholog_rows = [row for row in rows if row["Row Type"] == "ortholog"]
        skipped_rows = [row for row in rows if row["Row Type"] == "ortholog_retrieval" and row["Status"] == "skipped"]
        summaries = [row for row in rows if row["Row Type"] == "stage_status" and row["Stage"] == "blast"]
        self.assertEqual(code, 1)
        self.assertEqual([row["Ortholog Accession"] for row in ortholog_rows], ["NP_1", "NP_2", "NP_3"])
        self.assertEqual(len(ortholog_rows), 3)
        self.assertEqual(len(skipped_rows), 1)
        self.assertEqual(skipped_rows[0]["Ortholog Hit ID"], "native|hit4|")
        self.assertTrue(any("partial-retrieval-failure" == row["Status"] for row in summaries))
        metadata = json.loads((destination / "metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(metadata["summary"]["genes_failed"], 1)

    def test_multiple_motifs_each_have_one_correctly_supported_row(self):
        self.write_fasta("motifs.fa")
        database = StubMotifDatabase([
            make_motif("MA0001.1", "motif one"),
            make_motif("MA0002.1", "motif two"),
        ])
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[make_ortholog("NP_1", "AAAA")], database=database,
        )
        destination = self.root / "motifs"

        code, _ = self.run_cli(pipeline, *self.arguments("motifs"))

        motifs = [row for row in self.read_rows(destination) if row["Row Type"] == "motif"]
        self.assertEqual(code, 0)
        self.assertEqual(len(motifs), 2)
        self.assertEqual({row["Matrix ID"] for row in motifs}, {"MA0001.1", "MA0002.1"})
        self.assertTrue(all(row["Input Gene"] == "motifs" for row in motifs))
        self.assertTrue(all(row["Ortholog Accession"] == "NP_1" for row in motifs))
        self.assertTrue(all(row["Query Domain ID"] == "PF00046.1" for row in motifs))
        self.assertTrue(all(row["Ortholog Domain ID"] == "PF00046.1" for row in motifs))
        self.assertEqual(json.loads(motifs[0]["PFM JSON"]), {"A": [1, 2], "C": [3, 4]})

    def test_later_ortholog_failure_preserves_completed_first_ortholog_evidence(self):
        self.write_fasta("later_ortholog.fa")
        first = make_ortholog("NP_first", "AAAA")
        second = make_ortholog("NP_second", "CCCC")
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[first, second], scanner_fail_on="CCCC",
        )
        destination = self.root / "ortholog-failure"

        code, _ = self.run_cli(pipeline, *self.arguments("ortholog-failure"))

        rows = self.read_rows(destination)
        self.assertEqual(code, 1)
        self.assertTrue(any(row["Row Type"] == "ortholog_domain" and row["Ortholog Accession"] == "NP_first" for row in rows))
        self.assertTrue(any(row["Row Type"] == "gene_error" and row["Stage"] == "ortholog_pfam" for row in rows))

    def test_jaspar_failure_preserves_prior_motif_and_alignment_evidence(self):
        self.write_fasta("jaspar_failure.fa")
        sequence = "AAAAAAAA"
        database = StubMotifDatabase([make_motif("MA0001.1", "first")], fail_on_call=2)
        pipeline = self.make_pipeline(
            orfs={"ACGT": [sequence]},
            domains={sequence: [
                DomainRecord("First", "PF00046.1", 0, 4),
                DomainRecord("Second", "PF00046.2", 4, 8),
            ]},
            orthologs=[make_ortholog("NP_1", sequence)], database=database,
        )
        destination = self.root / "jaspar-failure"

        code, _ = self.run_cli(pipeline, *self.arguments("jaspar-failure"))

        rows = self.read_rows(destination)
        self.assertEqual(code, 1)
        self.assertEqual(len([row for row in rows if row["Row Type"] == "domain_comparison"]), 2)
        self.assertEqual(len([row for row in rows if row["Row Type"] == "motif"]), 1)
        self.assertTrue(any(row["Row Type"] == "gene_error" and row["Stage"] == "motif" for row in rows))

    def test_cli_continues_after_gene_execution_failure(self):
        self.write_fasta("first.fa")
        self.write_fasta("second.fa", "TGCA")
        first_orf = "AAAA"
        second_orf = "CCCC"
        pipeline = self.make_pipeline(
            orfs={"ACGT": [first_orf], "TGCA": [second_orf]},
            domains={first_orf: [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[], scanner_fail_on=first_orf,
        )
        destination = self.root / "continue"

        code, stderr = self.run_cli(pipeline, *self.arguments("continue"))

        rows = self.read_rows(destination)
        self.assertEqual(code, 1)
        self.assertIn("gene first", stderr)
        self.assertTrue(any(row["Input Gene"] == "first" and row["Row Type"] == "gene_error" for row in rows))
        self.assertTrue(any(row["Input Gene"] == "second" and row["Row Type"] == "input" for row in rows))
        self.assertEqual(json.loads((destination / "metadata.json").read_text())["summary"]["genes_failed"], 1)
        self.assertEqual(json.loads((destination / "metadata.json").read_text())["summary"]["genes_completed"], 1)

    def test_stop_after_alignment_never_calls_jaspar(self):
        self.write_fasta("alignment.fa")
        database = StubMotifDatabase([make_motif("MA0001.1", "unused")])
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[make_ortholog()], database=database,
        )
        destination = self.root / "alignment"

        code, _ = self.run_cli(pipeline, *self.arguments("alignment", "--stop-after", "alignment"))

        rows = self.read_rows(destination)
        metadata = json.loads((destination / "metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(code, 0)
        self.assertEqual(database.calls, 0)
        self.assertEqual(metadata["completion"], "intentionally-stopped")
        self.assertTrue(any(row["Stage"] == "alignment" and row["Row Type"] == "domain_comparison" for row in rows))
        self.assertTrue(any(row["Stage"] == "motif" and row["Status"] == "not-reached" for row in rows))

    def test_invalid_inputs_and_existing_output_exit_two(self):
        self.write_fasta("valid.fa")
        existing = self.root / "already-there"
        existing.mkdir()
        stderr = io.StringIO()
        with patch("lverage.cli._build_pipeline") as build:
            with redirect_stderr(stderr):
                code = main(self.arguments("already-there"))
        self.assertEqual(code, 2)
        self.assertIn("already exists", stderr.getvalue())
        build.assert_not_called()

        self.input_dir.joinpath("valid.fa").unlink()
        with redirect_stderr(io.StringIO()):
            self.assertEqual(main(self.arguments("missing-input")), 2)

    def test_parser_requires_output_dir_and_supports_stop_after(self):
        args = self.arguments("parsed", "--stop-after", "alignment")
        parsed = build_parser().parse_args(args)
        self.assertEqual(parsed.output_dir, str(self.root / "parsed"))
        self.assertEqual(parsed.stop_after, "alignment")
        with self.assertRaises(SystemExit):
            build_parser().parse_args([
                "--input-dir", str(self.input_dir), "--blast-db", "db",
                "--species-map", str(self.species_map), "--pfam-db", str(self.pfam_dir),
            ])

    def test_validate_input_discovery_preserves_order_and_fragment_lists(self):
        self.write_fasta("zeta.FNA")
        self.write_fasta("Alpha.fa", "TGCA")
        parsed = build_parser().parse_args(self.arguments("validate"))
        inputs, species_map, species_path = _validate_inputs(parsed)
        self.assertEqual([item["gene"] for item in inputs], ["Alpha", "zeta"])
        self.assertEqual(inputs[0]["fragments"], ["TGCA"])
        self.assertEqual(species_map, {"Homo sapiens": 9606})
        self.assertEqual(species_path, self.species_map)

    def test_build_pipeline_passes_existing_options(self):
        parsed = build_parser().parse_args(self.arguments(
            "configured", "--blastp", "/tools/blastp", "--blastdbcmd", "/tools/blastdbcmd",
            "--pfamscan", "/tools/pfamscan", "--identity-threshold", "0.83",
            "--valid-pfam", "PF00046", "--jaspar-hits", "7", "--jaspar-evalue", "2e-5",
        ))
        with patch("lverage.blast.LocalBlastSearcher") as blast, \
                patch("lverage.pfam.LocalPfamScanner") as pfam, \
                patch("lverage.orf_searcher.OrffinderOrfSearcher") as orf, \
                patch("lverage.jaspar.Jaspar2024MotifDB") as jaspar, \
                patch("lverage.pipeline.Lverage") as lverage:
            _build_pipeline(parsed, {"Homo sapiens": 9606})
        blast.assert_called_once_with(
            database_path=str(self.root / "blast"), species_map={"Homo sapiens": 9606},
            blastp_path="/tools/blastp", blastdbcmd_path="/tools/blastdbcmd",
        )
        pfam.assert_called_once_with(database_path=str(self.pfam_dir), pfamscan_path="/tools/pfamscan")
        jaspar.assert_called_once_with(n_hits=7, escore_threshold=2e-5)
        lverage.assert_called_once_with(
            motif_database_list=[jaspar.return_value], orf_searcher=orf.return_value,
            domain_scanner=pfam.return_value, ortholog_searcher=blast.return_value,
            valid_pfam_list=["PF00046"], dbd_identity_thresh=0.83,
        )


if __name__ == "__main__":
    unittest.main()
