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
            hit_group=4, description="Skipped protein [Homo sapiens]",
            species_name="Homo sapiens", species_tax_id=9606, evalue=1e-20,
            identity_ratio=0.8, identity_percent=80.0,
            query_coverage_ratio=0.9, query_coverage_percent=90.0,
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
        matrix_id, name, {"A": [1, 2], "C": [3, 4], "G": [4, 3], "T": [2, 1]},
        f"https://example.org/matrix/{matrix_id}/", "Homeobox", 1e-8,
        inference_dbd_identity=0.875,
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
        with (output_dir / "results.csv").open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter=","))

    def read_errors(self, output_dir):
        with (output_dir / "errors.csv").open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle))

    def matched_pipeline(self, database=None):
        return self.make_pipeline(
            orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)],
                     "CCCC": [DomainRecord("Homeobox", "PF00046.2", 0, 4)]},
            orthologs=[make_ortholog(sequence="CCCC")], database=database,
        )

    def test_exact_ten_column_contract_and_separate_error_file(self):
        self.write_fasta("gene.fa")
        pipeline = self.matched_pipeline(StubMotifDatabase([make_motif("MA0001.1", "Example")]))
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.9):
            code, _ = self.run_cli(pipeline, *self.arguments())
        destination = self.root / "bundle"
        self.assertEqual(code, 0)
        self.assertEqual({p.name for p in destination.iterdir()},
                         {"results.csv", "errors.csv", "orfs.fasta", "metadata.json"})
        rows = self.read_rows(destination)
        self.assertEqual(list(rows[0]), [
            "Gene", "Ortholog Species", "BLAST Description", "Query Pfam Name",
            "Hit Pfam Name", "Global Alignment Identity Percent", "Motif Name",
            "JASPAR Tool DBD Identity Ratio", "Motif PFM JSON", "ORF ID",
        ])
        self.assertEqual(rows[0]["Global Alignment Identity Percent"], "90.0")
        self.assertEqual(rows[0]["Motif Name"], "Example")
        self.assertEqual(rows[0]["JASPAR Tool DBD Identity Ratio"], "0.875")
        self.assertEqual(json.loads(rows[0]["Motif PFM JSON"])["A"], [1, 2])
        self.assertEqual(self.read_errors(destination), [])
        metadata = json.loads((destination / "metadata.json").read_text())
        self.assertEqual(metadata["schema_version"], 3)
        self.assertEqual(metadata["results_format"], "compact-csv")
        self.assertEqual(metadata["summary"]["genes_completed"], 1)

    def test_candidates_are_longest_first_and_shared_fasta_is_unique(self):
        self.write_fasta("alpha.fa")
        self.write_fasta("beta.fa", "TGCA")
        selected = "AAAAAA"
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["CCC", selected, "MMMMMMMM", selected], "TGCA": [selected]},
            domains={selected: [DomainRecord("Homeobox", "PF00046.1", 0, 4)]})
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        alpha = [r for r in rows if r["Gene"] == "alpha"]
        self.assertEqual([r["ORF ID"] for r in alpha],
                         ["orf_" + hashlib.sha256(s.encode()).hexdigest()
                          for s in ("MMMMMMMM", selected, "CCC")])
        self.assertEqual(alpha[0]["Query Pfam Name"], "")
        self.assertEqual(alpha[1]["Query Pfam Name"], "Homeobox")
        self.assertEqual(alpha[2]["Query Pfam Name"], "")
        self.assertEqual((self.root / "bundle/orfs.fasta").read_text().count(">"), 3)
        self.assertEqual(pipeline.domain_scanner.scanned, ["MMMMMMMM", selected, selected])

    def test_no_orfs_and_no_blast_hits_are_successful_with_blank_downstream_cells(self):
        self.write_fasta("empty.fa")
        code, _ = self.run_cli(self.make_pipeline(), *self.arguments("no-orfs"))
        rows = self.read_rows(self.root / "no-orfs")
        self.assertEqual(code, 0)
        self.assertEqual(rows, [{h: "empty" if h == "Gene" else "" for h in RESULT_HEADERS}])
        self.assertEqual(self.read_errors(self.root / "no-orfs"), [])
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]})
        code, _ = self.run_cli(pipeline, *self.arguments("no-hits"))
        row = self.read_rows(self.root / "no-hits")[0]
        self.assertEqual(code, 0)
        self.assertEqual(row["Query Pfam Name"], "Homeobox")
        self.assertEqual(row["BLAST Description"], "")
        self.assertEqual(self.read_errors(self.root / "no-hits"), [])

    def test_below_threshold_keeps_alignment_and_does_not_call_motif_search(self):
        self.write_fasta("gene.fa")
        database = StubMotifDatabase()
        pipeline = self.matched_pipeline(database)
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.4):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Global Alignment Identity Percent"], "40.0")
        self.assertEqual(rows[0]["Motif Name"], "")
        self.assertEqual(database.calls, 0)
        self.assertEqual(self.read_errors(self.root / "bundle"), [])

    def test_multiple_motifs_repeat_actual_parent_fields_without_deduplication(self):
        self.write_fasta("gene.fa")
        database = StubMotifDatabase([make_motif("MA0001.1", "Same name"),
                                      make_motif("MA0001.2", "Same name")])
        pipeline = self.matched_pipeline(database)
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.9):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0], rows[1])
        self.assertTrue(all(r["Query Pfam Name"] == r["Hit Pfam Name"] == "Homeobox" for r in rows))

    def test_same_named_domain_occurrences_keep_their_own_alignment_scores(self):
        self.write_fasta("gene.fa")
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["AAAAAAAA"]},
            domains={
                "AAAAAAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4),
                             DomainRecord("Homeobox", "PF00046.1", 4, 8)],
                "CCCCCCCC": [DomainRecord("Homeobox", "PF00046.2", 0, 8)],
            },
            orthologs=[make_ortholog(sequence="CCCCCCCC")],
            database=StubMotifDatabase([make_motif("MA0001.1", "Example")]),
        )
        with patch.object(pipeline, "_calculate_domain_identity", side_effect=[0.9, 0.8]):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual(len(rows), 2)
        self.assertEqual([r["Global Alignment Identity Percent"] for r in rows], ["90.0", "80.0"])
        self.assertTrue(all(r["Motif Name"] == "Example" for r in rows))

    def test_csv_quotes_commas_newlines_and_matrix_json(self):
        self.write_fasta("gene.fa")
        description = 'Protein, "example"\n[Homo sapiens]'
        pipeline = self.matched_pipeline(StubMotifDatabase([make_motif("MA0001.1", "TF, example")]))
        pipeline.ortholog_searcher.orthologs[0].description = description
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.9):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["BLAST Description"], description)
        self.assertEqual(rows[0]["Motif Name"], "TF, example")
        self.assertEqual(set(json.loads(rows[0]["Motif PFM JSON"])), set("ACGT"))
        self.assertEqual(set(rows[0]), set(RESULT_HEADERS))

    def test_gene_motifs_rank_across_hits_with_zero_before_blank_and_genes_together(self):
        self.write_fasta("alpha.fa")
        self.write_fasta("beta.fa", "TGCA")
        database = StubMotifDatabase()
        def predictions(request):
            scores = [0.8, 0.0] if request.ortholog_sequence == "CCCC" else [1.0, 0.9]
            records = [make_motif(f"MA000{i+1}.1", f"motif {i}") for i in range(2)]
            for record, score in zip(records, scores):
                record.inference_dbd_identity = score
            return records
        database.search = predictions
        domain = DomainRecord("Homeobox", "PF00046.1", 0, 4)
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["MMMMMMMM", "AAAA"], "TGCA": ["MMMMMMMM", "AAAA"]},
            domains={"AAAA": [domain], "CCCC": [domain], "GGGG": [domain]},
            orthologs=[make_ortholog("NP_1", "CCCC", "First hit"),
                       make_ortholog("NP_2", "GGGG", "Second hit")],
            database=database,
        )
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.9):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual([r["Gene"] for r in rows], ["alpha"] * 5 + ["beta"] * 5)
        for gene in ("alpha", "beta"):
            group = [r for r in rows if r["Gene"] == gene]
            self.assertEqual([r["JASPAR Tool DBD Identity Ratio"] for r in group],
                             ["1.0", "0.9", "0.8", "0.0", ""])
            self.assertEqual(group[0]["BLAST Description"], "Second hit")
            self.assertEqual(group[-1]["Query Pfam Name"], "")
            self.assertEqual(group[-1]["ORF ID"],
                             "orf_" + hashlib.sha256(b"MMMMMMMM").hexdigest())

    def test_equal_identity_preserves_inference_order(self):
        self.write_fasta("gene.fa")
        records = [make_motif("MA0002.1", "Second alphabetically"),
                   make_motif("MA0001.1", "First alphabetically")]
        pipeline = self.matched_pipeline(StubMotifDatabase(records))
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.9):
            code, _ = self.run_cli(pipeline, *self.arguments())
        self.assertEqual(code, 0)
        self.assertEqual([r["Motif Name"] for r in self.read_rows(self.root / "bundle")],
                         [record.motif_name for record in records])

    def test_unmatched_hit_domain_has_no_invented_query_or_alignment(self):
        self.write_fasta("gene.fa")
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]}, domains={
            "AAAA": [DomainRecord("Homeobox", "PF00046", 0, 4)],
            "CCCC": [DomainRecord("Other", "PF99999", 0, 4)]},
            orthologs=[make_ortholog(sequence="CCCC")])
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        hit = next(r for r in rows if r["Hit Pfam Name"])
        self.assertEqual(code, 0)
        self.assertEqual(hit["Hit Pfam Name"], "Other")
        self.assertEqual(hit["Query Pfam Name"], "")
        self.assertEqual(hit["Global Alignment Identity Percent"], "")

    def test_alignment_failure_keeps_domain_evidence_and_logs_error_separately(self):
        self.write_fasta("gene.fa")
        pipeline = self.matched_pipeline()
        with patch.object(pipeline, "_calculate_domain_identity", side_effect=RuntimeError("alignment failed")):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Hit Pfam Name"], "Homeobox")
        self.assertEqual(rows[0]["Global Alignment Identity Percent"], "")
        self.assertEqual(rows[0]["Motif Name"], "")
        error = self.read_errors(self.root / "bundle")[0]
        self.assertEqual(error["Stage"], "alignment")
        self.assertEqual(error["Error"], "alignment failed")
        self.assertEqual(error["ORF ID"], rows[0]["ORF ID"])

    def test_partial_motif_generator_keeps_first_prediction_without_error_result_row(self):
        self.write_fasta("gene.fa")
        database = StubMotifDatabase()
        def partial(request):
            yield make_motif("MA0001.1", "first")
            raise RuntimeError("later motif failed")
        database.search = partial
        pipeline = self.matched_pipeline(database)
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.9):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Motif Name"], "first")
        self.assertEqual(self.read_errors(self.root / "bundle")[0]["Stage"], "motif")

    def test_later_hit_pfam_failure_preserves_prior_motifs_and_failed_hit(self):
        self.write_fasta("gene.fa")
        pipeline = self.matched_pipeline(StubMotifDatabase([make_motif("MA0001.1", "first")]))
        pipeline.ortholog_searcher.orthologs.append(make_ortholog("NP_2", "GGGG", "Second hit"))
        pipeline.domain_scanner.fail_on = "GGGG"
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.9):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["Motif Name"], "first")
        self.assertEqual(rows[1]["BLAST Description"], "Second hit")
        self.assertEqual(rows[1]["Hit Pfam Name"], "")
        self.assertEqual(rows[1]["Motif Name"], "")
        self.assertEqual(self.read_errors(self.root / "bundle")[0]["BLAST Accession"], "NP_2")

    def test_skipped_retrieval_keeps_available_hit_details_and_separate_diagnostic(self):
        self.write_fasta("gene.fa")
        pipeline = self.matched_pipeline()
        pipeline.ortholog_searcher = IncrementalPartialSearcher([make_ortholog(sequence="CCCC")])
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.9):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        skipped = [r for r in rows if r["BLAST Description"] == "Skipped protein [Homo sapiens]"]
        self.assertEqual(len(skipped), 1)
        self.assertEqual(skipped[0]["Hit Pfam Name"], "")
        self.assertEqual(self.read_errors(self.root / "bundle")[0]["BLAST Accession"], "NP_4")

    def test_cli_continues_after_gene_failure_and_finalizes_both_genes(self):
        self.write_fasta("alpha.fa")
        self.write_fasta("beta.fa", "TGCA")
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"], "TGCA": ["CCCC"]},
            domains={"CCCC": [DomainRecord("Homeobox", "PF00046", 0, 4)]},
            scanner_fail_on="AAAA")
        code, _ = self.run_cli(pipeline, *self.arguments())
        self.assertEqual(code, 1)
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual([r["Gene"] for r in rows], ["alpha", "beta"])
        self.assertEqual(self.read_errors(self.root / "bundle")[0]["Gene"], "alpha")

    def test_stop_after_alignment_does_not_call_jaspar(self):
        self.write_fasta("gene.fa")
        database = StubMotifDatabase(fail_on_call=1)
        pipeline = self.matched_pipeline(database)
        with patch.object(pipeline, "_calculate_domain_identity", return_value=0.9):
            code, _ = self.run_cli(pipeline, *self.arguments("stopped", "--stop-after", "alignment"))
        self.assertEqual(code, 0)
        self.assertEqual(database.calls, 0)
        rows = self.read_rows(self.root / "stopped")
        self.assertEqual(rows[0]["Global Alignment Identity Percent"], "90.0")
        self.assertEqual(rows[0]["Motif Name"], "")
        self.assertEqual(self.read_errors(self.root / "stopped"), [])

    def test_initialization_failure_writes_gene_rows_and_separate_errors(self):
        self.write_fasta("alpha.fa")
        self.write_fasta("beta.fa")
        with patch("lverage.cli._build_pipeline", side_effect=RuntimeError("setup failed")):
            with redirect_stderr(io.StringIO()):
                code = main(self.arguments())
        rows = self.read_rows(self.root / "bundle")
        errors = self.read_errors(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertEqual([r["Gene"] for r in rows], ["alpha", "beta"])
        self.assertTrue(all(not value for r in rows for key, value in r.items() if key != "Gene"))
        self.assertEqual([r["Stage"] for r in errors], ["initialization", "initialization"])
        self.assertTrue(all(r["Error"] == "setup failed" for r in errors))

    def test_csv_is_not_written_on_events_and_is_published_only_at_finish(self):
        from lverage.cli import _Bundle
        destination = self.root / "deferred"
        destination.mkdir()
        bundle = _Bundle(destination, {})
        bundle.add_event("gene", {"event": "gene_start", "stage": "input"})
        self.assertFalse((destination / "results.csv").exists())
        self.assertFalse((destination / "errors.csv").exists())
        bundle.finish()
        self.assertTrue((destination / "results.csv").exists())
        self.assertTrue((destination / "errors.csv").exists())

    def test_formatter_failure_preserves_previous_csv(self):
        from lverage.cli import _Bundle
        destination = self.root / "atomic"
        destination.mkdir()
        bundle = _Bundle(destination, {})
        bundle.add_event("gene", {"event": "gene_start", "stage": "input"})
        bundle.finish()
        before = (destination / "results.csv").read_bytes()
        with patch("lverage.cli._grouped_result_rows", side_effect=RuntimeError("format failed")):
            with self.assertRaisesRegex(RuntimeError, "format failed"):
                bundle.finish()
        self.assertEqual((destination / "results.csv").read_bytes(), before)

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
            "--pfamscan", "/tools/pfamscan", "--hmmscan", "/tools/hmmscan",
            "--hmmalign", "/tools/hmmalign", "--identity-threshold", "0.83",
            "--valid-pfam", "PF00046", "--jaspar-hits", "7", "--jaspar-evalue", "2e-5",
        ))
        with patch("lverage.blast.LocalBlastSearcher") as blast, \
                patch("lverage.pfam.LocalPfamScanner") as pfam, \
                patch("lverage.orf_searcher.OrffinderOrfSearcher") as orf, \
                patch("lverage.jaspar.Jaspar2026MotifDB") as jaspar, \
                patch("lverage.pipeline.Lverage") as lverage:
            _build_pipeline(parsed, {"Homo sapiens": 9606})
        blast.assert_called_once_with(
            database_path=str(self.root / "blast"), species_map={"Homo sapiens": 9606},
            blastp_path="/tools/blastp", blastdbcmd_path="/tools/blastdbcmd",
        )
        pfam.assert_called_once_with(database_path=str(self.pfam_dir), pfamscan_path="/tools/pfamscan")
        jaspar.assert_called_once_with(
            n_hits=7, escore_threshold=2e-5, data_dir=None, sqlite_db=None,
            blastp_path="/tools/blastp", hmmscan_path="/tools/hmmscan",
            hmmalign_path="/tools/hmmalign",
        )
        lverage.assert_called_once_with(
            motif_database_list=[jaspar.return_value], orf_searcher=orf.return_value,
            domain_scanner=pfam.return_value, ortholog_searcher=blast.return_value,
            valid_pfam_list=["PF00046"], dbd_identity_thresh=0.83,
        )


if __name__ == "__main__":
    unittest.main()
