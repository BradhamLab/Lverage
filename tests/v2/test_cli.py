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
        with (output_dir / "results.csv").open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter=","))

    def test_bundle_contains_all_three_files_and_unique_orfs_with_gene_associations(self):
        self.write_fasta("alpha.fa")
        self.write_fasta("beta.fa", "TGCA")
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["AAAA", "CCCC"], "TGCA": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]})
        pipeline.ortholog_searcher.excluded_terms = ["configured description term"]
        code, _ = self.run_cli(pipeline, *self.arguments())
        destination = self.root / "bundle"
        self.assertEqual(code, 0)
        self.assertEqual({p.name for p in destination.iterdir()}, {"results.csv", "orfs.fasta", "metadata.json"})
        fasta = (destination / "orfs.fasta").read_text()
        shared_id = "orf_" + hashlib.sha256(b"AAAA").hexdigest()
        self.assertEqual(fasta.count(">"), 2)
        self.assertIn(f">{shared_id}\nAAAA\n", fasta)
        self.assertEqual({r["Gene"] for r in self.read_rows(destination) if r["ORF ID"] == shared_id}, {"alpha", "beta"})
        metadata = json.loads((destination / "metadata.json").read_text())
        self.assertEqual(metadata["summary"]["genes_completed"], 2)
        self.assertEqual(metadata["schema_version"], 2)
        self.assertEqual(metadata["results_format"], "grouped-csv")
        self.assertEqual(
            metadata["settings"]["blast_description_excluded_terms"],
            ["configured description term"],
        )
        with (destination / "results.csv").open() as f:
            lines = list(csv.reader(f, delimiter=","))
        self.assertEqual(lines[0], RESULT_HEADERS)
        self.assertTrue(all(len(line) == len(RESULT_HEADERS) for line in lines))

    def test_candidate_order_selection_and_blank_unexamined_stages(self):
        self.write_fasta("candidates.fa")
        selected = "AAAAAA"
        pipeline = self.make_pipeline(
            orfs={"ACGT": ["CCC", selected, "MMMMMMMM", selected]},
            domains={selected: [DomainRecord("Homeobox", "PF00046.1", 0, 4)]})
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual([r["ORF Length (aa)"] for r in rows], ["8", "6", "3"])
        self.assertEqual(pipeline.domain_scanner.scanned, ["MMMMMMMM", selected])
        self.assertEqual(rows[0]["Query Pfam Status"], "no accepted domain")
        self.assertEqual(rows[1]["ORF Selection Status"], "selected")
        self.assertEqual(rows[2]["ORF Examination Status"], "")
        self.assertEqual(rows[2]["ORF Selection Status"], "")
        self.assertTrue(all(not rows[2][h] for h in RESULT_HEADERS[10:]))
        self.assertTrue(all(r["Gene"] == "candidates" for r in rows))

    def test_no_orfs_and_no_blast_hits_are_biological_no_hits(self):
        self.write_fasta("empty.fa")
        code, _ = self.run_cli(self.make_pipeline(), *self.arguments("no-orfs"))
        rows = self.read_rows(self.root / "no-orfs")
        self.assertEqual(code, 0)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Gene Status"], "biological-no-hit")
        self.assertEqual(rows[0]["ORF ID"], "")
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]})
        code, _ = self.run_cli(pipeline, *self.arguments("no-hits"))
        rows = self.read_rows(self.root / "no-hits")
        self.assertEqual(code, 0)
        self.assertEqual(rows[0]["BLAST Status"], "biological-no-hit")
        self.assertEqual(rows[0]["BLAST Accession"], "")

    def test_below_threshold_comparison_is_exported_with_distinct_status(self):
        self.write_fasta("below.fa")
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)],
                     "TTTT": [DomainRecord("Homeobox", "PF00046.2", 0, 4)]},
            orthologs=[make_ortholog("NP_low", "TTTT")])
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Global Alignment Status"], "below-threshold")
        self.assertEqual(rows[0]["Global Alignment Threshold Passed"], "false")
        self.assertEqual(rows[0]["Global Alignment Identity Ratio"], "0.0")
        self.assertEqual(rows[0]["Motif Status"], "")

    def test_hits_one_to_three_are_exported_once_when_hit_four_is_skipped(self):
        self.write_fasta("partial.fa")
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]})
        pipeline.ortholog_searcher = IncrementalPartialSearcher([make_ortholog(f"NP_{i}") for i in range(1, 4)])
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertEqual([r["BLAST Accession"] for r in rows], ["NP_1", "NP_2", "NP_3", "NP_4"])
        self.assertEqual(rows[-1]["BLAST Status"], "skipped")
        self.assertEqual(rows[-1]["BLAST Hit ID"], "native|hit4|")
        self.assertEqual(rows[-1]["BLAST Hit Group"], "4")
        self.assertEqual(rows[-1]["BLAST Description"], "Skipped protein [Homo sapiens]")
        self.assertEqual(rows[-1]["BLAST Identity Percent"], "80.0")
        self.assertEqual(rows[-1]["Query Coverage Percent"], "90.0")
        self.assertIn("partial-retrieval-failure", rows[0]["BLAST Outcome"])
        self.assertEqual(rows[0]["Gene Status"], "failed")

    def test_multiple_motifs_each_have_one_correctly_supported_row(self):
        self.write_fasta("motifs.fa")
        database = StubMotifDatabase([make_motif("MA0001.1", "one"), make_motif("MA0002.1", "two")])
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[make_ortholog()], database=database)
        code, _ = self.run_cli(pipeline, *self.arguments())
        raw = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual(len(raw), 2)  # No extra ORF/hit/alignment-only rows.
        self.assertEqual([r["Motif Matrix ID"] for r in raw], ["MA0001.1", "MA0002.1"])
        repeated_fields = [
            "Gene", "Gene Status", "ORF Candidate Order", "ORF ID",
            "ORF Selection Status", "Query Pfam Occurrence", "Query Pfam ID",
            "BLAST Hit Group", "BLAST Hit ID", "BLAST Accession",
            "Hit Pfam Occurrence", "Hit Pfam ID", "Global Alignment ID",
            "Global Alignment Identity Ratio",
        ]
        self.assertTrue(all(all(row[field] for field in repeated_fields) for row in raw))
        self.assertEqual({row["Motif Matrix ID"] for row in raw}, {"MA0001.1", "MA0002.1"})
        self.assertTrue(all(row["Gene"] == "motifs" and row["BLAST Accession"] == "NP_1" for row in raw))
        self.assertTrue(all(row["Query Pfam ID"] == "PF00046.1" and row["Hit Pfam ID"] == "PF00046.1" for row in raw))
        self.assertEqual(raw[0]["Global Alignment ID"], raw[1]["Global Alignment ID"])
        self.assertEqual(json.loads(raw[0]["Motif PFM JSON"]), {"A": [1, 2], "C": [3, 4]})

    def test_unpaired_domain_is_visible_once_without_alignment(self):
        self.write_fasta("domains.fa")
        query = "AAAAAAAA"
        hit = "AAAAAAAT"
        pipeline = self.make_pipeline(orfs={"ACGT": [query]}, domains={
            query: [DomainRecord("Q1", "PF00046.1", 0, 4), DomainRecord("Q2", "PF00047.1", 4, 8)],
            hit: [DomainRecord("H1", "PF00046.2", 0, 4), DomainRecord("H2", "PF00047.2", 4, 8),
                  DomainRecord("Other", "PF99999.1", 0, 3)]},
            orthologs=[make_ortholog(sequence=hit)])
        code, _ = self.run_cli(pipeline, *self.arguments("domains", "--stop-after", "alignment"))
        rows = self.read_rows(self.root / "domains")
        self.assertEqual(code, 0)
        self.assertEqual(len(rows), 3)
        unmatched = [r for r in rows if r["Hit Pfam ID"] == "PF99999.1"]
        self.assertEqual(len(unmatched), 1)
        self.assertEqual(unmatched[0]["Global Alignment ID"], "")
        matched = [r for r in rows if r.get("Global Alignment ID")]
        self.assertEqual([(r["Query Pfam ID"].split('.')[0], r["Hit Pfam ID"].split('.')[0]) for r in matched],
                         [("PF00046", "PF00046"), ("PF00047", "PF00047")])

    def test_identical_domain_annotations_keep_occurrence_associations(self):
        self.write_fasta("occurrences.fa")
        sequence = "AAAAAAAA"
        domain1 = DomainRecord("Same", "PF00046.1", 0, 4)
        domain2 = DomainRecord("Same", "PF00046.1", 4, 8)
        pipeline = self.make_pipeline(orfs={"ACGT": [sequence]},
            domains={sequence: [domain1, domain2]}, orthologs=[make_ortholog(sequence=sequence)],
            database=StubMotifDatabase([make_motif("MA0001.1", "one")]))
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual(len(rows), 4)
        self.assertEqual({(r["Query Pfam Start (0-based aa)"], r["Hit Pfam Start (0-based aa)"]) for r in rows},
                         {("0", "0"), ("0", "4"), ("4", "0"), ("4", "4")})
        self.assertTrue(all(r["Motif Matrix ID"] == "MA0001.1" for r in rows))

    def test_later_ortholog_failure_preserves_completed_first_ortholog_evidence(self):
        self.write_fasta("later.fa")
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[make_ortholog("NP_first"), make_ortholog("NP_second", "CCCC")], scanner_fail_on="CCCC")
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertTrue(any(r["BLAST Accession"] == "NP_first" and r.get("Global Alignment ID") for r in rows))
        self.assertTrue(any(r["BLAST Accession"] == "NP_second" and r["Hit Pfam Status"] == "failed" for r in rows))

    def test_jaspar_failure_preserves_prior_motif_and_alignment_evidence(self):
        self.write_fasta("jaspar.fa")
        sequence = "AAAAAAAA"
        database = StubMotifDatabase([make_motif("MA0001.1", "first")], fail_on_call=2)
        pipeline = self.make_pipeline(orfs={"ACGT": [sequence]}, domains={sequence: [
            DomainRecord("First", "PF00046.1", 0, 4), DomainRecord("Second", "PF00046.2", 4, 8)]},
            orthologs=[make_ortholog(sequence=sequence)], database=database)
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertEqual(sum(bool(r["Global Alignment ID"]) for r in rows), 2)
        self.assertEqual(sum(bool(r["Motif Matrix ID"]) for r in rows), 1)
        failed = [r for r in rows if r["Motif Status"] == "failed"]
        self.assertEqual(len(failed), 1)
        self.assertIn("simulated JASPAR failure", failed[0]["Motif Error"])

    def test_cli_continues_after_gene_execution_failure(self):
        self.write_fasta("first.fa")
        self.write_fasta("second.fa", "TGCA")
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"], "TGCA": ["CCCC"]}, scanner_fail_on="AAAA")
        code, stderr = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertIn("gene first", stderr)
        self.assertEqual([r["Gene"] for r in rows], ["first", "second"])
        self.assertEqual(rows[0]["Query Pfam Status"], "failed")
        self.assertEqual(rows[1]["Query Pfam Status"], "no accepted domain")
        metadata = json.loads((self.root / "bundle/metadata.json").read_text())
        self.assertEqual(metadata["summary"]["genes_failed"], 1)
        self.assertEqual(metadata["summary"]["genes_completed"], 1)

    def test_stop_after_alignment_never_calls_jaspar(self):
        self.write_fasta("alignment.fa")
        database = StubMotifDatabase([make_motif("MA0001.1", "unused")])
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[make_ortholog()], database=database)
        code, _ = self.run_cli(pipeline, *self.arguments("alignment", "--stop-after", "alignment"))
        rows = self.read_rows(self.root / "alignment")
        self.assertEqual(code, 0)
        self.assertEqual(database.calls, 0)
        self.assertEqual(len(rows), 1)
        self.assertTrue(rows[0]["Global Alignment ID"])
        self.assertTrue(all(not rows[0][h] for h in RESULT_HEADERS if h.startswith("Motif ")))
        self.assertIn("intentional-stop", rows[0]["Gene Outcome"])
        self.assertEqual({r["Query Pfam Occurrence"] for r in rows}, {"1"})
        self.assertEqual({r["Hit Pfam Occurrence"] for r in rows}, {"1"})

    def test_alignment_failure_does_not_erase_domain_evidence(self):
        self.write_fasta("alignment_error.fa")
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]}, domains={
            "AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)],
            "CCCC": [DomainRecord("Homeobox", "PF00046.2", 0, 4)]},
            orthologs=[make_ortholog(sequence="CCCC")])
        with patch.object(pipeline, "_calculate_domain_identity", side_effect=RuntimeError("alignment failed")):
            code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertEqual(rows[0]["Hit Pfam ID"], "PF00046.2")
        self.assertEqual(rows[0]["Global Alignment Status"], "failed")
        self.assertEqual(rows[0]["Global Alignment Identity Ratio"], "")
        self.assertEqual(rows[0]["Motif Matrix ID"], "")

    def test_motif_failure_after_partial_generator_preserves_first_motif(self):
        self.write_fasta("partial_motifs.fa")
        database = StubMotifDatabase()
        def partial(request):
            yield make_motif("MA0001.1", "first")
            raise RuntimeError("later motif failed")
        database.search = partial
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAA"]},
            domains={"AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)]},
            orthologs=[make_ortholog()], database=database)
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 1)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["Motif Matrix ID"], "MA0001.1")
        self.assertEqual(rows[1]["Motif Status"], "failed")
        self.assertEqual(rows[1]["Global Alignment ID"], rows[0]["Global Alignment ID"])
        self.assertEqual(rows[1]["BLAST Accession"], rows[0]["BLAST Accession"])

    def test_filtered_query_domain_is_not_lost_when_another_is_accepted(self):
        self.write_fasta("filtered.fa")
        pipeline = self.make_pipeline(orfs={"ACGT": ["AAAAAAAA"]}, domains={"AAAAAAAA": [
            DomainRecord("Other", "PF99999.1", 0, 4), DomainRecord("Homeobox", "PF00046.1", 4, 8)]})
        pipeline.valid_pfam_list = ["PF00046"]
        code, _ = self.run_cli(pipeline, *self.arguments())
        rows = self.read_rows(self.root / "bundle")
        self.assertEqual(code, 0)
        self.assertEqual([r["Query Pfam ID"] for r in rows], ["PF99999.1", "PF00046.1"])
        self.assertEqual(rows[0]["Query Pfam Status"], "not-selected")
        self.assertEqual(rows[0]["BLAST Accession"], "")

    def test_atomic_result_write_preserves_previous_file_on_formatter_error(self):
        from lverage.cli import _Bundle
        destination = self.root / "atomic"
        destination.mkdir()
        bundle = _Bundle(destination, {})
        before = (destination / "results.csv").read_bytes()
        with patch("lverage.cli._grouped_result_rows", side_effect=RuntimeError("format failed")):
            with self.assertRaises(RuntimeError):
                bundle.add_event("gene", {"event": "gene_start", "stage": "input"})
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
