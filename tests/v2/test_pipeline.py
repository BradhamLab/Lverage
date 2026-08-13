import unittest

from lverage.domain_scanner import DomainRecord, DomainScannerTemplate
from lverage.motif_database import MotifDBRecordTemplate, MotifDBTemplate
from lverage.orf_searcher import OrfSearcherTemplate
from lverage.ortholog_searcher import OrthologRecord, OrthologSearcherTemplate
from lverage.pipeline import Lverage, LverageCode


class ExampleMotifRecord(MotifDBRecordTemplate):

    headers = ("Motif",)

    def __init__(self, motif):
        self.motif = motif

    def get_values(self):
        return [self.motif]


class StubMotifDB(MotifDBTemplate):

    def __init__(self, name="ExampleDB", motifs=None, species_available=True):
        self.name = name
        self.motifs = [] if motifs is None else motifs
        self.species_available = species_available
        self.requests = []

    def search(self, request):
        self.requests.append(request)
        return list(self.motifs)

    def check_species_validity(self, species_tax_id : int):
        return self.species_available


class StubOrfSearcher(OrfSearcherTemplate):

    def __init__(self, orfs):
        self.orfs = orfs

    def get_orfs(self, sequence : str):
        if isinstance(self.orfs, dict):
            return list(self.orfs.get(sequence, []))
        return list(self.orfs)


class StubDomainScanner(DomainScannerTemplate):

    def __init__(self, domains):
        self.domains = domains
        self.scanned = []

    def get_domains(self, sequence : str):
        self.scanned.append(sequence)
        if isinstance(self.domains, dict):
            return list(self.domains.get(sequence, []))
        return list(self.domains)


class StubOrthologSearcher(OrthologSearcherTemplate):

    def __init__(self, orthologs):
        self.orthologs = orthologs

    def get_orthologs(self, sequence : str):
        if isinstance(self.orthologs, Exception):
            raise self.orthologs
        return list(self.orthologs)


def make_ortholog(sequence="AAAA", species_tax_id=9606):
    return OrthologRecord(
        "NP_000001.1",
        "Protein [Homo sapiens]",
        "Homo sapiens",
        species_tax_id,
        sequence,
        1e-20,
        1.0,
        1.0,
    )


class LveragePipelineTests(unittest.TestCase):

    def make_pipeline(self,
                      orfs=None,
                      domains=None,
                      orthologs=None,
                      databases=None,
                      valid_pfam_list=None,
                      threshold=0.7):
        if orfs is None:
            orfs = ["AAAA"]
        if domains is None:
            domains = {
                "AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)],
            }
        if orthologs is None:
            orthologs = [make_ortholog()]
        if databases is None:
            databases = [StubMotifDB(motifs=[ExampleMotifRecord("MA0001.1")])]
        return Lverage(
            databases,
            StubOrfSearcher(orfs),
            StubDomainScanner(domains),
            StubOrthologSearcher(orthologs),
            valid_pfam_list=valid_pfam_list,
            dbd_identity_thresh=threshold,
        )

    def test_no_orf_code(self):
        pipeline = self.make_pipeline(orfs=[])

        self.assertEqual(pipeline.run("QUERY"), [])
        self.assertEqual(pipeline.lverage_code, LverageCode.NO_ORF)

    def test_no_valid_domain_code(self):
        pipeline = self.make_pipeline(domains={"AAAA": []})

        self.assertEqual(pipeline.run("QUERY"), [])
        self.assertEqual(pipeline.lverage_code, LverageCode.NO_VALID_DOMAIN)

    def test_no_orthologs_code(self):
        pipeline = self.make_pipeline(orthologs=[])

        self.assertEqual(pipeline.run("QUERY"), [])
        self.assertEqual(pipeline.lverage_code, LverageCode.NO_ORTHOLOGS)

    def test_no_valid_domain_pair_code(self):
        pipeline = self.make_pipeline(domains={
            "AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)],
        }, orthologs=[make_ortholog("TTTT")])
        pipeline.domain_scanner.domains["TTTT"] = [DomainRecord("Homeobox", "PF00046.2", 0, 4)]

        self.assertEqual(pipeline.run("QUERY"), [])
        self.assertEqual(pipeline.lverage_code, LverageCode.NO_VALID_DBD)

    def test_no_motif_code(self):
        pipeline = self.make_pipeline(databases=[StubMotifDB(motifs=[])])

        self.assertEqual(pipeline.run("QUERY"), [])
        self.assertEqual(pipeline.lverage_code, LverageCode.NO_MOTIF)

    def test_success_returns_nested_records(self):
        pipeline = self.make_pipeline()

        records = pipeline.run("QUERY")

        self.assertEqual(pipeline.lverage_code, LverageCode.SUCCESS)
        self.assertEqual(records[0].motif_database_name, "ExampleDB")
        self.assertEqual(records[0].domain_identity, 1.0)

    def test_identity_threshold_is_inclusive(self):
        pipeline = self.make_pipeline(threshold=1.0)

        self.assertEqual(len(pipeline.run("QUERY")), 1)

    def test_longest_distinct_orf_with_an_accepted_domain_is_selected(self):
        pipeline = self.make_pipeline(
            orfs={"first": ["AAAA", "AAAAAAAA", "AAAAAAAA"], "second": ["AAAAAA"]},
            domains={
                "AAAAAAAA": [],
                "AAAAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)],
                "AAAA": [DomainRecord("Homeobox", "PF00046.1", 0, 4)],
            },
            orthologs=[make_ortholog()],
        )

        pipeline.run(["first", "second"])

        self.assertEqual(pipeline.domain_scanner.scanned[:2], ["AAAAAAAA", "AAAAAA"])
        self.assertEqual(pipeline.motif_database_list[0].requests[0].query_sequence, "AAAAAA")

    def test_domains_match_by_base_pfam_accession(self):
        pipeline = self.make_pipeline(domains={
            "AAAA": [DomainRecord("Homeobox", "PF00046.32", 0, 4)],
            "AAAAT": [DomainRecord("Homeobox", "PF00046.2", 0, 4)],
        }, orthologs=[make_ortholog("AAAAT")])

        self.assertEqual(len(pipeline.run("QUERY")), 1)

    def test_multiple_matching_domain_pairs_are_evaluated(self):
        query_sequence = "AAAAAAAA"
        ortholog = make_ortholog(query_sequence)
        domains = [
            DomainRecord("First", "PF00046.1", 0, 4),
            DomainRecord("Second", "PF00046.2", 4, 8),
        ]
        pipeline = self.make_pipeline(orfs=[query_sequence], domains={query_sequence: domains}, orthologs=[ortholog])

        self.assertEqual(len(pipeline.run("QUERY")), 4)

    def test_unavailable_species_is_not_searched(self):
        database = StubMotifDB(motifs=[ExampleMotifRecord("MA0001.1")], species_available=False)
        pipeline = self.make_pipeline(databases=[database])

        self.assertEqual(pipeline.run("QUERY"), [])
        self.assertEqual(database.requests, [])
        self.assertEqual(pipeline.lverage_code, LverageCode.NO_MOTIF)

    def test_multiple_databases_produce_separate_records(self):
        databases = [
            StubMotifDB("First", [ExampleMotifRecord("MA0001.1")]),
            StubMotifDB("Second", [ExampleMotifRecord("MA0002.1")]),
        ]
        pipeline = self.make_pipeline(databases=databases)

        records = pipeline.run("QUERY")

        self.assertEqual([record.motif_database_name for record in records], ["First", "Second"])

    def test_repeated_runs_reset_the_code_and_keep_no_transient_evidence(self):
        pipeline = self.make_pipeline()
        pipeline.run("QUERY")
        pipeline.orf_searcher.orfs = []

        self.assertEqual(pipeline.run("QUERY"), [])
        self.assertEqual(pipeline.lverage_code, LverageCode.NO_ORF)
        self.assertFalse(hasattr(pipeline, "orf"))
        self.assertFalse(hasattr(pipeline, "orthologs"))

    def test_adapter_errors_propagate(self):
        pipeline = self.make_pipeline(orthologs=RuntimeError("service failure"))

        with self.assertRaises(RuntimeError):
            pipeline.run("QUERY")
        self.assertEqual(pipeline.lverage_code, LverageCode.NOT_SET)

    def test_run_requires_nonempty_sequences(self):
        pipeline = self.make_pipeline()

        for invalid in ["", "   ", []]:
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    pipeline.run(invalid)
        with self.assertRaises(TypeError):
            pipeline.run(["QUERY", 1])
        with self.assertRaises(TypeError):
            pipeline.run(1)


if __name__ == "__main__":
    unittest.main()
