import unittest
import math
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from lverage.jaspar_inference import (
    _rost_identity_cutoff,
    blast_search,
    infer_sequence,
    non_overlapping_domains,
    positional_alignment_score,
)


class JasparInferenceAlgorithmTests(unittest.TestCase):

    def test_upstream_non_overlap_rule_keeps_lowest_conditional_evalue(self):
        domains = [
            ("Homeobox", 0, 20, 1e-4),
            ("Forkhead", 10, 25, 1e-6),
            ("bZIP", 30, 45, 1e-3),
        ]
        self.assertEqual(
            non_overlapping_domains(domains),
            [("Forkhead", 10, 25, 1e-6), ("bZIP", 30, 45, 1e-3)],
        )

    def test_positional_features_match_upstream_identity_and_blosum_rules(self):
        aligned = ["A-C"]
        np.testing.assert_array_equal(
            positional_alignment_score(aligned, aligned, "identity"),
            np.array([1.0, 1.0, 1.0]),
        )
        np.testing.assert_array_equal(
            positional_alignment_score(aligned, aligned, "blosum62"),
            np.array([4.0, 1.0, 9.0]),
        )

    def test_rost_cutoff_matches_source_formula(self):
        length = 200
        expected = 5 + 480 * length ** (-0.32 * (1 + math.exp(-length / 1000)))
        self.assertAlmostEqual(_rost_identity_cutoff(length), expected)

    def test_blast_search_uses_local_db_and_upstream_outfmt_as_argv(self):
        record = SeqRecord(Seq("A" * 500), id="query")
        row = "target\t100\t500\t1\t500\t1\t500\t1e-20\t1000\t100\t500\t500\n"
        with patch(
            "lverage.jaspar_inference._checked_run",
            return_value=SimpleNamespace(stdout=row),
        ) as run:
            hits = blast_search(
                record, "/local/jaspar/files", "/tmp/query.fasta",
                blastp_path="/tools/blastp", taxons=("vertebrates",),
            )
        command = run.call_args.args[0]
        self.assertEqual(command[0], "/tools/blastp")
        self.assertIn("/local/jaspar/files/vertebrates.fa", command)
        self.assertEqual(command[command.index("-query") + 1], "/tmp/query.fasta")
        self.assertEqual(command[command.index("-outfmt") + 1],
                         "6 sseqid pident length qstart qend sstart send evalue bitscore ppos qlen slen")
        self.assertEqual(len(hits), 1)

    def test_sequence_inference_preserves_upstream_hsim_gate_and_versions(self):
        query_domains = [("Homeobox", "AC-G", 1, 4, 1e-9)]
        blast_hits = [("lverage_query", "target", "3-9", "2-8", 1e-12, 100.0,
                       90.0, 7, 95.0, 90.0)]
        data = {
            "target": {
                "pfam": [("Homeobox", "AC-G", 1, 4, 1e-9)],
                "profiles": [["MA0001.1", "TF"], ["MA0001.2", "TF"]],
            },
        }
        models = {None: {"fallback": True}}
        with patch("lverage.jaspar_inference.scan_and_align_domains",
                   return_value=query_domains), \
                patch("lverage.jaspar_inference.blast_search",
                      return_value=blast_hits), \
                patch("lverage.jaspar_inference.score_alignment_result",
                      return_value=(0.91, "HSim")) as score:
            predictions = infer_sequence(
                "ACGT", data_dir="/data", cisbp_models=models,
                jaspar_data=data, work_dir="/tmp/work", blastp_path="blastp",
                hmmscan_path="hmmscan", hmmalign_path="hmmalign",
            )
        self.assertEqual([row["matrix_id"] for row in predictions],
                         ["MA0001.1", "MA0001.2"])
        self.assertEqual(predictions[0]["evalue"], 1e-12)
        self.assertEqual(predictions[0]["dbd_identity"], 1.0)
        self.assertEqual(predictions[0]["query"], "lverage_query")
        self.assertEqual(score.call_args.args[1], models[None])


if __name__ == "__main__":
    unittest.main()
