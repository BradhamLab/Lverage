import unittest

from lverage.alignment import BiopythonGlobalAlignment


class BiopythonGlobalAlignmentTests(unittest.TestCase):

    def test_exact_identity_counts_mismatches_and_gap_columns(self):
        aligner = BiopythonGlobalAlignment()

        self.assertEqual(aligner.calc_identities("AAAA", ["AAAT"]), [0.75])
        self.assertEqual(aligner.calc_identities("AAAA", ["AAA"]), [0.75])

    def test_internal_and_terminal_gaps_count_in_denominator(self):
        aligner = BiopythonGlobalAlignment()

        self.assertEqual(aligner.calc_identities("ACGT", ["ACT"]), [0.75])
        self.assertEqual(aligner.calc_identities("ACGT", ["CGT"]), [0.75])

    def test_duplicate_sequences_return_distinct_ordered_results(self):
        identities = BiopythonGlobalAlignment().calc_identities(
            "AAAA",
            ["AAAT", "AAAT", "AAAA"],
        )

        self.assertEqual(identities, [0.75, 0.75, 1.0])

    def test_default_scores_select_affine_global_algorithm(self):
        aligner = BiopythonGlobalAlignment()
        configured = aligner._make_aligner()

        self.assertEqual(configured.mode, "global")
        self.assertEqual(configured.match_score, 2)
        self.assertEqual(configured.mismatch_score, -1)
        self.assertEqual(configured.open_gap_score, -2)
        self.assertEqual(configured.extend_gap_score, -0.5)
        self.assertEqual(configured.algorithm, "Gotoh global alignment algorithm")

    def test_linear_gap_configuration_selects_needleman_wunsch(self):
        aligner = BiopythonGlobalAlignment(
            gap_model="linear",
            open_gap_score=-2,
            extend_gap_score=-2,
        )

        self.assertEqual(aligner._make_aligner().algorithm, "Needleman-Wunsch")

    def test_linear_gap_configuration_requires_equal_gap_scores(self):
        with self.assertRaises(ValueError):
            BiopythonGlobalAlignment(gap_model="linear")

    def test_nonfinite_and_non_numeric_scores_are_rejected(self):
        with self.assertRaises(ValueError):
            BiopythonGlobalAlignment(match_score=float("inf"))
        with self.assertRaises(TypeError):
            BiopythonGlobalAlignment(match_score=True)

    def test_empty_sequences_are_rejected(self):
        aligner = BiopythonGlobalAlignment()

        with self.assertRaises(ValueError):
            aligner.calc_identities("", ["AAAA"])
        with self.assertRaises(ValueError):
            aligner.calc_identities("AAAA", [""])


if __name__ == "__main__":
    unittest.main()