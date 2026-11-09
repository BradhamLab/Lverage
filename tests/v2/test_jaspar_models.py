import unittest

import numpy as np

from lverage.jaspar_models import score_alignment_result


class JasparSimilarityRegressionTests(unittest.TestCase):

    def test_upstream_identical_domain_rule_accepts_exact_identity(self):
        model = {
            "Model.Class": "SequenceIdentity",
            "Model.Name": "PctID_L",
            "Threshold.HSim": 0.9,
            "Threshold.Dis": 0.5,
        }
        score, classification = score_alignment_result({"PctID_L": 1.0}, model)
        self.assertEqual(score, 1.0)
        self.assertEqual(classification, "HSim")

    def test_upstream_sequence_identity_threshold_is_inclusive(self):
        model = {
            "Model.Class": "SequenceIdentity",
            "Model.Name": "PctID_L",
            "Threshold.HSim": 0.8,
            "Threshold.Dis": 0.5,
        }
        result = {"PctID_L": 0.8}
        self.assertEqual(score_alignment_result(result, model), (0.8, "HSim"))
        score, classification = score_alignment_result({"PctID_L": 0.79}, model)
        self.assertTrue(np.isnan(score))
        self.assertTrue(np.isnan(classification))

    def test_upstream_similarity_regression_baseline_then_scaled_score(self):
        model = {
            "Model.Class": "SimilarityRegression",
            "Baseline": {"Name": "PctID_L", "Threshold.HSim": 1.1},
            "SR.Features": "PctID",
            "SR.FeatureScales.mean": np.asarray([0.0, 0.0]),
            "SR.FeatureScales.sd": np.asarray([1.0, 1.0]),
            "SR.Weights": np.asarray([0.25, 0.25]),
            "SR.Intercept": 0.0,
            "SR.LogisticTransform": False,
            "Threshold.HSim": 0.4,
            "Threshold.Dis": None,
        }
        result = {"PctID_L": 0.5, "ByPos.PctID": np.asarray([1.0, 1.0])}
        self.assertEqual(score_alignment_result(result, model), (0.5, "HSim"))


if __name__ == "__main__":
    unittest.main()
