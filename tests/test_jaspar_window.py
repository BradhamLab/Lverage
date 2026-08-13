import unittest
from unittest.mock import patch
import os
import sys
import types

TEST_DIR = os.path.dirname(__file__)
SRC_DIR = os.path.join(os.path.dirname(TEST_DIR), "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

if "ete3" not in sys.modules:
    ete3_stub = types.ModuleType("ete3")

    class _DummyNCBITaxa:
        pass

    ete3_stub.NCBITaxa = _DummyNCBITaxa
    sys.modules["ete3"] = ete3_stub

if "Bio" not in sys.modules:
    bio_stub = types.ModuleType("Bio")
    bio_align_stub = types.ModuleType("Bio.Align")

    class _DummyPairwiseAligner:
        def __init__(self):
            self.mode = None
            self.match_score = None
            self.mismatch_score = None

    bio_align_stub.PairwiseAligner = _DummyPairwiseAligner
    bio_stub.Align = bio_align_stub
    sys.modules["Bio"] = bio_stub
    sys.modules["Bio.Align"] = bio_align_stub

from DBDScanner import DBD
from MotifDB import JasparDB


class FakeResponse:
    def __init__(self, ok=True, payload=None):
        self.ok = ok
        self._payload = payload if payload is not None else {"results": []}

    def json(self):
        return self._payload


class JasparWindowTests(unittest.TestCase):
    def test_window_short_sequence_unchanged(self):
        dbd = DBD("Homeobox", 10, 20, "PF00046.1")
        seq = "A" * 1500
        self.assertEqual(JasparDB._window_ortholog_sequence(seq, dbd), seq)

    def test_window_centered_on_ortholog_dbd(self):
        seq = "A" * 5000
        dbd = DBD("Homeobox", 2400, 2450, "PF00046.1")
        window = JasparDB._window_ortholog_sequence(seq, dbd)

        self.assertLessEqual(len(window), 2000)
        self.assertEqual(len(window), 2000)
        self.assertEqual(window, seq[1425:3425])
        self.assertIn(seq[2400:2450], window)

    def test_window_handles_n_terminal_boundary(self):
        seq = "A" * 5000
        dbd = DBD("Homeobox", 50, 120, "PF00046.1")
        window = JasparDB._window_ortholog_sequence(seq, dbd)

        self.assertEqual(window, seq[:2000])
        self.assertIn(seq[50:120], window)

    def test_window_handles_c_terminal_boundary(self):
        seq = "A" * 5000
        dbd = DBD("Homeobox", 4900, 4980, "PF00046.1")
        window = JasparDB._window_ortholog_sequence(seq, dbd)

        self.assertEqual(window, seq[-2000:])
        self.assertIn(seq[4900:4980], window)

    @patch("MotifDB.requests.get")
    def test_search_sends_ortholog_sequence_not_gene_sequence(self, mock_get):
        observed_urls = []

        def fake_get(url, *args, **kwargs):
            observed_urls.append(url)
            return FakeResponse(ok=True, payload={"results": []})

        mock_get.side_effect = fake_get

        ortholog_seq = "ORTHOSEQ"
        gene_seq = "GENESEQ"
        gene_dbd = DBD("Homeobox", 1, 5, "PF00046.1")
        ortholog_dbd = DBD("Homeobox", 2, 6, "PF00046.1")

        mdb = JasparDB(dbd_scanner=object())
        mdb.search(ortholog_seq, "9606", gene_seq, gene_dbd, ortholog_dbd)

        self.assertTrue(observed_urls)
        self.assertEqual(
            observed_urls[0],
            f"{mdb.jaspar_rest_url}/infer/{ortholog_seq}"
        )

    @patch("MotifDB.requests.get")
    def test_search_uses_ortholog_dbd_for_window_not_gene_dbd(self, mock_get):
        observed_urls = []

        def fake_get(url, *args, **kwargs):
            observed_urls.append(url)
            return FakeResponse(ok=True, payload={"results": []})

        mock_get.side_effect = fake_get

        ortholog_seq = "A" * 5000
        gene_seq = "G" * 500
        gene_dbd = DBD("Homeobox", 100, 150, "PF00046.1")
        ortholog_dbd = DBD("Homeobox", 3200, 3250, "PF00046.1")

        expected_window = JasparDB._window_ortholog_sequence(ortholog_seq, ortholog_dbd)

        mdb = JasparDB(dbd_scanner=object())
        mdb.search(ortholog_seq, "9606", gene_seq, gene_dbd, ortholog_dbd)

        self.assertTrue(observed_urls)
        self.assertEqual(
            observed_urls[0],
            f"{mdb.jaspar_rest_url}/infer/{expected_window}"
        )


if __name__ == "__main__":
    unittest.main()
