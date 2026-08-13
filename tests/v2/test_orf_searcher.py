import unittest

from lverage.orf_searcher import OrfSearcherTemplate


class IncompleteOrfSearcher(OrfSearcherTemplate):
    pass


class ConcreteOrfSearcher(OrfSearcherTemplate):

    def get_orfs(self, sequence : str) -> list[str]:
        return [sequence]


class OrfSearcherTemplateTests(unittest.TestCase):

    def test_template_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            OrfSearcherTemplate()

    def test_incomplete_subclass_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            IncompleteOrfSearcher()

    def test_concrete_subclass_returns_orfs(self):
        searcher = ConcreteOrfSearcher()

        self.assertEqual(searcher.get_orfs("ATGAAATAG"), ["ATGAAATAG"])


if __name__ == "__main__":
    unittest.main()
