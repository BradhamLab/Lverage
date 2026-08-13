import unittest

from lverage.ortholog_searcher import OrthologRecord, OrthologSearcherTemplate


class IncompleteOrthologSearcher(OrthologSearcherTemplate):
    pass


class ConcreteOrthologSearcher(OrthologSearcherTemplate):

    def get_orthologs(self, sequence : str) -> list[OrthologRecord]:
        return []


class OrthologSearcherTemplateTests(unittest.TestCase):

    def test_template_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            OrthologSearcherTemplate()

    def test_incomplete_subclass_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            IncompleteOrthologSearcher()

    def test_concrete_subclass_returns_orthologs(self):
        self.assertEqual(ConcreteOrthologSearcher().get_orthologs("PROTEIN"), [])


if __name__ == "__main__":
    unittest.main()
