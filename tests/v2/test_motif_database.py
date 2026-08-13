import unittest

from lverage.motif_database import MotifDBTemplate


class IncompleteMotifDB(MotifDBTemplate):
    pass


class ConcreteMotifDB(MotifDBTemplate):

    def search(self, **kwargs):
        return []

    def check_species_validity(self, species_tax_id : int):
        return species_tax_id == 9606


class MotifDBTemplateTests(unittest.TestCase):

    def test_template_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            MotifDBTemplate()

    def test_incomplete_subclass_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            IncompleteMotifDB()

    def test_concrete_subclass_implements_database_interface(self):
        motif_database = ConcreteMotifDB()

        self.assertEqual(motif_database.search(sequence="MPEPTIDE"), [])
        self.assertTrue(motif_database.check_species_validity(9606))
        self.assertFalse(motif_database.check_species_validity(10090))


if __name__ == "__main__":
    unittest.main()
