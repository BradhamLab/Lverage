import unittest

from lverage.motif_database import MotifDBRecordTemplate


class IncompleteMotifRecord(MotifDBRecordTemplate):
    pass


class ConcreteMotifRecord(MotifDBRecordTemplate):

    def get_values(self):
        return ["MA0001.1"]


class MotifDBRecordTemplateTests(unittest.TestCase):

    def test_template_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            MotifDBRecordTemplate()

    def test_incomplete_subclass_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            IncompleteMotifRecord()

    def test_concrete_subclass_serializes_values(self):
        record = ConcreteMotifRecord()

        self.assertEqual(record.get_values(), ["MA0001.1"])


if __name__ == "__main__":
    unittest.main()
