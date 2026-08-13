import unittest

from lverage.motif_database import MotifDBRecordTemplate


class IncompleteMotifRecord(MotifDBRecordTemplate):
    pass


class ConcreteMotifRecord(MotifDBRecordTemplate):

    headers = ("Matrix ID",)

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

    def test_headers_are_available_from_record_class(self):
        self.assertEqual(ConcreteMotifRecord.get_headers(), ["Matrix ID"])

    def test_get_headers_returns_a_list_copy(self):
        headers = ConcreteMotifRecord.get_headers()
        headers.append("Name")

        self.assertEqual(ConcreteMotifRecord.headers, ("Matrix ID",))
        self.assertEqual(ConcreteMotifRecord.get_headers(), ["Matrix ID"])


if __name__ == "__main__":
    unittest.main()
