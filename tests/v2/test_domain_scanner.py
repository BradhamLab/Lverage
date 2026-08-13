import unittest

from lverage.domain_scanner import DomainRecord, DomainScannerTemplate


class IncompleteDomainScanner(DomainScannerTemplate):
    pass


class ConcreteDomainScanner(DomainScannerTemplate):

    def get_domains(self, sequence : str) -> list[DomainRecord]:
        return [DomainRecord("Homeobox", "PF00046.1", 0, len(sequence))]


class DomainScannerTemplateTests(unittest.TestCase):

    def test_template_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            DomainScannerTemplate()

    def test_incomplete_subclass_cannot_be_instantiated(self):
        with self.assertRaises(TypeError):
            IncompleteDomainScanner()

    def test_concrete_subclass_returns_domains(self):
        scanner = ConcreteDomainScanner()

        self.assertEqual(
            scanner.get_domains("MPEPTIDE"),
            [DomainRecord("Homeobox", "PF00046.1", 0, 8)]
        )


if __name__ == "__main__":
    unittest.main()
