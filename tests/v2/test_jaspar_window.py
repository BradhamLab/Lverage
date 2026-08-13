import unittest

from lverage.domain_scanner import DomainRecord
from lverage.jaspar import Jaspar2024MotifDB


class JasparWindowTests(unittest.TestCase):

    def test_short_sequence_is_unchanged(self):
        sequence = "A" * 1500
        domain = DomainRecord("Homeobox", "PF00046.1", 10, 20)

        self.assertEqual(Jaspar2024MotifDB._window_ortholog_sequence(sequence, domain), sequence)

    def test_window_is_centered_on_domain(self):
        sequence = "A" * 5000
        domain = DomainRecord("Homeobox", "PF00046.1", 2400, 2450)

        window = Jaspar2024MotifDB._window_ortholog_sequence(sequence, domain)

        self.assertEqual(window, sequence[1425:3425])
        self.assertEqual(len(window), 2000)

    def test_window_handles_n_terminal_domain(self):
        sequence = "A" * 5000
        domain = DomainRecord("Homeobox", "PF00046.1", 50, 120)

        self.assertEqual(
            Jaspar2024MotifDB._window_ortholog_sequence(sequence, domain),
            sequence[:2000],
        )

    def test_window_handles_c_terminal_domain(self):
        sequence = "A" * 5000
        domain = DomainRecord("Homeobox", "PF00046.1", 4900, 4980)

        self.assertEqual(
            Jaspar2024MotifDB._window_ortholog_sequence(sequence, domain),
            sequence[-2000:],
        )

    def test_invalid_bounds_raise_value_error(self):
        sequence = "A" * 5000
        invalid_domains = [
            DomainRecord("Homeobox", "PF00046.1", -1, 10),
            DomainRecord("Homeobox", "PF00046.1", 10, 10),
            DomainRecord("Homeobox", "PF00046.1", 0, 5001),
            DomainRecord("Homeobox", "PF00046.1", 0, 2001),
        ]

        for domain in invalid_domains:
            with self.subTest(domain=domain):
                with self.assertRaises(ValueError):
                    Jaspar2024MotifDB._window_ortholog_sequence(sequence, domain)


if __name__ == "__main__":
    unittest.main()
