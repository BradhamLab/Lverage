import unittest

from lverage.domain_scanner import DomainRecord
from lverage.pipeline import Lverage


class DomainIdentityTests(unittest.TestCase):

    def test_identical_domains_have_complete_identity(self):
        identity = Lverage._calculate_domain_identity(
            "XXAAAAZZ",
            DomainRecord("Domain", "PF00001", 2, 6),
            "YYAAAAQQ",
            DomainRecord("Domain", "PF00001", 2, 6),
        )

        self.assertEqual(identity, 1.0)

    def test_mismatches_count_as_alignment_columns(self):
        identity = Lverage._calculate_domain_identity(
            "AAAA",
            DomainRecord("Domain", "PF00001", 0, 4),
            "AAAT",
            DomainRecord("Domain", "PF00001", 0, 4),
        )

        self.assertEqual(identity, 0.75)

    def test_gaps_count_as_alignment_columns(self):
        identity = Lverage._calculate_domain_identity(
            "AAAA",
            DomainRecord("Domain", "PF00001", 0, 4),
            "AAA",
            DomainRecord("Domain", "PF00001", 0, 3),
        )

        self.assertEqual(identity, 0.75)

    def test_invalid_domain_bounds_raise_value_error(self):
        invalid_domains = [
            DomainRecord("Domain", "PF00001", -1, 2),
            DomainRecord("Domain", "PF00001", 2, 2),
            DomainRecord("Domain", "PF00001", 0, 5),
        ]

        for domain in invalid_domains:
            with self.subTest(domain=domain):
                with self.assertRaises(ValueError):
                    Lverage._calculate_domain_identity(
                        "AAAA",
                        domain,
                        "AAAA",
                        DomainRecord("Domain", "PF00001", 0, 4),
                    )


if __name__ == "__main__":
    unittest.main()
