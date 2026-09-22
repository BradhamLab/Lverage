import subprocess
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from lverage.domain_scanner import DomainRecord
from lverage.pfam import LocalPfamScanner


class LocalPfamScannerTests(unittest.TestCase):

    def setUp(self):
        self.which_patcher = patch(
            "lverage.pfam.shutil.which",
            side_effect=lambda path: f"/resolved/{path}",
        )
        self.run_patcher = patch("lverage.pfam.subprocess.run")
        self.which_patcher.start()
        self.mock_run = self.run_patcher.start()
        self.addCleanup(self.which_patcher.stop)
        self.addCleanup(self.run_patcher.stop)

    def test_constructor_resolves_executable_and_copies_database(self):
        scanner = LocalPfamScanner("pfamdb", "custom-pfam_scan.pl")

        self.assertEqual(scanner.database_path, "pfamdb")
        self.assertEqual(scanner.pfamscan_path, "/resolved/custom-pfam_scan.pl")
        self.assertIsInstance(scanner, LocalPfamScanner)

    def test_scan_builds_command_and_writes_temporary_fasta(self):
        scanner = LocalPfamScanner("pfamdb")
        captured = {}

        def run_scan(command, **kwargs):
            captured["path"] = command[2]
            captured["contents"] = Path(command[2]).read_text()
            return SimpleNamespace(stdout=(
                "# comment\n"
                "query 1 8 2 7 PF00046.1 Homeobox Domain 1 6 60\n"
            ))

        self.mock_run.side_effect = run_scan

        domains = scanner.get_domains("MPEPTIDE")

        self.assertEqual(domains, [DomainRecord("Homeobox", "PF00046.1", 1, 7)])
        self.assertEqual(captured["contents"], ">query\nMPEPTIDE\n")
        self.assertFalse(Path(captured["path"]).exists())
        command = self.mock_run.call_args.args[0]
        self.assertEqual(command[0], "/resolved/pfam_scan.pl")
        self.assertEqual(command[1], "-fasta")
        self.assertTrue(command[2].endswith(".fasta"))
        self.assertEqual(command[3:], ["-dir", "pfamdb"])

    def test_empty_scan_returns_no_domains(self):
        scanner = LocalPfamScanner("pfamdb")
        self.mock_run.return_value = SimpleNamespace(stdout="# no hits\n")

        self.assertEqual(scanner.get_domains("QUERY"), [])

    def test_scan_process_failure_propagates(self):
        scanner = LocalPfamScanner("pfamdb")
        self.mock_run.side_effect = subprocess.CalledProcessError(1, "pfam_scan.pl")

        with self.assertRaises(subprocess.CalledProcessError):
            scanner.get_domains("QUERY")

    def test_empty_sequence_is_rejected(self):
        scanner = LocalPfamScanner("pfamdb")

        with self.assertRaises(ValueError):
            scanner.get_domains(" ")


if __name__ == "__main__":
    unittest.main()