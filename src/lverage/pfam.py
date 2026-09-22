"""Pfam-based domain scanners."""

import logging
import shutil
import subprocess
import tempfile
import time

import requests

from .domain_scanner import DomainRecord, DomainScannerTemplate


LOGGER = logging.getLogger(__name__)
PFAMSCAN_URL = "https://www.ebi.ac.uk/Tools/services/rest/pfamscan/"
PFAMSCAN_RUN_URL = PFAMSCAN_URL + "run/"
PFAMSCAN_STATUS_URL = PFAMSCAN_URL + "status/"
PFAMSCAN_RESULT_TYPES_URL = PFAMSCAN_URL + "resulttypes/"
PFAMSCAN_RESULT_URL = PFAMSCAN_URL + "result/"
FAILED_STATUSES = {"ERROR", "FAILURE", "NOT_FOUND"}


def _resolve_executable(executable : str) -> str:
    resolved_executable = shutil.which(executable)
    if resolved_executable is None:
        raise FileNotFoundError(f"PfamScan executable not found: {executable}")
    return resolved_executable


def _validate_sequence(sequence : str) -> None:
    if not isinstance(sequence, str) or not sequence.strip():
        raise ValueError("sequence must be a nonempty string")


def _make_domain(result) -> DomainRecord:
    try:
        name = result["name"]
        accession = result["acc"]
        start = int(result["env"]["from"]) - 1
        end = int(result["env"]["to"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("invalid PfamScan domain result") from error
    if not isinstance(name, str) or not name:
        raise ValueError("invalid PfamScan domain name")
    if not isinstance(accession, str) or not accession:
        raise ValueError("invalid PfamScan domain accession")
    if start < 0 or end <= start:
        raise ValueError("invalid PfamScan domain coordinates")
    return DomainRecord(name, accession, start, end)


def _parse_online_results(results) -> list[DomainRecord]:
    if not isinstance(results, list):
        raise ValueError("PfamScan result must be a list")
    domains = []
    for result in results:
        try:
            domains.append(_make_domain(result))
        except ValueError:
            LOGGER.warning("Skipping malformed PfamScan domain", exc_info=True)
    return domains


def _parse_local_results(output) -> list[DomainRecord]:
    domains = []
    for line in output.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 7:
            LOGGER.warning("Skipping malformed local PfamScan result: %s", line)
            continue
        try:
            domains.append(_make_domain({
                "name": fields[6],
                "acc": fields[5],
                "env": {"from": fields[3], "to": fields[4]},
            }))
        except ValueError:
            LOGGER.warning("Skipping malformed local PfamScan result: %s", line)
    return domains


class LocalPfamScanner(DomainScannerTemplate):
    """
    Scan protein sequences with a local PfamScan executable.

    Parameters
    ----------
    database_path : str
        PfamScan database directory
    pfamscan_path : str, optional
        Executable name or path for ``pfam_scan.pl``
    """

    def __init__(self, database_path, pfamscan_path="pfam_scan.pl"):
        if not isinstance(database_path, str) or not database_path.strip():
            raise ValueError("database_path must be a nonempty string")

        self.database_path = database_path
        self.pfamscan_path = _resolve_executable(pfamscan_path)

    def get_domains(self, sequence : str) -> list[DomainRecord]:
        """
        Scan a protein sequence with the local PfamScan database.

        Parameters
        ----------
        sequence : str
            Protein sequence to scan

        Returns
        -------
        list
            Domain records parsed from PfamScan output
        """

        _validate_sequence(sequence)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta") as query_file:
            query_file.write(f">query\n{sequence}\n")
            query_file.flush()
            result = subprocess.run(
                [
                    self.pfamscan_path,
                    "-fasta", query_file.name,
                    "-dir", self.database_path,
                ],
                capture_output=True,
                check=True,
                text=True,
            )
        return _parse_local_results(result.stdout)


class OnlinePfamScanner(DomainScannerTemplate):
    """
    Scan protein sequences with the EMBL-EBI online PfamScan service.

    Parameters
    ----------
    email : str
        Contact email supplied to the service
    try_count : int, optional
        Maximum number of status requests
    time_interval : float, optional
        Seconds to wait between status requests
    request_timeout_seconds : float, optional
        Maximum duration of each HTTP request
    """

    def __init__(self,
                 email,
                 try_count=100,
                 time_interval=5,
                 request_timeout_seconds=30):
        if not isinstance(email, str) or not email.strip() or "@" not in email:
            raise ValueError("email must be a nonempty email address")
        if isinstance(try_count, bool) or not isinstance(try_count, int) or try_count <= 0:
            raise ValueError("try_count must be a positive integer")
        if not isinstance(time_interval, (int, float)) or time_interval < 0:
            raise ValueError("time_interval must be nonnegative")
        if not isinstance(request_timeout_seconds, (int, float)) or request_timeout_seconds <= 0:
            raise ValueError("request_timeout_seconds must be positive")

        self.email = email
        self.try_count = try_count
        self.time_interval = time_interval
        self.request_timeout_seconds = request_timeout_seconds

    def get_domains(self, sequence : str) -> list[DomainRecord]:
        """
        Submit a protein sequence to the online PfamScan service.

        Parameters
        ----------
        sequence : str
            Protein sequence to scan

        Returns
        -------
        list
            Domain records returned by PfamScan
        """

        _validate_sequence(sequence)
        response = requests.post(
            PFAMSCAN_RUN_URL,
            data={
                "email": self.email,
                "sequence": sequence,
                "format": "json",
            },
            timeout=self.request_timeout_seconds,
        )
        response.raise_for_status()
        job_id = response.text.strip()
        if not job_id:
            raise ValueError("PfamScan response did not contain a job identifier")

        for attempt in range(self.try_count):
            status_response = requests.get(
                PFAMSCAN_STATUS_URL + job_id,
                timeout=self.request_timeout_seconds,
            )
            status_response.raise_for_status()
            status = status_response.text.strip().upper()
            if status == "QUEUED" or status == "RUNNING":
                if attempt + 1 < self.try_count:
                    time.sleep(self.time_interval)
                continue
            if status == "FINISHED":
                break
            if status in FAILED_STATUSES:
                raise RuntimeError(f"PfamScan job failed with status: {status}")
            raise RuntimeError(f"Unknown PfamScan job status: {status}")
        else:
            raise TimeoutError("PfamScan job did not finish within try_count")

        result_types_response = requests.get(
            PFAMSCAN_RESULT_TYPES_URL + job_id,
            timeout=self.request_timeout_seconds,
        )
        result_types_response.raise_for_status()
        result_types = result_types_response.json()
        json_result_type = next(
            (
                result_type for result_type in result_types
                if result_type.get("mediaType", "").casefold() == "application/json"
            ),
            None,
        )
        if json_result_type is None:
            raise RuntimeError("PfamScan did not provide a JSON result type")

        result_response = requests.get(
            PFAMSCAN_RESULT_URL + job_id + "/" + json_result_type["identifier"],
            timeout=self.request_timeout_seconds,
        )
        result_response.raise_for_status()
        return _parse_online_results(result_response.json())