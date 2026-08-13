"""BLAST-based ortholog searchers."""

from collections.abc import Mapping, Sequence
from io import StringIO
import logging
import re
import shutil
import subprocess

from Bio.Blast import NCBIXML
from Bio import SeqIO
import requests

from .ortholog_searcher import OrthologRecord, OrthologSearcherTemplate


LOGGER = logging.getLogger(__name__)
DEFAULT_EXCLUDED_TERMS = (
    "hypothetical",
    "unnamed",
    "uncharacterized",
    "unknown",
    "partial",
    "isoform",
)
NCBI_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def _resolve_executable(executable : str) -> str:
    resolved_executable = shutil.which(executable)
    if resolved_executable is None:
        raise FileNotFoundError(f"BLAST+ executable not found: {executable}")
    return resolved_executable


def _copy_species_map(species_map : Mapping[str, int]) -> dict[str, int]:
    if not isinstance(species_map, Mapping) or not species_map:
        raise ValueError("species_map must be a nonempty mapping")

    copied_map = {}
    for species_name, species_tax_id in species_map.items():
        if not isinstance(species_name, str) or not species_name.strip():
            raise ValueError("species_map names must be nonempty strings")
        if isinstance(species_tax_id, bool) or not isinstance(species_tax_id, int) or species_tax_id <= 0:
            raise ValueError("species_map taxonomic identifiers must be positive integers")
        copied_map[species_name] = species_tax_id
    return copied_map


def _copy_excluded_terms(excluded_terms : Sequence[str] | None) -> list[str]:
    terms = DEFAULT_EXCLUDED_TERMS if excluded_terms is None else excluded_terms
    if isinstance(terms, str) or not isinstance(terms, Sequence):
        raise TypeError("excluded_terms must be a sequence of strings")
    if any(not isinstance(term, str) for term in terms):
        raise TypeError("excluded_terms must contain only strings")
    return list(terms)


def _parse_blast_hits(blast_xml, query_length, species_lookup, excluded_terms):
    blast_records = list(NCBIXML.parse(StringIO(blast_xml)))
    if len(blast_records) != 1:
        raise ValueError("BLAST output did not contain exactly one query result")

    hits = []
    for alignment in blast_records[0].alignments:
        try:
            description = alignment.hit_def
            if any(term.casefold() in description.casefold() for term in excluded_terms):
                LOGGER.warning("Skipping excluded BLAST hit: %s", description)
                continue

            species_match = re.search(r"\[([^\]]+)\]", description)
            if species_match is None:
                LOGGER.warning("Skipping BLAST hit without a bracketed species: %s", description)
                continue
            species = species_lookup.get(species_match.group(1).strip().casefold())
            if species is None:
                LOGGER.warning("Skipping BLAST hit with an unresolved species: %s", description)
                continue

            best_hsp = min(alignment.hsps, key=lambda hsp: hsp.expect)
            hits.append({
                "accession": alignment.accession,
                "description": description,
                "hit_id": alignment.hit_id,
                "species_name": species[0],
                "species_tax_id": species[1],
                "evalue": best_hsp.expect,
                "identity": best_hsp.identities / best_hsp.align_length,
                "query_coverage": best_hsp.align_length / query_length,
            })
        except (AttributeError, TypeError, ValueError, ZeroDivisionError):
            LOGGER.warning("Skipping malformed BLAST hit", exc_info=True)
    return hits


def _make_ortholog(hit, sequence):
    return OrthologRecord(
        accession=hit["accession"],
        description=hit["description"],
        species_name=hit["species_name"],
        species_tax_id=hit["species_tax_id"],
        sequence=sequence,
        evalue=hit["evalue"],
        identity=hit["identity"],
        query_coverage=hit["query_coverage"],
    )


class LocalBlastSearcher(OrthologSearcherTemplate):
    """
    Search a local BLAST protein database for orthologs.

    Parameters
    ----------
    database_path : str
        BLAST protein database prefix
    species_map : mapping
        Scientific species names mapped to positive NCBI taxonomic identifiers
    blastp_path : str, optional
        Executable name or path for ``blastp``
    blastdbcmd_path : str, optional
        Executable name or path for ``blastdbcmd``
    evalue_threshold : float, optional
        Maximum BLAST expectation value
    top_n : int, optional
        Maximum number of BLAST target sequences
    excluded_terms : sequence, optional
        Case-insensitive description terms to exclude; an empty sequence disables filtering
    """

    def __init__(self,
                 database_path,
                 species_map,
                 blastp_path="blastp",
                 blastdbcmd_path="blastdbcmd",
                 evalue_threshold=1e-6,
                 top_n=20,
                 excluded_terms=None):
        if not isinstance(database_path, str) or not database_path.strip():
            raise ValueError("database_path must be a nonempty string")
        if not isinstance(evalue_threshold, (int, float)) or evalue_threshold <= 0:
            raise ValueError("evalue_threshold must be positive")
        if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n <= 0:
            raise ValueError("top_n must be a positive integer")

        self.database_path = database_path
        self.species_map = _copy_species_map(species_map)
        self.blastp_path = _resolve_executable(blastp_path)
        self.blastdbcmd_path = _resolve_executable(blastdbcmd_path)
        self.evalue_threshold = evalue_threshold
        self.top_n = top_n
        self.excluded_terms = _copy_excluded_terms(excluded_terms)
        self.__species_lookup = {
            species_name.casefold(): (species_name, species_tax_id)
            for species_name, species_tax_id in self.species_map.items()
        }

        subprocess.run(
            [self.blastdbcmd_path, "-db", self.database_path, "-info"],
            capture_output=True,
            check=True,
            text=True,
        )

    def get_orthologs(self, sequence : str) -> list[OrthologRecord]:
        """
        Search the configured local database for complete ortholog proteins.

        Parameters
        ----------
        sequence : str
            Protein sequence used as the BLAST query

        Returns
        -------
        list
            Resolved ortholog records in BLAST result order
        """

        if not isinstance(sequence, str) or not sequence.strip():
            raise ValueError("sequence must be a nonempty string")

        result = subprocess.run(
            [
                self.blastp_path,
                "-db", self.database_path,
                "-outfmt", "5",
                "-evalue", str(self.evalue_threshold),
                "-max_target_seqs", str(self.top_n),
            ],
            input=f">query\n{sequence}\n",
            capture_output=True,
            check=True,
            text=True,
        )
        orthologs = []
        hits = _parse_blast_hits(
            result.stdout,
            len(sequence),
            self.__species_lookup,
            self.excluded_terms,
        )
        for hit in hits:
            try:
                complete_sequence = self.__get_complete_sequence(hit["hit_id"])
            except subprocess.CalledProcessError:
                LOGGER.warning(
                    "Skipping BLAST hit whose complete protein could not be retrieved: %s",
                    hit["hit_id"],
                )
                continue
            if not complete_sequence:
                LOGGER.warning("Skipping BLAST hit with no retrievable protein: %s", hit["hit_id"])
                continue
            orthologs.append(_make_ortholog(hit, complete_sequence))
        return orthologs

    def __get_complete_sequence(self, entry):
        result = subprocess.run(
            [
                self.blastdbcmd_path,
                "-db", self.database_path,
                "-entry", entry,
                "-outfmt", "%s",
            ],
            capture_output=True,
            check=True,
            text=True,
        )
        return "".join(result.stdout.split())


class RemoteBlastSearcher(OrthologSearcherTemplate):
    """
    Search NCBI's remote BLAST protein service for orthologs.

    Parameters
    ----------
    species_map : mapping
        Scientific species names mapped to positive NCBI taxonomic identifiers
    email : str
        Contact email supplied to NCBI E-utilities
    blastp_path : str, optional
        Executable name or path for ``blastp``
    database : str, optional
        Remote BLAST protein database name
    evalue_threshold : float, optional
        Maximum BLAST expectation value
    top_n : int, optional
        Maximum number of BLAST target sequences
    excluded_terms : sequence, optional
        Case-insensitive description terms to exclude; an empty sequence disables filtering
    timeout_seconds : float, optional
        Maximum number of seconds allowed for the remote BLAST+ process
    request_timeout_seconds : float, optional
        Maximum number of seconds allowed for each E-utilities request

    Notes
    -----
    NCBI BLAST and E-utilities are shared services. Avoid parallel searches and
    prefer off-peak hours for substantial workloads.
    """

    def __init__(self,
                 species_map,
                 email,
                 blastp_path="blastp",
                 database="nr",
                 evalue_threshold=1e-6,
                 top_n=20,
                 excluded_terms=None,
                 timeout_seconds=3600,
                 request_timeout_seconds=30):
        if not isinstance(email, str) or not email.strip() or "@" not in email:
            raise ValueError("email must be a nonempty email address")
        if not isinstance(database, str) or not database.strip():
            raise ValueError("database must be a nonempty string")
        if not isinstance(evalue_threshold, (int, float)) or evalue_threshold <= 0:
            raise ValueError("evalue_threshold must be positive")
        if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n <= 0:
            raise ValueError("top_n must be a positive integer")
        if not isinstance(timeout_seconds, (int, float)) or timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be positive")
        if not isinstance(request_timeout_seconds, (int, float)) or request_timeout_seconds <= 0:
            raise ValueError("request_timeout_seconds must be positive")

        self.species_map = _copy_species_map(species_map)
        self.email = email
        self.blastp_path = _resolve_executable(blastp_path)
        self.database = database
        self.evalue_threshold = evalue_threshold
        self.top_n = top_n
        self.excluded_terms = _copy_excluded_terms(excluded_terms)
        self.timeout_seconds = timeout_seconds
        self.request_timeout_seconds = request_timeout_seconds
        self.__species_lookup = {
            species_name.casefold(): (species_name, species_tax_id)
            for species_name, species_tax_id in self.species_map.items()
        }

    def get_orthologs(self, sequence : str) -> list[OrthologRecord]:
        """
        Search remote BLAST and retrieve complete proteins from NCBI.

        Parameters
        ----------
        sequence : str
            Protein sequence used as the BLAST query

        Returns
        -------
        list
            Resolved ortholog records in BLAST result order
        """

        if not isinstance(sequence, str) or not sequence.strip():
            raise ValueError("sequence must be a nonempty string")

        taxonomic_ids = list(dict.fromkeys(self.species_map.values()))
        entrez_query = "(" + " OR ".join(f"txid{tax_id}[ORGN]" for tax_id in taxonomic_ids) + ")"
        result = subprocess.run(
            [
                self.blastp_path,
                "-remote",
                "-db", self.database,
                "-entrez_query", entrez_query,
                "-outfmt", "5",
                "-evalue", str(self.evalue_threshold),
                "-max_target_seqs", str(self.top_n),
            ],
            input=f">query\n{sequence}\n",
            capture_output=True,
            check=True,
            text=True,
            timeout=self.timeout_seconds,
        )
        hits = _parse_blast_hits(
            result.stdout,
            len(sequence),
            self.__species_lookup,
            self.excluded_terms,
        )
        if not hits:
            return []

        proteins = self.__get_complete_sequences([hit["accession"] for hit in hits])
        orthologs = []
        for hit in hits:
            complete_sequence = proteins.get(hit["accession"])
            if complete_sequence is None:
                complete_sequence = proteins.get(hit["accession"].split(".", 1)[0])
            if complete_sequence is None:
                LOGGER.warning("Skipping remote BLAST hit with a missing protein: %s", hit["accession"])
                continue
            orthologs.append(_make_ortholog(hit, complete_sequence))
        return orthologs

    def __get_complete_sequences(self, accessions):
        response = requests.get(
            NCBI_EFETCH_URL,
            params={
                "db": "protein",
                "id": ",".join(accessions),
                "rettype": "fasta",
                "retmode": "text",
                "email": self.email,
                "tool": "Lverage",
            },
            timeout=self.request_timeout_seconds,
        )
        response.raise_for_status()

        proteins = {}
        for record in SeqIO.parse(StringIO(response.text), "fasta"):
            sequence = str(record.seq)
            proteins[record.id] = sequence
            proteins[record.id.split(".", 1)[0]] = sequence
        return proteins
