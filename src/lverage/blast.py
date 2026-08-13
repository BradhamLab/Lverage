"""BLAST-based ortholog searchers."""

from collections.abc import Mapping, Sequence
from io import StringIO
import logging
import re
import shutil
import subprocess

from Bio.Blast import NCBIXML

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
        blast_records = list(NCBIXML.parse(StringIO(result.stdout)))
        if len(blast_records) != 1:
            raise ValueError("BLAST output did not contain exactly one query result")

        orthologs = []
        for alignment in blast_records[0].alignments:
            ortholog = self.__create_ortholog(alignment, len(sequence))
            if ortholog is not None:
                orthologs.append(ortholog)
        return orthologs

    def __create_ortholog(self, alignment, query_length):
        try:
            description = alignment.hit_def
            if any(term.casefold() in description.casefold() for term in self.excluded_terms):
                LOGGER.warning("Skipping excluded BLAST hit: %s", description)
                return None

            species_match = re.search(r"\[([^\]]+)\]", description)
            if species_match is None:
                LOGGER.warning("Skipping BLAST hit without a bracketed species: %s", description)
                return None
            species = self.__species_lookup.get(species_match.group(1).strip().casefold())
            if species is None:
                LOGGER.warning("Skipping BLAST hit with an unresolved species: %s", description)
                return None

            best_hsp = min(alignment.hsps, key=lambda hsp: hsp.expect)
            complete_sequence = self.__get_complete_sequence(alignment.hit_id)
            if not complete_sequence:
                LOGGER.warning("Skipping BLAST hit with no retrievable protein: %s", alignment.hit_id)
                return None

            return OrthologRecord(
                accession=alignment.accession,
                description=description,
                species_name=species[0],
                species_tax_id=species[1],
                sequence=complete_sequence,
                evalue=best_hsp.expect,
                identity=best_hsp.identities / best_hsp.align_length,
                query_coverage=best_hsp.align_length / query_length,
            )
        except (AttributeError, TypeError, ValueError, ZeroDivisionError):
            LOGGER.warning("Skipping malformed BLAST hit", exc_info=True)
            return None
        except subprocess.CalledProcessError:
            LOGGER.warning("Skipping BLAST hit whose complete protein could not be retrieved: %s", alignment.hit_id)
            return None

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
