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


def _parse_blast_hits(blast_xml, query_length, species_lookup, excluded_terms,
                      evidence_collector=None):
    blast_records = list(NCBIXML.parse(StringIO(blast_xml)))
    if len(blast_records) != 1:
        raise ValueError("BLAST output did not contain exactly one query result")

    raw_alignments = blast_records[0].alignments
    hits = []
    for alignment in raw_alignments:
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
                "hit_id": getattr(alignment, "hit_id", None),
                "species_name": species[0],
                "species_tax_id": species[1],
                "evalue": best_hsp.expect,
                "identity": best_hsp.identities / best_hsp.align_length,
                "query_coverage": (
                    best_hsp.query_end - best_hsp.query_start + 1
                ) / query_length,
            })
        except (AttributeError, TypeError, ValueError, ZeroDivisionError):
            LOGGER.warning("Skipping malformed BLAST hit", exc_info=True)
    if evidence_collector is not None:
        evidence_collector.emit(
            "blast_search",
            status=("no-hits" if not raw_alignments else "hits-found" if hits else "no-retained-hits"),
            hit_count=len(raw_alignments), retained_hit_count=len(hits),
            retrieved_count=0, skipped_count=0, error="",
        )
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
        hit_id=hit.get("hit_id"),
    )


def _safe_retrieval_reason(error):
    """Return a diagnostic without echoing command arguments or credentials."""
    detail = getattr(error, "stderr", None)
    if isinstance(detail, bytes):
        detail = detail.decode("utf-8", errors="replace")
    if not isinstance(detail, str) or not detail.strip():
        detail = f"{type(error).__name__}: sequence retrieval failed"
    detail = re.sub(r"https?://[^\s]+", "<url-redacted>", detail)
    detail = re.sub(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", "<email-redacted>", detail)
    detail = re.sub(r"(?i)(token|password|secret|api[_-]?key)=([^\s&]+)", r"\1=<redacted>", detail)
    return detail.strip()


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

    supports_evidence_collector = True

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

    def get_orthologs(self, sequence : str, *, evidence_collector=None) -> list[OrthologRecord]:
        """
        Search the configured local database for complete ortholog proteins.

        Parameters
        ----------
        sequence : str
            Protein sequence used as the BLAST query
        evidence_collector : object, optional
            Evidence collector, defaulting to ``None``. When provided, BLAST
            search, retrieval attempts, and retained ortholog evidence are
            reported incrementally as the search progresses.

        Returns
        -------
        list
            Resolved ortholog records in BLAST result order
        """

        if not isinstance(sequence, str) or not sequence.strip():
            raise ValueError("sequence must be a nonempty string")

        if evidence_collector is not None:
            evidence_collector.set_stage("blast")

        result = subprocess.run(
            [
                self.blastp_path,
                "-db", self.database_path,
                                "-query", "-",
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
            evidence_collector=evidence_collector,
        )

        skipped_count = 0
        retrieval_exception = None
        record_exception = None
        evidence_exception = None

        def emit_evidence(event, **payload):
            nonlocal evidence_exception
            if evidence_collector is None:
                return
            try:
                evidence_collector.emit(event, **payload)
            except Exception as error:
                evidence_exception = error
                raise

        try:
            for hit in hits:
                complete_sequence = ""
                entries = [hit["hit_id"]]
                if hit["accession"] and hit["accession"] != hit["hit_id"]:
                    entries.append(hit["accession"])
                retrieval_errors = []
                for entry in entries:
                    if not entry:
                        continue
                    try:
                        complete_sequence = self.__get_complete_sequence(entry)
                    except subprocess.CalledProcessError as error:
                        LOGGER.warning("Failed to retrieve BLAST hit: %s", entry)
                        reason = _safe_retrieval_reason(error)
                        retrieval_errors.append(f"{entry}: {reason}")
                        emit_evidence(
                            "ortholog_retrieval", status="attempt-failed",
                            hit_id=hit.get("hit_id"), accession=hit.get("accession"),
                            attempted_identifier=entry, error=reason,
                        )
                        continue
                    except Exception as error:
                        retrieval_exception = error
                        try:
                            emit_evidence(
                                "ortholog_retrieval", status="failed",
                                hit_id=hit.get("hit_id"), accession=hit.get("accession"),
                                attempted_identifier=entry,
                                error=_safe_retrieval_reason(error),
                            )
                        except Exception:
                            LOGGER.exception("Failed to emit BLAST retrieval diagnostic")
                        raise
                    if complete_sequence:
                        break
                    retrieval_errors.append(f"{entry}: empty sequence returned")
                    emit_evidence(
                        "ortholog_retrieval", status="attempt-failed",
                        hit_id=hit.get("hit_id"), accession=hit.get("accession"),
                        attempted_identifier=entry,
                        error="empty sequence returned",
                    )
                if not complete_sequence:
                    LOGGER.warning(
                        "Skipping BLAST hit whose complete protein could not be retrieved: %s",
                        hit.get("hit_id") or hit.get("accession"),
                    )
                    skipped_count += 1
                    reason = "; ".join(retrieval_errors) or "no usable hit identifier was available"
                    emit_evidence(
                        "ortholog_retrieval", status="skipped",
                        hit_id=hit.get("hit_id"), accession=hit.get("accession"),
                        attempted_identifier="", error=reason,
                    )
                    continue
                try:
                    ortholog = _make_ortholog(hit, complete_sequence)
                except Exception as error:
                    record_exception = error
                    raise
                orthologs.append(ortholog)
                emit_evidence("ortholog", ortholog=ortholog, status="retained", error="")
        except Exception:
            if evidence_collector is not None:
                try:
                    if retrieval_exception is not None:
                        failure_status = "retrieval-failed"
                        failure_detail = _safe_retrieval_reason(retrieval_exception)
                    elif record_exception is not None:
                        failure_status = "record-construction-failed"
                        failure_detail = _safe_retrieval_reason(record_exception)
                    elif evidence_exception is not None:
                        failure_status = "evidence-write-failed"
                        failure_detail = "evidence callback failed while reporting BLAST results"
                    else:
                        failure_status = "processing-failed"
                        failure_detail = "unexpected BLAST hit-processing failure"
                    evidence_collector.emit(
                        "blast_search_summary",
                        status=failure_status,
                        hit_count=len(hits), retrieved_count=len(orthologs),
                        skipped_count=skipped_count,
                        error=failure_detail,
                    )
                except Exception:
                    LOGGER.exception("Failed to emit BLAST failure summary")
            raise

        if evidence_collector is not None and hits:
            retrieved_count = len(orthologs)
            search_status = (
                "retrieval-failed" if retrieved_count == 0
                else "partial-retrieval-failure" if skipped_count
                else "completed"
            )
            evidence_collector.emit(
                "blast_search_summary", status=search_status,
                hit_count=len(hits), retrieved_count=retrieved_count,
                skipped_count=skipped_count, error="",
            )
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
