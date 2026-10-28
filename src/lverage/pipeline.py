"""Core orchestration for motif inference."""

from enum import Enum
import logging

from Bio.Align import PairwiseAligner

from .domain_scanner import DomainRecord
from .domain_scanner import DomainScannerTemplate
from .motif_database import MotifDBTemplate, MotifSearchRequest
from .orf_searcher import OrfSearcherTemplate
from .ortholog_searcher import OrthologSearcherTemplate
from .records import LverageRecord


LOGGER = logging.getLogger(__name__)


class LverageCode(Enum):
    """
    Outcome code for the most recent pipeline run.

    Attributes
    ----------
    NOT_SET : int
        Lverage has not been run or the current run has not completed
    SUCCESS : int
        At least one motif was found
    NO_ORF : int
        No open reading frame was found
    NO_VALID_DOMAIN : int
        No accepted domain was found in an open reading frame
    NO_ORTHOLOGS : int
        No ortholog was found
    NO_VALID_DBD : int
        No matching domain pair met the identity threshold
    NO_MOTIF : int
        Domain evidence was found but no motif was returned
    """

    NOT_SET = -1
    SUCCESS = 0
    NO_ORF = 1
    NO_VALID_DOMAIN = 2
    NO_ORTHOLOGS = 3
    NO_VALID_DBD = 4
    NO_MOTIF = 5
    STOPPED_AFTER_ALIGNMENT = 6


class Lverage:
    """
    Infer transcription factor binding motifs using injected components.

    Parameters
    ----------
    motif_database_list : list
        Concrete motif databases to search
    orf_searcher : OrfSearcherTemplate
        Component that produces protein open reading frames
    domain_scanner : DomainScannerTemplate
        Component that detects protein domains
    ortholog_searcher : OrthologSearcherTemplate
        Component that finds complete orthologous proteins
    valid_pfam_list : list, optional
        Accepted PFAM accessions; an empty list accepts every detected domain
    dbd_identity_thresh : float, optional
        Minimum domain identity ratio from zero to one
    """

    lverage_code_info = {
        LverageCode.NOT_SET: "Lverage has not been run yet",
        LverageCode.SUCCESS: "Motif searching was successful",
        LverageCode.NO_ORF: "No Open Reading Frame was found in the transcription factor sequence",
        LverageCode.NO_VALID_DOMAIN: "No valid domain was found in the Open Reading Frame",
        LverageCode.NO_ORTHOLOGS: "No orthologs were found",
        LverageCode.NO_VALID_DBD: "No valid DNA-binding domain was found in the orthologs",
        LverageCode.NO_MOTIF: "No motif was found in the orthologs",
        LverageCode.STOPPED_AFTER_ALIGNMENT: "Intentionally stopped after alignment",
    }

    def __init__(self,
                 motif_database_list,
                 orf_searcher,
                 domain_scanner,
                 ortholog_searcher,
                 valid_pfam_list=None,
                 dbd_identity_thresh=0.7):
        if valid_pfam_list is None:
            valid_pfam_list = []

        self.motif_database_list = motif_database_list.copy() if isinstance(motif_database_list, list) else motif_database_list
        self.orf_searcher = orf_searcher
        self.domain_scanner = domain_scanner
        self.ortholog_searcher = ortholog_searcher
        self.valid_pfam_list = valid_pfam_list.copy() if isinstance(valid_pfam_list, list) else valid_pfam_list
        self.dbd_identity_thresh = dbd_identity_thresh
        self.lverage_code = LverageCode.NOT_SET

        self.validate_arguments()

    def validate_arguments(self) -> None:
        """Validate configuration supplied to the constructor."""

        if not isinstance(self.motif_database_list, list):
            raise TypeError("motif_database_list must be a list")
        if not self.motif_database_list:
            raise ValueError("motif_database_list must contain at least one motif database")
        if any(not isinstance(database, MotifDBTemplate) for database in self.motif_database_list):
            raise TypeError("motif_database_list must contain only MotifDBTemplate instances")

        if not isinstance(self.orf_searcher, OrfSearcherTemplate):
            raise TypeError("orf_searcher must be an OrfSearcherTemplate instance")
        if not isinstance(self.domain_scanner, DomainScannerTemplate):
            raise TypeError("domain_scanner must be a DomainScannerTemplate instance")
        if not isinstance(self.ortholog_searcher, OrthologSearcherTemplate):
            raise TypeError("ortholog_searcher must be an OrthologSearcherTemplate instance")

        if not isinstance(self.valid_pfam_list, list):
            raise TypeError("valid_pfam_list must be a list")
        if any(not isinstance(accession, str) for accession in self.valid_pfam_list):
            raise TypeError("valid_pfam_list must contain only strings")

        if isinstance(self.dbd_identity_thresh, bool) or not isinstance(self.dbd_identity_thresh, (int, float)):
            raise TypeError("dbd_identity_thresh must be a number")
        if self.dbd_identity_thresh < 0 or self.dbd_identity_thresh > 1:
            raise ValueError("dbd_identity_thresh must be between zero and one")

    def _select_query_orf(self, tf_sequences):
        return self._select_orf_from_candidates(self._get_sorted_orfs(tf_sequences))

    def _get_sorted_orfs(self, tf_sequences, evidence_collector=None):
        orf_list = []
        seen_orfs = set()
        for sequence in tf_sequences:
            for orf in self.orf_searcher.get_orfs(sequence):
                if orf in seen_orfs:
                    continue
                seen_orfs.add(orf)
                orf_list.append(orf)
                if evidence_collector is not None:
                    evidence_collector.emit(
                        "orf_candidate", orf_sequence=orf,
                        orf_index=len(orf_list), status="candidate", error="",
                    )

        orf_list.sort(key=len, reverse=True)
        return orf_list

    def _select_orf_from_candidates(self, orf_list, evidence_collector=None):
        for orf_index, orf in enumerate(orf_list, start=1):
            domains = self.domain_scanner.get_domains(orf)
            detected_domains = list(domains)
            if self.valid_pfam_list:
                domains = [
                    domain for domain in domains
                    if any(
                        ("." in valid_pfam and domain.accession == valid_pfam)
                        or ("." not in valid_pfam and domain.accession.split(".", 1)[0] == valid_pfam)
                        for valid_pfam in self.valid_pfam_list
                    )
                ]
            if evidence_collector is not None:
                for domain in detected_domains:
                    accepted = domain in domains
                    evidence_collector.emit(
                        "query_domain", orf_sequence=orf, orf_index=orf_index,
                        domain=domain, status="accepted" if accepted else "not-selected",
                        error="",
                    )
            if domains:
                if evidence_collector is not None:
                    evidence_collector.emit(
                        "orf_selection", orf_sequence=orf, orf_index=orf_index,
                        status="selected", error="",
                    )
                return orf, domains
            if evidence_collector is not None:
                evidence_collector.emit(
                    "orf_selection", orf_sequence=orf, orf_index=orf_index,
                        status="biological-no-hit", detail="no accepted query domain",
                )
        return None, []

    @staticmethod
    def _calculate_domain_identity(query_sequence : str,
                                   query_domain : DomainRecord,
                                   ortholog_sequence : str,
                                   ortholog_domain : DomainRecord) -> float:
        """
        Calculate exact identity across a global domain alignment.

        Parameters
        ----------
        query_sequence : str
            Complete query protein sequence
        query_domain : DomainRecord
            Zero-based, half-open query domain bounds
        ortholog_sequence : str
            Complete ortholog protein sequence
        ortholog_domain : DomainRecord
            Zero-based, half-open ortholog domain bounds

        Returns
        -------
        float
            Exact-match ratio over every alignment column, including gaps

        Raises
        ------
        ValueError
            A domain has invalid bounds or produces an empty slice
        """

        query_slice = Lverage._get_domain_slice(query_sequence, query_domain)
        ortholog_slice = Lverage._get_domain_slice(ortholog_sequence, ortholog_domain)

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        alignment = aligner.align(query_slice, ortholog_slice)[0]

        exact_matches = 0
        for query_index, ortholog_index in zip(alignment.indices[0], alignment.indices[1]):
            if query_index != -1 and ortholog_index != -1:
                if query_slice[query_index] == ortholog_slice[ortholog_index]:
                    exact_matches += 1
        return exact_matches / alignment.length

    @staticmethod
    def _get_domain_slice(sequence, domain):
        try:
            start = domain.start
            end = domain.end
        except AttributeError as error:
            raise ValueError("domain records must provide start and end bounds") from error

        if not isinstance(sequence, str):
            raise ValueError("domain sequences must be strings")
        if isinstance(start, bool) or isinstance(end, bool):
            raise ValueError("domain bounds must be integers")
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("domain bounds must be integers")
        if start < 0 or end > len(sequence) or start >= end:
            raise ValueError("domain bounds must define a nonempty slice within the sequence")
        return sequence[start:end]

    def run(self, tf_sequence : str | list[str], *, evidence_collector=None,
            stop_after=None) -> list[LverageRecord]:
        """
        Run motif inference for one transcription factor.

        Parameters
        ----------
        tf_sequence : str or list of str
            One sequence or multiple sequence fragments for one factor
        evidence_collector : object, optional
            Evidence collector, defaulting to ``None``. When provided,
            pipeline evidence is reported incrementally as stages complete.
        stop_after : {None, 'alignment'}, optional
            ``None`` (the default) runs the full pipeline, including motif
            lookup. ``'alignment'`` intentionally stops after domain
            alignment and before motif lookup.

        Returns
        -------
        list
            Motif result records. Returns an empty list when no motif records
            are found, including when stopped after alignment.

        Raises
        ------
        TypeError
            The input is not a string or a list containing only strings
        ValueError
            The input is empty, contains an empty sequence, or ``stop_after``
            is neither ``None`` nor ``'alignment'``
        Exception
            An exception raised by an injected pipeline component or evidence
            collector; such exceptions are propagated to the caller
        """

        self.lverage_code = LverageCode.NOT_SET
        if stop_after not in (None, "alignment"):
            raise ValueError("stop_after must be None or 'alignment'")
        tf_sequences = self._validate_run_input(tf_sequence)

        if evidence_collector is not None:
            evidence_collector.set_stage("orf")
        orf_list = self._get_sorted_orfs(tf_sequences, evidence_collector)
        if not orf_list:
            self.lverage_code = LverageCode.NO_ORF
            if evidence_collector is not None:
                evidence_collector.emit("gene_outcome", status="biological-no-hit", detail="no ORF")
            return []

        if evidence_collector is not None:
            evidence_collector.set_stage("query_pfam")
        query_orf, query_domains = self._select_orf_from_candidates(orf_list, evidence_collector)
        if query_orf is None:
            self.lverage_code = LverageCode.NO_VALID_DOMAIN
            if evidence_collector is not None:
                evidence_collector.emit("gene_outcome", status="biological-no-hit", detail="no accepted query domain")
            return []

        if evidence_collector is not None:
            evidence_collector.set_stage("blast")
            if hasattr(evidence_collector, "set_context"):
                evidence_collector.set_context(orf_sequence=query_orf)
        supports_evidence = getattr(
            self.ortholog_searcher, "supports_evidence_collector", False
        )
        if evidence_collector is not None and supports_evidence:
            orthologs = self.ortholog_searcher.get_orthologs(
                query_orf, evidence_collector=evidence_collector
            )
        else:
            orthologs = self.ortholog_searcher.get_orthologs(query_orf)
        if evidence_collector is not None and not supports_evidence:
            for ortholog_index, ortholog in enumerate(orthologs, start=1):
                evidence_collector.emit(
                    "ortholog", orf_sequence=query_orf, ortholog=ortholog,
                    ortholog_index=ortholog_index, status="retained", error="",
                )
        if not orthologs:
            self.lverage_code = LverageCode.NO_ORTHOLOGS
            if evidence_collector is not None:
                retrieval_failed = any(
                    event.get("event") == "blast_search_summary"
                    and event.get("status") == "retrieval-failed"
                    for event in getattr(evidence_collector, "events", [])
                )
                evidence_collector.emit(
                    "gene_outcome",
                    status="failed" if retrieval_failed else "biological-no-hit",
                    detail=("all BLAST hits failed sequence retrieval" if retrieval_failed
                            else "no retained ortholog"),
                )
            return []

        records = []
        has_valid_domain_pair = False
        for ortholog in orthologs:
            if evidence_collector is not None:
                evidence_collector.set_stage("ortholog_pfam")
            ortholog_domains = self.domain_scanner.get_domains(ortholog.sequence)
            if evidence_collector is not None:
                for domain_index, domain in enumerate(ortholog_domains, start=1):
                    evidence_collector.emit(
                        "ortholog_domain", orf_sequence=query_orf, ortholog=ortholog,
                        domain=domain, domain_index=domain_index,
                        status="retained", error="",
                    )
            if evidence_collector is not None:
                evidence_collector.set_stage("alignment")
            for query_domain in query_domains:
                for ortholog_domain in ortholog_domains:
                    if self._base_accession(query_domain.accession) != self._base_accession(ortholog_domain.accession):
                        continue
                    domain_identity = self._calculate_domain_identity(
                        query_orf,
                        query_domain,
                        ortholog.sequence,
                        ortholog_domain,
                    )
                    passes_threshold = domain_identity >= self.dbd_identity_thresh
                    if evidence_collector is not None:
                        evidence_collector.set_stage("alignment")
                        evidence_collector.emit(
                            "domain_comparison", orf_sequence=query_orf, ortholog=ortholog,
                            query_domain=query_domain, ortholog_domain=ortholog_domain,
                            domain_identity=domain_identity,
                            threshold=self.dbd_identity_thresh,
                            status="accepted" if passes_threshold else "below-threshold",
                            error="",
                        )
                    if domain_identity < self.dbd_identity_thresh:
                        continue

                    has_valid_domain_pair = True
                    request = MotifSearchRequest(
                        query_sequence=query_orf,
                        query_domain=query_domain,
                        ortholog_sequence=ortholog.sequence,
                        ortholog_domain=ortholog_domain,
                        ortholog_species_tax_id=ortholog.species_tax_id,
                    )
                    if stop_after == "alignment":
                        continue
                    if evidence_collector is not None:
                        evidence_collector.set_stage("motif")
                    for motif_database in self.motif_database_list:
                        if not motif_database.check_species_validity(ortholog.species_tax_id):
                            LOGGER.info(
                                "Skipping %s for unavailable species %s",
                                motif_database.name,
                                ortholog.species_tax_id,
                            )
                            continue
                        for motif_record in motif_database.search(request):
                            records.append(LverageRecord(
                                query_domain=query_domain,
                                ortholog=ortholog,
                                ortholog_domain=ortholog_domain,
                                domain_identity=domain_identity,
                                motif_database_name=motif_database.name,
                                motif_record=motif_record,
                            ))
                            if evidence_collector is not None:
                                evidence_collector.emit(
                                    "motif", orf_sequence=query_orf, ortholog=ortholog,
                                    query_domain=query_domain, ortholog_domain=ortholog_domain,
                                    domain_identity=domain_identity,
                                    threshold=self.dbd_identity_thresh,
                                    motif_database_name=motif_database.name,
                                    motif_record=motif_record, status="retained", error="",
                                )

        if stop_after == "alignment":
            self.lverage_code = LverageCode.STOPPED_AFTER_ALIGNMENT
            if evidence_collector is not None:
                evidence_collector.set_stage("alignment")
                evidence_collector.emit("gene_outcome", status="intentional-stop", error="")
            return []

        if records:
            self.lverage_code = LverageCode.SUCCESS
        elif has_valid_domain_pair:
            self.lverage_code = LverageCode.NO_MOTIF
        else:
            self.lverage_code = LverageCode.NO_VALID_DBD
        if evidence_collector is not None and self.lverage_code in (LverageCode.NO_MOTIF, LverageCode.NO_VALID_DBD):
            evidence_collector.emit(
                "gene_outcome", status="biological-no-hit",
                detail=self.lverage_code.name,
            )
        return records

    @staticmethod
    def _base_accession(accession):
        if not isinstance(accession, str) or not accession:
            raise ValueError("domain accessions must be nonempty strings")
        return accession.split(".", 1)[0]

    @staticmethod
    def _validate_run_input(tf_sequence):
        if isinstance(tf_sequence, str):
            if not tf_sequence.strip():
                raise ValueError("tf_sequence must be nonempty")
            return [tf_sequence]
        if not isinstance(tf_sequence, list):
            raise TypeError("tf_sequence must be a string or list of strings")
        if not tf_sequence:
            raise ValueError("tf_sequence must contain at least one sequence")
        if any(not isinstance(sequence, str) for sequence in tf_sequence):
            raise TypeError("tf_sequence must contain only strings")
        if any(not sequence.strip() for sequence in tf_sequence):
            raise ValueError("tf_sequence must contain only nonempty strings")
        return list(tf_sequence)
