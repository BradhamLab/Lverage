"""Core orchestration for motif inference."""

from enum import Enum

from .domain_scanner import DomainScannerTemplate
from .motif_database import MotifDBTemplate
from .orf_searcher import OrfSearcherTemplate
from .ortholog_searcher import OrthologSearcherTemplate


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
        orf_list = []
        for sequence in tf_sequences:
            orf_list.extend(self.orf_searcher.get_orfs(sequence))

        orf_list = list(dict.fromkeys(orf_list))
        orf_list.sort(key=len, reverse=True)

        for orf in orf_list:
            domains = self.domain_scanner.get_domains(orf)
            if self.valid_pfam_list:
                domains = [
                    domain for domain in domains
                    if any(
                        ("." in valid_pfam and domain.accession == valid_pfam)
                        or ("." not in valid_pfam and domain.accession.split(".", 1)[0] == valid_pfam)
                        for valid_pfam in self.valid_pfam_list
                    )
                ]
            if domains:
                return orf, domains
        return None, []

    def run(self, tf_sequence : str | list[str]):
        """
        Run motif inference for one transcription factor.

        Parameters
        ----------
        tf_sequence : str or list of str
            One sequence or multiple sequence fragments for one factor

        Raises
        ------
        NotImplementedError
            Pipeline orchestration is added in a later core commit
        """

        raise NotImplementedError
