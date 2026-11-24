"""Interfaces and Biopython implementations for domain alignment."""

from abc import ABC, abstractmethod
from collections.abc import Sequence
import math
from numbers import Real

from Bio.Align import PairwiseAligner


class AlignmentTemplate(ABC):
    """Abstract interface for calculating query-to-ortholog domain identities."""

    @abstractmethod
    def calc_identities(self, query_dbd: str, ortholog_dbds: Sequence[str]) -> list[float]:
        """Calculate one identity ratio for each ortholog domain, in input order.

        Parameters
        ----------
        query_dbd : str
            Nonempty query protein-domain sequence.
        ortholog_dbds : sequence of str
            Ordered, nonempty ortholog protein-domain sequences. Duplicate
            sequences represent distinct occurrences and must be retained.

        Returns
        -------
        list of float
            One ratio in the range zero to one for each input occurrence, in
            the same order.
        """

        raise NotImplementedError


class BiopythonGlobalAlignment(AlignmentTemplate):
    """Calculate exact identity from Biopython global pairwise alignments.

    The default ``gap_model="affine"`` retains the V2 scores: match 2,
    mismatch -1, gap open -2, and gap extension -0.5. This selects Gotoh's
    affine-gap global algorithm. Set ``gap_model="linear"`` and use equal
    ``open_gap_score`` and ``extend_gap_score`` values (for example, both -2)
    to select Biopython's Needleman-Wunsch linear-gap algorithm. The latter is
    explicit and does not alter the default scoring behavior.

    Identity is the number of exact, case-sensitive residue matches divided
    by the number of alignment columns; mismatches and columns containing a
    gap remain in the denominator.

    Parameters
    ----------
    gap_model : {"affine", "linear"}, optional
        Gap scoring model. Defaults to the existing affine-gap model.
    match_score : real number, optional
        Score for an exact residue match. Defaults to 2.
    mismatch_score : real number, optional
        Score for a residue mismatch. Defaults to -1.
    open_gap_score : real number, optional
        Gap-open score; defaults to -2.
    extend_gap_score : real number, optional
        Gap-extension score; defaults to -0.5. For linear scoring, this must
        equal ``open_gap_score``.

    Raises
    ------
    TypeError
        If the gap model or a score has an unsupported type.
    ValueError
        If a score is non-finite or linear gap scores are unequal.
    """

    def __init__(self,
                 gap_model: str = "affine",
                 match_score: float = 2,
                 mismatch_score: float = -1,
                 open_gap_score: float = -2,
                 extend_gap_score: float = -0.5):
        if gap_model not in ("affine", "linear"):
            raise ValueError("gap_model must be 'affine' or 'linear'")

        scores = {
            "match_score": match_score,
            "mismatch_score": mismatch_score,
            "open_gap_score": open_gap_score,
            "extend_gap_score": extend_gap_score,
        }
        for name, score in scores.items():
            if isinstance(score, bool) or not isinstance(score, Real):
                raise TypeError(f"{name} must be a real number")
            if not math.isfinite(score):
                raise ValueError(f"{name} must be finite")
        if gap_model == "linear" and open_gap_score != extend_gap_score:
            raise ValueError("linear gap scoring requires equal open and extend scores")

        self.gap_model = gap_model
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.open_gap_score = open_gap_score
        self.extend_gap_score = extend_gap_score

    def calc_identities(self, query_dbd: str, ortholog_dbds: Sequence[str]) -> list[float]:
        """Return exact identity ratios for the ordered ortholog domains.

        Parameters
        ----------
        query_dbd : str
            Nonempty query protein-domain sequence.
        ortholog_dbds : sequence of str
            Ordered, nonempty ortholog domain sequences.

        Returns
        -------
        list of float
            One exact-identity ratio per input sequence, preserving duplicates.

        Raises
        ------
        TypeError
            If inputs are not strings and an ordered sequence of strings.
        ValueError
            If any input sequence is empty.
        """

        if not isinstance(query_dbd, str):
            raise TypeError("query_dbd must be a string")
        if not query_dbd:
            raise ValueError("query_dbd must be nonempty")
        if isinstance(ortholog_dbds, (str, bytes)) or not isinstance(ortholog_dbds, Sequence):
            raise TypeError("ortholog_dbds must be an ordered sequence of strings")
        if any(not isinstance(sequence, str) for sequence in ortholog_dbds):
            raise TypeError("ortholog_dbds must contain only strings")
        if any(not sequence for sequence in ortholog_dbds):
            raise ValueError("ortholog_dbds must contain only nonempty strings")

        aligner = self._make_aligner()
        identities = []
        for ortholog_dbd in ortholog_dbds:
            alignment = aligner.align(query_dbd, ortholog_dbd)[0]
            exact_matches = sum(
                query_index != -1
                and ortholog_index != -1
                and query_dbd[query_index] == ortholog_dbd[ortholog_index]
                for query_index, ortholog_index in zip(alignment.indices[0], alignment.indices[1])
            )
            identities.append(exact_matches / alignment.length)
        return identities

    def _make_aligner(self) -> PairwiseAligner:
        """Build the configured Biopython aligner used for each comparison."""

        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = self.match_score
        aligner.mismatch_score = self.mismatch_score
        aligner.open_gap_score = self.open_gap_score
        aligner.extend_gap_score = self.extend_gap_score
        return aligner