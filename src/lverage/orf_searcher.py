"""
This file contains the template for all OrfSearcher classes. 
Any OrfSearcher tool that is utilized with Lverage MUST inherit from this class!

Copyright (C) <RELEASE_YEAR_HERE> Bradham Lab

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

Correspondence: 
    - Cynthia A. Bradham - cbradham@bu.edu - *
    - Anthony B. Garza   - abgarza@bu.edu  - **
    - Stephanie P. Hao   - sphao@bu.edu    - **
    - Yeting Li          - yetingli@bu.edu - **
    - Thomas Shin        - thshin@bu.edu   - **

    \* Principle Investigator, ** Software Developers
"""

from abc import ABC, abstractmethod
from collections.abc import Sequence

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from orffinder import orffinder

class OrfSearcherTemplate(ABC):
    """
    Abstract class for OrfSearcher classes, containing all necessary functionalities
    This should not be instantiated!
    """

    @abstractmethod
    def get_orfs(self, sequence : str) -> list[str]:
        """
        Method for searching a DNA sequence for open-reading-frames
        This method should be overridden by a concrete class

        Parameters
        ----------
        sequence : str
            DNA sequence to search for open-reading-frames

        Returns
        -------
        list
            List of str objects representing the open-reading-frames found in the sequence
        
        Raises
        ------
        NotImplementedError
        """

        raise NotImplementedError


class OrffinderOrfSearcher(OrfSearcherTemplate):
    """Find protein-producing ORFs using ``orffinder``.

    Parameters
    ----------
    minimum_length : int, optional
        Minimum ORF length in nucleotides. Complete ORFs include the stop
        codon in this count. The default of 75 matches ``orffinder`` and the
        legacy workflow.
    start_codons : sequence of str, optional
        Start codons to recognize. Installed ``orffinder`` 1.8 silently
        ignores custom start codons, so only ATG is supported.
    remove_nested : bool, optional
        Whether to ask ``orffinder`` to remove fully nested candidates.
        Defaults to False, which disables its nested-candidate filter.
    trim_trailing : bool, optional
        Whether to discard candidates that reach a sequence end without a
        stop codon. Defaults to False, matching the legacy workflow.

    Notes
    -----
    Input is a nonempty, case-insensitive DNA sequence using the IUPAC
    alphabet ``ACGTRYSWKMBDHVN``. Ambiguous codons may translate to ``X``.
    Whitespace, gaps, RNA ``U``, and other symbols are rejected rather than
    removed or converted. Sequences are passed without padding; the finder
    scans all six reading frames itself.

    ORFs are translated with the standard genetic code, excluding a complete
    ORF's terminal stop codon. No translation-table option is exposed because
    ``orffinder`` 1.8 accepts but ignores that argument. Its protein-producing
    function can truncate complete codons for incomplete terminal ORFs; this
    adapter instead extracts codons from the original sequence using the
    reported strand and start coordinate. Minimum length is checked against
    the actual nucleotide span, including any incomplete trailing bases.
    """

    _DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")

    def __init__(self,
                 minimum_length: int = 75,
                 start_codons: Sequence[str] = ("ATG",),
                 remove_nested: bool = False,
                 trim_trailing: bool = False):
        if isinstance(minimum_length, bool) or not isinstance(minimum_length, int):
            raise TypeError("minimum_length must be an integer")
        if minimum_length <= 0:
            raise ValueError("minimum_length must be greater than zero")

        if isinstance(start_codons, (str, bytes)) or not isinstance(start_codons, Sequence):
            raise TypeError("start_codons must be a sequence of codon strings")
        if not start_codons:
            raise ValueError("start_codons must contain ATG")
        if any(not isinstance(codon, str) for codon in start_codons):
            raise TypeError("start_codons must contain only strings")
        normalized_codons = tuple(codon.upper() for codon in start_codons)
        if any(len(codon) != 3 or set(codon) - set("ACGT") for codon in normalized_codons):
            raise ValueError("start codons must be three-base DNA codons")
        if set(normalized_codons) != {"ATG"}:
            raise ValueError("orffinder 1.8 supports only the ATG start codon")

        if not isinstance(remove_nested, bool):
            raise TypeError("remove_nested must be a boolean")
        if not isinstance(trim_trailing, bool):
            raise TypeError("trim_trailing must be a boolean")

        self.minimum_length = minimum_length
        self.start_codons = ("ATG",)
        self.remove_nested = remove_nested
        self.trim_trailing = trim_trailing

    def get_orfs(self, sequence: str) -> list[str]:
        """Return translated ORF candidates from a nucleotide sequence.

        Parameters
        ----------
        sequence : str
            Nonempty DNA sequence containing only IUPAC symbols
            ``ACGTRYSWKMBDHVN`` (case-insensitive).

        Returns
        -------
        list of str
            Translated protein candidates without terminal stop symbols or
            empty proteins. Returns an empty list when no usable candidate is
            found.

        Raises
        ------
        TypeError
            If ``sequence`` is not a string.
        ValueError
            If ``sequence`` is empty or contains unsupported symbols.
        """

        if not isinstance(sequence, str):
            raise TypeError("sequence must be a string")
        if not sequence:
            raise ValueError("sequence must not be empty")
        invalid_symbols = sorted(set(sequence.upper()) - self._DNA_ALPHABET)
        if invalid_symbols:
            raise ValueError(
                "sequence contains unsupported symbols: " + "".join(invalid_symbols)
            )

        nucleotide_sequence = Seq(sequence.upper())
        sequence_length = len(nucleotide_sequence)
        # orffinder 1.8 counts incomplete terminal candidates one nucleotide
        # short, so scan one nt below our threshold and enforce it accurately
        # against the extracted span below.
        finder_minimum_length = max(1, self.minimum_length - 1)
        loci = orffinder.getORFs(
            SeqRecord(nucleotide_sequence),
            minimum_length=finder_minimum_length,
            remove_nested=self.remove_nested,
            trim_trailing=self.trim_trailing,
        )

        forward = str(nucleotide_sequence)
        reverse = str(nucleotide_sequence.reverse_complement())
        proteins = []
        for locus in loci:
            is_trailing = locus["trailing"]
            strand_sequence = forward if locus["sense"] == "+" else reverse
            start_index = (
                locus["start"] - 1
                if locus["sense"] == "+"
                else sequence_length - locus["start"] + 1
            )

            if is_trailing:
                nucleotide_length = sequence_length - start_index
                if nucleotide_length < self.minimum_length:
                    continue
                coding_length = nucleotide_length - nucleotide_length % 3
            else:
                nucleotide_length = locus["length"]
                if nucleotide_length < self.minimum_length:
                    continue
                coding_length = nucleotide_length - 3  # Exclude terminal stop.

            if coding_length <= 0:
                continue
            coding_sequence = strand_sequence[start_index:start_index + coding_length]
            protein = str(Seq(coding_sequence).translate()).rstrip("*")
            if protein:
                proteins.append(protein)

        return proteins
