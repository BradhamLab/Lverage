"""
This file contains the DomainScannerTemplate class which details a template for scanning domains from a given sequence

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

    \* Principle Investigator, ** Software Developers
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class DomainRecord:
    """
    Record describing a domain found in a protein sequence.

    Attributes
    ----------
    name : str
        Name of the domain
    accession : str
        Accession identifier of the domain
    start : int
        Zero-based, inclusive start position of the domain
    end : int
        Zero-based, exclusive end position of the domain
    """

    name : str
    accession : str
    start : int
    end : int

class DomainScannerTemplate(ABC):

    @abstractmethod
    def get_domains(self, sequence : str) -> list[DomainRecord]:
        """
        Method for scanning a protein sequence for domains.

        Parameters
        ----------
        sequence : str
            Protein sequence to scan for domains

        Returns
        -------
        list
            List of DomainRecord objects representing the domains found

        Raises
        ------
        NotImplementedError
        """

        raise NotImplementedError
