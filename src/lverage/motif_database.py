"""
This file contains the template for all MotifDB classes. 
Any MotifDB that is utilized with Lverage MUST inherit from this class!

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

#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# IMPORTS
from abc import ABC, abstractmethod
from dataclasses import dataclass

from .domain_scanner import DomainRecord


@dataclass
class MotifSearchRequest:
    """
    Record containing the evidence used to search for a motif.

    Attributes
    ----------
    query_sequence : str
        Protein sequence of the transcription factor being investigated
    query_domain : DomainRecord
        DNA-binding domain from the query protein sequence
    ortholog_sequence : str
        Protein sequence of the supporting ortholog
    ortholog_domain : DomainRecord
        Matching DNA-binding domain from the ortholog protein sequence
    ortholog_species_tax_id : int
        NCBI taxonomic identifier of the ortholog species
    """

    query_sequence : str
    query_domain : DomainRecord
    ortholog_sequence : str
    ortholog_domain : DomainRecord
    ortholog_species_tax_id : int


class MotifDBTemplate(ABC):
    """
    Abstract class for MotifDB classes, containing all necessary functionalities
    This should not be instantiated!

    Attributes
    ----------
    name: str
        Name of the motif database
    metainfo: dict
        Any extra information regarding this motif database
    """

    name = ""
    metainfo = {}

    @abstractmethod
    def search(self, request : MotifSearchRequest) -> list["MotifDBRecordTemplate"]:
        """
        Method for searching the database for a motif.
        This method should be overridden by a concrete class

        Parameters
        ----------
        request : MotifSearchRequest
            Query and ortholog evidence used for motif searching

        Returns
        -------
        list
            List of MotifDBRecordTemplate objects representing records in the database
        
        Raises
        ------
        NotImplementedError
        """

        raise NotImplementedError

    @abstractmethod
    def check_species_validity(self, species_tax_id : int):
        """
        Method to check if a species appears in the database.
        
        Parameters
        ----------
        species_tax_id : int
            Taxonomic identifier of the species
        
        Returns
        ----------
        is_present : bool
            If the species appears in the database
        """

        raise NotImplementedError
    
class MotifDBRecordTemplate(ABC):
    """
    Abstract class for MotifDB records, containing all necessary functionalities. 
    All MotifDB classes should have a record class that inherits from this class.
    Subclasses should be returned by the search method of the MotifDB class.
    """

    headers : tuple[str, ...] = () # Tuple of headers in the record; should be overridden by subclasses

    @classmethod
    def get_headers(cls) -> list[str]:
        """
        Method to get the headers of the record.
        
        Returns
        -------
        list
            List of headers (str) the record uses
        """

        return list(cls.headers)

    @abstractmethod
    def get_values(self):
        """
        Method to get the values of the record.
        
        Returns
        -------
        list
            List of values in the record
        """

        raise NotImplementedError
    




    
    
