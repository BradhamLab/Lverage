"""
This file contains the Jaspar2024 Motif Database class which can be found at https://jaspar.elixir.no/.
Citation:
    Rauluseviciute I, Riudavets-Puig R, Blanc-Mathieu R, Castro-Mondragon JA, Ferenc K, Kumar V, Lemma RB, 
    Lucas J, Chèneby J, Baranasic D, Khan A, Fornes O, Gundersen S, Johansen M, Hovig E, Lenhard B, 
    Sandelin A, Wasserman WW, Parcy F, Mathelier A JASPAR 2024: 20th anniversary of the open-access 
    database of transcription factor binding profiles Nucleic Acids Res. 2024 Jan 5;52(D1):D174-D182.; 
    doi: 10.1093/nar/gkad1059

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
# Imports
from .motif_database import MotifDBTemplate, MotifSearchRequest
import requests

class Jaspar2024MotifDB(MotifDBTemplate):
    """
    Class for the JASPAR2024 Motif Database, which can be found at https://jaspar.elixir.no/

    Parameters
    ----------
    n_hits : int
        Number of motif hits to return
    escore_threshold : float
        Minimum e-score threshold for motif hits

    Attributes
    ----------
    name: str
        Name of the motif database
    metainfo: dict
        Any extra information regarding this motif database
    jaspar_rest_url: str
        URL for the JASPAR REST API
    jaspar_logo_url: str
        URL for retrieving motif logos
    n_hits: int
        Number of motif hits to return
    escore_threshold: float
        Minimum e-score threshold for motif hits
    """

    name = "JASPAR2024"
    metainfo = {
                "citation":"Rauluseviciute I, Riudavets-Puig R, Blanc-Mathieu R, Castro-Mondragon JA, Ferenc K, Kumar V, Lemma RB, Lucas J, Chèneby J, Baranasic D, Khan A, Fornes O, Gundersen S, Johansen M, Hovig E, Lenhard B, Sandelin A, Wasserman WW, Parcy F, Mathelier A JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles Nucleic Acids Res. 2024 Jan 5;52(D1):D174-D182.; doi: 10.1093/nar/gkad1059",
                "version": "JASPAR2024",
                "collection":"JASPAR CORE",
                "url":"https://jaspar.elixir.no/"
            }
    
    # URL for the JASPAR REST API
    jaspar_rest_url = "https://jaspar.elixir.no/api/v1"

    # URL for retrieving species information
    jaspar_rest_species_url = "https://jaspar.elixir.no/api/v1/species"

    # URL for retrieving motif logos
    jaspar_logo_url = "https://jaspar2020.genereg.net/static/logos/all/"



    # Parameters for getting all species
    species_params = {"page":1,
                      "page_size":1000,
                      "release":"2024"
                      }

    def __init__(self, 
                 n_hits : int = 10, 
                 escore_threshold : float = 10**-6):
        """Constructor"""
        
        self.n_hits = n_hits
        self.escore_threshold = escore_threshold
        self.jaspar_species = None

    def search(self, request : MotifSearchRequest):
        """
        Search JASPAR for motif records.

        This method is not implemented yet.
        """

        raise NotImplementedError

    def check_species_validity(self, species_tax_id : int) -> bool:
        """
        Check whether a species appears in JASPAR.

        Parameters
        ----------
        species_tax_id : int
            Taxonomic identifier of the species

        Returns
        -------
        bool
            If the species appears in the database
        """

        if self.jaspar_species is None:
            species_result = requests.get(self.jaspar_rest_species_url, params=self.species_params).json()['results']
            self.jaspar_species = [species["tax_id"] for species in species_result]

        return species_tax_id in self.jaspar_species

if __name__ == "__main__":
    print(Jaspar2024MotifDB.jaspar_species)
