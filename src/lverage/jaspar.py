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
from dataclasses import dataclass

from .motif_database import MotifDBRecordTemplate, MotifDBTemplate, MotifSearchRequest
import requests


@dataclass
class JasparRecord(MotifDBRecordTemplate):
    """
    Motif record returned by the JASPAR adapter.

    Parameters
    ----------
    matrix_id : str
        JASPAR matrix identifier
    motif_name : str
        Name of the motif
    pfm : dict
        Position frequency matrix supplied by JASPAR
    motif_url : str
        URL for the motif page or logo
    motif_class : str
        JASPAR motif class
    inference_evalue : float
        E-value reported by sequence inference
    """

    matrix_id : str
    motif_name : str
    pfm : dict
    motif_url : str
    motif_class : str
    inference_evalue : float

    headers = (
        "Matrix ID",
        "Motif Name",
        "PFM",
        "Motif URL",
        "Motif Class",
        "Inference E-value",
    )

    def get_values(self) -> list:
        """
        Return values in the JASPAR record schema.

        Returns
        -------
        list
            Serialized JASPAR values
        """

        return [
            self.matrix_id,
            self.motif_name,
            self.pfm,
            self.motif_url,
            self.motif_class,
            self.inference_evalue,
        ]


class Jaspar2024MotifDB(MotifDBTemplate):
    """
    Class for the JASPAR2024 Motif Database, which can be found at https://jaspar.elixir.no/

    Parameters
    ----------
    n_hits : int
        Number of motif hits to return
    escore_threshold : float
        Maximum e-score threshold for motif hits
    request_timeout_seconds : float
        Maximum number of seconds allowed for each API request

    Attributes
    ----------
    name: str
        Name of the motif database
    metainfo: dict
        Any extra information regarding this motif database
    jaspar_rest_url: str
        URL for the JASPAR REST API
    n_hits: int
        Number of motif hits to return
    escore_threshold: float
        Minimum e-score threshold for motif hits
    request_timeout_seconds: float
        Maximum number of seconds allowed for each API request
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
    jaspar_rest_species_url = "https://jaspar.elixir.no/api/v1/species/"

    # Parameters for getting all species
    species_params = {"page":1,
                      "page_size":1000,
                      "release":"2024"
                      }

    def __init__(self, 
                 n_hits : int = 10, 
                 escore_threshold : float = 10**-6,
                 request_timeout_seconds : float = 30):
        """
        Initialize a JASPAR motif database adapter.

        Parameters
        ----------
        n_hits : int, optional
            Maximum number of accepted motif hits
        escore_threshold : float, optional
            Maximum inference E-value
        request_timeout_seconds : float, optional
            Maximum number of seconds allowed for each API request
        """

        if isinstance(n_hits, bool) or not isinstance(n_hits, int) or n_hits <= 0:
            raise ValueError("n_hits must be a positive integer")
        if not isinstance(escore_threshold, (int, float)) or escore_threshold <= 0:
            raise ValueError("escore_threshold must be positive")
        if not isinstance(request_timeout_seconds, (int, float)) or request_timeout_seconds <= 0:
            raise ValueError("request_timeout_seconds must be positive")
        
        self.n_hits = n_hits
        self.escore_threshold = escore_threshold
        self.request_timeout_seconds = request_timeout_seconds
        self.jaspar_species = None

    def search(self, request : MotifSearchRequest) -> list[JasparRecord]:
        """
        Search JASPAR for motif records.

        Parameters
        ----------
        request : MotifSearchRequest
            Query and ortholog evidence used for motif inference

        Returns
        -------
        list
            Accepted JASPAR motif records ordered by inference E-value
        """

        if not isinstance(request, MotifSearchRequest):
            raise TypeError("request must be a MotifSearchRequest")
        if isinstance(request.ortholog_species_tax_id, bool):
            raise ValueError("ortholog_species_tax_id must be a positive integer")
        if not isinstance(request.ortholog_species_tax_id, int) or request.ortholog_species_tax_id <= 0:
            raise ValueError("ortholog_species_tax_id must be a positive integer")

        sequence = self._window_ortholog_sequence(
            request.ortholog_sequence,
            request.ortholog_domain,
        )
        response = requests.get(
            f"{self.jaspar_rest_url}/infer/{sequence}/",
            timeout=self.request_timeout_seconds,
        )
        response.raise_for_status()
        payload = response.json()
        results = self.__get_results(payload, "JASPAR inference")

        try:
            results = sorted(results, key=lambda result: float(result["evalue"]))
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("JASPAR inference returned malformed results") from error

        records = []
        for result in results:
            try:
                inference_evalue = float(result["evalue"])
                motif_url = result["url"]
                if not isinstance(motif_url, str) or not motif_url:
                    raise ValueError
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError("JASPAR inference returned a malformed hit") from error
            if inference_evalue > self.escore_threshold:
                continue

            motif_response = requests.get(
                motif_url,
                timeout=self.request_timeout_seconds,
            )
            motif_response.raise_for_status()
            motif = motif_response.json()
            if not isinstance(motif, dict):
                raise ValueError("JASPAR returned malformed motif details")

            try:
                motif_species = motif["species"]
                species_ids = {int(species["tax_id"]) for species in motif_species}
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError("JASPAR motif details contain malformed species") from error
            if request.ortholog_species_tax_id not in species_ids:
                continue

            records.append(self.__create_record(motif, inference_evalue))
            if len(records) == self.n_hits:
                break
        return records

    def __create_record(self, motif, inference_evalue):
        try:
            matrix_id = motif["matrix_id"]
            motif_name = motif["name"]
            pfm = motif["pfm"]
            motif_class = motif["class"]
            if isinstance(motif_class, list):
                motif_class = motif_class[0]
            if not all(isinstance(value, str) and value for value in [matrix_id, motif_name, motif_class]):
                raise ValueError
            if not isinstance(pfm, dict):
                raise ValueError
        except (KeyError, IndexError, TypeError, ValueError) as error:
            raise ValueError("JASPAR returned malformed motif details") from error

        return JasparRecord(
            matrix_id=matrix_id,
            motif_name=motif_name,
            pfm=pfm,
            motif_url=f"https://jaspar.elixir.no/matrix/{matrix_id}/",
            motif_class=motif_class,
            inference_evalue=inference_evalue,
        )

    @staticmethod
    def __get_results(payload, service_name):
        if not isinstance(payload, dict) or not isinstance(payload.get("results"), list):
            raise ValueError(f"{service_name} returned a malformed response")
        return payload["results"]

    @staticmethod
    def _window_ortholog_sequence(sequence, ortholog_domain):
        """
        Limit an ortholog sequence while retaining its complete domain.

        Parameters
        ----------
        sequence : str
            Complete ortholog protein sequence
        ortholog_domain : DomainRecord
            Domain with zero-based, half-open bounds

        Returns
        -------
        str
            Sequence window of at most 2,000 residues

        Raises
        ------
        ValueError
            The sequence or domain bounds are invalid
        """

        try:
            start = ortholog_domain.start
            end = ortholog_domain.end
        except AttributeError as error:
            raise ValueError("ortholog domains must provide start and end bounds") from error

        if not isinstance(sequence, str) or not sequence:
            raise ValueError("ortholog sequences must be nonempty strings")
        if isinstance(start, bool) or isinstance(end, bool):
            raise ValueError("ortholog domain bounds must be integers")
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError("ortholog domain bounds must be integers")
        if start < 0 or end > len(sequence) or start >= end:
            raise ValueError("ortholog domain bounds must define a nonempty slice within the sequence")

        window_length = 2000
        if end - start > window_length:
            raise ValueError("ortholog domains must fit inside the JASPAR sequence window")
        if len(sequence) <= window_length:
            return sequence

        domain_center = (start + end) // 2
        window_start = domain_center - window_length // 2
        window_start = max(0, min(window_start, len(sequence) - window_length))
        return sequence[window_start:window_start + window_length]

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
            species_ids = set()
            species_url = self.jaspar_rest_species_url
            species_params = self.species_params.copy()
            while species_url is not None:
                response = requests.get(
                    species_url,
                    params=species_params,
                    timeout=self.request_timeout_seconds,
                )
                response.raise_for_status()
                payload = response.json()
                results = self.__get_results(payload, "JASPAR species")
                try:
                    species_ids.update(int(species["tax_id"]) for species in results)
                except (KeyError, TypeError, ValueError) as error:
                    raise ValueError("JASPAR species returned malformed results") from error
                next_url = payload.get("next")
                if next_url is not None and (not isinstance(next_url, str) or not next_url):
                    raise ValueError("JASPAR species returned malformed pagination")
                species_url = next_url
                species_params = None
            self.jaspar_species = species_ids

        return species_tax_id in self.jaspar_species
