"""Interfaces and records for ortholog searching."""

from dataclasses import dataclass


@dataclass
class OrthologRecord:
    """
    Record describing an orthologous protein found by a searcher.

    Parameters
    ----------
    accession : str
        Protein accession identifier
    description : str
        Description of the protein
    species_name : str
        Scientific name of the species
    species_tax_id : int
        NCBI taxonomic identifier of the species
    sequence : str
        Complete protein sequence
    evalue : float
        Best BLAST expectation value
    identity : float
        Exact identity ratio from zero to one
    query_coverage : float
        Query coverage ratio from zero to one
    """

    accession : str
    description : str
    species_name : str
    species_tax_id : int
    sequence : str
    evalue : float
    identity : float
    query_coverage : float
