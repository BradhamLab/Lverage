"""
This file contains functions for general use throughout the LvERAGE pipeline.

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
    - Nofal Ouardaoui    - naouarda@bu.edu - **

    \* Principle Investigator, ** Software Developers
"""

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Imports
import numpy as np
from ete3 import NCBITaxa
from DBDScanner import DBDScanner, DBD
from Bio.Align import PairwiseAligner

class PrependedOutput:

    def __init__(self, original_stdout, prepended_text):
        self.original_stdout = original_stdout
        self.prepended_text = prepended_text

    def write(self, s):
        if s.strip():  # avoid prepending to empty lines
            self.original_stdout.write(self.prepended_text + s)
        else:
            self.original_stdout.write(s)

    def flush(self):
        self.original_stdout.flush()

def align_dbds(sequence_one : str, dbd_one : DBD, sequence_two : str, dbd_scanner : DBDScanner, aligner : PairwiseAligner):
    """This function aligns two DBDs from two different sequences.

    Parameters
    ----------
    sequence_one : str
        The first sequence.
    dbd_one : DBD
        The first DBD.
    sequence_two : str
        The second sequence.
    dbd_scanner : DBDScanner
        The DBD scanner.
    aligner : PairwiseAligner
        For aligning DBDs.

    Returns
    -------
    list
        List of alignments.

    """

    # Get the second DBD
    dbd_two_list = dbd_scanner.find_dbds(sequence_two)

    # for each DBD in the second sequence, find the corresponding DBD
    dbd_two = None
    for dbd in dbd_two_list:
        if dbd.get_name() == dbd_one.get_name():
            dbd_two = dbd
            break

    # Only align if the dbd_two is found
    alignments = None
    if dbd_two is not None:

        alignments = aligner.align(sequence_one[dbd_one.get_start():dbd_one.get_end()], sequence_two[dbd_two.get_start():dbd_two.get_end()])

    return alignments


def calculate_alignment_similarity(align1, align2):
    """This function calculates the similarity between two alignments.

    Parameters
    ----------
    align1 : str
        The first alignment.
    align2 : str
        The second alignment.

    Returns
    -------
    float
        The similarity between the two alignments.

    Raises
    ------
    ValueError
        If the alignments are not the same length.
    """

    # Alignments MUST be of same length
    if len(align1) != len(align2):
        raise ValueError("Alignments are not the same length.")

    similarity = 0 # counter for matching characters

    # for every character in the alignments
    for i in range(len(align1)):
        
        # If the characters match across alignment and are not a gap, increase similarity counter
        if align1[i] == align1[i] and align1[i] != '-':
            similarity += 1

    # returning the similarity; ratio of similar characters over alignment length
    return similarity / len(align1)


def get_species_name(tax_id):
    """NCBI Taxonomy Database is queried using the taxonomic ID to get the species name.

    Parameters
    ----------
    tax_id : int
        The taxonomic ID.

    Returns
    -------
    str
        The species name.

    Raises
    ------
    ValueError
        1. If the tax ID is a string but is not numeric.
        2.If the tax ID is not a positive number.
    TypeError
        If the tax ID is not a number
    """
    

    if type(tax_id) == str:
        if not tax_id.isdigit():
            raise ValueError("Tax ID must be a number.")
        tax_id = int(tax_id)
    
    elif type(tax_id) != int:
        raise TypeError("Tax ID must be a number.")
    
    if tax_id <= 0:
        raise ValueError("Tax ID must be a positive number.")

    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(tax_id)
    names = ncbi.get_taxid_translator(lineage)

    return names.get(tax_id, None)


def get_tax_id(species_name):
    """NCBI Taxonomy Database is queried using the species name to get the taxonomic ID.

    Parameters
    ----------
    species_name : str
        The species name.

    Returns
    -------
    int
        The taxonomic ID.
    None
        If the species name is not found in the database.

    Raises
    ------
    TypeError
        If the species name is not a string.
    """


    if type(species_name) != str:
        raise TypeError("Species name must be a string.")

    ncbi = NCBITaxa()
    name_dict = ncbi.get_name_translator([species_name])
    name_list = name_dict.get(species_name, None)

    return name_list[0] if name_list is not None else None


def check_valid_species(species_name):
    """NCBI Taxonomy Database is queried to check if the species name is valid.

    Parameters
    ----------
    species_name : str
        The species name.

    Returns
    -------
    bool
        True if the species name is valid.
        False if the species name is invalid.
    """

    ncbi = NCBITaxa()

    r = False # result

    # If the species name is a number, then it is a taxonomic ID
    if species_name.isdigit():
        tax_id = int(species_name)
        r = ncbi.get_rank([tax_id]).get(tax_id, None) is not None
    
    # If the species name is a string, then it is a species name
    else:
        name_dict = ncbi.get_name_translator([species_name])
        r = name_dict.get(species_name, None) is not None

    return r




