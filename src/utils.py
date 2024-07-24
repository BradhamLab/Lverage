#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
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

Correspondence: Correspondence: 
    Cynthia A. Bradham - cbradham@bu.edu - *
    Anthony B. Garza   - abgarza@bu.edu  - **
    Stephanie P. Hao   - sphao@bu.edu    - **
    Yeting Li          - yetingli@bu.edu - **
    Nofal Ouardaoui    - naouarda@bu.edu - **

    * - Principle Investigator
    ** - Software Developers
'''
#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Imports
import numpy as np
from ete3 import NCBITaxa
from src.DBDScanner import DBDScanner, DBD
from Bio.Align import PairwiseAligner



#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Functions

def align_dbds(sequence_one : str, dbd_one : DBD, sequence_two : str, dbd_scanner : DBDScanner, aligner : PairwiseAligner):
    '''
    This function aligns the DBDs of two sequences globally.
    The DBD of the first sequence is already provided.
    The second DBD is found using a DBDScanner object.

    Arguments:
    sequence_one: first sequence
    dbd_one: DBD object (DBDScanner.DBD) of first sequence
    sequence_two: second sequence
    dbd_scanner: DBDScanner object for second sequence
    aligner: PairwiseAligner object for alignment

    Returns:
    Alignments of the two DBDs or None if the second DBD is not found.
    '''

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
    '''
    This function calculates the similarity between two alignments.

    Arguments:
    align1: first alignment
    align2: second alignment, of same length as align1

    Returns:
    float value of similarity percentage
    '''

    # Alignments MUST be of same length
    assert len(align1) == len(align2), "Alignments are not the same length."

    similarity = 0 # counter for matching characters

    # for every character in the alignments
    for i in range(len(align1)):
        
        # If the characters match across alignment and are not a gap, increase similarity counter
        if align1[i] == align1[i] and align1[i] != '-':
            similarity += 1

    # returning the similarity; ratio of similar characters over alignment length
    return similarity / len(align1)


def get_species_name(tax_id):
    '''
    This function returns the species name of a given taxonomic ID.
    
    Arguments:
    tax_id: taxonomic ID of species
    
    Returns:
    species name of taxonomic ID
    '''

    if type(tax_id) == str:
        assert tax_id.isdigit(), "Tax ID must be a number."
        tax_id = int(tax_id)
    
    elif type(tax_id) != int:
        assert False, "Tax ID must be a number."
    
    assert tax_id > 0, "Tax ID must be a positive number."

    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(tax_id)
    names = ncbi.get_taxid_translator(lineage)

    return names.get(tax_id, None)


def get_tax_id(species_name):
    '''
    This function returns the taxonomic ID of a given species name.
    
    Arguments:
    species_name: name of species
    
    Returns:
    taxonomic ID of species
    '''

    assert type(species_name) == str, "Species name must be a string."

    ncbi = NCBITaxa()
    name_dict = ncbi.get_name_translator([species_name])
    name_list = name_dict.get(species_name, None)

    return name_list[0] if name_list is not None else None


def check_valid_species(species_name):
    '''
    This function checks if a species name is valid.
    
    Arguments:
    species_name: name of species or taxonomic ID
    
    Returns:
    boolean value of whether species name is valid
    '''

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




