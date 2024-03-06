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

Correspondence: anthonygarza124@gmail.com OR ...cyndi's email here...
'''
#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Imports
import numpy as np


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Functions

def calculate_similarity(alignment_list, start = None, size = None):
    '''
    This function calculates the similarity between the first alignment and every other alignment
    If start and size are given, then calculate similarity in a region defined by where start and start + size in the original sequence of the first alignment are (mapped to alignments)

    Arguments:
    alignment_list: a list of alignments
    start: start index of region in original sequence of first alignment
    size: size of region in original sequence of first alignment

    Returns:
    numpy array of similarities, where each value is the similarity of the corresponding alignment to the first alignment.
    '''

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # validating arguments
    assert start + size <= len(alignment_list[0]), f"End index can not be larger than the first alignment: {size} vs. {len(alignment_list[0])}"

    region_list = []

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # if start and size are given (corresponding to first alignment), grab those regions across all alignments. Else, grab entire region
    if start is not None and size is not None:

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step one: map the start to the first alignment

        first_align = alignment_list[0] # First alignment we map start and start + size to

        count = 0
        align_start = -1
        for i, c in enumerate(first_align):
            if c != '-':
                if count == start:
                    align_start = i
                    break
                
                count += 1

        assert align_start != -1, "Start index not found."

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step two: find the end of the region
        count = 0
        align_end = -1
        for i, c in enumerate(first_align[align_start:]):
            if c != '-':
                if count == size:
                    align_end = i + align_start
                    break
                
                count += 1

        assert align_end != -1, "End index not found."

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step three: get the regions
        for align in alignment_list:
            region_list.append(align[align_start:align_end])

    
    else:
        region_list = alignment_list


    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step four: calculate similarity to first region
    similarity_array = np.zeros((len(region_list)), dtype=np.float32)

    similarity_array[0] = 1
    for i in range(1, len(region_list)):
        similarity_array[i] = _calculate_similarity(region_list[0], region_list[i])

    
    return similarity_array



def _calculate_similarity(align1, align2):
    '''
    This function calculates the similarity between two alignments.
    Hidden function; use calculate_similarity to calculate across list of alignments instead.

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