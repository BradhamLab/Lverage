#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
This file contains a class that uses the tool ORFfinder and parses its output to find the longest open reading frame (ORF) of a DNA sequence.

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
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime # for random file generation

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class ProteinFinder:
    '''
    This class searches for the longest open reading frame (ORF) of a DNA sequence.
    '''

    def __init__(self, table=1, start_codons=['ATG', 'CTG', 'TTG', 'GTG'], verbose = False):
        '''
        Arguments:
        table: translation table according to NCBI
        start_codons: start codons to consider
        verbose: if True, prints messages to the console
        '''
        self.table = table
        self.start_codons = start_codons
        self.verbose = verbose
        
        self.is_message_printed = False

        assert self.table >= 1, 'Invalid translation table'

    def find_protein(self, sequence):
        '''
        This function finds the protein sequence of a DNA sequence by searching for the longest open reading frame (ORF).

        Arguments:
        sequence: a DNA sequence
        '''

    def find_protein(self, sequence):
        '''
        This function finds the protein sequence of a DNA sequence by searching for the longest open reading frame (ORF).

        Arguments:
        sequence: a DNA sequence
        '''

        protein = "" # return value

        # if sequence is a string, convert it to a Seq object; if a Seq object, do nothing; otherwise, raise an error
        if isinstance(sequence, str):
            sequence = Seq(sequence)
            reverse_sequence = sequence.reverse_complement()
        elif isinstance(sequence, Seq):
            reverse_sequence = sequence.reverse_complement()
        else:
            raise ValueError('Invalid sequence type')
        
        
        # For both strands
        for strand, nuc in [(+1, sequence), (-1, reverse_sequence)]:

            # For each reading frame
            for frame in range(3):

                # Extend to the length of a complete codon frame
                length = 3 * ((len(nuc) - frame) // 3)

                # Search for a start codon
                for start in range(frame, length, 3):

                    # Check if the current codon is a start codon
                    if nuc[start:start+3] in self.start_codons:

                        # Translate the sequence from this point
                        orf = nuc[start:].translate(table=self.table, to_stop=True)
                        if len(orf) > len(protein):

                            if self.verbose:

                                if not self.is_message_printed:
                                    print('.')
                                    self.is_message_printed = True
                                else:
                                    print('.', end='')

                            protein = orf

        return protein
    
    
    def reset_message(self):
        '''
        When find_protein is called, it prints a message.
        Every time it is called, this message is extended.
        When this function is called, it prints a new message to the next line.
        '''
        self.is_message_printed = False

        

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# FOR TESTING PURPOSES   
if __name__ == "__main__":

    from Bio import SeqIO

    sequence = str(SeqIO.read("testcases/lverg.fa", "fasta").seq)

    pf = ProteinFinder()

    print(pf.find_protein(sequence))

