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

Correspondence: Correspondence: 
    Cynthia A. Bradham - cbradham@bu.edu - *
    Anthony B. Garza   - abgarza@bu.edu  - **
    Stephanie P. Hao   - sphao@bu.edu    - **
    Yeting Li          - yetingli@bu.edu - **
    Nofal Ouardaoui    - naouarda@bu.edu - **

    *  - Principle Investigator
    ** - Software Developers
'''
#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Imports
from Bio import SeqIO
from Bio.Seq import Seq

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class ProteinTranslator:
    '''
    This class searches for the longest open reading frame (ORF) of a DNA sequence, and translates the DNA to protein.
    '''

    def __init__(self, table=1, start_codons=['ATG', 'CTG', 'TTG', 'GTG']):
        '''
        Arguments:
        table: translation table according to NCBI
        start_codons: start codons to consider
        '''

        self.table = table
        self.start_codons = start_codons
        
        assert self.table >= 1, 'Invalid translation table'


    def translate(self, sequence):
        '''
        This method finds the protein sequence of a DNA sequence by searching for the longest open reading frame (ORF).

        Arguments:
        sequence: a DNA sequence

        Returns:
        protein: the protein sequence of the DNA sequence
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
                        nuc_to_translate = nuc[start:]

                        # pad the sequence with Ns to make it a multiple of 3
                        nuc_to_translate += 'N' * (3 - len(nuc_to_translate) % 3)

                        orf = nuc_to_translate.translate(table=self.table, to_stop=True)
                        if len(orf) > len(protein):

                            protein = orf

        return str(protein)

        

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# FOR TESTING PURPOSES   
if __name__ == "__main__":

    from Bio import SeqIO

    sequence = str((list(SeqIO.parse("../Data/TestCases/Lv-Alx1.fa", "fasta")))[0].seq)

    pf = ProteinTranslator()

    print(pf.translate(sequence))