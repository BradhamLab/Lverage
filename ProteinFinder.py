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
from datetime import datetime # for random file generation

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class ProteinFinder:
    '''
    This class uses the tool ORFfinder and parses its output to find the longest open reading frame (ORF) of a DNA sequence.s
    '''

    def __init__(self, orf_finder_path, verbose = False):
        '''
        Arguments:
        orf_finder_path: path to ORFfinder
        verbose: if True, prints messages to the console
        '''
        self.orf_finder_path = orf_finder_path
        self.verbose = verbose
        
        self.is_message_printed = False

        assert shutil.which(self.orf_finder_path) is not None or os.path.exists(self.orf_finder_path), "ORFfinder not found."

    def find_protein(self, sequence):
        '''
        This function finds the protein sequence of a DNA sequence using ORF Finder.
        It goes through all open reading frames and returns the longest one.

        Arguments:
        sequence: a DNA sequence
        '''

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step one: write sequence to a temporary file

        time_string = datetime.now().strftime('%Y%m%d%H%M%S')
        temp_in_name = f'temp{time_string}.txt'
        temp_out_name = f'temp_orf{time_string}.txt'

        with open(temp_in_name, "w") as f:
            f.write(sequence)

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step two: run ORF Finder

        if self.verbose:
            if not self.is_message_printed:
                print("\tRunning ORFfinder", end='')
                self.is_message_printed = True
            else:
                print('.', end='')
                
            
        os.system(f"{self.orf_finder_path} -in {temp_in_name} -out {temp_out_name} -outfmt 0 -s 1")

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step three: find longest ORF  
                  
        orf = ""
        for rec in SeqIO.parse(temp_out_name, 'fasta'):
            seq = str(rec.seq)
            orf = seq if len(seq) > len(orf) else orf
            
        os.remove(temp_in_name)
        os.remove(temp_out_name)

        return orf
    
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

    pf = ProteinFinder("./ORFfinder")

    print(pf.find_protein(sequence))

