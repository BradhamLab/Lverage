#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
This file contains a class that uses Clustal Omega to align multiple protein sequences.
Similarity between regions of the alignment can be calculated.

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
from Bio.Align.Applications import ClustalOmegaCommandline
import os
from Bio import SeqIO


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
class Aligner:
    '''
    Uses Clustal Omega to perform multi-alignment and has functionality to calculate similarity across these alignments in specific regions
    '''

    def __init__(self, clustal_omega_path, verbose = False):
        '''
        Arguments:
        clustal_omega_path: Path to clustal omega executable; empty if on PATH
        verbose: If steps taken should be printed
        '''
        self.clustal_omega_path = clustal_omega_path
        self.verbose = verbose

    def align(self, sequences):
        '''
        This function aligns multiple protein sequences using Clustal Omega.

        Arguments:
        sequences: list of sequences to perform multi-alignment on

        Returns:
        list of alignments
        '''

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step one: write sequences to a temporary file
        with open("temp.fa", "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq{i}\n")
                f.write(f"{seq}\n")

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step two: run Clustal Omega
        if self.verbose:
            print("\tRunning Clustal Omega.", flush=True)

        clustalomega_cline = ClustalOmegaCommandline(self.clustal_omega_path, infile="temp.fa", outfile="temp.aln", verbose=True, seqtype = "Protein", auto=True, force=True)
        clustalomega_cline()

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step three: read Clustal Omega output
        if self.verbose:
            print("\tReading Clustal Omega output.", flush=True)

        with open("temp.aln", "r") as f:
            content = f.read()

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step four: parse Clustal Omega output
        align_list = [str(x.seq) for x in list(SeqIO.parse("temp.aln", "fasta"))]


        os.remove("temp.fa")
        os.remove("temp.aln")

        return align_list


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# FOR TESTING PURPOSES
if __name__ == "__main__":
    aligner = Aligner("./clustalo")
    path = "testcases/proteins.fa"

    from Bio import SeqIO
    records = list(SeqIO.parse(path, "fasta"))
    sequences = [str(record.seq) for record in records]

    alignment_list = aligner.align(sequences)
