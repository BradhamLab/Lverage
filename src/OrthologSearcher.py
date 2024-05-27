#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
This file contains a class that is used to search for orthologous sequences of a given sequence using blastp

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
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from itertools import product
from Bio import Entrez, SeqIO

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class AlignmentRecord:
    '''
    BLAST alignment record, contains useful information about a BLAST hit
    '''

    # List of words that are uneeded in the name of the BLAST description; we can remove them
    uneeded_words_list = ['transcription factor']


    def __init__(self, alignment, query_coverage, sequence):        
        '''
        Arguments:
        alignment: a Bio.Blast.Record.Alignment object
        query_coverage: the percentage of the query sequence that is covered by the alignment
        sequence: the sequence of the query
        '''                        
        


        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Get title of alignment
        self.title = alignment.title

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Get name of species
        left_bracket = self.title.find("[")
        right_bracket = self.title.find("]")
        self.species = self.title[left_bracket + 1:right_bracket]

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Get name of sequence
        left_bar = self.title.find("|")
        right_bar = self.title.find("|", left_bar + 1)
        self.name = self.title[right_bar + 2: left_bracket - 1]
        
        # removing unwanted descriptions from name
        self.name.replace(', partial', '')
        self.name.replace('protein', '')    

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Get the first alignment sequence without gaps
        self.seq = sequence

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Getting other information
        self.accession = alignment.accession

        self.expect = alignment.hsps[0].expect

        hsp = alignment.hsps[0]

        self.percent_identity = (hsp.identities / hsp.align_length) * 100

        self.query_coverage = query_coverage

        
        

    def __str__(self):
        return f"{self.title} {self.name} {self.species} {self.seq}"
    
    def get_title(self):
        return self.title
    
    def get_name(self):
        return self.name
    
    def get_usable_name(self):
        '''
        This function returns the name of the sequence without any uneeded words
        Also removes parentheses if they are at the beginning and end of the name
        '''

        usable_name = self.name
        for uneeded in self.uneeded_words_list:
            usable_name = usable_name.replace(uneeded, '')
            
        usable_name = ' '.join(usable_name.split()) # Removing extra spaces

        # remove parentheses if they are at the beginning and end of the name
        if usable_name.startswith('(') and usable_name.endswith(')'):
            usable_name = usable_name[1:-1]

        return usable_name
    
    def get_usable_name_combos(self):
        '''

        This function returns all possible combinations of the usable name of the sequence
        '''
        usable_name = self.get_usable_name()
        words = usable_name.split()
        
        c = []
        for combination in product(*[(word, '') for word in words]):
            r = ' '.join(filter(None, combination))
            
            # if the combination is not empty and if the combination isn't just a number
            if r and not r.replace(' ', '').isdigit():
                c.append(r)

        c.sort(key=lambda x: len(x.split()), reverse=True)


        return c
    
    def get_species(self):
        return self.species
    
    def get_seq(self):
        return self.seq

    def get_accession(self):
        return self.accession
    
    def get_escore(self):
        return self.expect
    
    def get_percent_identity(self):
        return self.percent_identity
    
    def get_query_coverage(self):
        return self.query_coverage
    

class OrthologSearcher:
    '''
    Searches for orthologous sequences of a given sequence using blastp
    '''
    bad_words = ['hypothetical', 'unnamed', 'uncharacterized', 'unknown', 'partial', 'isoform'] # Words that we don't want in the name of the sequence; discard the alignment if it contains any of these words


    def __init__(self, species_list=['homo sapiens', 'xenopus laevis', 'drosophila melanogaster', 'mus musculus'], hitlist_size = 100, email = None, verbose = False):
        '''
        Arguments:
        species_list: list of orthologous species to look through
        hitlist_size: how many hits should we look through
        email: email to use for Entrez
        verbose: should we print out the current status and results?
        '''

        assert hitlist_size > 0, "hitlist_size must be greater than 0"
        assert email is not None, "Email must be provided!"


        self.species_list = [x.lower() for x in species_list]
        self.entrez_query = f'({" OR ".join([species.lower() + "[ORGN]" for species in species_list])})'
        self.hitlist_size = hitlist_size
        self.email = email
        self.verbose = verbose

        Entrez.email = self.email


    def search(self, sequence):
        '''
        This function searches for orthologous sequences of a given sequence using blastp.

        Arguments:
        sequence: a string representing a protein sequence
        '''

        assert len(sequence) > 0, "Sequence must be non-empty"

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step one: run blastp
        if self.verbose:
            print("\tRunning NCBI blastp.", flush=True)

        result_handle = NCBIWWW.qblast("blastp", "nr", sequence, entrez_query=self.entrez_query, hitlist_size=self.hitlist_size)

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step two: read blastp output
        if self.verbose:
            print("\tParsing NCBI blastp output", flush=True)

        self.blast_record = NCBIXML.read(result_handle)
        blast_record = self.blast_record

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step three: filter blastp output

        if self.verbose:
            print("\tFiltering hits.")
            
        filtered_hits = []

        # default sorted by e-score
        for alignment in blast_record.alignments:
 
            handle = Entrez.efetch(db="protein", id=alignment.accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()

            ortholog_sequence = str(record.seq)

            query_coverage = alignment.hsps[0].align_length / len(sequence)

            record = AlignmentRecord(alignment, query_coverage, ortholog_sequence)

            # E-score must be less than 10^-6. If not, we stop searching completely as the rest of the hits will be worse
            if record.expect > 10**-6:
                break

            # If the title contains any of the bad words, we skip this alignment
            if not all(sub not in record.name for sub in self.bad_words):
                continue
        
            # Check if species is in species list
            is_species = False
            for species in self.species_list:
                if species in record.get_species().lower():
                    is_species = True
                    break

            if is_species:

                filtered_hits.append(record)

                if self.verbose:
                    print(f"\tObtained ortholog from {record.get_species()}: {record.get_name()}, with e-score {record.expect} and similarity {record.get_percent_identity()}", flush=True)

        return filtered_hits

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# FOR TESTING PURPOSES
if __name__ == "__main__":
    import sys

    sequence ="LQRTTRSXLALLLCMWPVTRTNDRVLVWFQNRRAKWRKRERFQQFQNMRGLGPGSGYEMPIAPRPDAYSQVNSPGFHVLGDTHQPPAPVEGAMLRICRNLQNLRREFDSRKIGCHPSSSVGGTPGTSTTNESQDTSNHSSMIHQSSPWATAANMASPLASSMSPVGQQPQMPGQNPINSCMAPQSTLPSFMGVPAHQMNNTGVMNPMSNMTSMTSMPTSMPPSSGTAPVSSPSSNFMSSVGGLNAAAYNGQYTDMHPTVEGVGGVDRRTNSIAALRLRAKEHSSVMGMMNGYS"

    ortho = OrthologSearcher(verbose=True)

    hits = ortho.search(sequence)
    
