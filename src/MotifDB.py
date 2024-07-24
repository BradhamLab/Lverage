#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
This file contains an interface to Motif databases.

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
import requests
import src.utils as lvutils
from Bio.Align import PairwiseAligner
from src.DBDScanner import DBDScanner, DBD

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class MotifRecord:
    '''
    This abstract class represents a motif record
    This should not be instantiated
    '''

    header_list = ["Matrix ID"]

    def __init__(self, matrix_id):
        '''
        Arguments:
        matrix_id: a string representing the matrix ID of the motif
        '''

        self.matrix_id = matrix_id

    def get_matrix_id(self):
        return self.matrix_id
    
    def get_values(self):
        raise [self.matrix_id]
    
    @staticmethod
    def get_test_motif():
        return MotifRecord("test")

class JasparRecord(MotifRecord):
    '''
    This class represents a JASPAR motif record
    '''

    header_list = MotifRecord.header_list + ["Motif Name", "Motif PFM", "Motif Logo Link", "Motif Class", "UniProt Validation ID", "Gene-UniProt Percent Identity", "Gene_DBD-UniProt_DBD Percent Identity"]

    def __init__(self, matrix_id, name, pfm, logo_link, motif_class, uniprot_id, sequence_perc_id, dbd_perc_id):
        '''
        Arguments:
        matrix_id: a string representing the matrix ID of the motif
        name: a string representing the name of the motif
        pfm: a list of lists representing the position frequency matrix of the motif
        logo_link: a string representing the link to the motif's logo
        motif_class: a string representing the class of the motif
        uniprot_id: a string representing the UniProt ID of the motif
        sequence_perc_id: a float representing the global sequence percent identity of the gene of interest and the validation sequence from uiprot
        dbd_perc_id: a float representing the DNA-binding domain percent identity between the gene of interest and the validation sequence from uniprot
        '''
        
        super().__init__(matrix_id)
        self.name = name
        self.pfm = pfm
        self.logo_link = logo_link
        self.motif_class = motif_class
        self.uniprot_id = uniprot_id
        self.sequence_perc_id = sequence_perc_id
        self.dbd_perc_id = dbd_perc_id

    def get_values(self):
        r = [self.matrix_id, self.name, self.pfm, self.logo_link, self.motif_class, self.uniprot_id, self.sequence_perc_id, self.dbd_perc_id]
        assert len(r) == len(JasparRecord.header_list), "Length of values does not match length of header"
        return r
    
    @staticmethod
    def get_test_motif():
        import random
        l = range(1, 400)
        matrix_id = random.choice(l)

        # pad the matrix_id with 0 up to 4 digits
        matrix_id = str(matrix_id).zfill(4)

        return JasparRecord("Sample Matrix ID", "Sample Motif Name", [[1, 1, 1, 1], [1, 1, 1, 1]], f"https://jaspar2020.genereg.net/static/logos/all/MA{matrix_id}.1.png", "Sample Motif Class", "Sample UniPort Validation ID", 0.9, 0.9)

class MotifDBInterface:
    '''
    This class serves as an interface for motif databases.
    '''

    name = ""

    help_message = ""

    def __init__(self):
        pass

    def search(self, description, tax_id):
        '''
        This function searches the database for a given description and returns the corresponding motif

        Arguments:
        description: a string representing a description of a gene
        tax_id: a string representing the taxonomy ID of the species
        '''
        raise NotImplementedError
    

class JasparDB(MotifDBInterface):

    name = "JASPAR"

    jaspar_rest_url = "https://jaspar.elixir.no/api/v1"
    jaspar_logo_url = "https://jaspar2020.genereg.net/static/logos/all/"

    uniprot_rest_url = "https://rest.uniprot.org/uniprotkb/"

    help_message =  """
                    JASPAR is a database of transcription factor binding profiles.
                    The profile inference tool requires a protein sequence to be less than 2000 amino acids.
                    As such, if any protein sequence is longer, we create a window of 2000 amino acids including the DNA-binding domain.
                    """

    def __init__(self, n_hits = 10, dbd_threshold = 0.85, escore_threshold = 10**-6, dbd_scanner = None, email = None):
        '''
        Arguments:
        n_hits: an integer representing the number of hits to return
        dbd_threshold: a float representing the DNA-binding domain percent identity threshold
        escore_threshold: a float representing the e-score threshold
        dbd_scanner: a DBDScanner object to use for finding DNA-binding domains for validation. If not provided, a new DBDScanner object will be created.
        email: a string representing the email to use for the DBDScanner object if one needs to be created; if neither dbd_scanner nor email is provided, and error is raised.
        '''
        super().__init__()

        assert n_hits > 0, f"n_hits must be greater than 0, currently {n_hits} is given."
        assert 0 <= dbd_threshold <= 1, f"dbd_threshold must be between 0 and 1, currently {dbd_threshold} is given."

        self.n_hits = n_hits
        self.dbd_threshold = dbd_threshold
        self.escore_threshold = escore_threshold
        
        if dbd_scanner is None and email is None:
            raise ValueError("Either a DBDScanner object or an email must be provided.")
        
        self.dbd_scanner = dbd_scanner if dbd_scanner is not None else DBDScanner(email = email)

        self.global_aligner = PairwiseAligner()
        self.global_aligner.mode = 'global'
        self.global_aligner.match_score = 2
        self.global_aligner.mismatch_score = -1

        

    def search(self, ortholog_seq, ortholog_tax_id, protein_seq, dbd):
        '''
        This method searches the JASPAR database using its protein inference tool to find motifs given an ortholog's protein sequence.

        Arguments:
        ortholog_seq: the protein sequence of an ortholog
        ortholog_tax_id: the taxonomy ID of the species the ortholog belongs to
        protein_seq: the protein sequence of the gene the user is trying to find motifs for
        dbd: DNA-binding domain (DBD) object corresponding to protein_seq;
        '''

        assert len(ortholog_seq) > 0, "Ortholog sequence must be non-empty"
        assert ortholog_tax_id.isdigit(), "Taxonomy ID must be a number"
        assert len(protein_seq) > 0, "Protein sequence must be non-empty"
        assert isinstance(dbd, DBD), "DBD must be a DBD object"

        r = [] # list of motifs to return

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # If the ortholog sequence is longer than 2000 amino acids, we create a window of 2000 amino acids including the DNA-binding domain

        if len(ortholog_seq) > 2000:
            dbd_start = dbd.get_start()
            dbd_end = dbd.get_end()

            # Get the new start; if the dbd_start is less than 1000, set the new start to 0
            new_start = max(0, dbd_start - 1000)

            # Calculate how many amino acids are left to use out of 2000
            remaining = 2000 - dbd_end - new_start

            # Get the new end; if the dbd_end + remaining is greater than the length of the ortholog sequence, set the new end to the length of the ortholog sequence
            new_end = min(len(ortholog_seq), dbd_end + remaining)

            # Get remaining amino acids to use
            remaining = 2000 - (new_end - new_start)

            # If there are still remaining amino acids, add them to the start
            new_start = max(0, new_start - remaining)

            ortholog_seq = ortholog_seq[new_start:new_end]


        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # search for motif using ortholog sequence and species
        infer_query = f'{self.jaspar_rest_url}/infer/{ortholog_seq}'
        matrix_download = requests.get(infer_query)

        # If the matrix download fails, return an empty list
        if matrix_download.ok:
            data = matrix_download.json()['results']

            # Sort the data list by e-score
            data.sort(key = lambda x : x['evalue'])

            # For each hit motif
            for i, hit in enumerate(data):

                if i >= self.n_hits:
                    break

                # Check if the e-score is below the threshold
                if hit['evalue'] > self.escore_threshold:
                    print(hit['evalue'], 'evalue is bad')
                    continue

                # retrieve the motif information
                motif_query = hit['url']
                motif_download = requests.get(motif_query)

                # If the motif download fails, skip this motif
                if not motif_download.ok:
                    continue

                # get motif information
                motif = motif_download.json()

                # If the motif's tax_id does not match the ortholog's tax_id, skip this motif
                motif_tax_id = str(motif['species'][0]['tax_id'])
                if motif_tax_id != ortholog_tax_id:
                    continue

                matrix_id = motif['matrix_id']
                pfm = motif['pfm']
                motif_class = motif['class'][0]
                uniprot_ids = motif['uniprot_ids'] 
                motif_name = motif['name']
                logo = self.jaspar_logo_url + matrix_id + ".png"


                #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
                # For each uniprot id, VALIDATE THE MOTIF with a percent identity threshold
                best_validation = None
                for uniprot_id in uniprot_ids:

                    # obtain the uniprot sequence
                    uniprot_download = requests.get(self.uniprot_rest_url + uniprot_id + ".fasta")

                    if not uniprot_download.ok:
                        continue

                    uniprot_fasta = uniprot_download.text.split('>')[1] # get the first fasta entry
                    uniprot_seq = ''.join(uniprot_fasta.split('\n')[1:]) # remove the header and join the sequence

                    # globally align the uniprot sequence with the protein sequence
                    alignments = self.global_aligner.align(protein_seq, uniprot_seq)[0]
                    sequence_perc_id = lvutils.calculate_alignment_similarity(*alignments)

                    # globally align the DBDs of the two sequences.
                    alignments = lvutils.align_dbds(protein_seq, dbd, uniprot_seq, self.dbd_scanner, self.global_aligner)
                    dbd_perc_id = lvutils.calculate_alignment_similarity(*alignments[0]) if alignments is not None else 0.0 

                    # If the dbd percent identity is below the threshold, skip this uniprot id
                    if dbd_perc_id < self.dbd_threshold:
                        # continue
                        pass

                    # If the validation is the best so far, update the best validation
                    if best_validation is None or dbd_perc_id > best_validation[2]:
                        best_validation = (uniprot_id, sequence_perc_id, dbd_perc_id)

                if best_validation:
                    uniprot_id, sequence_perc_id, dbd_perc_id = best_validation
                    r.append(JasparRecord(matrix_id, motif_name, pfm, logo, motif_class, uniprot_id, sequence_perc_id, dbd_perc_id))

        return r


        

class MotifDBFactory:
    '''
    This class is a factory that returns a motif database based on the database name
    '''

    # Database name -> MotifDB class
    motif_db_dict = {
        "JASPAR": JasparDB
    }

    # Database name -> MotifRecord class
    motif_record_dict = {
        "JASPAR": JasparRecord
    }

    @staticmethod
    def get_motif_db(db_name, **kwargs):
        '''
        This function returns a motif database object based on the database name

        Arguments:
        db_name: a string representing the name of the database
        **kwargs: keyword arguments to pass to the motif database
        '''

        assert db_name in MotifDBFactory.motif_db_dict, "Database name not found"

        return MotifDBFactory.motif_db_dict[db_name](**kwargs)
    
    @staticmethod
    def get_motif_db_names():
        '''
        This function returns the names of all motif databases
        '''
        return list(MotifDBFactory.motif_db_dict.keys())
    
    @staticmethod
    def get_dbs():
        return list(MotifDBFactory.motif_db_dict.values())

    @staticmethod
    def has_db(db_name):
        return db_name in MotifDBFactory.motif_db_dict
    
    @staticmethod
    def get_db_record(db_name):
        '''
        This function returns the MotifRecord class for a given database name.
        Note that this is the CLASS and not an instance of the class.
        
        Arguments:
        db_name: a string representing the name of the database
        '''


        assert db_name in MotifDBFactory.motif_record_dict, "Database name not found"


        return MotifDBFactory.motif_record_dict[db_name]



#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# FOR TESTING PURPOSES
if __name__ == "__main__":
    from DBDScanner import DBD

    mdb = MotifDBFactory.get_motif_db(db_name="JASPAR", n_hits=10)
    l = mdb.search("homeodomain", "9606", "MKDDKKMYCYQCSTIHHAGSPAAAAHLSNGGVCPHACDDSPYSELSYGGDLDETFARRKQRRNRTTFTVQQLEELESAFAKTHYPDVFTREDLALRINLTEARVQVWFQNRRAKWRKAERTKQERGPSSTSSPENEDRLSSSEGAVGQSREELNMSPGEPSEERRRDKMTVDGEKSDSQNDDGSPLHSLDRAPSSGRLMPPTFANPSASTMLNPFYHPGGAARLLLASQPYESLRGHGSSARFPSLISPSYASQLMSFASARKEPVPTSGGS", DBD("Homeodomain", 60, 57))

    for rec in l:
        print(rec.get_matrix_id())
        print(rec.get_name())
        print(rec.get_uniprot_id())
        print(rec.get_global_percent_identity())
        print(rec.get_dbd_percent_identity())
        print("\n")



