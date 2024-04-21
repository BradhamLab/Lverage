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

Correspondence: anthonygarza124@gmail.com OR ...cyndi's email here...
'''
#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Imports
import requests
import src.utils
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class MotifRecord:
    '''
    This class represents a motif record
    '''

    def __init__(self, matrix_id, name, pfm, logo_link, motif_class, uniprot_id, global_perc_id, dbd_perc_id):
        self.matrix_id = matrix_id
        self.name = name
        self.pfm = pfm
        self.logo_link = logo_link
        self.motif_class = motif_class
        self.uniprot_id = uniprot_id
        self.global_perc_id = global_perc_id
        self.dbd_perc_id = dbd_perc_id

    def get_matrix_id(self):
        return self.matrix_id

    def get_name(self):
        return self.name
    
    def get_pfm(self):
        return self.pfm
    
    def get_logo_link(self):
        return self.logo_link
    
    def get_motif_class(self):
        return self.motif_class
    
    def get_uniprot_id(self):
        return self.uniprot_id
    
    def get_global_percent_identity(self):
        return self.global_perc_id
    
    def get_dbd_percent_identity(self):
        return self.dbd_perc_id
    



class MotifDBInterface:
    '''
    This class serves as an interface for motif databases.
    '''

    name = ""

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

    jaspar_rest_url = "https://jaspar2020.genereg.net/api/v1/"
    jaspar_logo_url = "https://jaspar2020.genereg.net/static/logos/all/"

    uniprot_rest_url = "https://rest.uniprot.org/uniprotkb/"

    def __init__(self, n_hits = 10, dbd_threshold = 0.85):
        '''
        
        '''
        super().__init__()

        assert n_hits > 0, f"n_hits must be greater than 0, currently {n_hits} is given."
        assert 0 <= dbd_threshold <= 1, f"dbd_threshold must be between 0 and 1, currently {dbd_threshold} is given."

        self.n_hits = n_hits
        self.dbd_threshold = dbd_threshold

        self.local_aligner = PairwiseAligner()
        self.local_aligner.mode = 'local'

        self.global_aligner = PairwiseAligner()
        self.global_aligner.mode = 'global'
        self.global_aligner.match_score = 2
        self.global_aligner.mismatch_score = -1



    def search(self, description, tax_id, protein_seq, dbd):
        '''
        This function searches the JASPAR database for a given description and gets the corresponding motif

        Arguments:
        description: a string representing a description of a gene
        tax_id: a string representing the taxonomy ID of the species
        protein_seq: the protein sequence of a gene from the species the user is trying to find motifs for
        dbd: DNA-binding domain object corresponding to the protein-seq; should be able to retrieve start and size of DBD

        Returns:
        A list of validated (by calculating identity) motifs
        '''

        assert len(description) > 0, "Description must be non-empty"
        assert tax_id.isdigit(), "Taxonomy ID must be a number"

        r = [] # list of motifs to return
        found_exact = False # condition to break the loop if a perfect match is found

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # search for motif using description and species
        query = f"matrix/?tax_id={tax_id}&search={description}&page_size={self.n_hits}"
        matrix_download = requests.get(self.jaspar_rest_url + query)


        if matrix_download.ok:
            data = matrix_download.json()['results']

            # Sort the data list based on similarity to the description variable
            data = sorted(data, key=lambda motif: motif['name'] != description)

            # For each hit motif
            for motif in data:
                matrix_id = motif['matrix_id']
                name = motif['name']

                logo = self.jaspar_logo_url + matrix_id + ".png"

                profile_download = requests.get(self.jaspar_rest_url + "matrix/" + matrix_id)

                if not profile_download.ok:
                    continue

                profile = profile_download.json()
                pfm = profile['pfm']
                motif_class = profile['class']
                uniprot_ids = profile['uniprot_ids']

                #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
                # For each uniprot id, VALIDATE THE MOTIF with a percent identity threshold
                for uniprot_id in uniprot_ids:

                    # obtain the uniprot sequence
                    uniprot_download = requests.get(self.uniprot_rest_url + uniprot_id + ".fasta")

                    if not uniprot_download.ok:
                        continue

                    uniprot_fasta = uniprot_download.text.split('>')[1] # get the first fasta entry
                    uniprot_seq = ''.join(uniprot_fasta.split('\n')[1:]) # remove the header and join the sequence

                    # globally align the uniprot sequence with the protein sequence
                    alignments = self.global_aligner.align(protein_seq, uniprot_seq)[0]
                    global_perc_id = utils._calculate_similarity(*alignments)

                    # locally align the uniprot sequence with the protein sequence
                    alignments = self.local_aligner.align(protein_seq, uniprot_seq)[0]

                    # If the dbd has low similarity with the unitprot sequence, skip this motif
                    dbd_perc_id = utils.calculate_similarity(alignments, dbd.get_start(), dbd.get_size())[1]

                    if dbd_perc_id < self.dbd_threshold:
                        continue

                    r.append(MotifRecord(matrix_id, name, pfm, logo, motif_class, uniprot_id, global_perc_id * 100, dbd_perc_id * 100))

                    # If the motif found has an exact match to the description, break the loop and return the motif
                    if name == description:
                        found_exact = True
                        break

                if found_exact:
                    break

        return r




        

class MotifDBFactory:
    '''
    This class is a factory that returns a motif database based on the database name
    '''

    # Database name -> MotifDB class
    motif_db_dict = {
        "JASPAR": JasparDB
    }

    @staticmethod
    def get_motif_db(db_name, **kwargs):
        '''
        This function returns a motif database based on the database name

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



