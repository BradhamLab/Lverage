#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
The purpose of this script is to build a directory of FASTA files for each gene in the LvEDGE database.
Multi-FASTA format files are used if a gene has multiple scaffolds.
The name of each FASTA file is: <gene name>.fa

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

# Imports
from src.LvEdgeScraper import LvEdgeScraper
import numpy as np
import os
import argparse
from Bio import SeqIO

# Classes
class SpeciesDBInterface:
    '''
    This class is an interface for species databases.
    '''
    name = ""

    def __init__(self):
        self._current_index = 0  # Used for iteration


    def __iter__(self):
        self._current_index = 0  # Reset the index each time iter is called
        return self

    def __next__(self):
        raise NotImplementedError

    

class GreenSeaUrchinDB(SpeciesDBInterface):
    '''
    This class is a database for the green sea urchin.
    '''

    name = "green sea urchin"

    def __init__(self, file_path, skip_list = []):
        '''
        Arguments:
        file_path: path to the file containing the LvEdge IDs
        skip_list: any gene ids to skip
        '''
        super().__init__()
        self.file_path = file_path
        self.skip_list = skip_list

        # read id database in; first row is header, column 2 and 3 are spu id and gene id
        self.id_db = np.loadtxt(self.file_path, dtype=str, skiprows=1, delimiter=',')

        # Initialize LvEdgeScraper to read database online
        self.scraper = LvEdgeScraper()

    def __next__(self):

        if self._current_index < len(self.id_db):

            #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
            # Get the spu id and gene id

            while True:
                _, spu_id, gene_id = self.id_db[self._current_index]
                spu_id, gene_id = spu_id.strip('"'), gene_id.strip('"')

                # Skip gene if in skip list
                if gene_id in self.skip_list:
                    self._current_index += 1

                    if self._current_index >= len(self.id_db):
                        raise StopIteration
                else:
                    break

            self._current_index += 1

            #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
            # Get the list of rows for the spu id
            for i in range(3):
                try:
                    matrix = self.scraper.search(spu_id)
                    break
                except Exception as e:
                    raise(e)

            # If no rows for the spu id was found, return an empty list   
            if not matrix:
                sequence_list = []

            else:
                #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
                # Get the first gene hit
                row = matrix[0]

                sequence_list = row[3]

            return gene_id, sequence_list
        
        else:
            raise StopIteration
        

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
        '''The purpose of this script is to build a directory of FASTA files for each gene in the LvEDGE database.
        Multi-FASTA format files are used if a gene has multiple scaffolds.
        The name of each FASTA file is: <gene name>.fa
        Copyright: 
        Lab: Bradham Lab at Boston University
        Correspondence: anthonygarza124@gmail.com''')

    parser.add_argument('-o', '--output', type=str, required=True, help='The directory to save the FASTA files. Recommended to be empty.')
    parser.add_argument('-db', '--database', type=str, required=True, help='The path to the LvEDGE IDs file')

    args = parser.parse_args()

    output_path = args.output
    db_path = args.database



    # Validate
    if not os.path.exists(output_path):
        raise FileNotFoundError(f'The output directory does not exist: {output_path}')

    if not os.path.isdir(output_path):
        raise NotADirectoryError(f'The output path is not a directory: {output_path}')

    if not os.path.exists(db_path):
        raise FileNotFoundError(f'The LvEDGE IDs file does not exist: {db_path}')

    # Getting list of genes that already have a FASTA file
    fasta_files = os.listdir(output_path)
    gene_ids = [file.split('.')[0] for file in fasta_files]

    # Scrape the LvEDGE database
    gsdb = GreenSeaUrchinDB(db_path, skip_list=gene_ids)

    for gene_id, record_list in gsdb:

        if len(record_list) == 0:
            print(gene_id)
            continue

        gene_id = gene_id.replace('/', '_')
        # Write the gene to a FASTA file
        output_file = os.path.join(output_path, f'{gene_id}.fa')
        with open(output_file, 'w') as f:
            SeqIO.write(record_list, f, 'fasta')
            