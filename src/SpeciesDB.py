#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
This file contains the database interface for different species and a factory that returns a species database based on the species name.

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
from src.LvEdgeScraper import LvEdgeScraper

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
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

    def __init__(self, file_path):
        '''
        Arguments:
        file_path: path to the file containing the LvEdge IDs
        '''
        super().__init__()
        self.file_path = file_path

        # read id database in; first row is header, column 2 and 3 are spu id and gene id
        self.id_db = np.loadtxt(self.file_path, dtype=str, skiprows=1, delimiter=',')

        # Initialize LvEdgeScraper to read database online
        self.scraper = LvEdgeScraper()

    def __next__(self):

        if self._current_index < len(self.id_db):

            #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
            # Get the spu id and gene id
            _, spu_id, gene_id = self.id_db[self._current_index]
            spu_id, gene_id = spu_id.strip('"'), gene_id.strip('"')

            self._current_index += 1

            #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
            # Get the list of rows for the spu id
            for i in range(3):
                try:
                    matrix = self.scraper.search(spu_id)
                    break
                except Exception as e:
                    raise(e)
            
            #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
            # Get the first gene hit
            row = matrix[0]

            sequence_list = row[3]

            #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
            # Get the largest sequence corresponding to the hit gene
            sequence_list = [str(x.seq) for x in sequence_list]

            return gene_id, sequence_list
        
        else:
            raise StopIteration



class HumanDB(SpeciesDBInterface):
    '''
    This class is a database for the human.
    '''

    name = "human"

    def __init__(self, **kwargs):
        super().__init__()
        




class SpeciesDBFactory:
    '''
    This class is a factory that returns a species database based on the species name.
    It also contains useful information on species names and taxon ids
    '''

    # Species taxon id -> SpeciesDB class
    species_db_dict = {
        "7654": GreenSeaUrchinDB,
        "9606": HumanDB
    }

    # Species can have different names; this maps all names to the taxon id of the species
    # Species name -> taxon id
    species_name_taxon_dict = {}
    
    # Taxon id or Alternative Species Names -> standard scientific name
    species_taxon_name_dict = {}

    # taxon id set
    taxon_id_set = set()

    # standard scientific name set
    standard_name_set = set()

    @staticmethod
    def load_species_names():
        with open("species.txt", "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    name_list = [name.strip() for name in line.split(',')]
                    taxon_id = name_list[0]      # Taxon id should always be first column.
                    standard_name = name_list[1] # standard scientific name should always be second column.

                    SpeciesDBFactory.taxon_id_set.add(taxon_id)
                    SpeciesDBFactory.standard_name_set.add(standard_name)

                    for name in name_list:
                        SpeciesDBFactory.species_name_taxon_dict[name] = taxon_id
                        SpeciesDBFactory.species_taxon_name_dict[name] = standard_name

    @staticmethod
    def get_species_db(species_type, **kwargs):

        # Check if species type is valid
        species_type = species_type.strip().lower()

        # check if species is listed in species.txt file
        in_species_name = species_type in SpeciesDBFactory.species_name_taxon_dict

        # check if species has a database if it is listed in species.txt file
        in_species_db = True if in_species_name and SpeciesDBFactory.species_name_taxon_dict[species_type] in SpeciesDBFactory.species_db_dict else False

        if not in_species_db:

            db_names = SpeciesDBFactory.get_db_names()
            if not in_species_name:

                raise ValueError(f"Species type {species_type} is unknown! Please check your spelling and refer to species.txt. Known species with database parsers are: {', '.join(db_names)}")
            

            raise ValueError(f"Species type {species_type} does not have an available database parser! The following species are supported: {', '.join(db_names)}")
        
        # Get taxon id of species and return the database mapped to that species
        taxon_id = SpeciesDBFactory.species_name_taxon_dict[species_type]
        return SpeciesDBFactory.species_db_dict[taxon_id](**kwargs)
    
    @staticmethod
    def get_species():
        ''' 
        Get only species with available database parsers
        '''
        return list(SpeciesDBFactory.species_db_dict.keys())
    
    @staticmethod
    def get_all_species():
        return list(SpeciesDBFactory.standard_name_set)
    
    @staticmethod
    def get_dbs():
        return list(SpeciesDBFactory.species_db_dict.values())
    
    @staticmethod
    def get_taxon_id(species):
        return SpeciesDBFactory.species_name_taxon_dict[species.lower()]
    
    @staticmethod
    def get_standard_name(species):
        return SpeciesDBFactory.species_taxon_name_dict[species.lower()]
    
    @staticmethod
    def has_species(species):
        return species.lower() in SpeciesDBFactory.species_taxon_name_dict
    
    @staticmethod
    def get_db_names():
        return [taxon_id + ": " + SpeciesDBFactory.species_db_dict[taxon_id].name for taxon_id in SpeciesDBFactory.species_db_dict]
        
    
SpeciesDBFactory.load_species_names()

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# FOR TESTING PURPOSES    
if __name__ == "__main__":
    a = GreenSeaUrchinDB('../Data/LvEdgeIDs.csv')

    for gene, sequence in a:
        print(gene, sequence)
        break