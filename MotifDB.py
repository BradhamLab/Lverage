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
from bs4 import BeautifulSoup

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class MotifRecord:
    '''
    This class represents a motif record
    '''

    def __init__(self, name, pfm, logo_link):
        self.name = name
        self.pfm = pfm
        self.logo_link = logo_link

    def get_name(self):
        return self.name
    
    def get_pfm(self):
        return self.pfm
    
    def get_logo_link(self):
        return self.logo_link


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

    search_url = "https://jaspar.elixir.no/search?q="
    rest_url = "https://jaspar2020.genereg.net/api/v1/"
    logo_url = "https://jaspar2020.genereg.net/static/logos/all/"

    def __init__(self):
        super().__init__()


    def search(self, description, tax_id):
        '''
        This function searches the JASPAR database for a given description and gets the corresponding motif

        Arguments:
        description: a string representing a description of a gene
        tax_id: a string representing the taxonomy ID of the species
        '''

        assert len(description) > 0, "Description must be non-empty"
        assert tax_id.isdigit(), "Taxonomy ID must be a number"


        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Get the matrix ID for the motif

        search = requests.get(self.search_url + description + "&tax_id=" + tax_id)

        # Search for matrix ID
        matrix_id = None
        if search.ok:
            table = BeautifulSoup(search.text, "html.parser").find("table")

            for row in table.find_all("tr"):
                a_tag = row.find('a', href=True)
                if a_tag:
                    matrix_id = a_tag['href'].split('/')[-1]
                    break

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Get motif name, pfm, and logo
                
        name = None
        pfm = None
        logo = None

        if matrix_id:
            motif_download = requests.get(self.rest_url + "matrix/" + matrix_id)

            if motif_download.ok:
                data = motif_download.json()
                name = data['name']
                pfm = data['pfm']
                logo = self.logo_url + matrix_id + ".png"

        return MotifRecord(name, pfm, logo) if name else None

class MotifDBFactory:
    '''
    This class is a factory that returns a motif database based on the database name
    '''

    # Database name -> MotifDB class
    motif_db_dict = {
        "JASPAR": JasparDB
    }

    def get_motif_db(self, db_name, **kwargs):
        '''
        This function returns a motif database based on the database name

        Arguments:
        db_name: a string representing the name of the database
        **kwargs: keyword arguments to pass to the motif database
        '''

        assert db_name in self.motif_db_dict, "Database name not found"

        return self.motif_db_dict[db_name](**kwargs)
    
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
    mdb = JasparDB()
    rec = mdb.search("ETS transcription factor ERG", "9606")
    print(rec.get_name())
    print(rec.get_pfm())
    print(rec.get_logo_link())



