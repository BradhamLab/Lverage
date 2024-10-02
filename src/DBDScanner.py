"""This module contains functionality for scanning protein sequences for their dna binding domain (DBD) 

Copyright (C) <RELEASE_YEAR_HERE> Bradham Lab

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

Correspondence: 
    - Cynthia A. Bradham - cbradham@bu.edu - *
    - Anthony B. Garza   - abgarza@bu.edu  - **
    - Stephanie P. Hao   - sphao@bu.edu    - **
    - Yeting Li          - yetingli@bu.edu - **
    - Nofal Ouardaoui    - naouarda@bu.edu - **

    \* Principle Investigator, ** Software Developers
"""

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Imports
import requests
import time

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class DBD:
    """DNA-Binding domain (DBD) record that contains the name, start (in originating sequence), end of the DBD, and accession identifier.

    The purpose of this class is to store the information of a DBD in a protein sequence.

    Attributes
    ----------
    name : str
        name of the DBD
    start : int
        start of the DBD in its originating sequence
    end : int
        end of the DBD in its originating sequence
    accession : str
        accession identifier of the DBD

    """

    def __init__(self, name : str, start : int, end : int, accession: str):
        """Constructor

        Parameters
        ----------
        name : str
            Name of the DBD.
        start : int
            Start position of the DBD in its originating sequence.
        end : int
            End position of the DBD in its originating sequence.
        accession : str
            Accession identifier of the DBD.

        Raises
        ------
        ValueError
            1. If the name is empty, start is negative, end is less than 1, or accession is empty.
            2. If the start is negative.
            3. If the end is less than 1.
            4. If the accession is empty.
        """

        if len(name) == 0:
            raise ValueError('Invalid name for DBD!')
        if start < 0:
            raise ValueError(f'Invalid start {start}')
        if end < 1:
            raise ValueError(f'Invalid end {end}')
        if len(accession) == 0:
            raise ValueError('Invalid accession for DBD!')

        self.name = name
        self.start = start
        self.end = end
        self.accession = accession


    def get_name(self):
        """
        Returns the name of the DBD.

        Returns
        -------
        str
            Name of the DBD.
        """
        return self.name

    def get_start(self):
        """
        Returns the start position of the DBD in its originating sequence.

        Returns
        -------
        int
            Start position of the DBD.
        """
        return self.start

    def get_end(self):
        """
        Returns the end position of the DBD in its originating sequence.

        Returns
        -------
        int
            End position of the DBD.
        """
        return self.end

    def get_size(self):
        """
        Returns the size of the DBD.

        Returns
        -------
        int
            Size of the DBD.
        """
        return self.end - self.start

    def get_accession(self):
        """
        Returns the accession identifier of the DBD.

        Returns
        -------
        str
            Accession identifier of the DBD.
        """
        return self.accession
    
    def __str__(self):
        return f"DBD {self.name} {self.start} {self.size} {self.accession}"
    

    
def build_dbd_from_dict(d: dict) -> DBD:
    """
    Builds a DBD object from a dictionary returned by the pfamscan service.

    Parameters
    ----------
    d : dict
        Dictionary containing the DBD information with keys 'name', 'type', 'env', and 'accession'.

    Returns
    -------
    DBD
        A DBD object initialized with the dictionary values.

    Raises
    ------
    ValueError
        1. If the dictionary is missing required fields: 'name', 'type', 'env', or 'accession'.
        2. If the 'env' sub-dictionary is missing 'from' or 'to' keys.
        3. If the 'from' or 'to' fields from the 'env' sub-dictionary are not integers.
    """

    # Validate the required keys exist and their types are correct
    required_keys = ['name', 'type', 'env', 'acc']
    for key in required_keys:
        if key not in d:
            raise ValueError(f"Invalid dictionary, missing '{key}' key!")

    # Validate the structure of the 'env' sub-dictionary
    if 'from' not in d['env'] or 'to' not in d['env']:
        raise ValueError("Invalid dictionary, 'env' must contain 'from' and 'to' keys!")

    # Extract values from the dictionary and cast appropriately
    name = d['name']
    accession = d['acc']

    try:
        start = int(d['env']['from']) - 1  # Convert to 0-based index
        end = int(d['env']['to'])
    except (ValueError, TypeError):
        raise ValueError("The 'from' and 'to' fields in 'env' must be integers.")

    # Create and return the DBD object
    return DBD(name, start, end, accession)



class DBDScanner:
    """
    DBDScanner is a class that runs pfamscan online to locate the DNA-Binding Domains (DBD) of a protein sequence.

    Attributes
    ----------
    email : str
        User email required by ebi services.
    verbose : bool
        Flag to print out status and results. Default is False.
    try_count : int
        Number of times to request the online tool before giving up. Default is 10.
    time_interval : int
        Number of seconds to wait between requests. Default is 5.
    url : str
        Base URL for the pfamscan service.
    run_url : str
        URL to submit a pfamscan job.
    status_url : str
        URL to check the status of a pfamscan job.
    result_url : str
        URL to retrieve the results of a pfamscan job.

    Raises
    ------
    ValueError
        1. If the email is empty
        2. If try_count is less than or equal to 0
        3. If time_interval is less than 0
    """


    url = "https://www.ebi.ac.uk/Tools/services/rest/pfamscan/"
    run_url = url + "run/"
    status_url = url + "status/"
    result_url = url + "result/"


    def __init__(self, email, verbose = False, try_count = 100, time_interval = 5):
        """Constructor
        
        Parameters
        ----------
        email : str
            User email required by ebi services.
        verbose : bool
            Flag to print out status and results. Default is False.
        try_count : int
            Number of times to request the online tool before giving up. Default is 10.
        time_interval : int
            Number of seconds to wait between requests. Default is 5.
            
        Raises
        ------
        ValueError
            1. If the email is empty
            2. If try_count is less than or equal to 0
            3. If time_interval is less than 0
        """

        if not email:
            raise ValueError('Invalid email!')
        if try_count <= 0:
            raise ValueError(f'Invalid try_count of {try_count}!')
        if time_interval < 0:
            raise ValueError(f'Invalid time_interval of {time_interval}!')

        if time_interval < 1:
            print("\tWARNING: time_interval is less than 1 second. This may cause issues with the online tool.")

        self.email = email
        self.verbose = verbose
        self.try_count = try_count
        self.time_interval = time_interval

        # self.url = "https://www.ebi.ac.uk/Tools/services/rest/pfamscan/"
        # self.run_url = self.url + "run/"
        # self.status_url = self.url + "status/"
        # self.result_url = self.url + "result/"

    def find_dbds(self, protein_seq):
        """
        Finds the DNA binding domains of a protein sequence.

        Interfaces with the pfamscan online tool to locate the DBDs in a provided protein sequence.
        If the tool fails to run, retrieve the results, or find any DBDs, an empty list is returned.

        Parameters
        ----------
        protein_seq : str
            Protein sequence to scan for DNA binding domains.
        
        Returns
        -------
        list
            List of DBD objects found in the protein sequence.
        """


        r = [] # list of DBD objects to return in the end

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step one: run pfamscan
        data = {
            'email' : self.email,
            'sequence': protein_seq
        }

        if self.verbose:
            print("Running pfamscan.", flush=True)

        run = requests.post(self.run_url, data=data)
        for i in range(self.try_count):
            time.sleep(self.time_interval)
            if run.ok:
                break
            else:
                run = requests.post(self.run_url, data=data)
        
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step two: get the job id and wait for it to finish
        if run.ok:

            job_id = run.text

            if self.verbose:
                print(f"\tpfamscan Job ID: {job_id}", flush=True)

            status = requests.get(self.status_url + job_id)
            if status.ok:

                while status.text != "FINISHED":
                    time.sleep(self.time_interval)
                    status = requests.get(self.status_url + job_id)    

                if self.verbose:
                    print("\tpfamscan finished, parsing results.", flush=True)

                #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
                # Step three: get the results
                result = requests.get(self.result_url + job_id + "/out")
                if result.ok:

                    dbd_dict_list = result.json()

                    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
                    # Step four: parse the results, create the DBD objects with the name, start, and size for each DBD.
                    if len(dbd_dict_list) > 0:
                        for dbd_dict in dbd_dict_list:

                            dbd = build_dbd_from_dict(dbd_dict)
                            r.append(dbd)

                            if self.verbose:
                                print(f"\t\tDBD with name {dbd.get_name()} and accession {dbd.get_accession()} was found at {dbd.get_start()} with size {dbd.get_size()}.", flush=True)


        else:
            print(f"Pfamscan failed to run. Perhaps check the status of the tool online?")

                    
        return r

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# FOR TESTING PURPOSES
if __name__ == "__main__":

    test_sequence = "ORF35MLPKGGSCNICIGRVLLSDCLLSKATSCTIGDYLFMILGVLTDRTSSNKVEALSVVSDDQTLFENTYRDVNKSNTIVSSPQSSNHILTTPVHKPKLSVVASPSNESHIWSSCVQDQRMKQEIEQPTAGAHDGLPPPQVSRGEPESPLDCSVAKPRQQPPPHGAPPGAINAEPTQAAPPYPSCPPQPYLGNTGSTPEERSASAGTTGRSPPPNVTTNEKRVIVPADPNMWTAEHVQQWVQWAVREYSLQDVQVARFNMDGKHLCKMTKDDFSRLTNIVNVDVLISHLNFLKQTPLPNLTSDDIDKALQPSPRNPPSSQAGYTPGAKSNLDDKYNAADHQTTLDTNPSGGSAFPYPTSTTTTVDSVHRLPRTAPSCDSLVRGRPNAWPTTVPSAVSKGYTQTSPTLPKPGIDSSATQIRPGGRPAYGGSEFDPYQVFGHTSRTLANPVIPSDWQNLQSVILALRHNKGSGQIQLWQFLLELLSDSSNANCITWEGTNGEFKMTDPDEVARRWGERKSKPNMNYDKLSRALRYYYDKNIMTKVHGKRYAYKFDFAGLAQAMQPVQADPSMYRYQSDLTYLPGYHPTKLNFVGTPINPSTNASLFSSHSSYWSSPTGANIYPSGHVTHPHASHMSSHIGTYYG"

    ds = DBDScanner("antiny124@gmail.com", verbose=True)
    print(ds.find_dbd(test_sequence))

