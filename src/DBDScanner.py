#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
This file contains a class that is used to scan a protein sequence for its dna binding domain (DBD) 

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
import time

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class DBD:
    '''
    DNA-Binding domain (DBD) record that contains the name, start (in originating sequence), and size of the DBD.
    '''

    def __init__(self, name, start, size):
        '''
        Arguments:
        name: name of the DBD
        start: start of the DBD in its originating sequence
        size: size of the DBD in its originating sequence
        '''

        assert len(name) > 0, 'Invalid name for DBD!'
        assert start >= 0, f'Invalid start {start}'
        assert size >= 1, f'Invalid size {size}'

        self.name = name
        self.start = start
        self.size = size

    def get_name(self):
        return self.name
    
    def get_start(self):
        return self.start
    
    def get_size(self):
        return self.size
    
    def get_end(self):
        return self.start + self.size
    
    def __str__(self):
        return f"DBD {self.name} {self.start} {self.size}"



class DBDScanner:
    '''
    Runs pfamscan online to locate the DNA-Binding Domains (DBD) of a protein sequence
    '''

    # URLs for the pfamscan service
    url = "https://www.ebi.ac.uk/Tools/services/rest/pfamscan/"
    run_url = url + "run/"
    status_url = url + "status/"
    result_url = url + "result/"   

    def __init__(self, email, verbose = False, try_count = 10, time_interval = 5):
        '''
        Arguments:
        email: User email. This is required by ebi services
        verbose: Should it print out status and results?
        try_count: how many times should we request the online tool before giving up?
        time_interval: how many seconds should we wait between requests?
        '''

        assert len(email) > 0, 'Invalid email!'
        assert try_count > 0, f'Invalid try_count of {try_count}!'
        assert time_interval >= 0, f'Invalid time_interval of {time_interval}!'

        if time_interval < 1:
            print("\WARNING: time_interval is less than 1 second. This may cause issues with the online tool.")

        self.email = email
        self.verbose = verbose
        self.try_count = try_count
        self.time_interval = time_interval

    def find_dbds(self, protein_seq):
        '''
        This function finds the DNA binding domains of a protein sequence.

        Arguments:
        protein_sequence: a string representing a protein sequence

        Returns:
        dbd objects in a list
        '''


        r = [] # list of DBD objects to return in the end

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step one: run pfamscan
        data = {
            'email' : self.email,
            'sequence': protein_seq
        }

        if self.verbose:
            print("\tRunning pfamscan.", flush=True)

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
                print("\tpfamscan Job ID:", job_id, flush=True)

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

                    dbd_list = result.json()

                    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
                    # Step four: parse the results, create the DBD objects with the name, start, and size for each DBD.
                    if len(dbd_list) > 0:
                        for dbd in dbd_list:
                            location_dict = dbd['env']
                            start = int(location_dict['from']) - 1
                            size = int(location_dict['to']) - start
                            name = dbd['name']

                            if self.verbose:
                                print(f"\tDBD with name {name} found at {start} with size {size}.", flush=True)

                            r.append(DBD(name, start, size))

        else:
            print(f"\tPfamscan failed to run. Perhaps check the status of the tool online?")

                    
        return r

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# FOR TESTING PURPOSES
if __name__ == "__main__":

    test_sequence = "ORF35MLPKGGSCNICIGRVLLSDCLLSKATSCTIGDYLFMILGVLTDRTSSNKVEALSVVSDDQTLFENTYRDVNKSNTIVSSPQSSNHILTTPVHKPKLSVVASPSNESHIWSSCVQDQRMKQEIEQPTAGAHDGLPPPQVSRGEPESPLDCSVAKPRQQPPPHGAPPGAINAEPTQAAPPYPSCPPQPYLGNTGSTPEERSASAGTTGRSPPPNVTTNEKRVIVPADPNMWTAEHVQQWVQWAVREYSLQDVQVARFNMDGKHLCKMTKDDFSRLTNIVNVDVLISHLNFLKQTPLPNLTSDDIDKALQPSPRNPPSSQAGYTPGAKSNLDDKYNAADHQTTLDTNPSGGSAFPYPTSTTTTVDSVHRLPRTAPSCDSLVRGRPNAWPTTVPSAVSKGYTQTSPTLPKPGIDSSATQIRPGGRPAYGGSEFDPYQVFGHTSRTLANPVIPSDWQNLQSVILALRHNKGSGQIQLWQFLLELLSDSSNANCITWEGTNGEFKMTDPDEVARRWGERKSKPNMNYDKLSRALRYYYDKNIMTKVHGKRYAYKFDFAGLAQAMQPVQADPSMYRYQSDLTYLPGYHPTKLNFVGTPINPSTNASLFSSHSSYWSSPTGANIYPSGHVTHPHASHMSSHIGTYYG"

    ds = DBDScanner("antiny124@gmail.com", verbose=True)
    print(ds.find_dbd(test_sequence))

