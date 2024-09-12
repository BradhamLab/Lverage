"""This file contains a class that is used to search for orthologous sequences of a given sequence using blastp

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
"""

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Imports
from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO
import io

import subprocess
import os

import requests
import time

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

class AlignmentRecord:
    """
    BLAST alignment record, contains useful information about a BLAST hit

    Attributes
    ----------
    title : str
        Title of the alignment
    species : str
        Name of the species
    name : str
        Name of the sequence
    seq : str
        Sequence of the alignment
    accession : str
        Accession number of the sequence
    expect : float
        E-score of the alignment
    percent_identity : float
        Percent identity of the alignment
    query_coverage : float
        Query coverage of the alignment

    Methods
    -------
    get_title()
        Get the title of the alignment
    get_name()
        Get the name of the sequence
    get_species()
        Get the name of the species
    get_seq()
        Get the sequence of the alignment
    get_accession()
        Get the accession number of the sequence
    get_escore()
        Get the E-score of the alignment
    get_percent_identity()
        Get the percent identity of the alignment
    get_query_coverage()
        Get the query coverage of the alignment
    """

    # List of words that are uneeded in the name of the BLAST description; we can remove them
    uneeded_words_list = ['transcription factor']


    def __init__(self, alignment, query_coverage, seq):        
        """Constructor
        """
        


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
        self.seq = seq

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
        """Get the title of the alignment

        Returns
        -------
        str
            Title of the alignment
        """
        return self.title
    
    def get_name(self):
        """Get the name of the sequence

        Returns
        -------
        str
            Name of the sequence
        """
        return self.name
    
    def get_species(self):
        """Get the name of the species

        Returns
        -------
        str
            Name of the species
        """
        return self.species
    
    def get_seq(self):
        """Get the sequence of the alignment

        Returns
        -------
        str
            Sequence of the alignment
        """
        return self.seq

    def get_accession(self):
        """Get the accession number of the sequence

        Returns
        -------
        str
            Accession number of the sequence
        """
        return self.accession
    
    def get_escore(self):
        """Get the E-score of the alignment

        Returns
        -------
        float
            E-score of the alignment
        """
        return self.expect
    
    def get_percent_identity(self):
        """Get the percent identity of the alignment

        Returns
        -------
        float
            Percent identity of the alignment
        """
        return self.percent_identity
    
    def get_query_coverage(self):
        """Get the query coverage of the alignment

        Returns
        -------
        float
            Query coverage of the alignment
        """
        return self.query_coverage
    

class OrthologSearcher:
    """Class to search for orthologous sequences of a given sequence using blastp

    Attributes
    ----------
    excluded_terms : list
        Words that we don't want in the name of the sequence. Discards an alignment from BLAST if it contains any of these words
    url : str
        URL for the BLAST API search
    species_list : list
        List of orthologous species to look through. Default contains human, frog, fly, and mouse
    hitlist_size : int
        How many hits should we look through from BLAST. Default is 100
    email : str
        Email to use for Entrez.
    escore_threshold : float
        E-score threshold for BLAST. Default is 10^-6
    verbose : bool
        Should we print out the current status and results? Default is False
    blastp_path : str
        Path to the blastp executable if being run locally
    database_path : str
        If BLAST is being ran locally, where is the database to query?
    """
    excluded_terms = ['hypothetical', 'unnamed', 'uncharacterized', 'unknown', 'partial', 'isoform'] # Words that we don't want in the name of the sequence; discard the alignment if it contains any of these words
    url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    def __init__(self, species_list=['homo sapiens', 'xenopus laevis', 'drosophila melanogaster', 'mus musculus'], 
                 hitlist_size = 100, email = None, escore_threshold = 10**-6, verbose = False,
                 blastp_path = None, database_path = None):
        """Constructor

        Raises
        ------
        TypeError
            The parameter 'hitlist_size' must receive an integer
        ValueError
            1. The parameter 'hitlist_size' must receive a value greater than 0
            2. Email must be provided if running remotely OR both the BLAST database and blastp executable paths must be provided
            3. E-score threshold must be greater than 0
        """

        if not isinstance(hitlist_size, int):
            raise TypeError("hitlist_size must be an integer!")
        if hitlist_size <= 0:
            raise ValueError("hitlist_size must be greater than 0")
        if not ((database_path and blastp_path) or email):
            raise ValueError("Email must be provided if running remotely OR both the BLAST database and blastp executable paths must be provided!")
        if escore_threshold <= 0:
            raise ValueError("E-score threshold must be greater than 0!")


        self.species_list = [x.lower() for x in species_list]
        self.entrez_query = f'({" OR ".join([species.lower() + "[ORGN]" for species in species_list])})'
        self.hitlist_size = hitlist_size
        self.email = email
        self.escore_threshold = escore_threshold
        self.verbose = verbose

        self.blastp_path = blastp_path
        self.database_path = database_path
        self.is_local = database_path and blastp_path

        Entrez.email = self.email

    def search_local(self, sequence):
        """This function searches for orthologous sequences of a given sequence using blastp locally.

        Parameters
        ----------
        sequence : str
            A string representing a protein sequence

        Returns
        -------
        blast_record : Bio.Blast.Record object

        Raises
        ------
        subprocess.CalledProcessError
            If the subprocess call to blastp fails
        """

        #@#@#@@
        # Step one: run blastp
        if self.verbose:
            print("\tRunning local blastp.", flush=True)

        # Write sequence to a file
        query_file = "query.fasta"
        with open(query_file, "w") as f:
            f.write(">query\n")
            f.write(sequence)

        # Run blastp
        output_file = "blastp_output.xml"
        blastp_cmd = [self.blastp_path, 
                  "-query", query_file, 
                  "-db", self.database_path,
                  "-out", output_file,
                  "-outfmt", "5", # XML format
                  "-evalue", str(self.escore_threshold), 
                  "-max_target_seqs", str(self.hitlist_size)]

        try:
            subprocess.run(blastp_cmd, check=True)
        
        except subprocess.CalledProcessError as e:
            print(f"Unable to run blastp: {e}")

        finally:
            # Delete query_file and output_file
            if os.path.exists(query_file):
                os.remove(query_file)
            if os.path.exists(output_file):
                os.remove(output_file)

        #@#@#@@
        # Step two: read blastp output
        if self.verbose:
            print("\tParsing local blastp output", flush=True)

        with open(output_file, "r") as f:
            blast_record = NCBIXML.read(f)

        # We have to adjust for the accession. If a custom database is used, the accession gets replaced with the local sequence ID.
        for alignment in blast_record.alignments:
            alignment.accession = alignment.hit_def.split()[0]

        return blast_record

    def search_remote(self, sequence):
        """This function searches for orthologous sequences of a given sequence using blastp remotely.

        Parameters
        ----------
        sequence : str
            A string representing a protein sequence

        Returns
        -------
        blast_record : Bio.Blast.Record object

        Raises
        ------
        requests.exceptions.HTTPError
            1. If the request to the NCBI server fails
            2. If the result of the BLAST job cannot be retrieved
        ValueError
            1. If the RID cannot be found in the response
            2. If the search fails
            3. If the status is unknown
        """
        
        #@#@#@@
        # Step one: run blastp
        if self.verbose:
            print("\tRunning NCBI blastp.", flush=True)

        # REST API parameters for blastp query
        params = {
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": "nr",
            "QUERY": sequence,
            "EXPECT": self.escore_threshold,
            "ENTREZ_QUERY": self.entrez_query,
            "EMAIL": self.email,
            "HITLIST_SIZE": self.hitlist_size,
            "FORMAT_TYPE": "XML",
        }

        # posting query to NCBI
        response = requests.post(self.url, data=params)
        response.raise_for_status()


        # rid is the request ID; lets us track the job and status
        try:
            rid = response.text.split("RID = ")[1].split("\n")[0]
        except IndexError:
            raise ValueError("Failed to retrieve RID from BLAST response")
        
        if self.verbose:
            print(f"\tUse the following link to check the status of the search: {self.url}?CMD=Get&RID={rid}", flush=True)

        # REST API parameters for checking status of job
        check_params = {
            "CMD": "Get",
            "FORMAT_TYPE": "HTML",
            "RID": rid
        }

        # Until the server is done processing the request, we keep checking the status
        while True:
            check_response = requests.post(self.url, data=check_params)

            if "Status=WAITING" in check_response.text:
                time.sleep(30)

            elif "Status=FAILED" in check_response.text:
                raise ValueError("Search failed")
            
            elif "Status=READY" in check_response.text:
                break

            else:
                raise ValueError("Unknown status")
            
        # REST API parameters for getting the results of the job
        result_params = {
            "CMD": "Get",
            "FORMAT_TYPE": "XML",
            "RID": rid
        }
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}

        result_handle = requests.get(self.url, params=result_params, headers=headers)
        result_handle.raise_for_status()
        
        #@#@#@@
        # Step two: read blastp output
        if self.verbose:
            print("\tParsing NCBI blastp output", flush=True)

        blast_record = NCBIXML.read(io.StringIO(result_handle.text))

        return blast_record


    def search(self, sequence):
        """This function searches for orthologous sequences of a given sequence using blastp.
        
        Parameters
        ----------
        sequence : str
            A string representing a protein sequence
            
        Returns
        -------
        hit_list : list
            List of AlignmentRecord objects representing the orthologous sequences found

        Raises
        ------
        ValueError
            Sequence must be non-empty
        """

        if len(sequence) == 0:
            raise ValueError("Sequence must be non-empty")

        hit_list = []

        if self.is_local:
            blast_record = self.search_local(sequence)
        
        else:
            blast_record = self.search_remote(sequence)


        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step three: filter blastp output

        if self.verbose:
            print("\tFiltering hits.")
            
        # default sorted by e-score
        for alignment in blast_record.alignments:

            # get the sequence of the ortholog
            try:
                handle = Entrez.efetch(db="protein", id=alignment.accession, rettype="fasta", retmode="text")
            except Exception as e:
                if self.verbose:
                    print(f"\tUnable to fetch sequence with accession {alignment.accession}: {e}")
                continue
            
            record = SeqIO.read(handle, "fasta")
            handle.close()
            time.sleep(1) # no more than 3 requests per second, so we need to sleep


            ortholog_sequence = str(record.seq)

            query_coverage = alignment.hsps[0].align_length / len(sequence)

            record = AlignmentRecord(alignment, query_coverage, ortholog_sequence)

            # If the title contains any of the excluded terms, we skip this alignment
            if not all(sub not in record.name for sub in self.excluded_terms):
                continue
        
            # Check if species is in species list
            is_species = False
            for species in self.species_list:
                if species in record.get_species().lower():
                    is_species = True
                    break

            if is_species:

                hit_list.append(record)

                if self.verbose:
                    print(f"\tObtained ortholog from {record.get_species()}: {record.get_name()}, with e-score {record.expect} and similarity {record.get_percent_identity()}", flush=True)

        return hit_list

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# FOR TESTING PURPOSES
if __name__ == "__main__":
    import sys

    sequence ="LQRTTRSXLALLLCMWPVTRTNDRVLVWFQNRRAKWRKRERFQQFQNMRGLGPGSGYEMPIAPRPDAYSQVNSPGFHVLGDTHQPPAPVEGAMLRICRNLQNLRREFDSRKIGCHPSSSVGGTPGTSTTNESQDTSNHSSMIHQSSPWATAANMASPLASSMSPVGQQPQMPGQNPINSCMAPQSTLPSFMGVPAHQMNNTGVMNPMSNMTSMTSMPTSMPPSSGTAPVSSPSSNFMSSVGGLNAAAYNGQYTDMHPTVEGVGGVDRRTNSIAALRLRAKEHSSVMGMMNGYS"

    ortho = OrthologSearcher(verbose=True, email = 'anthonygarza124@gmail.com', species_list=['homo sapiens', 'mus musculus'], hitlist_size=10)

    hits = ortho.search(sequence)
    
