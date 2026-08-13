"""
This file contains the Lverage class

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

    \* Principle Investigator, ** Software Developers
"""

#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# IMPORTS
from ete3 import NCBITaxa
from validate_email import validate_email
import os
import shutil
from enum import Enum

from MotifDBTemplate import MotifDBTemplate
from OrfSearcherTemplate import OrfSearcherTemplate
from DomainScannerTemplate import DomainScannerTemplate

class LverageCode(Enum):
    """
    Enum for Lverage result codes. These are used to check where and why motif searching failed.

    Attributes
    ----------
    NOT_SET : int
        Lverage has not been run yet
    SUCCESS : int
        Motif searching was successful
    NO_ORF : int
        No Open Reading Frame was found in the transcription factor sequence
    NO_VALID_DOMAIN : int
        No valid domain was found in the Open Reading Frame
    NO_ORTHOLOGS : int
        No orthologs were found with blastp
    NO_VALID_DBD : int
        No valid DNA-binding domain was found in the orthologs
    NO_MOTIF : int
        No motif was found in the orthologs
    """
    NOT_SET = -1
    SUCCESS = 0
    NO_ORF = 1
    NO_VALID_DOMAIN = 2
    NO_ORTHOLOGS = 3
    NO_VALID_DBD = 4
    NO_MOTIF = 5

class Lverage:
    """
    
    Attributes
    ----------
    motif_database_list : list
        MotifDB objects representing databases to search in with orthologous species. Each MotifDB must be a subclass of MotifDBTemplate
    orf_searcher : subclass of OrfSearcherTemplate
        Open-Reading-Frame searching object. Must be a subclass of OrfSearcherTemplate
    domain_scanner : subclass of DomainScannerTemplate
        Domain scanning object. Must be a subclass of DomainScannerTemplate
    ortholog_species_list : list
        Taxonomic IDs of orthologous species to use for searching in Motif Database
    valid_pfam_list : list
        List of valid PFAM domains (str) that are DNA-binding domains. By default or if none is passed, all PFAMs are considered valid.
    blastp_escore_thresh : float
        The minimum e-score threshold for blastp
    blastp_top_n : int
        The maximum number of top hits to consider for each query
    blastp_database_path : str
        Path to the blastp database
    blastp_path : str
        Path to the blastp executable
    dbd_identity_thresh : float
        The minimum identity threshold for DBD to ortholog DBD to be considered for motif searching
    email : str
        Email to use for any services requiring email to use, e.g., EMBL tools.
    verbosity : int
        How verbose should Lverage be? 0 is silent, 1 is minimal output (Status of current sequence), 2 is full output (every step)
    user_interrupt : bool
        Interrupt the program at specific stages to ask for user input. This automatically sets verbosity to 2 if True.
            1. When checking is species are available in the Motif Databases, ask user if they want to continue. Verbosity 2 will print which species are not available.
    
    Methods
    ----------
    __init__
        Constructor for Lverage class
    validate_arguments
        Validates class attributes that are passed to the constructor
    """

    lverage_code_info = {
        LverageCode.NOT_SET: "Lverage has not been run yet",
        LverageCode.SUCCESS: "Motif searching was successful",
        LverageCode.NO_ORF: "No Open Reading Frame was found in the transcription factor sequence",
        LverageCode.NO_VALID_DOMAIN: "No valid domain was found in the Open Reading Frame",
        LverageCode.NO_ORTHOLOGS: "No orthologs were found with blastp",
        LverageCode.NO_VALID_DBD: "No valid DNA-binding domain was found in the orthologs",
        LverageCode.NO_MOTIF: "No motif was found in the orthologs"
    }

    def __init__(self, 
            motif_database_list : list,
            orf_searcher : OrfSearcherTemplate,
            domain_scanner : DomainScannerTemplate,
            ortholog_species_list : list = [9606, 10090], 
            valid_pfam_list : list = [], 
            blastp_escore_thresh : float = 10**-6, 
            blastp_top_n : int = 20, 
            blastp_database_path : str = "",
            blastp_path : str = "",
            dbd_identity_thresh : float = 0.7,
            email : str = "", 
            verbosity : int = 0,
            user_interrupt : bool = False
            ):
        """Constructor"""

        # Assigning class attributes
        self.motif_database_list = motif_database_list
        self.orf_searcher = orf_searcher
        self.domain_scanner = domain_scanner
        self.ortholog_species_list = ortholog_species_list
        self.valid_pfam_list = valid_pfam_list
        self.blastp_escore_thresh = blastp_escore_thresh
        self.blastp_top_n = blastp_top_n
        self.blastp_database_path = blastp_database_path
        self.blastp_path = blastp_path
        self.dbd_identity_thresh = dbd_identity_thresh
        self.email = email
        self.verbosity = verbosity
        self.user_interrupt = user_interrupt

        self.is_blastp_local = True if self.blastp_database_path else False
        self.blastp_path = shutil.which("blastp") if not self.blastp_path and self.is_blastp_local else self.blastp_path

        # result code for checking where and why motif searching failed; by default, it's not set
        self.lverage_code = LverageCode.NOT_SET

        # Attributes set during run
        self.tf_sequences = None # list of sequences for a specific transcription factor
        self.orf = None # Open Reading Frame scanned from the transcription factor sequence
        self.valid_domains = None # list of valid domains (in valid_pfam_list) found in the ORF
        self.orthologs = None # list of orthologs found with blastp

        # Validating arguments
        self.validate_arguments()

        # modify verbosity if user wishes to interrupt so they can know whats happening
        if self.user_interrupt:
            self.verbosity = 2

        # Check if species are available in the motif databases; at least one is required
        self.motif_available_species = {mdb.name : [] for mdb in self.motif_database_list} # {motif database name: [present species]}
        is_present = False # has at least one species present in at least one motif database
        has_missing = False # has at least one species missing from at least one motif database
        for motif_database in self.motif_database_list:
            for ortholog_species in self.ortholog_species_list:
                if motif_database.check_species_validity(ortholog_species):
                    is_present = True
                    self.motif_available_species[motif_database.name].append(ortholog_species)
                else:
                    has_missing = True
                    if self.verbosity in [1, 2]:
                        print(f"Species {ortholog_species} is not available in {motif_database.name}")
        
        if has_missing and self.user_interrupt:
            if input("Currently there are species missing from the database. Do you want to continue with the remaining species? (y/n): ").lower() != "y":
                raise ValueError("Lverage constructor error! User chose to not continue with missing species")
            
        if not is_present:
            raise ValueError(f"Lverage constructor error! Argument 'ortholog_species_list' must contain at least one species present in at least one motif database!")

        self.ortholog_species_list = sorted(set(
            species for species_list in self.motif_available_species.values()
            for species in species_list
        ))

    def print_result(self) -> None:
        pass

    def validate_arguments(self) -> None:
        """Validates class attributes that are passed to the constructor"""

        # validating motif_database_list
        if not isinstance(self.motif_database_list, list):
            raise TypeError(f"Lverage constructor error! Argument 'motif_database_list' must be a list, currently {type(self.motif_database_list)}")
        if len(self.motif_database_list) == 0:
            raise ValueError("Lverage constructor error! Argument 'motif_database_list' must contain at least one motif database")
        for motif_database in self.motif_database_list:
            if motif_database.__class__ is MotifDBTemplate:
                raise TypeError("Lverage constructor error! Elements in 'motif_database_list' must not be instances of MotifDBTemplate, but subclasses of it")
            if not isinstance(motif_database, MotifDBTemplate):
                raise TypeError(
                    f"Lverage constructor error! All elements in 'motif_database_list' must be instances of MotifDBTemplate, "
                    f"currently {type(motif_database).__name__} (Base classes: {[cls.__name__ for cls in motif_database.__class__.__mro__[1:]]})"
                )


        # validating orf_searcher
        if self.orf_searcher.__class__ is OrfSearcherTemplate:
            raise TypeError("Lverage constructor error! Object 'orf_searcher' must not be an instance of OrfSearcherTemplate, but a subclass of it")
        if not isinstance(self.orf_searcher, OrfSearcherTemplate):
            raise TypeError(
                f"Lverage constructor error! Argument 'orf_searcher' must be an instance of OrfSearcherTemplate, "
                f"currently {type(self.orf_searcher).__name__} (Base classes: {[cls.__name__ for cls in self.orf_searcher.__class__.__mro__[1:]]})"
            )
        
        # validating domain_scanner
        if self.domain_scanner.__class__ is DomainScannerTemplate:
            raise TypeError("Lverage constructor error! Object 'domain_scanner' must not be an instance of DomainScannerTemplate, but a subclass of it")
        if not isinstance(self.domain_scanner, DomainScannerTemplate):
            raise TypeError(
                f"Lverage constructor error! Argument 'domain_scanner' must be an instance of DomainScannerTemplate, "
                f"currently {type(self.domain_scanner).__name__} (Base classes: {[cls.__name__ for cls in self.domain_scanner.__class__.__mro__[1:]]})"
            )


        # validating ortholog_species_list
        if not isinstance(self.ortholog_species_list, list):
            raise TypeError(f"Lverage constructor error! Argument 'ortholog_species_list' must be a list, currently {type(self.ortholog_species_list)}")
        if len(self.ortholog_species_list) == 0:
            raise ValueError("Lverage constructor error! Argument 'ortholog_species_list' must contain at least one orthologous species")
        for ortholog_species in self.ortholog_species_list:
            if not isinstance(ortholog_species, int):
                raise TypeError(f"Lverage constructor error! Argument 'ortholog_species_list' must contain integer taxonomic IDs, currently {ortholog_species} is {type(ortholog_species)}")
            ncbi = NCBITaxa()
            if ortholog_species not in ncbi.get_rank([ortholog_species]).keys():
                raise ValueError(f"Lverage constructor error! Argument 'ortholog_species_list' must contain valid taxonomic IDs, currently {ortholog_species} is invalid")
            
        is_present = False # if at least one species is present in at least one motif database
        for motif_database in self.motif_database_list:
            for ortholog_species in self.ortholog_species_list:
                if motif_database.check_species_validity(ortholog_species):
                    is_present = True
                    break
            if is_present:
                break
        if not is_present:
            raise ValueError(f"Lverage constructor error! Argument 'ortholog_species_list' must contain at least one species present in at least one motif database!")
            
        # validating valid_pfam_list
        if not isinstance(self.valid_pfam_list, list):
            raise TypeError(f"Lverage constructor error! Argument 'valid_pfam_list' must be a list, currently {type(self.valid_pfam_list)}")
        if len(self.valid_pfam_list) == 0:
            raise ValueError("Lverage constructor error! Argument 'valid_pfam_list' must contain at least one valid PFAM domain")
        for valid_pfam in self.valid_pfam_list:
            if not isinstance(valid_pfam, str):
                raise TypeError(f"Lverage constructor error! Argument 'valid_pfam_list' must contain strings, currently {valid_pfam} is {type(valid_pfam)}")
        
        # validating blastp_escore_thresh
        if not isinstance(self.blastp_escore_thresh, float):
            raise TypeError(f"Lverage constructor error! Argument 'blastp_escore_thresh' must be a float, currently {type(self.blastp_escore_thresh)}")
        if self.blastp_escore_thresh <= 0:
            raise ValueError(f"Lverage constructor error! Argument 'blastp_escore_thresh' must be greater than 0, currently {self.blastp_escore_thresh}")
        
        # validating blastp_top_n
        if not isinstance(self.blastp_top_n, int):
            raise TypeError(f"Lverage constructor error! Argument 'blastp_top_n' must be an integer, currently {type(self.blastp_top_n)}")
        if self.blastp_top_n <= 0:
            raise ValueError(f"Lverage constructor error! Argument 'blastp_top_n' must be greater than 0, currently {self.blastp_top_n}")
        
        # validating blastp_database_path
        if not isinstance(self.blastp_database_path, str):
            raise TypeError(f"Lverage constructor error! Argument 'blastp_database_path' must be a string, currently {type(self.blastp_database_path)}")
        if self.blastp_database_path:
            if not os.path.exists(self.blastp_database_path):
                raise FileNotFoundError(f"Lverage constructor error! Argument 'blastp_database_path' does not exist at {self.blastp_database_path}")

        # validating blastp_path
        if not isinstance(self.blastp_path, str):
            raise TypeError(f"Lverage constructor error! Argument 'blastp_path' must be a string, currently {type(self.blastp_path)}")
        if self.is_blastp_local:
            if not os.path.exists(self.blastp_path):
                raise FileNotFoundError(f"Lverage constructor error! Argument 'blastp_path' does not exist at {self.blastp_path}")
  
        # validating dbd_identity_thresh
        if not isinstance(self.dbd_identity_thresh, float):
            raise TypeError(f"Lverage constructor error! Argument 'dbd_identity_thresh' must be a float, currently {type(self.dbd_identity_thresh)}")
        if self.dbd_identity_thresh <= 0:
            raise ValueError(f"Lverage constructor error! Argument 'dbd_identity_thresh' must be greater than 0, currently {self.dbd_identity_thresh}")
        
        # validating email
        if not isinstance(self.email, str):
            raise TypeError(f"Lverage constructor error! Argument 'email' must be a string, currently {type(self.email)}")
        if not validate_email(self.email):
            raise ValueError(f"Lverage constructor error! Argument 'email' must be a valid email address, currently {self.email}")
        
        # validating verbosity
        if not isinstance(self.verbosity, int):
            raise TypeError(f"Lverage constructor error! Argument 'verbosity' must be an integer, currently {type(self.verbosity)}")
        if self.verbosity < 0 or self.verbosity > 2:
            raise ValueError(f"Lverage constructor error! Argument 'verbosity' must be 0, 1, or 2, currently {self.verbosity}")
        
        # validating user_interrupt
        if not isinstance(self.user_interrupt, bool):
            raise TypeError(f"Lverage constructor error! Argument 'user_interrupt' must be a boolean, currently {type(self.user_interrupt)}")
        
    def __search_orfs(self):
        """
        Search for Open Reading Frames (ORFs) in the transcription factor sequence.
        Sets class attribute self.orf to the first ORF found.
        Also sets class attribute self.valid_domains to the list of valid domains found in the ORF.
        """

        self.orf = None
        self.valid_domains = []

        # Getting the Open Reading Frames
        orf_list = []
        for sequence in self.tf_sequences:
            orf_list.extend(self.orf_searcher.get_orfs(sequence))     

        orf_list.sort(key=lambda x: len(x), reverse=True) # sort by length

        # Getting the first ORF that has a domain (in provided valid_pfam_list)
        for orf in orf_list:
            domains = self.domain_scanner.get_domains(orf)

            if not domains:
                continue

            if self.valid_pfam_list:
                for domain in domains:
                    if domain in self.valid_pfam_list:
                        self.orf = orf
                        self.valid_domains.append(domain)
            else:
                self.orf = orf
                self.valid_domains = domains

            if self.orf:
                break
        
    def __blastp(self):
        pass

    def __blastp_local(self):
        pass

    def __blastp_remote(self):
        pass
        
    def run(self, tf_sequence : str | list[str]) -> list[int]:
        """
        Run Lverage on a transcription factor sequence or a list of sequences.

        Parameters
        ----------
        tf_sequence : str | list[str]
            Sequence of a transcription factor to search for DNA-Binding Motifs with. If there are multiple sequences for this one, specific transcription factor (such as scaffolds), pass a list of sequences.
        """

        motif_records = []

        if not isinstance(tf_sequence, str) and not isinstance(tf_sequence, list):
            raise TypeError(f"Lverage run error! Argument 'tf_sequence' must be a string or a list of strings, currently {type(tf_sequence)}")
        if isinstance(tf_sequence, list):
            for seq in tf_sequence:
                if not isinstance(seq, str):
                    raise TypeError(f"Lverage run error! Argument 'tf_sequence' must be a string or a list of strings, currently {type(seq)} in the list")
                
        if isinstance(tf_sequence, str):
            tf_sequence = [tf_sequence]

        self.tf_sequences = tf_sequence

        # Searching for Open Reading Frames
        self.__search_orfs()

        if not self.orf:
            if self.verbosity in [1, 2]:
                print("No ORF found with a valid domain")
            return motif_records
        
        # Searching for orthologs with blastp
        self.orthologs = self.__blastp()

        # Get domains from orthologs
        self.ortholog_domains = {} # {ortholog: [domains]}
        for ortholog in self.orthologs:
            self.ortholog_domains[ortholog] = self.domain_scanner.get_domains(ortholog)

        





#@#@#@#@#@#@#@#@ FOR TESTING ONLY #@#@#@#@#@#@#@#@
if __name__ == "__main__":
    motif_database_list = [MotifDBTemplate()]
    orf_searcher = [OrfSearcherTemplate()]
    ortholog_species_list=[9606, 10090]
    valid_pfam_list=[]
    blastp_escore_thresh=10**-6
    blastp_top_n=20
    dbd_identity_thresh=0.7
    email="anthonygarza124@gmail.com"
    verbosity=2
    user_interrupt=False

    lverage = Lverage(
        motif_database_list,
        orf_searcher,
        ortholog_species_list,
        valid_pfam_list,
        blastp_escore_thresh,
        blastp_top_n,
        dbd_identity_thresh,
        email,
        verbosity,
        user_interrupt
    )
    print("Lverage object created successfully!")