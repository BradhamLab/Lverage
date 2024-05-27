#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
The purpose of this script is to predict DNA-binding domain motifs of a species using orthologous species.

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
# global variables
is_optimize = False


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Importing modules

import warnings
import os
import re
from Bio import BiopythonDeprecationWarning, SeqIO
from Bio.Align import PairwiseAligner
from dataclasses import dataclass
import time


from src.ProteinTranslator import ProteinTranslator
from src.DBDScanner import DBDScanner, DBD
from src.OrthologSearcher import OrthologSearcher
from src.MotifDB import MotifDBFactory
import src.utils as lvutils
import src.exceptions as lvexceptions

warnings.simplefilter('ignore', BiopythonDeprecationWarning)



#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

@dataclass
class DiagnosticRecord:
    '''
    This class contains the intermediate results of the Lverage pipeline

    Attributes:
    dbd_name: name of DNA-binding domain
    ortholog_description: BLAST description of ortholog species
    ortholog_species: BLAST's name for ortholog species
    ortholog_escore: BLAST E-value
    ortholog_percent_identity: BLAST percent identity
    ortholog_query_coverage: BLAST query coverage
    ortholog_dbd_percent_identity: percent identity of ortholog's DNA-binding domain with the main species' DNA-binding domain
    '''

    dbd_name: str
    ortholog_description: str
    ortholog_species: str
    ortholog_escore: float
    ortholog_percent_identity: float
    ortholog_query_coverage: float
    ortholog_dbd_percent_identity: float

    header_list = ["Gene DBD Name, Ortholog BLAST Description, Ortholog BLAST Species, Ortholog BLAST E-value, Ortholog BLAST Percent Identity, Ortholog BLAST Query Coverage, Gene_DBD-Ortholog_DBD Percent Identity"]

    def get_values(self):
        return [self.dbd_name, self.ortholog_description, self.ortholog_species, self.ortholog_escore, self.ortholog_percent_identity, self.ortholog_query_coverage, self.ortholog_dbd_percent_identity]

class Lverage:
    '''
    This class contains a pipeline that predicts DNA-binding domain motifs of a species using orthologous species.
    '''

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Class variables


    def __init__(self, motif_database : str, ortholog_name_list : list, email : str, identity_threshold : float, verbose : bool,
                        blast_hit_count : int = 20, motif_hit_count : int = 10, escore_threshold : float = 10**-6):
        '''
        Arguments:
        motif_database: name of motif database to use
        ortholog_name_list: list of ortholog species names to search for motif in
        email: email address for EBML tools
        identity_threshold: identity threshold for similarity between sequences (0 to 1)
        verbose: if True, prints messages to the console
        blast_hit_count: number of hits to use from BlastP
        motif_hit_count: number of hits to use from motif database
        escore_threshold: E-value threshold required for any alignment to be considered
        '''

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@# 
        # Validating arguments

        # Checking if motif database is available
        if not MotifDBFactory.has_db(motif_database):
            raise lvexceptions.MotifDatabaseError()
        
        # Checking if ortholog species are in NCBI database; if tax ID, converting to scientific name
        for i, ortholog_species in enumerate(ortholog_name_list):
            if not lvutils.check_valid_species(ortholog_species):
                raise lvexceptions.OrthologSpeciesError()

        
        # Checking if email is valid
        pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
        if not re.match(pattern, email):
            raise lvexceptions.EmailError()
        
        # Checking if identity threshold is between 0 and 1
        if not 0 <= identity_threshold <= 1:
            raise lvexceptions.IdentityThresholdError()
        
        # Checking if number of blast hits to search for is greater than 0
        if blast_hit_count <= 0:
            raise lvexceptions.BlastHitCountError()
        
        # Checking if number of motif hits to search for is greater than 0
        if motif_hit_count <= 0:
            raise lvexceptions.MotifHitCountError()
        
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Assigning arguments to class variables

        self.motif_database = motif_database
        self.ortholog_name_list = ortholog_name_list
        self.email = email
        self.identity_threshold = identity_threshold
        self.verbose = verbose
        self.blast_hit_count = blast_hit_count
        self.motif_hit_count = motif_hit_count


        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Setting up tools

        # Setting up ORFfinder
        self.pf = ProteinTranslator()

        # Setting up pfamscan
        self.ds = DBDScanner(email = self.email, verbose = self.verbose)

        # Setting up aligner
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1

        # Setting up BlastP
        self.ors = OrthologSearcher(hitlist_size = self.blast_hit_count, verbose = self.verbose, species_list = self.ortholog_name_list, email = self.email)

        # Setting up motif database
        self.mdb = MotifDBFactory.get_motif_db(db_name = self.motif_database, n_hits = self.motif_hit_count, dbd_scanner = self.ds)

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # setting up call variables

        # For each call, this list will contain a DiagnosticRecord object for each hit. This list is overwritten on each call
        self.diagnostic_list = []


    def find_protein(self, gene_sequence_list : list):
        '''
        Step 1 in the pipeline:
        Get the protein sequence of the DNA sequence (or scaffolds of a gene)

        Arguments:
        gene_sequence_list: list of gene DNA sequences

        Returns:
        protein: the protein sequence of the DNA sequence
        '''

        if self.verbose:
            print(f">Obtaining protein sequence.", flush=True)

        if is_optimize:
            start_time = time.time()

        # translate DNA to protein
        protein_sequence = ""

        # Get longest protein sequence
        for gene_sequence in gene_sequence_list:
            ps = self.pf.translate(gene_sequence)
            protein_sequence = ps if len(ps) > len(protein_sequence) else protein_sequence

        if is_optimize:
                print(f"Time to find protein sequence: {time.time() - start_time}", flush=True)    

        # Only print protein sequence if found
        if protein_sequence and self.verbose:
            print("\t" + protein_sequence, flush=True)

        elif self.verbose:
            print(f"\tWARNING: Unable to find protein sequence.", flush=True)


        return protein_sequence
    
    
    def find_orthologs(self, protein_sequence : str):
        '''
        Step 2 in the pipeline:
        Get ortholog sequences
        
        Arguments:
        protein_sequence: protein sequence
        
        Returns:
        List of Ortholog objects
        '''

        if self.verbose:
            print(f">Searching for orthologs.", flush=True)

        if is_optimize:
            start_time = time.time()

        # search for orthologs
        ortholog_hit_list = self.ors.search(protein_sequence)

        if is_optimize:
            print(f"Time to find orthologs: {time.time() - start_time}", flush=True)

        if not ortholog_hit_list and self.verbose:
            print(f"\tWARNING: Unable to find orthologs.", flush=True)

        return ortholog_hit_list
    

    def find_dbds(self, protein_sequence_list : list):
        '''
        Step 3 in the pipeline:
        Get the start and size of DNA-binding domain within protein sequences
        
        Arguments:
        protein_sequence_list: list of protein sequences
        
        Returns:
        List of list of DBD objects (an inner list for each protein sequence)
        '''

        if self.verbose:
            print(f">Searching for DBD.", flush=True)

        if is_optimize:
            start_time = time.time()
        
        # search for DBDs
        dbd_list = []
        for protein_sequence in protein_sequence_list:
            dbd_list.append(self.ds.find_dbds(protein_sequence))
            time.sleep(5) # to prevent overloading the server

        if is_optimize:
            print(f"Time to find DBD: {time.time() - start_time}", flush=True)

        if not dbd_list and self.verbose:
            print(f"\tWARNING: Unable to find DBDs.", flush=True)

        return dbd_list
    
    def find_alignments(self, protein_sequence : str, main_dbd : DBD, ortholog_hit_list : list, ortholog_dbd_list : list):
        '''
        Step 4a in the pipeline:
        Align main species' protein dbd with orthologs' protein dbd
        
        Arguments:
        protein_sequence: protein sequence
        main_dbd: DBD object of the main species' protein
        ortholog_hit_list: list of Ortholog objects
        ortholog_dbd_list: list of list of DBD objects (an inner list for each ortholog sequence)
        
        
        Returns:
        List of alignments (may contain None if alignment is not found)
        '''

        assert len(ortholog_hit_list) == len(ortholog_dbd_list), "Ortholog hit list and ortholog DBD list are not the same length"

        if self.verbose:
            print(f"\tAligning sequences...", flush=True)

        if is_optimize:
            start_time = time.time()

        # for each ortholog sequence
        alignment_list = []
        main_dbd_seq = protein_sequence[main_dbd.get_start():main_dbd.get_end()] # Getting just the main DBD sequence
        for ortholog_hit, ortholog_dbd_list in zip(ortholog_hit_list, ortholog_dbd_list):
            
            # get the correct DBD
            ortholog_dbd = None
            for dbd in ortholog_dbd_list:
                if main_dbd.get_name() == dbd.get_name():
                    ortholog_dbd = dbd
                    break

            # If corresponding DBD is not found in ortholog, skip
            if ortholog_dbd is None:
                alignment_list.append(None)
                if self.verbose:
                    print(f"\t\tWARNING: Unable to find matching DBD for {ortholog_hit.get_name()} for main DBD {main_dbd.get_name()}.", flush=True)
                continue

            # Getting just the ortholog DBD sequence
            ortholog_dbd_seq = ortholog_hit.get_seq()[ortholog_dbd.get_start():ortholog_dbd.get_end()]

            # Aligning DBDs
            alignments = self.aligner.align(main_dbd_seq, ortholog_dbd_seq)

            # getting best alignment
            alignment_list.append(alignments[0])

            
        if is_optimize:
            print(f"Time to align sequences: {time.time() - start_time}", flush=True)

        return alignment_list
    
    def find_motifs(self, dbd_list : list, ortholog_hit_list : list, protein_sequence : str):
        '''
        Step 4 in the pipeline:
        Find motifs from orthologous sequences that meet the identity threshold with the DBD in the provided
        
        Arguments:
        dbd_list: list of list of DBD objects
        ortholog_hit_list: list of Ortholog objects
        protein_sequence: protein sequence of the gene

        Returns:
        List of Motif objects
        '''

        assert len(dbd_list) - 1 == len(ortholog_hit_list), "DBD list should have one more element than ortholog list (including the main gene of interest)"

        # return value
        result_motif_list = []

        if is_optimize:
            start_time = time.time

        main_dbd_list = dbd_list[0]
        ortholog_dbd_list = dbd_list[1:]

        for dbd in main_dbd_list:

            if self.verbose:
                print(f">Finding conserved motifs of {gene_id} with dbd {dbd.get_name()}.", flush=True)

            # Step 4a: aligning DBD for similarity threshold check
            alignment_list = self.find_alignments(protein_sequence, dbd, ortholog_hit_list, ortholog_dbd_list)
            if not alignment_list:
                continue

            sim_array = []
            for alignment in alignment_list:
                if alignment:
                    sim_array.append(lvutils.calculate_alignment_similarity(*alignment))
                else:
                    sim_array.append(None)

            # for each ortholog sequence
            for i in range(len(sim_array)):

                ortholog_hit = ortholog_hit_list[i]
                ortholog_description = ortholog_hit.get_name()

                # If the DBD similarity is greater than or equal to the identity threshold, get the motif
                dbd_similarity = sim_array[i]

                if dbd_similarity is None:
                    if self.verbose:
                        print(f'\tWARNING: No similarity for {ortholog_description}! Skipping...', flush=True)
                    continue

                if dbd_similarity >= self.identity_threshold:

                    ortholog_species = ortholog_hit.get_species()
                    ortholog_taxon_id = str(lvutils.get_tax_id(ortholog_species))

                    if ortholog_taxon_id is None:
                        if self.verbose:
                            print(f'\tWARNING: Unable to find taxon ID for {ortholog_species}! Skipping...', flush=True)
                        continue

                    ortholog_escore = ortholog_hit.get_escore()
                    ortholog_percent_identity = ortholog_hit.get_percent_identity()
                    ortholog_query_coverage = ortholog_hit.get_query_coverage()
                    ortholog_sequence = ortholog_hit.get_seq()

                    if self.verbose:
                        print(f"\n\t-Looking for ({ortholog_description}) from {ortholog_hit.get_species()}", flush=True)
                        print("\t\t", end='', flush=True)

                    found_motif = False

                    try:
                        for motif in self.mdb.search(ortholog_sequence, ortholog_taxon_id, protein_sequence, dbd):
                            result_motif_list.append(motif)
                            found_motif = True

                            self.diagnostic_list.append(DiagnosticRecord(dbd.get_name(), ortholog_description, ortholog_species, ortholog_escore, 
                                                                         ortholog_percent_identity, ortholog_query_coverage, dbd_similarity))
                            
                    except Exception as e:
                        if self.verbose:
                            print(f"\tERROR: SKIPPING {ortholog_description} from {ortholog_hit.get_species()} - {e}", flush=True)
                        continue
                    

                    if not found_motif:
                        if self.verbose:
                            print(f'\t\tWARNING: Unable to find {ortholog_description} in the database', flush=True)

                else:
                    if self.verbose:
                        print(f'\tWARNING: {dbd.get_name()} is not conserved enough with the DBD in {ortholog_description}! {sim_array[i]} similarity! Skipping...', flush=True)

        if is_optimize:
            print(f"Time to find motifs: {time.time() - start_time}", flush=True)

        return result_motif_list


    def call(self, gene_sequence = None, gene_file = None):
        '''
        This function predicts DNA-binding domain motifs of a species using orthologous species.

        Arguments:
        gene_sequence: gene DNA sequence
        gene_file: path to FASTA file containing gene sequence (scaffolds allowed); if gene_sequence is provided, this is ignored

        Returns:
        List of MotifDB.Motif objects
        '''

        # Clearing diagnostic list
        self.diagnostic_list.clear()

        # Getting gene sequence
        gene_sequence_list = []
        if gene_sequence:
            gene_sequence_list.append(gene_sequence)
        elif gene_file:
            gene_sequence_list = [str(x.seq) for x in SeqIO.parse(gene_file, "fasta")]

        # return value
        motif_list = []

        if self.verbose:
            print('---------------------------------', flush=True)


        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step one: get the protein sequence of the gene
        
        protein_sequence = self.find_protein(gene_sequence_list)

        if not protein_sequence:
            return motif_list

            
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step two: get ortholog sequences

        ortholog_hit_list = self.find_orthologs(protein_sequence)

        # if failed to find orthologs, return
        if not ortholog_hit_list:
            return motif_list # BAD PRACTICE; find some way to not need multiple return statements
        
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step three: get start and size of DBD within protein

        dbd_list = self.find_dbds([protein_sequence] + [ortholog_hit.get_seq() for ortholog_hit in ortholog_hit_list])

        # if failed to find DBD for gene of interest or any DBD for the orthologous sequences, return
        if not dbd_list[0] or not dbd_list[1:]:
            return motif_list # BAD PRACTICE; find some way to not need multiple return statements

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step four: find conserved motifs from orthologous sequences that had 70% or more identity with the gene from given species
        
        motif_list = self.find_motifs(dbd_list, ortholog_hit_list, protein_sequence)


        return motif_list
        


if __name__ == "__main__":

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Getting user inputs
    import argparse

    parser = argparse.ArgumentParser(description=
        '''The purpose of this script is to predict DNA-binding domain motifs of a species using orthologous species.
        Copyright: 
        Lab: Bradham Lab at Boston University
        Correspondence: anthonygarza124@gmail.com''')

    # Pipeline arguments
    parser.add_argument('-f', '--fasta', help='Path to folder of fasta files. Each fasta file should be for a singular gene. A fasta file may contain multiple scaffolds of this gene in multi-FASTA format.')

    parser.add_argument('-mdb', '--motif_database', help='Name of motif database to use. Default is JASPAR. Use -lmdb to see list of available motif databases.', default='JASPAR')
    parser.add_argument('-or', '--orthologs', nargs='*', help='Ortholog species to search for motif in. Each species should be enclosed in quotes and separated by spaces. Default is only Homo Sapiens') 
    parser.add_argument('-o', '--output', help='Output file path; provide a path to a file or a directory where output.tsv will be made.', default='output.tsv')

    # Pipeline TOOL Arguments
    parser.add_argument('-e', '--email', help='Email address for EBML tools')

    # Pipeline parameters
    parser.add_argument('-it', '--identity_threshold', help='Identity threshold for similarity between sequences. Default is 0.7', default=0.7, type=float)

    # User info requests
    parser.add_argument('-lmdb', '--list_motif_database', help='Prints list of motif databases available', action='store_true')

    # options
    parser.add_argument('-v', '--verbose', help='Prints out more information', action='store_true', default=False)


    args = parser.parse_args()

    fasta_path = args.fasta
    motif_database = args.motif_database
    ortholog_name_list = args.orthologs
    output_path = args.output

    email = args.email

    identity_threshold = args.identity_threshold

    list_motif_db = args.list_motif_database

    verbose = args.verbose



    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Printing list if asked for by user

    should_exit = False

    # Printing list of motif databases available
    if list_motif_db:
        print("List of motif databases available:", flush=True)
        print(', '.join(MotifDBFactory.get_motif_db_names()), flush=True)
        print()
        should_exit = True

    if should_exit:
        exit(0)


    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Validating arguments

    # Checking if required arguments are provided
    if not fasta_path:
        raise RuntimeError("Please provide a path to a folder of fasta files using -f or --fasta")

    if not os.path.isdir(fasta_path):
        raise RuntimeError("Fasta folder does not exist or is not a folder")

    if not email:
        raise RuntimeError("Please provide an email address using -e or --email")

    # Check if ortholog species passed; if not, use homo sapiens
    if not ortholog_name_list:
        ortholog_name_list = ['Homo Sapiens']
        

    # Checking if output path exists
    if os.path.isdir(output_path):
        output_path = os.path.join(output_path, "output.tsv")
    else:
        output_dir = os.path.dirname(output_path)
        if not os.path.exists(output_dir):
            raise RuntimeError("Output directory does not exist")


    # Getting all fasta files in the directory
    gene_file_list = [os.path.join(fasta_path, x) for x in os.listdir(fasta_path) if x.endswith('.fasta') or x.endswith('.fa') or x.endswith('.fna')]
    assert len(gene_file_list) > 0, "No fasta files found in the directory. Make sure the files end with .fasta, .fa, or .fna"
    gene_file_list.sort(key=lambda x: os.path.basename(x))

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Setting up Lverage

    lverage = Lverage(motif_database = motif_database, ortholog_name_list = ortholog_name_list, email = email, 
                      identity_threshold = identity_threshold, verbose = verbose)



    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # For every gene sequence in the fasta file directory, search for the motifs

    if verbose:
        print(f"Utilizing the following orthologs: {', '.join(ortholog_name_list)}", flush = True)
        print(f"Using an identity threshold of {identity_threshold * 100}%", flush = True)

    with open(output_path, 'w') as f:

        headers = ["Gene ID"] + DiagnosticRecord.header_list + MotifDBFactory.get_db_record(motif_database).header_list
        f.write('\t'.join(headers) + '\n')

        for gene_file in gene_file_list:

            gene_id = os.path.basename(gene_file).split('.')[0]

            if verbose:
                print(flush=True) # newline

            motif_list = lverage.call(gene_file = gene_file)
            diagnostic_list = lverage.diagnostic_list

            assert len(motif_list) == len(diagnostic_list), "Motif list and diagnostic list are not the same length"

            for i in range(len(motif_list)):
                diagnostic = diagnostic_list[i]
                motif = motif_list[i]

                diagnostic_values = diagnostic.get_values()
                motif_values = motif.get_values()

                # convert values to strings
                diagnostic_values = [str(x) for x in diagnostic_values]
                motif_values = [str(x) for x in motif_values]

                # write to file as tab-separated values
                s = '\t'.join([gene_id] + diagnostic_values + motif_values) + '\n'
                f.write(s)
                f.flush()





