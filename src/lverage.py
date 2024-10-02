#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
"""
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

Correspondence: 
    - Cynthia A. Bradham - cbradham@bu.edu - *
    - Anthony B. Garza   - abgarza@bu.edu  - **
    - Stephanie P. Hao   - sphao@bu.edu    - **
    - Yeting Li          - yetingli@bu.edu - **
    - Nofal Ouardaoui    - naouarda@bu.edu - **

    \* Principle Investigator, ** Software Developers
"""
#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Importing modules

import warnings
import os
import sys
from Bio import BiopythonDeprecationWarning, SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from orffinder import orffinder
from dataclasses import dataclass
import time
import shutil
from validate_email import validate_email


from DBDScanner import DBDScanner, DBD
from OrthologSearcher import OrthologSearcher
from MotifDB import MotifRecord
from MotifDB import MotifDBFactory
import lvutils

warnings.simplefilter('ignore', BiopythonDeprecationWarning)


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Classes

@dataclass
class LverageRecord:
    """
    """

    motif : MotifRecord # pass in a concrete class that inherits from abstract class MotifRecord
    dbd_name: str
    dbd_accession: str
    ortholog_description: str
    ortholog_species: str
    ortholog_escore: float
    ortholog_percent_identity: float
    ortholog_query_coverage: float
    ortholog_dbd_percent_identity: float

    header_list = ["Gene DBD Name", "Gene DBD Accession", "Ortholog BLAST Description", "Ortholog BLAST Species", "Ortholog BLAST E-value", "Ortholog BLAST Percent Identity", "Ortholog BLAST Query Coverage", "Gene DBD vs. Ortholog DBD Percent Identity"]

    def get_values(self):
        """

        """

        return [self.dbd_name, self.dbd_accession, self.ortholog_description, self.ortholog_species, 
                self.ortholog_escore, self.ortholog_percent_identity, self.ortholog_query_coverage, 
                self.ortholog_dbd_percent_identity] + self.motif.get_values()

class Lverage:
    """
    """


    def __init__(self, motif_database : str, ortholog_name_list : list, email : str, valid_dbd_list : list, 
                    identity_threshold : float = 0.7, verbose : bool = False, blast_hit_count : int = 20, 
                    motif_hit_count : int = 10, escore_threshold : float = 10**-6, blast_db_path : str = None,
                    blastp_path : str = None, start_codons : list = ['ATG'], orf_hit_count : int = 10):
        """Constructor
        """

        # Validating arguments

        # Checking if motif database is available
        if not MotifDBFactory.has_db(motif_database):
            raise ValueError("Motif database not found")
        
        # Checking if ortholog species are in NCBI database;
        for i, ortholog_species in enumerate(ortholog_name_list):
            if not lvutils.check_valid_species(ortholog_species):
                raise ValueError(f"Ortholog species ({ortholog_species}) not found in NCBI database")

        # Checking for at least one ortholog species
        if len(ortholog_name_list) == 0:
            raise ValueError("No ortholog species provided")

        # Checking if email is valid
        if not validate_email(email):
            raise ValueError("Email address is not valid")
        
        # Checking if identity threshold is between 0 and 1
        if not 0 <= identity_threshold <= 1:
            raise ValueError("Identity threshold must be between 0 and 1")
        
        # Checking if number of blast hits to search for is greater than 0
        if blast_hit_count <= 0:
            raise ValueError("Blast hit count must be greater than 0")
        
        # Checking if number of motif hits to search for is greater than 0
        if motif_hit_count <= 0:
            raise ValueError("Motif hit count must be greater than 0")
        
        # Checking if e-score threshold is greater than 0
        if escore_threshold <= 0:
            raise ValueError("E-score threshold must be greater than 0")
        
        # If using local BLAST
        if blast_db_path:
            if not os.path.exists(blast_db_path + ".pot"):
                raise FileNotFoundError(f"BLAST Database at {blast_db_path} doesn't exist!")
            
            # Checking if BLASTP executable is available
            if not os.path.exists(blastp_path):
                raise FileNotFoundError(f"BLASTP executable does not exist in {blastp_path}")
            
        # Checking if at least one start codon is provided
        if len(start_codons) == 0:
            raise ValueError("At least one start codon must be provided")
        
        # Checking if valid DBD list is not empty
        if len(valid_dbd_list) == 0:
            raise ValueError("Valid DBD list is empty")
        
        # Checking if ORF hit count is greater than 0
        if orf_hit_count <= 0:
            raise ValueError("ORF hit count must be greater than 0")
        
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Assigning arguments to class variables

        self.motif_database = motif_database
        self.ortholog_name_list = ortholog_name_list
        self.email = email
        self.identity_threshold = identity_threshold
        self.verbose = verbose
        self.blast_hit_count = blast_hit_count
        self.motif_hit_count = motif_hit_count
        self.escore_threshold = escore_threshold
        self.blast_db_path = blast_db_path
        self.blastp_path = blastp_path
        self.start_codons = start_codons
        self.valid_dbd_list = valid_dbd_list
        self.orf_hit_count = orf_hit_count


        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Setting up tools

        # Setting up pfamscan
        self.ds = DBDScanner(email = self.email, verbose = self.verbose)

        # Setting up aligner
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1

        # Setting up BlastP
        self.ors = OrthologSearcher(hitlist_size = self.blast_hit_count, verbose = self.verbose, 
                                    species_list = self.ortholog_name_list, email = self.email, 
                                    escore_threshold = self.escore_threshold, database_path = self.blast_db_path, 
                                    blastp_path = self.blastp_path)

        # Setting up motif database
        self.mdb = MotifDBFactory.get_motif_db(db_name = self.motif_database, n_hits = self.motif_hit_count, dbd_scanner = self.ds)


    
    @classmethod
    def from_dict(cls, dict_obj):
        """
        """
        return cls(
            motif_database = dict_obj['motif_database'],
            ortholog_name_list = dict_obj['ortholog_name_list'],
            email = dict_obj['email'],
            identity_threshold = dict_obj['identity_threshold'],
            verbose = dict_obj['verbose'],
            blast_hit_count = dict_obj['blast_hit_count'],
            motif_hit_count = dict_obj['motif_hit_count'],
            escore_threshold = dict_obj['escore_threshold']
        )

    def call(self, gene_sequences = None, gene_file = None):
        """
        """

        # Getting gene sequence
        gene_sequence_list = []
        if gene_sequences:
            gene_sequence_list.extend(gene_sequences)
        elif gene_file:
            gene_sequence_list = [str(x.seq) for x in SeqIO.parse(gene_file, "fasta")]


        # Getting protein sequence of gene
        protein_sequence = None
        gene_dbd_list = [] # list of valid DNA-binding domains of the gene

        if self.verbose:
            print(f"Searching for best Open Reading Frame for gene sequence", flush=True)

        orf_list = []
        for gene_sequence in gene_sequence_list:
            # Padding sequence with Ns if it's not a multiple of 3
            sequence = 'N' * (3 - len(gene_sequence) % 3) + gene_sequence if len(gene_sequence) % 3 != 0 else gene_sequence

            ps_list = orffinder.getORFProteins(SeqRecord(Seq(sequence)), start_codons=self.start_codons)
            if ps_list:
                orf_list.extend([str(ps.rstrip('*')) for ps in ps_list])

        orf_list = list(set(orf_list))
        orf_list.sort(key = lambda x: len(x), reverse = True)
        for orf in orf_list[:self.orf_hit_count]:

            if self.verbose:
                print(f"\tSearching for DNA-binding domains in ORF", flush=True)

            sys.stdout = lvutils.PrependedOutput(sys.stdout, "\t\t")
            dbd_list = self.ds.find_dbds(orf)
            sys.stdout = sys.stdout.original_stdout

            if not dbd_list:
                if self.verbose:
                    print(f"\t\tCould not find any DNA-binding domains, trying next ORF", flush=True)
                time.sleep(5)
                continue

            # Find all valid DBDs
            for dbd in dbd_list:
                accession = dbd.get_accession()
                base_accession = accession.split('.')[0]
                
                for valid_dbd in self.valid_dbd_list:

                    # If version number is in the valid DBD list, match accession and version. Otherwise, match just the accession
                    if ('.' in valid_dbd and accession == valid_dbd) or base_accession == valid_dbd:
                        protein_sequence = orf
                        gene_dbd_list.append(dbd)
                        break

            if gene_dbd_list:
                break

            if self.verbose:
                print(f"\t\tCould not find valid DNA-binding domains, trying next ORF", flush=True)

            time.sleep(5)

        if not protein_sequence:
            if self.verbose:
                print(f"\tCould not find any valid DNA-binding domains in any ORF", flush=True)
            return
        
        if self.verbose:
            print(f"Found protein sequence", flush=True)
            print(protein_sequence, flush=True)

        # Getting orthologs
        if self.verbose:
            print(f"Searching for orthologs of gene sequence", flush=True)

        sys.stdout = lvutils.PrependedOutput(sys.stdout, "\t")
        ortholog_hit_list = self.ors.search(protein_sequence)
        sys.stdout = sys.stdout.original_stdout

        if self.verbose:
            print(f"Removing invalid orthologs", flush=True)

        _ortholog_hit_list = []
        for ortholog_hit in ortholog_hit_list:
            if str(lvutils.get_tax_id(ortholog_hit.get_species())) is not None:
                _ortholog_hit_list.append(ortholog_hit)

            else:
                if self.verbose:
                    print(f"\t{ortholog_hit.get_title()} has a species ({ortholog_hit.get_species()}) not found in NCBI database", flush=True)
        ortholog_hit_list = _ortholog_hit_list

        if not ortholog_hit_list:
            if self.verbose:
                print(f"Could not find any valid orthologs", flush=True)
            return
        
        # Getting DNA-binding domains from orthologs that are shared with the gene's DNA-binding domains
        ortholog_dbd_list = []
        valid_ortholog_hits = []
        accession_order = [dbd.get_accession() for dbd in gene_dbd_list]

        if self.verbose:
            print(f"Searching for DNA-binding domains in orthologs and mapping to gene's DBDs", flush=True)

        for ortholog_hit in ortholog_hit_list:

            sys.stdout = lvutils.PrependedOutput(sys.stdout, "\t")
            dbd_list = self.ds.find_dbds(ortholog_hit.get_seq()) 
            sys.stdout = sys.stdout.original_stdout

            dbd_dict = {dbd.get_accession(): dbd for dbd in dbd_list} 

            ordered_dbd_list = []
            for accession in accession_order:
                if accession in dbd_dict:
                    ordered_dbd_list.append(dbd_dict[accession])
                # If the gene of interest has a DBD that the ortholog doesn't have
                else:
                    ordered_dbd_list.append(None) 

            # Only append if there's at least one valid DBD match
            if ordered_dbd_list:
                ortholog_dbd_list.append(ordered_dbd_list)
                valid_ortholog_hits.append(ortholog_hit)
            elif self.verbose:
                print(f"\t{ortholog_hit.get_accession()} does not have any DNA-binding domains that match the gene's DNA-binding domains. Removing ortholog", flush=True)

        ortholog_hit_list = valid_ortholog_hits

        if not ortholog_dbd_list:
            return
        
        # Aligning DBDs for similarity threshold check
        alignment_matrix = [] # [gene_dbd][ortholog]
        for i, dbd in enumerate(gene_dbd_list):
            gene_dbd_seq = protein_sequence[dbd.get_start():dbd.get_end()]
            across_ortholog_dbds = [dbd_list[i] for dbd_list in ortholog_dbd_list]
            alignment_list = []

            # for each ortholog
            for ortholog_hit, ortholog_dbd in zip(ortholog_hit_list, across_ortholog_dbds):

                # if a corresponding DBD from the gene exists for the ortholog, if not assume 0 similarity
                if ortholog_dbd:
                    ortholog_seq = ortholog_hit.get_seq()
                    ortholog_dbd_seq = ortholog_seq[ortholog_dbd.get_start():ortholog_dbd.get_end()]
                    alignment = self.aligner.align(gene_dbd_seq, ortholog_dbd_seq)[0]
                    similarity = lvutils.calculate_alignment_similarity(*alignment)
                    alignment_list.append(similarity)
                else:
                    alignment_list.append(0.0)

            alignment_matrix.append(alignment_list)

        # Finding motifs
        if self.verbose:
            print(f"Searching for motifs in gene sequence", flush=True)

        for i, dbd in enumerate(gene_dbd_list):

            if self.verbose:
                print(f"\tSearching with DBD {dbd.get_name()} {dbd.get_accession()}", flush=True)
                
            for j, ortholog_hit in enumerate(ortholog_hit_list):
                similarity = alignment_matrix[i][j]
                if similarity >= self.identity_threshold:

                    if self.verbose:
                        print(f"\t\tOrtholog {ortholog_hit.get_accession()} passed with a similarity of {similarity}", flush=True)

                    ortholog_description = ortholog_hit.get_name()
                    ortholog_species = ortholog_hit.get_species()
                    ortholog_taxon_id = str(lvutils.get_tax_id(ortholog_species))
                    ortholog_escore = ortholog_hit.get_escore()
                    ortholog_percent_identity = ortholog_hit.get_percent_identity()
                    ortholog_query_coverage = ortholog_hit.get_query_coverage()
                    ortholog_seq = ortholog_hit.get_seq()

                    sys.stdout = lvutils.PrependedOutput(sys.stdout, "\t\t\t")
                    try:
                        motif_list = self.mdb.search(protein_sequence, ortholog_taxon_id, ortholog_seq, dbd)

                    
                    except Exception as e:
                        if self.verbose:
                            print(f"\t\t\tAn error occurred while searching for motifs and the ortholog will be skipped. Error: {e}", flush=True)
                        continue

                    finally:
                        sys.stdout = sys.stdout.original_stdout

                    for motif in motif_list:
                        if self.verbose:
                            print(f"\t\t\tFound motif {motif.name} {motif.matrix_id}", flush=True)
                        yield LverageRecord(motif, dbd.get_name(), dbd.get_accession(), ortholog_description, ortholog_species, ortholog_escore, ortholog_percent_identity, ortholog_query_coverage, similarity)
                    
                    if not motif_list and self.verbose:
                        print(f"\t\t\tNo motifs found with the ortholog's DBD.", flush=True)

                else:
                    if self.verbose:
                        print(f"\t\tOrtholog {ortholog_hit.get_accession()} failed with a similarity of {similarity}", flush=True)


    
    
    def to_dict(self):
        """
        Returns a dictionary representation of the parameters used to create the object

        Returns
        -------
        dict
            dictionary representation of the object
        """

        return {
            'motif_database': self.motif_database,
            'ortholog_name_list': self.ortholog_name_list,
            'email': self.email,
            'identity_threshold': self.identity_threshold,
            'verbose': self.verbose,
            'blast_hit_count': self.blast_hit_count,
            'motif_hit_count': self.motif_hit_count,
            'escore_threshold': self.escore_threshold
        }
        


if __name__ == "__main__":

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Getting user inputs
    import argparse

    parser = argparse.ArgumentParser(description=
'''The purpose of this script is to predict DNA-binding domain motifs of a species using orthologous species.''')

    # Pipeline arguments
    parser.add_argument("-f", "--fasta", help="REQUIRED. Path to folder of fasta files. Each fasta file should be for a singular gene. A fasta file may contain multiple scaffolds of this gene in multi-FASTA format.")

    parser.add_argument("-mdb", "--motif_database", help="Name of motif database to use. Default is JASPAR. Use -lmdb to see list of available motif databases.", default="JASPAR")
    parser.add_argument("-or", "--orthologs", nargs="*", help="Ortholog species to search for motif in. Each species should be enclosed in quotes and separated by spaces. Default is only Homo Sapiens",
                        default=["Homo Sapiens"]) 
    parser.add_argument("-o", "--output", help="REQUIRED. Output file path; provide a path to a file or a directory where the Comma-Separated-Values (CSV) output file will be made.")
    parser.add_argument("-vd", "--valid_dbd", help="Valid DBDs file path; provide a path to a file containing valid DNA-Binding Domains. Each line except the first should contain a valid PFAM ID. The first row should contain the header 'PFAM ID'. If the user wishes to add extra columns detailing information about the domain, ensure the first column remains the PFAM ID and the other columns are separated by tabs.")

    # Protein translation arguments
    parser.add_argument("-sc", "--start_codons", nargs="*", help="When translating a gene DNA sequence to Amino Acids, which start codons should be used? Provide a space separated list of codons, e.g., \"ATG\" \"TTG\" \"GTG\". Default is ATG only.",
                        default=["ATG"])

    # NCBI Arguments
    parser.add_argument("-e", "--email", help="Email address for EBML tools")
    parser.add_argument("-b", "--blast", help="Path to the BLAST tools. If not provided, assumed to be on PATH.", default="")
    parser.add_argument("-bdb", "--blastdb", help="Path to the protein database to use with BLAST if being done locally. Ensure this database has the ortholog species of interest. If this argument is provided, BLAST WILL be run locally!", default = None)


    # Pipeline parameters
    parser.add_argument("-it", "--identity_threshold", help="Identity threshold for similarity between sequences. Default is 0.7", default=0.7, type=float)

    # User info requests
    parser.add_argument("-lmdb", "--list_motif_database", help="Prints list of motif databases available", action="store_true")

    # options
    parser.add_argument("-v", "--verbose", help="Prints out more information", action="store_true", default=False)


    args = parser.parse_args()

    fasta_path = args.fasta
    motif_database = args.motif_database
    ortholog_name_list = args.orthologs
    output_path = args.output
    valid_dbd_path = args.valid_dbd

    email = args.email
    blast_path = args.blast
    blast_db_path = args.blastdb

    identity_threshold = args.identity_threshold

    list_motif_db = args.list_motif_database

    verbose = args.verbose



    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Printing list of motif databases available if asked then exiting
    if list_motif_db:
        print("List of motif databases available:", flush=True)
        print(', '.join(MotifDBFactory.get_motif_db_names()), flush=True)
        print()
        exit(0)


    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Validating arguments

    def check_required(condition, error_message):
        if not condition:
            parser.print_help()
            parser.error(error_message)

    def validate_argument(condition, error_message, exception):
        if not condition:
            parser.print_help()
            raise exception(error_message)

    check_required(fasta_path, "Please provide a path to a folder of fasta files using -f or --fasta",)
    check_required(valid_dbd_path, "Please provide a path to a valid DBD file using -vd or --valid_dbd",)
    check_required(email, "Please provide an email address using -e or --email")

    validate_argument(os.path.isdir(fasta_path), "Fasta folder does not exist or is not a folder", FileNotFoundError)
    validate_argument(os.path.exists(valid_dbd_path), "Valid DBD file does not exist", FileNotFoundError)
    validate_argument(os.path.exists(fasta_path), "Fasta folder does not exist", FileNotFoundError)
    validate_argument(validate_email(email), "Email address is not valid", ValueError)

    if os.path.isdir(output_path):
        output_path = os.path.join(output_path, "output.tsv")
    else:
        output_dir = os.path.dirname(output_path)
        validate_argument(os.path.isdir(output_dir), f"Output directory {output_dir} does not exist", FileNotFoundError)

    # Getting valid DBDs
    with open(valid_dbd_path, 'r') as f:
        lines = [line.strip() for line in f.readlines()]
        if len(lines) <= 1:
            raise ValueError("Valid DBDs file is empty")

        if len(lines[0].split('\t')) > 1:
            valid_dbd_list = [line.split('\t')[0] for line in lines[1:]]
        else:
            valid_dbd_list = lines[1:]


    # Checking if using local BLAST and if so, if the blastp executable is available and if the database is avaialble
    blastp_path = None
    if blast_db_path:

        validate_argument(os.path.exists(blast_db_path + ".pot"), f"BLAST Database at {blast_db_path} doesn't exist!", FileNotFoundError)

        # checking if blastp executable is available
        if blast_path:
            blastp_path = os.path.join(blast_path, "blastp")
            validate_argument(os.path.exists(blastp_path), f"BLASTP executable does not exist in {blast_path}", FileNotFoundError)
        else:
            blastp_path = shutil.which("blastp")
            validate_argument(blastp_path, "blastp executable is not in PATH", FileNotFoundError)



    # Getting all fasta files in the directory
    gene_file_list = [os.path.join(fasta_path, x) for x in os.listdir(fasta_path)]
    gene_file_list.sort(key=lambda x: os.path.basename(x))
    validate_argument(len(gene_file_list) > 0, "No fasta files found in the directory. Make sure the files end with .fasta, .fa, or .fna", FileNotFoundError)
    
    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Setting up Lverage

    lverage = Lverage(motif_database = motif_database, ortholog_name_list = ortholog_name_list, email = email, 
                      identity_threshold = identity_threshold, verbose = verbose, blast_db_path = blast_db_path, blastp_path = blastp_path,
                      valid_dbd_list = valid_dbd_list)

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # For every gene sequence in the fasta file directory, search for the motifs

    if verbose:
        print(f"Utilizing the following orthologs: {', '.join(ortholog_name_list)}", flush = True)
        print(f"Using an identity threshold of {identity_threshold * 100}%", flush = True)

    with open(output_path, 'w') as f:

        headers = ["Gene ID"] + LverageRecord.header_list + MotifDBFactory.get_db_record(motif_database).header_list
        f.write('\t'.join(headers) + '\n')
        f.flush()

        for gene_file in gene_file_list:

            gene_id = os.path.basename(gene_file).split('.')[0]

            if verbose:
                print(flush=True)
                print(f"Processing gene {gene_id}", flush=True)

            for record in lverage.call(gene_file = gene_file):

                record_values = [str(x) for x in record.get_values()]

                # write to file as comma-separated values
                s = '\t'.join([gene_id] + record_values) + '\n'
                f.write(s)
                f.flush()





