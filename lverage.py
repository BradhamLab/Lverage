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
# global variables
is_optimize = True

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Getting user inputs before importing other modules
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
parser.add_argument('-c', '--clustalo', help='Path to Clustalo executable. If not provided, assumed to be in PATH')
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
ortholog_list = args.orthologs
output_path = args.output

clustalo_path = args.clustalo
email = args.email

identity_threshold = args.identity_threshold

list_motif_db = args.list_motif_database

verbose = args.verbose




#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Importing modules

import warnings
import os
import shutil
from Bio import BiopythonDeprecationWarning, SeqIO
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

from src.ProteinFinder import ProteinFinder
from src.DBDScanner import DBDScanner
from src.OrthologSearcher import OrthologSearcher
from src.Aligner import Aligner
from src.MotifDB import MotifDBFactory
import src.utils as lvutils

if is_optimize:
    import time # for optimization only

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Printing lists if asked for by user

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

# Check if motif database available
if not MotifDBFactory.has_db(motif_database):
    raise RuntimeError(f"{motif_database} not found in motif databases")

# Check if ortholog species passed; if not, use homo sapiens
if not ortholog_list:
    ortholog_list = ['Homo Sapiens']
    

# Checking if output path exists
if os.path.isdir(output_path):
    output_path = os.path.join(output_path, "output.tsv")
else:
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        raise RuntimeError("Output directory does not exist")


# Checking if Clustalo is available
if not clustalo_path:
    if shutil.which('clustalo') is None:
        raise RuntimeError("clustalo not found. Please provide path to clustalo executable using -c or --clustalo")
    clustalo_path = "clustalo"
else:
    assert os.path.exists(clustalo_path), "clustalo not ofund. Please provide path to clustalo executable using -c or --clustalo or ensure that clustal is in PATH"
    

# Checking if ortholog species are in NCBI database; if tax ID, converting to scientific name
for i, ortholog_species in enumerate(ortholog_list):
    if not lvutils.check_valid_species(ortholog_species):
        raise RuntimeError(f"{ortholog_species} not found in the NCBI database")

    # if tax ID, convert to species name
    if ortholog_species.isdigit():
        ortholog_list[i] = lvutils.get_species_name(ortholog_species)


# Getting all fasta files in the directory
gene_file_list = [os.path.join(fasta_path, x) for x in os.listdir(fasta_path) if x.endswith('.fasta') or x.endswith('.fa') or x.endswith('.fna')]
assert len(gene_file_list) > 0, "No fasta files found in the directory. Make sure the files end with .fasta, .fa, or .fna"
gene_file_list.sort(key=lambda x: os.path.basename(x))

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Setting up tools

# Setting up ORFfinder
pf = ProteinFinder(verbose=verbose)

# Setting up pfamscan
ds = DBDScanner(email=email, verbose=verbose)

# Setting up BlastP
ors = OrthologSearcher(hitlist_size=20, verbose=verbose, species_list=ortholog_list)

# Setting up Clustal Omega
aligner = Aligner(clustalo_path, verbose=verbose)

# Setting up motif database
mdb = MotifDBFactory.get_motif_db(db_name=motif_database, n_hits=10)

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Setting up output file

f = open(output_path, "w")

headers = ["ID", "DBD Name", "BLAST Hit Description", "BLAST Hit Species", "BLAST Hit E-value", "BLAST Hit Percent Identity (%)", "BLAST Hit Query Coverage",  "Motif Database Query", "Motif ID", "Motif Name", "Motif Class", "Motif Logo", "Motif PWM", "Uniprot ID", "Lv-Uniprot Identity (%)", "LvDBD-UniprotDBD Identity (%)"]
f.write('\t'.join(headers) + '\n')
f.flush()

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# For every gene sequence in the fasta file directory, search for the motifs

if verbose:
    print(f"Utilizing the following orthologs: {', '.join(ortholog_list)}", flush = True)
    print(f"Using an identity threshold of {identity_threshold * 100}%", flush = True)

for gene_file in gene_file_list:

    gene_id = os.path.basename(gene_file).split('.')[0]
    gene_sequence_list = [str(x.seq) for x in SeqIO.parse(gene_file, "fasta")]

    if verbose:
        print(flush=True) # newline

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step one: get the protein sequence of the gene
    if is_optimize:
        start_time = time.time()

    if verbose:
        print(f">Obtaining protein sequence of {gene_id}.", flush=True)

    # Get longest Open Reading Frame
    try:
        protein_sequence = ""
        for gene_sequence in gene_sequence_list:
            ps = pf.find_protein(gene_sequence)
            protein_sequence = ps if len(ps) > len(protein_sequence) else protein_sequence
            
    except Exception as e:
        print(f"\tERROR: SKIPPING {gene_id} - {e}")
        continue

    pf.reset_message()
    print(flush=True)
        
    if len(protein_sequence) == 0:
        if verbose:
            print(f"\tWARNING:Unable to obtain a protein sequence for {gene_id}", flush=True)
            continue
            
    if verbose:
        print("\t" + protein_sequence, flush=True)

    if is_optimize:
        print(f"Time to find protein sequence: {time.time() - start_time}", flush=True)


    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step two: get start and size of DBD within protein
    if is_optimize:
        start_time = time.time()

    if verbose:
        print(f">Searching for DBD in {gene_id}.", flush=True)

    try:
        dbd_list = ds.find_dbds(protein_sequence)
    except Exception as e:
        print(f"\tERROR: SKIPPING {gene_id} - {e}")
        continue
    
    # if failed to find DBD
    if len(dbd_list) == 0:

        if verbose:
            print(f"\tWARNING: Unable to find DBDs in {gene_id}, skipping.", flush=True)

        continue

    if is_optimize:
        print(f"Time to find DBD: {time.time() - start_time}", flush=True)
        
        
    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step three: get ortholog sequences
    if is_optimize:
        start_time = time.time()

    if verbose:
        print(f">Searching for orthologs of {gene_id}.", flush=True)

    try:
        ortholog_hit_list = ors.search(protein_sequence)
    except Exception as e:
        print(f"\tERROR: SKIPPING {gene_id} - {e}")
        continue

    # if failed to find orthologs
    if len(ortholog_hit_list) == 0:

        if verbose:
            print(f"\tWARNING:Unable to find orthologs of {gene_id}, skipping.", flush=True)

        continue

    if is_optimize:
        print(f"Time to find orthologs: {time.time() - start_time}", flush=True)
        
    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step four: multiple sequence alignment of protein sequences between species of interest and orthologs
    if is_optimize:
        start_time = time.time()

    if verbose:
        print(f">Aligning sequences of {gene_id}...", flush=True)

    try:
        alignment_list = aligner.align([protein_sequence] + [x.get_seq() for x in ortholog_hit_list])
    except Exception as e:
        print(f"\tERROR: SKIPPING {gene_id} - {e}")
        continue

    if is_optimize:
        print(f"Time to align sequences: {time.time() - start_time}", flush=True)

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step five: find conserved motifs from orthologous sequences that had 70% or more identity with the gene from given species
    start_time = time.time()

    for dbd in dbd_list:

        # getting similarity between alignments
        try:
            sim_array = lvutils.calculate_similarity(alignment_list, dbd.get_start(), dbd.get_size()) # similarity of all sequences to the first sequence at the dbd's location
        except Exception as e:
            print(f"\tERROR: SKIPPING DBD {dbd.get_name()} - {e}")
            continue

        if verbose:
            print(f">Finding conserved motifs of {gene_id} with dbd {dbd.get_name()}.", flush=True)
        

        
        # for each ortholog sequence (sim_array at 0 is species of interest, 1 and after are orthologs)
        for i in range(1, len(sim_array)):

            ortholog_hit = ortholog_hit_list[i-1] # -1 because ortholog_hit_list does not include species of interest
            ortholog_description = ortholog_hit.get_name()
            ortholog_reduced_description = ortholog_hit.get_usable_name()
            

            # If the similarity is greater than or equal to the identity threshold, get the motif
            if sim_array[i] >= identity_threshold:
                ortholog_species = ortholog_hit.get_species()
                ortholog_taxon_id = str(lvutils.get_tax_id(ortholog_species))

                if ortholog_taxon_id is None:
                    if verbose:
                        print(f'\tWARNING:Unable to find taxon ID for {ortholog_species}! Skipping...', flush=True)
                    continue

                ortholog_escore = ortholog_hit.get_escore()
                ortholog_percent_identity = ortholog_hit.get_percent_identity()
                ortholog_query_coverage = ortholog_hit.get_query_coverage()

                #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
                # Step six: for each orthologous sequence with a conserved dbd, search motif database for motif

                if verbose:
                    print(f"\n\t-Looking for ({ortholog_description}) with combinations of ({ortholog_reduced_description})", flush=True)
                    print("\t\t", end='', flush=True)

                found_motif = False # if a motif is found at least once, set this to true so user can later know if no motif was found
                
                # Use combinations of words from the ortholog's hit name; sometimes the full name doesn't work but a reduced form of it does
                for j, mdb_query in enumerate(ortholog_hit.get_usable_name_combos()):

                    try:
                        motif_list = mdb.search(mdb_query, ortholog_taxon_id, protein_sequence, dbd)
                    except Exception as e:
                        print(f"\tERROR: SKIPPING QUERY {mdb_query} - {e}")
                        continue

                    # if motif was found, add to output
                    if motif_list:

                        for motif in motif_list:


                            # This is same format as headers variable above the outermost for loop; keep the same!!!
                            output = [gene_id,
                                      dbd.get_name(),
                                      ortholog_description,
                                      ortholog_species,
                                      ortholog_escore,
                                      ortholog_percent_identity,
                                      ortholog_query_coverage,
                                      mdb_query,
                                      motif.get_matrix_id(),
                                      motif.get_name(),
                                      motif.get_motif_class(),
                                      motif.get_logo_link(),
                                      motif.get_pfm(),
                                      motif.get_uniprot_id(),
                                      motif.get_global_percent_identity(),
                                      motif.get_dbd_percent_identity()
                            ]
                            output = [str(x) for x in output]
                                      

                            f.write('\t'.join(output) + '\n')
                            f.flush()
                            found_motif = True

                            if verbose:
                                print('.', end='', flush=True)
                            #     print(f"\t\tFound logo for {mdb_query}", flush=True)

                        # If exact match found, break the loop early    
                        if j == 0:
                            break


                if not found_motif:
                    if verbose:
                        print(f'\t\tWARNING: Unable to find {ortholog_description} in the database', flush=True)
                        
            else:
                if verbose:
                    print(f'\tWARNING:{dbd.get_name()} is not conserved enough with dbd in {ortholog_description}! {sim_array[i]} similarity', flush=True)

    if is_optimize:
        print(f"Time to find motifs: {time.time() - start_time}", flush=True)
        
    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#

f.close()
