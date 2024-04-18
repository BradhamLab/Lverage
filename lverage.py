#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
The purpose of this script is to predict DNA-binding domain motifs of a species using orthologous species.
This program works similar to a pipeline such that modules can be replaced with other modules to perform the same task but possibly on different species.

# ADD REQUIREMENTS HERE ANTHONY JESUS CHRIST DON'T FORGET

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
# Getting user inputs before importing other modules
import argparse

parser = argparse.ArgumentParser(description=
    '''The purpose of this script is to predict DNA-binding domain motifs of a species using orthologous species.
    This program works similar to a pipeline such that modules can be replaced with other modules to perform the same task but possibly on different species.
    Copyright: 
    Lab: Bradham Lab at Boston University
    Correspondence: anthonygarza124@gmail.com''')

# Pipeline arguments
parser.add_argument('-db', '--database', help='Path to database file for species of interest if required')
parser.add_argument('-s', '--species', help='Species name associated with the database and its parser. Use -lsdb to see list of available species.')
parser.add_argument('-mdb', '--motif_database', help='Name of motif database to use. Default is JASPAR. Use -lmdb to see list of available motif databases.', default='JASPAR')
parser.add_argument('-or', '--orthologs', nargs='*', help='Ortholog species to search for motif in. Each species should be enclosed in quotes and separated by spaces. Default is only Homo Sapiens') 
parser.add_argument('-o', '--output', help='Output file path; provide a path to a file or a directory where output.tsv will be made.', default='output.tsv')

# Pipeline TOOL Arguments
parser.add_argument('-c', '--clustalo', help='Path to Clustalo executable. If not provided, assumed to be in PATH')
parser.add_argument('-e', '--email', help='Email address for EBML tools')

# Pipeline parameters
parser.add_argument('-it', '--identity_threshold', help='Identity threshold for similarity between sequences. Default is 0.7', default=0.7, type=float)

# User info requests
parser.add_argument('-lsdb', '--list_species_database', help='Prints list of species who have database parsers available', action='store_true')
parser.add_argument('-lmdb', '--list_motif_database', help='Prints list of motif databases available', action='store_true')
parser.add_argument('-ls', '--list_species', help='Prints list of all species available as orthologs', action='store_true')

parser.add_argument('-v', '--verbose', help='Prints out more information', action='store_true', default=False)




args = parser.parse_args()

database_path = args.database
species = args.species
motif_database = args.motif_database
ortholog_list = args.orthologs
output_path = args.output

clustalo_path = args.clustalo
email = args.email

identity_threshold = args.identity_threshold

list_species_db = args.list_species_database
list_motif_db = args.list_motif_database
list_species = args.list_species

verbose = args.verbose




#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Importing modules

import warnings
from Bio import BiopythonDeprecationWarning
import shutil
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

from SpeciesDB import SpeciesDBFactory
from ProteinFinder import ProteinFinder
from DBDScanner import DBDScanner
from OrthologSearcher import OrthologSearcher
from Aligner import Aligner
from MotifDB import MotifDBFactory
import utils
import os

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Printing lists if asked for by user

should_exit = False
# Printing list of species available to use
if list_species_db:
    print("List of species with available database parsers:", flush=True)
    print(', '.join(SpeciesDBFactory.get_species()), flush=True)
    print()
    should_exit = True

# Printing list of motif databases available
if list_motif_db:
    print("List of motif databases available:", flush=True)
    print(', '.join(MotifDBFactory.get_motif_db_names()), flush=True)
    print()
    should_exit = True

# Printing list of species available as orthologs
if list_species:
    print("List of species available as orthologs:", flush=True)
    print(', '.join(SpeciesDBFactory.get_all_species()), flush=True)
    should_exit = True

if should_exit:
    exit(0)


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Validating arguments

# Checking if required arguments are provided
if not database_path or not species or not email:
    parser.print_help()
    raise RuntimeError("Please provide a database path, a species name, and an email")
    
# Check if ortholog species passed; if not, use homo sapiens
if not ortholog_list:
    ortholog_list = ['Homo Sapiens']
    
# Checking if ortholog species available in species.txt; if so, convert to standard name
for i, ortholog_species in enumerate(ortholog_list):
    if not SpeciesDBFactory.has_species(ortholog_species):
        raise RuntimeError(f"{ortholog_species} not found in species.txt")

    ortholog_list[i] = SpeciesDBFactory.get_standard_name(ortholog_species)


# Checking if database path exists
if not os.path.exists(database_path):
    raise RuntimeError("Database file does not exist")

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
    
#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Setting up tools

# Getting species database
sdb = SpeciesDBFactory.get_species_db(species, file_path=database_path)

# Setting up ORFfinder
pf = ProteinFinder(verbose=verbose)

# Setting up pfamscan
ds = DBDScanner(email=email, verbose=verbose)

# Setting up BlastP
os = OrthologSearcher(hitlist_size=20, verbose=verbose, species_list=ortholog_list)

# Setting up Clustal Omega
aligner = Aligner(clustalo_path, verbose=verbose)

# Setting up motif database
mdb = MotifDBFactory.get_motif_db(db_name=motif_database, n_hits=10)

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Setting up output file

f = open(output_path, "w")

headers = ["ID", "DBD Name", "BLAST Hit Description", "BLAST Hit Species", "BLAST Hit E-value", "BLAST Hit Percent Identity (%)", "BLAST Hit Query Coverage",  "Motif Database Query", "Motif ID", "Motif Name", "Motif Class", "Motif Logo", "Motif PWM", "Uniprot ID", "Lv-Uniprot Identity (%)", "LvDBD-UniprotDBD Identity (%)"]
f.write(species + '\n')
f.write('\t'.join(headers) + '\n')
f.flush()

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# For every gene sequence in the database, search for the motifs

if verbose:
    print(f"Searching the {species} database.", flush=True)
    print(f"Utilizing the following orthologs: {', '.join(ortholog_list)}", flush = True)
    print(f"Using an identity threshold of {identity_threshold * 100}%", flush = True)

for gene_id, gene_sequence_list in sdb:

    if verbose:
        print(flush=True) # newline

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step one: get the protein sequence of the gene
        
    if verbose:
        print(f">Obtaining protein sequence of {gene_id}.", flush=True)

    # Get longest Open Reading Frame
    protein_sequence = ""
    for gene_sequence in gene_sequence_list:
        ps = pf.find_protein(gene_sequence)
        protein_sequence = ps if len(ps) > len(protein_sequence) else protein_sequence
        
    pf.reset_message()
    print(flush=True)
        
    if len(protein_sequence) == 0:
        if verbose:
            print(f"\tWARNING:Unable to obtain a protein sequence for {gene_id}", flush=True)
            continue
            
    if verbose:
        print(protein_sequence, flush=True)

    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step two: get start and size of DBD within protein

    if verbose:
        print(f">Searching for DBD in {gene_id}.", flush=True)

    dbd_list = ds.find_dbds(protein_sequence)
    
    # if failed to find DBD
    if len(dbd_list) == 0:

        if verbose:
            print(f"\tWARNING: Unable to find DBDs in {gene_id}, skipping.", flush=True)

        continue
        
    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step three: get ortholog sequences

    if verbose:
        print(f">Searching for orthologs of {gene_id}.", flush=True)

    ortholog_hit_list = os.search(protein_sequence)

    # if failed to find orthologs
    if len(ortholog_hit_list) == 0:

        if verbose:
            print(f"\tWARNING:Unable to find orthologs of {gene_id}, skipping.", flush=True)

        continue
        
    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step four: multiple sequence alignment of protein sequences between species of interest and orthologs

    if verbose:
        print(f">Aligning sequences of {gene_id}...", flush=True)

    alignment_list = aligner.align([protein_sequence] + [x.get_seq() for x in ortholog_hit_list])
    
    for dbd in dbd_list:

        # getting similarity between alignments
        sim_array = utils.calculate_similarity(alignment_list, dbd.get_start(), dbd.get_size()) # similarity of all sequences to the first sequence at the dbd's location

        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
        # Step five: find conserved motifs from orthologous sequences that had 70% or more identity with the gene from given species

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
                ortholog_taxon_id = SpeciesDBFactory.get_taxon_id(ortholog_species)

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
                    motif_list = mdb.search(mdb_query, ortholog_taxon_id, protein_sequence, dbd)

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


    #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
    # Step seven: profit (forever out of our reach)

    # Step eight: ???

    # Step nine: cry
    

f.close()
