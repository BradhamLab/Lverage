#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
'''
This file contains an interactive script that allows users to search and order the results of LvERAGE.

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

parser = argparse.ArgumentParser(description='Interactive script for searching and ordering LvERAGE results')
parser.add_argument('-i', '--input', type=str, help='Path to the LvERAGE output file', required=True)
parser.add_argument('-g', '--gui', action='store_true', help='If True, use the GUI interface', required=False, default=False)

args = parser.parse_args()
input_path = args.input
use_gui = args.gui


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Importing Modules
import os
import tkinter as tk
from tkinter import ttk
import tkinter.filedialog as filedialog


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Validating Arguments

assert os.path.exists(input_path), f"Input file does not exist: {input_path}"

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Defining Functions
    
def query_similarity(record):
    '''
    This function calculates the similarity of a query to the description it came from
    '''
    
    description = record.get_value("BLAST Hit Description")
    query = record.get_value("Motif Database Query")

    description = description.split()
    query = query.split()

            
    return len(query) / len(description)

def sort_records(records):
    '''
    This function sorts records by the following criteria:
    1. BLAST Hit Percent Identity
    2. Query Similarity (how similar the query is to the description it came from)
    3. LvDBD-UniprotDBD Identity
    '''
    
    records.sort(key=lambda x: (float(x.get_value('BLAST Hit Percent Identity (%)')), query_similarity(x), float(x.get_value('LvDBD-UniprotDBD Identity (%)'))), reverse=True)
    
    return records

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Defining Classes

class MotifRecord:
    '''
    This class represents a motif record, i.e., a line from the LvERAGE output file.
    '''

    record_headers = {x: i for i, x in enumerate(["ID", "DBD Name", "BLAST Hit Description", "BLAST Hit Species", "BLAST Hit E-value", "BLAST Hit Percent Identity (%)", "BLAST Hit Query Coverage",  "Motif Database Query", "Motif ID", "Motif Name", "Motif Class", "Motif Logo", "Motif PWM", "Uniprot ID", "Lv-Uniprot Identity (%)", "LvDBD-UniprotDBD Identity (%)"])}

    def __init__(self, data):
        '''
        Arguments:
        data: a list of strings representing the data in the motif record
        '''
        self.data = data

    def get_value(self, header):
        '''
        Arguments:
        header: the header of the value to get

        Returns:
        the value of the record corresponding to the header
        '''
        return self.data[self.record_headers[header]]
    
    def __str__(self):
        return '\n'.join([f"{x}: {self.data[self.record_headers[x]]}" for x in self.record_headers])
    

class LverageBrowser(tk.Tk):
    '''
    This class is a window that allows users to browse the results of LvERAGE
    '''

    def __init__(self, gene_dbd_dict):
        '''
        gene_dbd_dict: a dictionary that maps gene names to dictionaries that map DNA-Binding Domain names to dictionaries that map BLAST hit descriptions to lists of records
        '''

        super().__init__()
        self.gene_dbd_dict = gene_dbd_dict
        self.title("LvERAGE Browser")

        self.gene = ""
        self.dbd = ""
        self.blast_hit = ""

        self.create_widgets()

    def create_widgets(self):
        '''
        This function creates the widgets for the window and sets up the events
        '''

        # Create a frame for the listboxes
        self.listbox_frame = tk.Frame(self)
        self.listbox_frame.pack(side="top", fill="both", expand=True)
        
        # Create listboxes within the frame
        self.gene_listbox = tk.Listbox(self.listbox_frame)
        self.gene_listbox.bind('<<ListboxSelect>>', self.update_dbd_list)
        self.dbd_listbox = tk.Listbox(self.listbox_frame)
        self.dbd_listbox.bind('<<ListboxSelect>>', self.update_blast_hit_list)
        self.blast_hit_listbox = tk.Listbox(self.listbox_frame)
        self.blast_hit_listbox.bind('<<ListboxSelect>>', self.update_record_list)
        
        # Pack the listboxes side by side
        self.gene_listbox.pack(side="left", fill="both", expand=True)
        self.dbd_listbox.pack(side="left", fill="both", expand=True)
        self.blast_hit_listbox.pack(side="left", fill="both", expand=True)
        
        # Create and pack the records_text widget separately at the bottom
        self.records_text = tk.Text(self, height=10)
        self.records_text.pack(side="top", fill="both", expand=True)

        self.download_button = tk.Button(self, text="Download Records", command=self.download_records)
        self.download_button.pack(side="bottom", pady=10)

        self.populate_gene_list() # populating the gene listbox from the get-go

    def populate_gene_list(self):
        '''
        This function populates the gene listbox with the gene names
        '''

        for gene in self.gene_dbd_dict.keys():
            self.gene_listbox.insert(tk.END, gene)


    def update_dbd_list(self, event):
        '''
        This function updates the DNA-Binding Domain listbox when a gene is selected
        '''

        gene_selection = self.gene_listbox.curselection()
        if gene_selection:
            self.gene = self.gene_listbox.get(gene_selection[0])
            self.dbd_listbox.delete(0, tk.END) # clearing the dbd listbox
            self.blast_hit_listbox.delete(0, tk.END) # clearing the blast hit listbox
            self.records_text.delete('1.0', tk.END)

            # populating the dbd listbox given the gene
            for dbd in self.gene_dbd_dict[self.gene].keys():
                self.dbd_listbox.insert(tk.END, dbd)


    def update_blast_hit_list(self, event):
        '''
        This function updates the BLAST hit listbox when a DNA-Binding Domain is selected
        '''

        dbd_selection = self.dbd_listbox.curselection()
        if dbd_selection:
            self.dbd = self.dbd_listbox.get(dbd_selection[0])
            self.blast_hit_listbox.delete(0, tk.END) # clearing the blast hit listbox
            self.records_text.delete('1.0', tk.END)

            # populating the blast hit listbox given the gene and dbd
            for blast_hit in self.gene_dbd_dict[self.gene][self.dbd]:
                self.blast_hit_listbox.insert(tk.END, blast_hit)

    def update_record_list(self, event):
        '''
        This function updates the records text box when a BLAST hit is selected
        '''

        blast_hit_selection = self.blast_hit_listbox.curselection()
        if blast_hit_selection:
            self.blast_hit = self.blast_hit_listbox.get(blast_hit_selection[0])
            self.records_text.delete('1.0', tk.END)

            for record in self.gene_dbd_dict[self.gene][self.dbd][self.blast_hit]:
                self.records_text.insert(tk.END, f"{str(record)}\n")


    def download_records(self):
        '''
        This function allows users to download the records in the records text box
        '''
        save_path = filedialog.asksaveasfilename(defaultextension=".txt")
        if save_path:
            with open(save_path, 'w') as f:
                f.write(self.records_text.get('1.0', tk.END))

            print(f"Records downloaded to {save_path}")


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# Loading the LvERAGE output file

c = MotifRecord.record_headers

gene_dbd_dict = {} # gene -> dbd -> blast hit description -> list of records
with open(input_path, 'r') as f:

    # Skipping headers
    species = f.readline().replace("\n", "")
    f.readline()

    # Reading data
    line = f.readline()
    while line:
        data = line.split("\t")

        gene = data[c["ID"]]
        if gene not in gene_dbd_dict:
            gene_dbd_dict[gene] = {}

        dbd = data[c["DBD Name"]]
        if dbd not in gene_dbd_dict[gene]:
            gene_dbd_dict[gene][dbd] = {}

        blast_hit_description = data[c["BLAST Hit Description"]]
        if blast_hit_description not in gene_dbd_dict[gene][dbd]:
            gene_dbd_dict[gene][dbd][blast_hit_description] = []

        gene_dbd_dict[gene][dbd][blast_hit_description].append(MotifRecord(data))

        line = f.readline()

# Sorting records
for gene in gene_dbd_dict:
    for dbd in gene_dbd_dict[gene]:
        for blast_hit_description in gene_dbd_dict[gene][dbd]:
            gene_dbd_dict[gene][dbd][blast_hit_description] = sort_records(gene_dbd_dict[gene][dbd][blast_hit_description])

#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# If GUI interface
if use_gui:
    app = LverageBrowser (gene_dbd_dict)
    app.mainloop()


#@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#
# If Terminal Interactive
else:

    # User Querying

    print(f"Looking through results for {species} with {len(gene_dbd_dict)} genes", flush=True)
    print(f"To restart the interactive queries, type 'retry' at any of the inputs.", flush=True)
    print(f"To exit the interactive queries, type 'exit' at any of the inputs.", flush=True)
    print(flush=True)

    while True:
        usr_retry = False # if user wishes to restart the interactive input
        usr_exit = False # if user wishes to exit the interactive input
        
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@##
        
        # Step 1; ask usr for gene id
        # allow user to ask for gene ids to be displayed
        usr_gene = ""
        while not usr_gene:
            usr_gene = input(f"Please enter Gene ID. To display all Gene IDs ({len(gene_dbd_dict)}), press enter without providing any input: ")
            
            # if no input was provided, list gene ids
            if not usr_gene:
                for gene_id in gene_dbd_dict.keys():
                    print(gene_id, flush=True)
                print()
                continue
                
            # restart 
            if usr_gene == "retry":
                usr_retry = True
                continue
                
            if usr_gene == "exit":
                usr_exit = True
                break
                
            # validate user input
            if usr_gene not in gene_dbd_dict:
                usr_gene = ""
                print("Invalid Gene ID, try again.", flush=True)
                
            print(flush=True)
                
        if usr_retry:
            continue
            
        elif usr_exit:
            break
            
        dbd_dict = gene_dbd_dict[usr_gene]
                
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@##
        
        # Step 2: ask for DNA-Binding Domain name
        # allow user to ask for DBDs to be displayed
        usr_dbd = ""
        while not usr_dbd:
            usr_dbd = input(f"Please enter the DNA-Binding Domain to search through. To display all DNA-Binding Domains ({len(dbd_dict)}), press enter without providing any input: ")
            
            # if no input was provided, list DBD names
            if not usr_dbd:
                for dbd_name in dbd_dict.keys():
                    print(dbd_name, flush=True)
                print()
                continue
                
            # restart
            if usr_dbd == "retry":
                usr_retry = True
                continue
                
            if usr_dbd == "exit":
                usr_exit = True
                break
                
            # validate user input
            if usr_dbd not in dbd_dict:
                usr_dbd = ""
                print("Invalid DBD ID, try again.", flush=True)
                
            print(flush=True)
                
                
        if usr_retry:
            continue
            
        elif usr_exit:
            break
            
        description_dict = dbd_dict[usr_dbd]
        
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@##
        
        # Step 3: ask for which description to use
        # Allow user to ask for descriptions to be displayed

        usr_description = ""
        while not usr_description:
            usr_description = input(f"Please enter the ortholog description to search through. To display all descriptions ({len(description_dict)}), press enter without providing any input: ")
            
            # if no input was provided, list DBD names
            if not usr_description:
                for desc in description_dict.keys():
                    print(desc, flush=True)
                print()
                continue
                
            # restart
            if usr_description == "retry":
                usr_retry = True
                continue
                
            if usr_description == "exit":
                usr_exit = True
                break
                
            # validate user input
            if usr_description not in description_dict:
                usr_description = ""
                print("Invalid Description ID, try again.", flush=True)
                
            print(flush=True)
                
                
        if usr_retry:
            continue
            
        elif usr_exit:
            break
            
        description = usr_description.split() # converting usr inputted description into a list of words for later comparison
        motif_list = description_dict[usr_description]
        
        #@#@#@@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@##
        
        # Step 4: printing out motif records
        
        for record in motif_list:
            print(str(record), flush=True)
    

# done :(