# Lverage
Motif Finder pipeline - searching for motifs through orthologous species.

Many species such as homo sapien and mus musculus have been well studied throughout the years with motifs well-characterized. However, this focus on specific model species has left a dearth in motifs for other interesting species. Motivated by [Insert Paper Here about 70%], we developed Lverage to locate candidate motifs in the genes of a specific species, leveraging information about motifs from orthologous species to do so.

[Insert graphic design of flow]

Given a gene DNA sequence, Lverage will first convert this sequence to the protein sequence, using the largest open-reading frame. Pfamscan from [ecbi.uk link here] is then used to obtain the DNA-binding domains (DBD). This is stored for later. BLASTP is ran on the protein sequence of the gene with the standard database to find orthologous hits. We then perform multiple-alignment [with Clustalo? Ask team. Maybe do multiple pair-wise alignments?] with the protein sequence and the BLASTP hits. If the DBD is highly conserved between species, this alignment should match along it. Using the DBD we located earlier, we can calculate the similarity of the DBD with every hit's DBD. Based on [paper about 70%], we use a threshold of 70% similarity as any DBD that meets this value implies that a motif is conserved between the species. If a hit meets this criterion, we use the JASPAR motif database to obtain the orthologous motif. We then validate this motif by utilizing the corresponding uniprot sequence delineated by JASPAR and ask for an 85% match (through local alignment) with our gene's protein sequence.

test

Note this is a readme draft. Imma just be jotting down things that are needed and I'll beautify this later. One day. Eventually.

I'll get to it when I get to it.


## Requirements
  1. Clustalo http://www.clustal.org/omega/
  3. Python 3.8 or above

Clone the directory and navigate inside of it.
  ```
  git clone ...
  cd LvERAGE
  ```
Please install the required python modules. @@@ TO ME ADD A REQUIREMENTS.TXT FILE PLEASE AND THANK YOU :) @@@
  ```
  pip install -r requirements.txt
  ```

## How To
  Not set in stone yet, please hold

  For the GUI, call...
   ```
   python3 reader.py -i lverage_output.tsv -g
   ```
  Where lverage_output.tsv is the output file of lverage.py

## Cite Us
  Maybe
