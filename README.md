![Logo](static/logo.jpg "Lverage Logo")

# Lverage (WIP)
Motif Finder pipeline - searching for motifs through orthologous species.

![Workflow of Lverage. Obtaining and comparing DBDs between transcription factors to obtain motifs.](static/LVerage_v5.png "Lverage Workflow")

## Install Requirements

### Preparing the Directory
Prepare a directory for Lverage. Then either download the ZIP file or clone the repository to this directory. Once downloaded, navigate inside.
```
git clone https://github.com/BradhamLab/Lverage.git
cd Lverage
```

### Python
Lverage requires Python3.8 at minimum to run. The tool was developed using Python3.10.12. To obtain python and the required libraries, we advise using Conda.

#### Downloading Conda
To install Conda, follow the instructions on their [website user guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html "https://conda.io/projects/conda/en/latest/user-guide/install/index.html"). We suggest downloading Miniconda.

#### Installing Python and Packages
The following command installs our configuration of Python.
```
conda env create -f conda_lverage.yml
conda activate lverage
```

#### Alternative
The user can also download Python and the required libraries themselves. Follow the instructions at the [Python website](www.python.org "www.python.org") to download Python3.10.12 (or any other alternative with a minimum of 3.8). 

Then the user must install the required libraries. Located in the repository is a requirements.txt file which contains the required libraries. Assuming that calling *python3* calls the user-downloaded python, please call the following command.

```
python3 -m pip install -r requirements.txt
```

## Local BLAST Database Setup for V2

Lverage V2 can search a local protein database with BLAST+. BLAST+ is an
external requirement and is not installed with the Python package. Confirm that
both required programs are available:

```bash
blastp -version
blastdbcmd -version
```

The source FASTA must contain complete protein sequences. Each record must have
a unique sequence identifier immediately after `>`, followed by a description
containing the scientific name in square brackets. For example:

```text
>NP_000001.1 Homeobox protein [Homo sapiens]
MPEPTIDESEQUENCE
>NP_000002.1 Homeobox protein [Mus musculus]
MSEQUENCE
```

Build the protein database with `-parse_seqids`. This option is required because
Lverage retrieves the complete sequence of each BLAST hit with `blastdbcmd`.
The value passed to `-out` is a database prefix, not a directory or an
individual database file. See the
[NCBI BLAST+ documentation](https://www.ncbi.nlm.nih.gov/books/NBK52637/table/blast_setup_pc.T.programs_and_utilities/)
for the specific-retrieval requirement.

```bash
mkdir -p databases

makeblastdb \
    -in proteins.fasta \
    -dbtype prot \
    -parse_seqids \
    -out databases/lverage_proteins
```

Do not omit `-parse_seqids`. A database without parsed sequence identifiers can
be searched by `blastp`, but Lverage will be unable to retrieve complete subject
proteins from it.

Validate both the database and accession retrieval before using it:

```bash
blastdbcmd -db databases/lverage_proteins -info

blastdbcmd \
    -db databases/lverage_proteins \
    -entry NP_000001.1 \
    -outfmt '%f'
```

The second command should print the expected FASTA record. If it reports that
the database contains no accession information, rebuild the database from its
source FASTA with `-parse_seqids`. Use a new output prefix rather than
overwriting a database that may be in use.

If a taxonomy-enabled source database such as NCBI `nr` is already installed,
a smaller FASTA can first be extracted. The following example selects human and
mouse proteins:

```bash
blastdbcmd \
    -db /path/to/nr \
    -taxids 9606,10090 \
    -target_only \
    -outfmt '%f' \
    -out human_mouse.fasta

makeblastdb \
    -in human_mouse.fasta \
    -dbtype prot \
    -parse_seqids \
    -out databases/human_mouse_nr
```

Taxonomic extraction requires the taxonomy files associated with the source
BLAST database. Check the resulting FASTA headers and sequence count before
building the new database.

Pass the database prefix and a mapping of the scientific names used in its
FASTA headers to `LocalBlastSearcher`:

```python
from lverage.blast import LocalBlastSearcher


ortholog_searcher = LocalBlastSearcher(
    database_path="databases/human_mouse_nr",
    species_map={
        "Homo sapiens": 9606,
        "Mus musculus": 10090,
    },
)
```

## How to Use
We warn against moving any file within the directory anywhere else as this will create errors. If you wish to access from other places, we suggest appending the directory to your PATH environment variable, creating an alias, or creating a shortcut.

To call Lverage, ensure that the requirements above are all met. Within the directory is a file called *lverage.py*. All calls should be made with this file.

The following table shows all arguments for Lverage.


|Argument|Description|
|---|---|
|-h/--help| Provides a description of the tool and arguments |
|-f/--fasta| Path to folder of fasta files. Each fasta file should be for a singular gene. A fasta file may contain multiple scaffolds of this gene in multi-FASTA format. The name of the file should be the gene's name. |
|-mdb/--motif_database| Motif database to search; currently only JASPAR which is default|
|-or/--orthologs|Ortholog species to search through. Povide the NCBI Tax IDs or scientific names, each one enclosed in quotes and separated by spaces|
|-o/--output|Output file path; if a directory is provided, output.tsv will be made there|
|-e/--email|Email address for EMBL Tools|
|-it/--identity_threshold| According to <insert paper here>, 70% similarity with an ortholog sequence means that the motif is conserved in the ortholog. This parameter asks that any ortholog sequence must be 70% similar to a provided gene's sequence.
|-v/--verbose|If provided, will print out every step along the way as well as intermittent reuslts|



### Examples

#### Calling Lverage on Green Sea Urchin

```
python3 lverage.py -f Data/LvGenes/ -e useremail@mail.com -o ../Output/lvedge_output.tsv -v
```
Here we provide a directory of fasta files with -f. The fasta files used in this example were gathered from LvEDGE and are provided in this repository. In -e, we provide an email for any EMBL tools. We provide a path for an output file we wish to be created with -o. Finally, we ask that it prints out each step with -v.


# Contributions

Lverage is an open-source tool made for the community. As such, we are welcoming of any contributions to the project!

If a bug/error is found, we suggest adding it as an issue in the GitHub Repository.

For directly contributing (fixing the error or adding new funcionality), we suggest that the contributor follows the (https://github.com/firstcontributions/first-contributions)[standard process]. Begin by forking the repository, making local changes, and then submit a pull request to the master branch. Be sure to name your branch a meaningful name that represents what changes you have added! We suggest that if a contributor wishes to add/fix multiple parts that they create a fork for each part.

Please note that this project uses the GNU AFFERO GENERAL PUBLIC LICENSE Version 3. Any contributions made will fall under this license.

# Contact Us

The team members are available to be contacted for any queries relating to **Lverage** usage and issues. We suggest first and foremost that any issue be posted to the issue board.

To contact us, please use the following contact information.


| Name               | Email               |
|--------------------|---------------------|
| Cynthia A. Bradham | cbradham@bu.edu     |
| Anthony B. Garza   | abgarza@bu.edu      |
| Stephanie P. Hao   | sphao@bu.edu        |
| Yeting Li          | yetingli@bu.edu     |
| Nofal Ouardaoui    | naouarda@bu.edu     |
| Thomas Shin        | thshin@bu.edu       |
