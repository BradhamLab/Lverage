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
The legacy workflow requires Python 3.8 at minimum and was developed with
Python 3.10.12. The V2 package requires Python 3.11 or newer; see the V2 CLI
section below. To obtain Python and the required libraries, we advise using
Conda.

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

Lverage finds orthologs by running BLAST against a local protein database that
you build once. The database holds the proteins of the species you want
orthologs from. If you study species A and want orthologs from B, C, and D,
the database contains the proteomes of B, C, and D; A is the query, not part
of the database.

### 1. Install the tools

BLAST+ and the NCBI datasets CLI are external requirements and are not
installed with the Python package. With conda:

```bash
conda install -c conda-forge -c bioconda blast ncbi-datasets-cli unzip
```

### 2. Build the database

Look up the NCBI taxonomy ID of each ortholog species at
<https://www.ncbi.nlm.nih.gov/taxonomy> and pass them to the script. This
example builds a human, mouse, and fruit fly database:

```bash
scripts/make_blast_db.sh -t 9606,10090,7227 -o databases/orthologs
```

The script downloads each species' RefSeq reference proteome, combines them,
builds the database with `-parse_seqids` (Lverage needs this to retrieve each
hit's complete protein), and checks that retrieval works. It finishes by
printing the species names found in the headers.

If you already have a protein FASTA, build from it with `-i proteins.fasta`
instead. Headers must carry a unique identifier and the species in brackets:

```text
>NP_000001.1 Homeobox protein [Homo sapiens]
MPEPTIDESEQUENCE
```

If you have a taxonomy-enabled BLAST database such as NCBI `nr` on disk, add
`-n /path/to/nr` to extract the species from it rather than downloading. Run
`scripts/make_blast_db.sh -h` for all options.

### 3. Keep the output for Lverage

Note the database prefix you chose and the species names the script prints
at the end. Lverage needs both: the prefix says where to search, and the
species names, mapped to their taxonomy IDs, tell it which hits count as
orthologs. The names must match the headers exactly, including any strain
suffix such as `Saccharomyces cerevisiae S288C`. Bacterial proteomes label
many records with only a genus or group name, such as `Bacillus`, so map
every printed name you want to count, not just the species. Two species of
the same genus cannot be told apart by such a label, so expect their shared
proteins to be attributed to whichever taxonomy ID you map that name to. See
the
[API documentation](docs/api/index.rst) for how to pass them to the V2
pipeline.

## V2 Command-Line Interface

The V2 CLI is a new interface to the modular V2 package. It does not reuse or
extend the legacy CLI below. From the repository root, install the package in
your Python environment with `python -m pip install .` to get the `lverage`
command, or run it directly from a source checkout with
`PYTHONPATH=src python -m lverage.cli`. The CLI's `--help` works without
BLAST+, PfamScan, databases, or network access.

### Prerequisites

- Python 3.11+ and the package dependencies from this repository.
- Local BLAST+ `blastp` and `blastdbcmd` executables.
- A local protein BLAST database suitable for complete sequence retrieval.
- A local PfamScan installation and database.
- Network access when running motif inference: JASPAR2024 is queried online.

The local BLAST database setup above describes how to build a usable protein
database. Build and configure PfamScan separately according to its
distribution's instructions.

### Input mapping and example

Create a JSON file mapping bracketed scientific species names in BLAST hit
descriptions to positive NCBI taxonomy IDs, for example
`config/species-taxids.json`:

```json
{
	"Homo sapiens": 9606,
	"Mus musculus": 10090,
	"Drosophila melanogaster": 7227
}
```

Put one nucleotide FASTA file per gene/TF under `inputs/`. Supported suffixes
are `.fa`, `.fasta`, and `.fna`, matched case-insensitively. Multiple records
in one file are passed as independent fragments for that filename's gene;
they are never assembled, and headers are not used to verify that records
belong to the same gene. See the [V2 FASTA and CLI guide](docs/user_guide.rst)
for details.

```bash
lverage --input-dir ./inputs --blast-db ./db/proteins \
	--species-map ./config/species-taxids.json --pfam-db ./db/pfam \
	--output-dir ./results
```

The required options are `--input-dir`, `--blast-db`, `--species-map`,
`--pfam-db`, and `--output-dir`. The output path must not already exist; the
CLI creates it and writes `results.tsv`, `orfs.fasta`, and `metadata.json`.
Optional executable overrides are
`--blastp` (default `blastp`), `--blastdbcmd` (default `blastdbcmd`), and
`--pfamscan` (default `pfam_scan.pl`). `--identity-threshold` defaults to
`0.7`. Repeat `--valid-pfam ACCESSION` to limit accepted domains; when omitted,
any detected domain is accepted. JASPAR defaults to 10 hits and E-value
threshold `1e-6`; the existing ORF and BLAST component defaults are retained
(75-nt minimum ORF; BLAST E-value `1e-6`, top 20 hits, and standard description
filters).

The TSV contains explicit stage/status/error/outcome-detail fields, one motif
per motif row plus typed evidence rows for inputs, ORFs, query domains,
orthologs, ortholog domains, comparisons, stage statuses, and gene errors.
Earlier-stage rows are not motif or
ortholog hits; stage-status rows distinguish biological no-hits, below-
threshold completion, failure, and not-reached stages. Missing fields are blank. Rows retain BLAST accession,
description and native hit ID when returned, E-value, identity and query
coverage ratios and percentages, Pfam names/IDs, and each computed domain-pair
identity whether or not it meets the threshold. Candidate ORFs are uniquely
stored once in the FASTA by full-sequence SHA-256 ID; the TSV retains per-gene
associations and selected/not-selected/not-reached state. The existing
longest-first ORF-selection algorithm remains unchanged; downstream stages
run only on the selected ORF. `metadata.json` records settings, input
checksums, timestamps, tool versions when established (otherwise `unknown`),
errors, and completion counts. Bundle files are flushed as evidence arrives;
this is not checkpoint/resume or arbitrary-crash recovery.

Use `--stop-after alignment` for pre-JASPAR validation. This intentional stop
does not contact JASPAR and exits 0. Exit 0 also means full completion,
including biological no-hits; exit 2 means preflight validation failure; exit
1 means one or more execution failures. Per-gene errors are recorded and
processing continues to other genes, preserving completed partial evidence.
Diagnostics go to stderr. The output directory is not removed after execution
starts, including when a downstream operation fails.

## Legacy Command-Line Workflow

This older workflow is separate from the V2 package and CLI. The options and
examples in this section refer only to the legacy implementation.

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

### Input FASTA files

Lverage takes a folder of FASTA files, one file per gene or transcription factor. Name each file after its gene, since that name is used to label the output. Do not put different genes in the same file.

A file normally holds one record. If a gene's sequence is only available as fragments (for example several scaffolds or contigs), put all of the fragments in that gene's file as a multi-FASTA. Every record in the file is run through the pipeline for that gene: each fragment is searched for open reading frames on its own, the candidates from all fragments are pooled, and the best one is used. Fragments are not assembled or joined together.

```text
>geneX_fragment_1
ATGGCTAGCTAGCTAGGCTAGCTAGC...
>geneX_fragment_2
GGCTAGCATCGATCGTAGCTAGCTAG...
```

These are nucleotide sequences. They are separate from the protein FASTA used to build the local BLAST database described above. For more detail, see the [FASTA input guide](docs/user_guide.rst).



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
