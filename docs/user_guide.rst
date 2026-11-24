Gene and transcription-factor FASTA inputs
===========================================

Use nucleotide FASTA input to provide the gene or transcription-factor (TF)
sequence from which Lverage infers a protein and its DNA-binding domain.
Prepare one FASTA file per gene/TF and name the file after the gene, for
example ``geneX.fa``. Keep records for different genes in separate files. The
V2 command-line interface derives the gene identifier from each filename
stem; for the Python API, the caller reads files and supplies any identifier.

V2 command-line interface
-------------------------

The V2 CLI uses the current modular pipeline and is distinct from the legacy
command-line workflow described at the end of this guide. It can be invoked as
an installed ``lverage`` command or, from a source checkout, as
``PYTHONPATH=src python -m lverage.cli``. ``--help`` only displays arguments;
it does not construct adapters, inspect external tools, or access a network.

Prerequisites are the Python package dependencies, local BLAST+ executables
(``blastp`` and ``blastdbcmd``), a usable local protein database, and a local
PfamScan installation/database. Local motif inference uses the JASPAR profile
inference calculations adapted into Lverage from
``amathelier/JASPAR-inference-tool`` and local JASPAR 2026 reference data; it
makes no online requests. Users do not need that source checkout, its Git
history, an upstream entry point, or a separate Python environment. A run
stopped after alignment does not require or inspect JASPAR data. The local BLAST
database must support complete protein retrieval with ``blastdbcmd``. See the
README's V2 local BLAST database instructions for building one.

Configure the external inference data directory and local
``JASPAR2026.sqlite`` using CLI options or the ``LVERAGE_JASPAR_DATA_DIR`` and
``LVERAGE_JASPAR_SQLITE_DB`` environment variables. The implementation records
upstream source revision ``ac4024148f1f140c7370d7b13a70556b4d13f7a3`` and
version ``2026.1.0`` as algorithm provenance; neither is a runtime checkout
requirement. Keep reference data, databases, and temporary files outside the
repository. Configure BLASTP through ``--blastp`` and HMMER through ``--hmmscan``
and ``--hmmalign`` (or their executable names on ``PATH``); those external
programs are invoked directly with argument lists. The Python inference
implementation imports Biopython and NumPy from the normal Lverage environment.
The configured data directory is a standalone data tree, not the upstream
source checkout. If reusing its downloaded reference files, stage the data
outside the checkout; Lverage does not read the upstream source tree at runtime.

The adapted routines correspond to the upstream sequence-level inference flow:
load the per-taxon Pfam/profile/UniProt and Cis-BP models; detect and align
query Pfam domains with HMMER; search the local taxon BLAST databases and apply
the upstream Rost identity curve; remove HMM alignment insertions; calculate
the upstream positional identity and BLOSUM62 features; and apply the Cis-BP
identity or similarity-regression ``HSim`` decision before returning associated
JASPAR profiles. The upstream source revision is provenance only. The external
data directory contains the precomputed models, domain mappings, and BLAST
databases generated for the JASPAR 2026 profile set; the separate SQLite file
supplies exact-version matrices and annotations.

Create a JSON species-name-to-NCBI-tax-ID map. Species names must match the
bracketed names in BLAST descriptions (case-insensitive), and tax IDs must be
positive integers:

.. code-block:: json

   {
     "Homo sapiens": 9606,
     "Mus musculus": 10090
   }

Example using portable placeholder paths:

.. code-block:: console

      lverage --input-dir ./inputs --blast-db /path/to/local-data/blast/proteins \
         --species-map ./config/species-taxids.json \
         --pfam-db /path/to/local-data/pfam --output-dir /path/outside/repository/results \
         --jaspar-data-dir /path/to/local-data/jaspar-inference/files \
         --jaspar-sqlite-db /path/to/local-data/jaspar/JASPAR2026.sqlite \
         --blastp /path/to/blastp --hmmscan /path/to/hmmscan \
         --hmmalign /path/to/hmmalign

For a local pre-JASPAR run that stops after domain alignment:

.. code-block:: console

      lverage --input-dir ./inputs --blast-db /path/to/local-data/blast/proteins \
         --species-map ./config/species-taxids.json \
         --pfam-db /path/to/local-data/pfam \
         --output-dir /path/outside/repository/alignment-results \
         --stop-after alignment

Required options are ``--input-dir``, ``--blast-db``, ``--species-map``,
``--pfam-db``, and ``--output-dir``. The output directory must be new. It
contains four files: ``results.csv`` contains the scientific results,
``errors.csv`` contains execution diagnostics, ``orfs.fasta`` stores each
unique candidate protein sequence once, and ``metadata.json`` records run
settings and summaries. Optional executable
overrides are ``--blastp`` (``blastp``),
``--blastdbcmd`` (``blastdbcmd``), and ``--pfamscan`` (``pfam_scan.pl``).
The adapted motif inference also accepts ``--hmmscan`` and ``--hmmalign``
overrides.
``--identity-threshold`` defaults to ``0.7``. Repeat ``--valid-pfam`` to
restrict accepted domains; if it is omitted, every detected Pfam domain is
accepted. The pipeline selects at most one query ORF per gene: candidates from
the input fragments are deduplicated and sorted longest-first, then the first
candidate with an accepted domain is selected for downstream analysis.
JASPAR options default to 10 hits and a maximum inference-tool BLAST E-value of
``1e-6``; this cutoff applies to the BLAST E-value reported by the local
inference tool. Equivalence between that field and the former REST service's
``evalue`` field has not been established, so equivalent cutoff behavior or
identical motif results must not be assumed. The adapter retains all exact
versioned matrix IDs returned by the inference before applying the species
filter, this cutoff, and the configured hit limit. Set data paths using
``--jaspar-data-dir`` and ``--jaspar-sqlite-db``, or the corresponding
``LVERAGE_JASPAR_DATA_DIR`` and ``LVERAGE_JASPAR_SQLITE_DB`` environment
variables. Configure ``blastp``, ``hmmscan`` and ``hmmalign`` through their
executable options. The data directory must contain the locally stored JASPAR
profile-inference files, Pfam HMMs and pressed indexes, Cis-BP models, and the
per-taxon JSON and BLAST-database files used by the adapted algorithm. The
SQLite database can be the local JASPAR 2026 database distributed with
pyJASPAR; it supplies exact-version PFMs, TF classes, and species taxonomy
IDs. Other component options retain their V2 defaults: ORF minimum length
75 nucleotides; BLAST E-value ``1e-6``, at most 20 hits, and the standard
description filters. By default, descriptions containing ``hypothetical``,
``unnamed``, ``uncharacterized``, ``unknown``, or ``partial`` are excluded
case-insensitively. An ``isoform`` label alone does not exclude a protein:
isoform-labelled proteins are retained unless another exclusion applies.

The CLI discovers ``.fa``, ``.fasta``, and ``.fna`` files nonrecursively,
case-insensitively, in filename order. The filename stem labels results and
must be unique ignoring case. Every FASTA record in a file is passed as an
independent sequence fragment to the same gene run. The CLI does not assemble
fragments, and it cannot determine from record headers whether they really
belong to the same gene. Empty files/records and unsupported DNA symbols are
rejected before external-tool adapters are constructed.

``results.csv`` is a comma-separated UTF-8 table with exactly ten columns,
in this order:

.. list-table:: Results columns
   :header-rows: 1
   :widths: 35 65

   * - Column
     - Meaning
   * - Gene
     - Input FASTA filename without its extension.
   * - Ortholog Species
     - Species assigned to the retained BLAST hit, or parsed from a hit whose
       protein retrieval failed.
   * - BLAST Description
     - Description of the local protein BLAST hit.
   * - Query Pfam Name
     - Pfam domain detected in the candidate ORF protein.
   * - Hit Pfam Name
     - Pfam domain detected in the BLAST-hit protein.
   * - Global Alignment Identity Percent
     - Exact amino-acid matches divided by all columns, including gaps, in
       the global alignment of the paired Pfam domains, multiplied by 100.
       Values below the configured threshold are also retained.
   * - Motif Name
     - Name of the retained JASPAR profile resolved from the local database.
   * - JASPAR Tool DBD Identity Ratio
     - Domain identity reported by the adapted JASPAR inference algorithm,
       on a 0 to 1 scale. This is separate from Lverage's global-alignment
       identity; it is not a calibrated confidence or probability.
   * - Motif PFM JSON
     - Position frequency matrix serialized as a JSON object with A, C, G,
       and T arrays. Each array gives that nucleotide's count at successive
       motif positions. Standard CSV quoting keeps the JSON in one cell.
   * - ORF ID
     - ``orf_`` followed by the full SHA-256 of the candidate protein sequence.
       Use this identifier to retrieve its sequence from ``orfs.fasta``.

Genes are listed in input-filename order (alphabetical, case-insensitive).
Within each gene, rows are sorted by ``JASPAR Tool DBD Identity Ratio`` from
highest to lowest, across all supporting ORFs, BLAST hits, and domain
comparisons. The highest-identity retained motif appears first for that gene.
Rows with blank ratios follow scored rows; zero is a valid score and precedes
blanks. Equal scores preserve the original evidence order: longest-first ORF,
query-domain occurrence, BLAST hit, hit-domain occurrence, then inference
result order. A gene's rows remain together, but its domain comparisons may be
interleaved by identity rank. This export sort does not change which ORF is
selected, which motifs pass filtering, or the configured motif hit limit.
There is one row per retained motif hit per supporting domain comparison.
Every applicable supporting field is repeated on each motif row. Associations
are taken from the actual query and hit domains, not inferred from matching
names or from preceding spreadsheet rows.

Completed earlier evidence remains visible when there is no downstream result:
for example, an alignment below threshold retains its identity percentage and
has blank motif fields. A failed hit-domain scan keeps the available BLAST
fields and leaves the hit-domain and later fields blank. A hit domain without
a matching query domain has a blank query-domain field. ORFs without further
evidence still have a row with their gene and ORF ID; a gene with no ORFs has
only its gene filled in. Candidate export does not run downstream analysis on
unselected candidates.

There are no status, error, or outcome columns in this table. A blank cell
means that value was unavailable or that stage did not produce a result;
blanks alone do not distinguish an unexamined candidate, biological no-hit,
intentional stop, and execution failure. Do not forward-fill them. Consult
``errors.csv`` for execution failures and ``metadata.json`` for run completion
and settings. Biological no-hits are not execution errors.

The compact table intentionally omits matrix IDs, domain coordinates, candidate
order numbers, and selection flags. Different exact matrix versions or domain
occurrences can therefore produce identical-looking rows. They are preserved
as distinct results, not deduplicated. Versioned matrix IDs are still resolved
exactly during inference, with no fallback to another version, but cannot be
reconstructed from the exported motif name alone.

``errors.csv`` contains only events with an execution error, with these columns:
``Gene``, ``ORF ID``, ``Stage``, ``BLAST Hit ID``, ``BLAST Accession``,
``Query Pfam Name``, ``Hit Pfam Name``, ``Attempted Identifier``, and ``Error``.
Identifying context is included when available. Retrieval diagnostics include
failed attempts even when a subsequent attempt recovers the protein. A run
without recorded errors produces a header-only error file. Error messages do
not create additional diagnostic-only rows in ``results.csv``.

``orfs.fasta`` stores each unique discovered candidate protein once under its
full-sequence SHA-256 identifier. Identical proteins share a FASTA record;
the results table retains their per-gene associations.

The input to local motif inference remains the ortholog-protein window:
at most 2,000 amino acids centered on the supporting ortholog domain, retaining
the complete domain. Species are filtered against exact local NCBI taxonomy
IDs. The inference BLAST E-value cutoff still applies, although that metric is
not exported in the compact table. Its equivalence to the former REST E-value
has not been established; identical results between those services are not
claimed. The aligned domain sequences are not exported.

Metadata schema version 3 identifies this layout with
``results_format="compact-csv"`` and lists both CSV column sets. This replaces
the wider grouped CSV layout. Fields containing commas, quotation marks, or
line breaks follow standard CSV quoting rules.

``metadata.json`` records run timestamps, settings, input paths/checksums,
sequence-ID conventions, known tool versions (unknown versions are reported as
``unknown``), and completion/failure summaries. Progress metadata is updated
during execution. Scientific evidence and diagnostics are accumulated in memory,
ordered, and written to ``results.csv`` and ``errors.csv`` after processing all
genes. ``orfs.fasta`` is also written at that point. Each CSV is published via
a temporary file and atomic replacement. Execution failures are caught per
gene; processing continues to the other genes and the final files include
completed evidence from failed genes. This is not checkpoint/resume or crash
recovery: forcibly terminating the process may leave progress metadata without
final CSV or FASTA files.

Use ``--stop-after alignment`` for pre-JASPAR validation. This intentional stop
does not call JASPAR and does not require the JASPAR options or data; it is not
a failure. Exit status ``0`` means completed or
intentionally stopped, including biological no-hits; ``2`` means preflight
validation failure; ``1`` means one or more execution failures despite
preserved partial results. Diagnostics go to standard error.

FASTA format
------------

A FASTA record starts with a header line whose first character is ``>``.
Sequence data follows on one or more lines, up to the next header or end of
file. Sequence lines may be wrapped across lines for readability; wrapping does
not create additional records.

A single-record DNA file:

.. code-block:: text

   >geneX
   ATGCGTACGTAGCTAGCTACGATCG

A file may instead contain multiple records when they are fragments of the
same gene/TF, for example several scaffolds or contigs:

.. code-block:: text

   >geneX_fragment_1
   ATGCGTACGTAGCTAG
   CTACGATCG
   >geneX_fragment_2
   GCTAGCATCGATCGTA

Treat all records in the file as intended fragments of the same gene. The
headers themselves are not validated as evidence that fragments belong to
that gene. Keep unrelated genes or TFs in separate files.

The V2 ORF adapter accepts nonempty, case-insensitive DNA strings using the
IUPAC symbols ``A C G T R Y S W K M B D H V N``. It does not accept RNA ``U``,
gaps, or other characters. Wrapped sequence lines are fine; a DNA string passed
directly to the adapter must not contain whitespace.

Fragments are independent inputs
--------------------------------

For a multi-record gene file, read each record as a separate sequence. Lverage
searches each fragment independently for ORFs. It does **not** concatenate,
assemble, or splice fragments. Candidate protein sequences found across the
fragments are pooled, deduplicated, and sorted from longest to shortest; the
pipeline then scans candidates in that order and selects the first one with an
accepted domain. Thus, fragments are treated as alternative pieces of input for
the same gene, not as pieces that Lverage joins into a reconstructed gene.

V2 Python API: read records, pass sequence strings
---------------------------------------------------

The V2 ``Lverage.run()`` method accepts either one DNA string or a list of DNA
strings for one gene/TF. It does not accept a FASTA filename. Read the records
from the file with Biopython, then pass the record sequences to an already
configured V2 pipeline:

.. code-block:: python

   from Bio import SeqIO

   # ``pipeline`` is an already configured lverage.pipeline.Lverage instance.
   with open("inputs/geneX_homeobox.fa", encoding="utf-8") as fasta_handle:
       fragments = [
           str(record.seq)
           for record in SeqIO.parse(fasta_handle, "fasta")
       ]

   results = pipeline.run(fragments)

For a file with one record, this list contains one sequence. Passing the list
also works for that case. The configured pipeline supplies the ORF searcher,
domain scanner, ortholog searcher, and motif database; FASTA parsing is the
caller's responsibility.

Nucleotide FASTA and protein FASTA are different inputs
-------------------------------------------------------

The per-gene FASTA described above contains nucleotide DNA for motif inference.
It is not the same as a protein FASTA used by protein-analysis tools. PfamScan
(or another protein-domain scanner) operates on translated protein sequences;
a BLAST protein database is built from protein FASTA records and searched with
BLASTP. Do not use those protein records as the gene DNA input to motif
inference. The README's V2 local BLAST database section describes the protein
FASTA used to build that database.

Legacy command-line workflow and V2 Python API
----------------------------------------------

The legacy command-line workflow is separate from the V2 Python API. Its
``-f/--fasta`` option takes a directory containing one file per gene/TF; a
file may contain one record or multiple records for fragments of that same
gene. The legacy implementation reads all FASTA records in each file.

The V2 CLI is documented above. For the V2 Python API, configure the pipeline
and call ``Lverage.run()`` with a DNA string or list of DNA strings as shown
above.
