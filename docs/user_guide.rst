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
contains three files: ``results.csv`` contains typed pipeline evidence,
``orfs.fasta`` stores each unique candidate protein sequence once, and
``metadata.json`` records run settings and summaries. Optional executable
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

``results.csv`` is a comma-separated UTF-8 table ordered by gene, ORF candidate,
query Pfam, BLAST hit, hit Pfam, global alignment, and motif hit. Each row is a
self-contained evidence record: every applicable gene, ORF, query-Pfam,
BLAST-hit, hit-Pfam, and alignment field is repeated on each descendant result
row. Columns are named for their stage (for example, ``Query Pfam ID``,
``BLAST Accession``, and ``Global Alignment Identity Ratio``). Candidates are
listed in longest-first examination order. The ORF examination and selection
columns identify examined candidates, candidates with no accepted domain, and
the selected candidate. Unexamined candidates still show their order, ID, and
length; their examination, selection, and downstream cells remain blank. All
unique ORFs are represented.

Within the selected ORF, each query-domain occurrence has its own occurrence
number and retains its Pfam ID, name, and zero-based, end-exclusive amino-acid
coordinates. BLAST hits are grouped beneath query Pfams and retain accession,
native hit ID, description, species, E-value, identity, and query coverage.
Occurrence and hit-group numbers, BLAST descriptive fields, and supporting
alignment fields are repeated on every associated motif row. When a hit is
evaluated against multiple query-Pfam occurrences, each row identifies the
specific query/hit-domain occurrence pair. Hit Pfam occurrences remain visible
even when no query Pfam matches them. Domain comparisons are associated by the
actual query-domain and hit-domain occurrences, not by protein name. Each
computed global domain alignment has an identity ratio, threshold, and pass/fail
decision. There is one row per motif hit per supporting domain comparison;
multiple motifs and multiple comparisons therefore remain distinct rows.

Rows remain grouped contiguously for convenient reading, but their meaning does
not depend on carrying values forward from preceding rows. Do not forward-fill
blank cells: a field is blank only when the information is unavailable or its
stage was not reached. No-hit outcomes and execution errors are recorded in the
relevant stage status, outcome, and error columns. Sorting or filtering may
change presentation order, but each row retains its own gene, candidate, hit,
domain-occurrence, alignment, and motif associations.
BLAST no hits and failed sequence retrieval are distinct outcomes. An empty
Pfam result, a domain pair below threshold, no matching Pfam, and a motif
search with no returned motifs are also distinguishable.

``orfs.fasta`` stores each unique discovered candidate protein once under a
deterministic full-sequence SHA-256 identifier referenced by the table; gene
associations and selection status are retained in the grouped rows. Candidate export does not run downstream analysis on
unselected candidates.

BLAST identity is the identity ratio from the HSP with the lowest E-value and
is distinct from domain-pair identity, which is the exact-match ratio across
every column of a global amino-acid alignment of matching Pfam domains
(including gaps).
The CSV reports BLAST identity and query coverage as ratios and percentages,
and records every computed domain-pair identity, including values below the
threshold. Domain start and end coordinates in exported rows are zero-based,
end-exclusive amino-acid positions. The domain-alignment identity score is
exported as a ratio and percentage, but the aligned sequences themselves are
not currently exported.

Motif rows separately record the upstream inference tool's BLAST E-value,
query/TF regions, TF name, and upstream DBD identity ratio. This upstream DBD
identity is not the Lverage global-alignment identity and is not substituted
for it. The local BLAST E-value is not established to be equivalent to the old
REST inference E-value; identical motif results between those services are not
claimed. The input sent to local inference is the same ortholog-protein window
previously sent to the REST endpoint: at most 2,000 amino acids centered on the
supporting ortholog domain, retaining the complete domain. It is not replaced
by the query protein or by a domain-only sequence. Matrix IDs must be
versioned, and each ID is resolved exactly in the configured local SQLite
database; there is no version fallback. Motif species are filtered against
exact local NCBI taxonomy IDs, not broad taxonomic groups.

Metadata schema version 2 identifies this layout with
``results_format="grouped-csv"``. It replaces the former event-row TSV schema. Fields containing commas,
quotation marks, or line breaks are quoted using standard CSV rules.

``metadata.json`` records run timestamps, settings, input paths/checksums,
sequence-ID conventions, known tool versions (unknown versions are reported as
``unknown``), and completion/failure summaries. Results and metadata are
flushed as pipeline evidence completes. This is incremental evidence retention,
not checkpoint/resume or recovery from arbitrary crashes. After execution
begins, the bundle is preserved on failure. Execution failures are caught per
gene and processing continues to other genes, retaining completed evidence if
a later stage fails.

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
