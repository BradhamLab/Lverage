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
PfamScan installation/database. JASPAR2024 is queried through its online
service, so each actual run also requires network access. The local BLAST
database must support complete protein retrieval with ``blastdbcmd``. See the
README's V2 local BLAST database instructions for building one.

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

   lverage --input-dir ./inputs --blast-db ./db/proteins \
       --species-map ./config/species-taxids.json --pfam-db ./db/pfam \
       --output-dir ./results

Required options are ``--input-dir``, ``--blast-db``, ``--species-map``,
``--pfam-db``, and ``--output-dir``. The output directory must be new. It
contains ``results.tsv``, ``orfs.fasta``, and ``metadata.json``. Optional
executable overrides are ``--blastp`` (``blastp``),
``--blastdbcmd`` (``blastdbcmd``), and ``--pfamscan`` (``pfam_scan.pl``).
``--identity-threshold`` defaults to ``0.7``. Repeat ``--valid-pfam`` to
restrict accepted domains; if it is omitted, any detected domain is accepted.
JASPAR options default to 10 hits and an inference E-value threshold of
``1e-6``. Other component options retain their V2 defaults: ORF minimum length
75 nucleotides; BLAST E-value ``1e-6``, at most 20 hits, and the standard
description filters.

The CLI discovers ``.fa``, ``.fasta``, and ``.fna`` files nonrecursively,
case-insensitively, in filename order. The filename stem labels results and
must be unique ignoring case. Every FASTA record in a file is passed as an
independent sequence fragment to the same gene run. The CLI does not assemble
fragments, and it cannot determine from record headers whether they really
belong to the same gene. Empty files/records and unsupported DNA symbols are
rejected before external-tool adapters are constructed.

The TSV has explicit ``Row Type``, ``Stage``, ``Status``, ``Error``, and
``Outcome Detail`` fields. ``Error`` is reserved for execution failures;
biological reasons and not-reached explanations use ``Outcome Detail``.
Row types include ``input``, ``orf``, ``query_domain``, ``ortholog``,
``ortholog_retrieval``, ``ortholog_domain``, ``domain_comparison``, ``motif``,
and ``gene_error``.
``stage_status`` rows record biological no-hits, below-threshold completion,
intentional stops, and stages not reached.
``ortholog_retrieval`` rows retain failed retrieval attempts and skipped hits,
including available native hit IDs/accessions, the attempted identifier, and a
diagnostic reason. A BLAST no-hit status is distinct from retrieval failure.
Earlier-stage rows are not ortholog or motif hits; missing values are blank.
Each motif has one row and identifies its input gene, query ORF ID, ortholog,
and matching query/ortholog domains. Multiple hits, domain occurrences, and
motifs remain distinct. ORF rows retain candidates even when no accepted
domain exists, and an input row records a gene with no ORF. ``orfs.fasta``
stores each unique discovered candidate protein once under a deterministic
full-sequence SHA-256 identifier referenced by the TSV; gene associations and
selection status remain in the TSV. Candidate export does not change the
longest-first selection algorithm or run downstream analysis on unselected
candidates. BLAST fields include accession, native hit ID when returned, full
description, E-value, identity, and query coverage. Ratios and percentages are
separately labeled. Pfam IDs/names and every computed domain-pair identity are
kept, including values below the configured threshold.

``metadata.json`` records run timestamps, settings, input paths/checksums,
sequence-ID conventions, known tool versions (unknown versions are reported as
``unknown``), and completion/failure summaries. Results and metadata are
flushed as pipeline evidence completes. This is incremental evidence retention,
not checkpoint/resume or recovery from arbitrary crashes. After execution
begins, the bundle is preserved on failure. Execution failures are caught per
gene and processing continues to other genes, retaining completed evidence if
a later stage fails.

Use ``--stop-after alignment`` for pre-JASPAR validation. This intentional stop
does not call JASPAR and is not a failure. Exit status ``0`` means completed or
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
