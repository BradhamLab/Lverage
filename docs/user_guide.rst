Gene and transcription-factor FASTA inputs
===========================================

Use nucleotide FASTA input to provide the gene or transcription-factor (TF)
sequence from which Lverage infers a protein and its DNA-binding domain.
Prepare one FASTA file per gene/TF and name the file after the gene, for
example ``geneX.fa``. Keep records for different genes in separate files. The
command-line workflow requires this layout; for the Python API it is a
convention, since the caller reads the file.

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

The record headers identify fragments of the same gene in that one file. Do
not place unrelated genes or TFs together in a multi-record input file: use a
separate FASTA file for each one.

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

Legacy command-line workflow versus V2 API
------------------------------------------

The legacy command-line workflow is separate from the V2 Python API. Its
``-f/--fasta`` option takes a directory containing one file per gene/TF; a
file may contain one record or multiple records for fragments of that same
gene. The legacy implementation reads all FASTA records in each file.

For V2, configure the Python pipeline and call ``Lverage.run()`` with a DNA
string or list of DNA strings as shown above.
