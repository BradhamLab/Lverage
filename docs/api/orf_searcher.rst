ORF searcher
============

The V2 package provides an adapter for the tested ``orffinder==1.8`` release.
It receives nucleotide DNA and returns translated protein candidates; it does
not return nucleotide ORFs. The pipeline remains responsible for sorting and
deduplicating those candidates.

The default minimum length is 75 **nucleotides**, matching the package and the
legacy workflow. For a complete ORF, this includes the terminal stop codon: for
example, ``ATGAAATAG`` has a 9-nt ORF and yields the 2-aa protein ``MK`` after
the terminal stop is removed. This is not a protein-length threshold.

By default, only ATG start codons are supported. In the installed 1.8 release,
custom ``start_codons`` are silently ignored, so the adapter rejects custom
codons rather than suggesting that they work. Translation uses Biopython's
standard genetic code; the release accepts but ignores its ``translation_table``
argument, so the adapter does not expose a translation-table option.

The adapter does not use ``getORFProteins`` for translation: in 1.8 that
function can truncate a full codon for candidates ending at the sequence
boundary (and can return an empty protein for ``ATG``). Instead it uses
``getORFs`` to locate candidates and extracts complete codons from the original
input in the reported strand direction. It excludes the stop codon for complete
ORFs, and trims only the 1 or 2 leftover nucleotides from incomplete terminal
ORFs. No artificial bases are added.

``remove_nested=False`` disables the package's nested-candidate filter. It does
not guarantee that every overlapping ORF is enumerated; that depends on the
finder's candidate-generation logic. ``trim_trailing=False`` (the default)
allows ORFs without a terminal stop at a sequence boundary. The installed
release internally reports incomplete terminal candidates one nucleotide
shorter than their actual nucleotide span; the adapter compensates for this
and enforces ``minimum_length`` against the actual span. For example,
``ATGAAATG`` has an 8-nt terminal candidate and yields ``MK`` at minimums of 7
or 8 nt, but not 9 nt. Setting ``trim_trailing=True`` excludes that candidate.

Inputs must be nonempty and contain only case-insensitive IUPAC DNA symbols
``ACGTRYSWKMBDHVN``. Ambiguous codons can translate to ``X``. Whitespace,
gaps, RNA ``U``, and other symbols are rejected; input is not padded or
silently cleaned. ORFFinder scans all six frames itself.

.. autoclass:: lverage.orf_searcher.OrfSearcherTemplate
   :members:
   :show-inheritance:

.. autoclass:: lverage.orf_searcher.OrffinderOrfSearcher
   :members:
   :show-inheritance:
