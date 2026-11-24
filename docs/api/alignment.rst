Alignment
=========

The alignment component receives one query domain and an ordered sequence of
matching ortholog-domain strings. Its ``calc_identities`` method returns one
identity ratio per input occurrence, in the same order. Duplicate sequences
must not be collapsed: the pipeline retains the corresponding ortholog and
domain records by their candidate positions.

``BiopythonGlobalAlignment`` is the default V2 implementation. Its default
affine-gap scoring preserves the existing global alignment configuration:
match 2, mismatch -1, gap open -2, and gap extension -0.5. That configuration
uses Gotoh's affine-gap global algorithm; it is not linear-gap
Needleman-Wunsch. To explicitly select Biopython's Needleman-Wunsch
implementation, set ``gap_model="linear"`` and give equal open and extension
scores, for example both -2. This option does not change the default.

Identity is exact, case-sensitive residue matches divided by all alignment
columns. Mismatches and gap columns remain in the denominator. The supplied
implementation uses pairwise alignments; the collection-based interface allows
a future multiple-sequence implementation without changing the pipeline
contract.

Pass a custom implementation as ``alignment_component`` when constructing
``Lverage``. It must inherit from ``AlignmentTemplate`` and return ordered
ratios in the range zero to one:

.. code-block:: python

   from lverage.alignment import BiopythonGlobalAlignment
   from lverage.pipeline import Lverage

   pipeline = Lverage(
       motif_database_list,
       orf_searcher,
       domain_scanner,
       ortholog_searcher,
       alignment_component=BiopythonGlobalAlignment(),
   )

.. autoclass:: lverage.alignment.AlignmentTemplate
   :members:
   :show-inheritance:

.. autoclass:: lverage.alignment.BiopythonGlobalAlignment
   :members:
   :show-inheritance: