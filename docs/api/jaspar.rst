JASPAR
======

The V2 adapter implements the JASPAR profile-inference calculations within
Lverage, adapted from source revision
``ac4024148f1f140c7370d7b13a70556b4d13f7a3`` of
``amathelier/JASPAR-inference-tool``. It uses normal package imports and local
JASPAR reference data; no upstream checkout, Git history, CLI entry point, or
separate inference interpreter is required. Runtime inference makes no online
requests. Exact-version motif matrices and annotations are resolved from a local
``JASPAR2026.sqlite`` database. The algorithm applies the upstream HMMER domain
thresholds and overlap selection, BLASTP/Rost filter, insertion removal,
positional identity/BLOSUM62 scoring, and Cis-BP HSim sequence-identity or
similarity-regression decisions.

The adapter filters by the upstream BLAST E-value, exact matrix species
taxonomy IDs, and then the configured hit limit.
The configured E-value cutoff applies specifically to the local inference
tool's BLAST E-value. Its equivalence to the former REST endpoint's ``evalue``
field has not been established; identical motif results between the services
are not claimed. Exact matrix versions returned by the fork are retained,
subject to species filtering, the E-value cutoff, and the hit limit.

The input is the unchanged ortholog-protein window of at most 2,000 residues
centered on the complete supporting domain. The inference BLAST E-value and DBD
identity are retained independently from Lverage BLAST metrics and Lverage's
global-alignment identity. Exact versioned matrix identifiers are required;
missing identifiers or local annotations are errors rather than an implicit
version substitution.

.. autoclass:: lverage.jaspar.JasparRecord
   :members:
   :show-inheritance:

.. autoclass:: lverage.jaspar.Jaspar2026MotifDB
   :members:
   :show-inheritance:
