BLAST searchers
===============

BLAST+ is an external prerequisite and is not installed as a Python package.
Configure the executable and database locations for the environment where
Lverage is run.

Local databases must be built with ``makeblastdb -parse_seqids`` because
specific complete-sequence retrieval through ``blastdbcmd`` requires a parsed
sequence-ID index. Rebuild databases that lack this index before using them
with ``LocalBlastSearcher``.

``LocalBlastSearcher`` returns :class:`lverage.ortholog_searcher.OrthologRecord`
objects.
Identity and query coverage are ratios from zero to one. Individual malformed,
excluded, unresolved, or unretrievable hits are skipped and logged; executable,
database, subprocess, timeout, and whole-result parsing failures are raised to
the caller.

The existing ``get_orthologs(sequence)`` call remains unchanged. The shared
pipeline may additionally pass an evidence collector; in that mode each
retrieved ortholog and each failed/skipped sequence-retrieval attempt is
reported as it completes. This lets the CLI retain earlier hits if a later
retrieval fails and distinguishes a BLAST result with no alignments, alignments
filtered from retrieval, and retained hits whose sequences could not be
retrieved. The optional reporting mode does not rerun BLAST.

``LocalBlastSearcher`` validates a database prefix with ``blastdbcmd -info``
and retrieves complete subject proteins from that database.

.. autoclass:: lverage.blast.LocalBlastSearcher
   :members:
   :show-inheritance:
