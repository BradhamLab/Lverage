BLAST searchers
===============

BLAST+ is an external prerequisite and is not installed as a Python package.
On the SCC, ``blast+/2.12.0`` and the ``human_mouse_nr`` database are useful
development examples; applications must still pass their own executable and
database locations.

Local databases must be built with ``makeblastdb -parse_seqids`` because
specific complete-sequence retrieval through ``blastdbcmd`` requires a parsed
sequence-ID index. The existing SCC ``human_mouse_nr`` database was built
without that index and must be rebuilt before it can support an end-to-end
``LocalBlastSearcher`` run.

``LocalBlastSearcher`` returns :class:`lverage.ortholog_searcher.OrthologRecord`
objects.
Identity and query coverage are ratios from zero to one. Individual malformed,
excluded, unresolved, or unretrievable hits are skipped and logged; executable,
database, subprocess, timeout, and whole-result parsing failures are raised to
the caller.

``LocalBlastSearcher`` validates a database prefix with ``blastdbcmd -info``
and retrieves complete subject proteins from that database.

.. autoclass:: lverage.blast.LocalBlastSearcher
   :members:
   :show-inheritance:
