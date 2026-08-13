BLAST searchers
===============

BLAST+ is an external prerequisite and is not installed as a Python package.
On the SCC, ``blast+/2.12.0`` and the ``human_mouse_nr`` database are useful
development examples; applications must still pass their own executable and
database locations.

Both adapters return :class:`lverage.ortholog_searcher.OrthologRecord` objects.
Identity and query coverage are ratios from zero to one. Individual malformed,
excluded, unresolved, or unretrievable hits are skipped and logged; executable,
database, subprocess, HTTP, timeout, and whole-result parsing failures are
raised to the caller.

``LocalBlastSearcher`` validates a database prefix with ``blastdbcmd -info``
and retrieves complete subject proteins from that database. The remote adapter
runs BLAST+ with ``-remote`` and batch-fetches complete proteins through NCBI
E-utilities. NCBI services are shared: do not parallelize remote searches, and
prefer off-peak hours for larger workloads.

.. autoclass:: lverage.blast.LocalBlastSearcher
   :members:
   :show-inheritance:

.. autoclass:: lverage.blast.RemoteBlastSearcher
   :members:
   :show-inheritance:
