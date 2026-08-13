Pipeline
========

``Lverage`` uses direct constructor injection for ORF searching, domain
scanning, ortholog searching, and motif databases. Constructors validate only
local configuration and do not contact remote services.

Each run deduplicates ORFs, selects the first longest ORF with an accepted
domain, compares matching base PFAM accessions, and searches motif databases
for domain pairs meeting the identity threshold. Intermediate evidence remains
local to the run; the instance retains its copied configuration and the latest
:class:`lverage.pipeline.LverageCode`.

Valid biological no-hit outcomes return an empty list and set a deterministic
code. Adapter and service failures remain exceptions.

.. autoclass:: lverage.pipeline.LverageCode
   :members:

.. autoclass:: lverage.pipeline.Lverage
   :members:
   :show-inheritance:
