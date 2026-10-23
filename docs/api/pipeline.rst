Pipeline
========

``Lverage`` uses direct constructor injection for ORF searching, domain
scanning, ortholog searching, and motif databases. Constructors validate only
local configuration and do not contact remote services.

Each run deduplicates ORFs, selects the first longest ORF with an accepted
domain, compares matching base PFAM accessions, and searches motif databases
for domain pairs meeting the identity threshold. Intermediate evidence is not
retained by default. Callers may pass an optional ``evidence_collector`` to
``Lverage.run()`` to receive completed stage events as they occur; this supports
incremental reporting but is not checkpoint or crash recovery. The instance
retains its copied configuration and latest :class:`lverage.pipeline.LverageCode`.

``stop_after="alignment"`` deliberately returns before motif-database searches
while preserving full-pipeline behavior when the option is omitted.

Valid biological no-hit outcomes return an empty list and set a deterministic
code. Adapter and service failures remain exceptions.

.. autoclass:: lverage.pipeline.LverageCode
   :members:

.. autoclass:: lverage.pipeline.Lverage
   :members:
   :show-inheritance:
