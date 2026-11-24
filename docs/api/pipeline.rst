Pipeline
========

``Lverage`` uses direct constructor injection for ORF searching, domain
scanning, ortholog searching, motif databases, and domain alignment.
Constructors validate only local configuration and do not contact remote
services. See :doc:`alignment` for the identity contract, scoring options,
and how to inject a different alignment implementation.

Each run deduplicates ORFs, selects the first longest ORF with an accepted
domain, compares matching base PFAM accessions, and searches motif databases
for domain pairs meeting the identity threshold. Intermediate evidence is not
retained by default. Callers may pass an optional ``evidence_collector`` to
``Lverage.run()`` to receive completed stage events as they occur; this supports
incremental reporting but is not checkpoint or crash recovery. The instance
retains its copied configuration and latest :class:`lverage.pipeline.LverageCode`.

``stop_after="alignment"`` deliberately returns before motif-database searches
while preserving full-pipeline behavior when the option is omitted.

For each query-domain occurrence, the pipeline gathers matching ortholog
domains in ortholog/domain order and calls the alignment component once. It
keeps score-to-domain associations by position, including duplicate domain
sequences, and emits accepted results in the existing ortholog/query-domain/
ortholog-domain order.

Valid biological no-hit outcomes return an empty list and set a deterministic
code. Adapter and service failures remain exceptions.

.. autoclass:: lverage.pipeline.LverageCode
   :members:

.. autoclass:: lverage.pipeline.Lverage
   :members:
   :show-inheritance:
