Records
========

.. autoclass:: lverage.records.LverageRecord
   :members:

``LverageRecord.get_headers()`` and ``get_values()`` flatten nested evidence
using the schema of that record's motif database. Group tabular output by motif
database when database-specific schemas differ.

The other version 2 interfaces define
:class:`lverage.domain_scanner.DomainRecord` and
:class:`lverage.motif_database.MotifSearchRequest` in the modules where they
are used.
