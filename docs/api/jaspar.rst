JASPAR
======

The adapter uses the current ``jaspar.elixir.no`` REST API. Network requests
use finite timeouts, species availability is loaded lazily, and valid no-hit
responses return an empty list. Service, protocol, and malformed-response
failures are raised to the caller.

Inference uses a window of at most 2,000 ortholog residues while retaining the
complete ortholog domain. Accepted hits are sorted by inference E-value, checked
against every species listed in their motif details, and limited only after
species and threshold filtering.

.. autoclass:: lverage.jaspar.JasparRecord
   :members:
   :show-inheritance:

.. autoclass:: lverage.jaspar.Jaspar2024MotifDB
   :members:
   :show-inheritance:
