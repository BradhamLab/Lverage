JASPAR
======

The adapter uses the current ``jaspar.elixir.no`` REST API. Network requests
use finite timeouts, species availability is loaded lazily, and valid no-hit
responses return an empty list. Service, protocol, and malformed-response
failures are raised to the caller.

.. autoclass:: lverage.jaspar.JasparRecord
   :members:
   :show-inheritance:

.. autoclass:: lverage.jaspar.Jaspar2024MotifDB
   :members:
   :show-inheritance:
