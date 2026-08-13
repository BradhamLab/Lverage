BLAST searchers
===============

BLAST+ is an external prerequisite and is not installed as a Python package.
On the SCC, ``blast+/2.12.0`` and the ``human_mouse_nr`` database are useful
development examples; applications must still pass their own executable and
database locations.

.. autoclass:: lverage.blast.LocalBlastSearcher
   :members:
   :show-inheritance:
