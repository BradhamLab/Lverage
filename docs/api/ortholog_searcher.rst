Ortholog searcher
=================

Ortholog searchers isolate sequence-similarity services from pipeline
orchestration. Result identity and coverage fields are stored as ratios from
``0.0`` to ``1.0`` and should be converted to percentages only for display.

.. autoclass:: lverage.ortholog_searcher.OrthologRecord
   :members:

.. autoclass:: lverage.ortholog_searcher.OrthologSearcherTemplate
   :members:
   :show-inheritance:
