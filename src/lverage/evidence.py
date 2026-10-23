"""Small event collector for incremental pipeline evidence."""


class EvidenceCollector:
    """Forward completed pipeline evidence to a consumer as it is produced.

    The collector is intentionally in-memory and per invocation. It is not a
    checkpoint or crash-recovery mechanism.
    """

    def __init__(self, callback=None):
        self.callback = callback
        self.events = []
        self.stage = "preflight"
        self.context = {}

    def set_stage(self, stage):
        self.stage = stage

    def set_context(self, **context):
        """Attach identifying context to subsequently reported events."""
        self.context.update(context)

    def emit(self, event, **payload):
        item = {"event": event, "stage": self.stage, **self.context, **payload}
        self.events.append(item)
        if self.callback is not None:
            self.callback(item)
        return item
