"""Single logging entry point for eHT-PCR.

Replaces the four divergent ``logPrint`` helpers copy-pasted across the legacy
files with one stdlib ``logging`` configuration.
"""
from __future__ import annotations

import logging
import sys

_CONFIGURED = False


class _StderrHandler(logging.StreamHandler):
    """A handler bound to the *current* ``sys.stderr`` at emit time.

    A plain ``StreamHandler`` pins the stream at construction. typer's
    ``CliRunner`` swaps ``sys.stderr`` for a buffer it may later close, so a
    pinned handler can raise ``ValueError: I/O operation on closed file`` on a
    subsequent log. Resolving ``sys.stderr`` lazily (as the stdlib's own
    last-resort handler does) avoids that.
    """

    def __init__(self) -> None:
        logging.Handler.__init__(self)

    @property
    def stream(self):
        return sys.stderr

    @stream.setter
    def stream(self, value):   # StreamHandler.__init__/setStream try to set it; ignore
        pass


def configure(level: str | int = "INFO") -> None:
    """Configure the root ``ehtpcr`` logger once (idempotent)."""
    global _CONFIGURED
    logger = logging.getLogger("ehtpcr")
    if not _CONFIGURED:
        handler = _StderrHandler()
        handler.setFormatter(
            logging.Formatter("[%(asctime)s] %(levelname)s %(name)s - %(message)s",
                              datefmt="%Y-%m-%d %H:%M:%S")
        )
        logger.addHandler(handler)
        _CONFIGURED = True
    logger.setLevel(level)


def get_logger(name: str = "ehtpcr") -> logging.Logger:
    """Return a child logger under the ``ehtpcr`` namespace."""
    if name == "ehtpcr" or name.startswith("ehtpcr."):
        return logging.getLogger(name)
    return logging.getLogger(f"ehtpcr.{name}")
