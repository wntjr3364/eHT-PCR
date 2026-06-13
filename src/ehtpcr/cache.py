"""Shared resolver for eHT-PCR's on-disk index/DB caches.

Every cached artifact — bwa indexes, BLAST databases, miniprot indexes, the
pyfaidx FASTA ``.fai``, and the gffutils GFF SQLite DB — lives under one root, so
a single writable location controls them all. Resolution order (first that is
set wins):

    explicit arg  >  $EHTPCR_CACHE_DIR  >  $XDG_CACHE_HOME/ehtpcr  >  ~/.cache/ehtpcr

The ``--cache-dir`` flag sets ``$EHTPCR_CACHE_DIR`` for the whole run (parent and
spawned workers), so it applies uniformly to the specificity index *and* every
target_fasta locator cache.
"""
from __future__ import annotations

import os
from pathlib import Path


def cache_root(override: str | None = None) -> Path:
    """The eHT-PCR cache root (not created here; see :func:`ensure_writable`)."""
    if override:
        return Path(override)
    if env := os.environ.get("EHTPCR_CACHE_DIR"):
        return Path(env)
    if xdg := os.environ.get("XDG_CACHE_HOME"):
        return Path(xdg) / "ehtpcr"
    return Path(os.path.expanduser("~")) / ".cache" / "ehtpcr"


def ensure_writable(path: Path) -> Path:
    """``mkdir -p`` *path*, raising a clear, actionable error when the cache
    location is read-only.

    The default cache lives under ``~/.cache``; on containers / shared HPC nodes
    that is often read-only, and the bare failure is an opaque ``Errno 30`` from
    deep inside an index build. Surface it here with the fix instead.
    """
    try:
        path.mkdir(parents=True, exist_ok=True)
        writable = os.access(path, os.W_OK)
    except OSError as e:
        writable = False
        reason = f" ({e.strerror})"
    else:
        reason = ""
    if not writable:
        raise RuntimeError(
            f"cache directory is not writable: {path}{reason}. "
            "Point eHT-PCR at a writable location, e.g. "
            "`--cache-dir /tmp/ehtpcr-cache` or `export EHTPCR_CACHE_DIR=/tmp/ehtpcr-cache` "
            "(or set XDG_CACHE_HOME)."
        )
    return path
