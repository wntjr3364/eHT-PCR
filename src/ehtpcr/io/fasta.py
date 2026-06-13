"""FASTA access via pyfaidx (replaces the hand-rolled byte-offset Fasta class).

The legacy ``Sequence.py`` reimplemented ``.fai`` random access by hand and
force-``.upper()``-ed everything (destroying soft-masking). pyfaidx is
battle-tested and preserves case. All coordinates here are 0-based half-open,
matching :class:`ehtpcr.records.Region` (pyfaidx ``[]`` slicing is 0-based).
"""
from __future__ import annotations

import fcntl
import hashlib
import os

from pyfaidx import Fasta

from ..cache import cache_root, ensure_writable
from ..records import Region
from ..sequtils import reverse_complement


def _fai_cache_path(path: str) -> str:
    d = ensure_writable(cache_root() / "fai")
    # key on content (size+mtime), like every other on-disk cache, so an in-place
    # edit re-keys and two distinct references at the same abspath don't collide.
    st = os.stat(path)
    key = hashlib.sha1(
        f"{os.path.abspath(path)}:{st.st_size}:{st.st_mtime_ns}".encode()
    ).hexdigest()[:16]
    return str(d / f"{key}.fai")


def _build_fasta_locked(path: str, indexname: str) -> Fasta:
    """Open a pyfaidx Fasta with the .fai redirected to *indexname*, serialized by a
    lock — pyfaidx writes the .fai non-atomically, so concurrent workers (jobs>1 over
    a read-only reference) could otherwise tear it."""
    with open(indexname + ".lock", "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        return Fasta(path, indexname=indexname)


class FastaRef:
    """Thin wrapper over ``pyfaidx.Fasta`` providing 0-based half-open fetch."""

    def __init__(self, path: str) -> None:
        self.path = path
        # pyfaidx writes a sibling .fai. If the reference directory is read-only
        # (shared db/ of reference genomes) and no .fai exists yet, redirect the
        # index to a writable per-user cache instead of crashing.
        indexname = None
        fai = os.path.abspath(path) + ".fai"
        ref_dir = os.path.dirname(os.path.abspath(path)) or "."
        if not os.path.exists(fai) and not os.access(ref_dir, os.W_OK):
            indexname = _fai_cache_path(path)
        try:
            self._fa = (Fasta(path) if indexname is None
                        else _build_fasta_locked(path, indexname))
        except ValueError as e:
            # two headers share a first whitespace token -> pyfaidx "Duplicate key"
            raise ValueError(
                f"{path}: cannot index FASTA ({e}). Header names (first whitespace "
                f"token) must be unique."
            ) from e
        except (OSError, PermissionError):
            # .fai write failed (read-only dir we couldn't pre-detect) -> cache it
            self._fa = _build_fasta_locked(path, _fai_cache_path(path))

    def names(self) -> list[str]:
        return list(self._fa.keys())

    def find_names(self, substring: str) -> list[str]:
        """Header names containing *substring* (legacy ``locus --contain``), sorted."""
        return sorted(n for n in self._fa.keys() if substring in n)

    def length(self, seqid: str) -> int:
        return len(self._fa[seqid])

    def fetch(self, region: Region) -> str:
        """Sequence for *region* (0-based half-open), case preserved; revcomp on '-'.

        Raises ``KeyError`` for an unknown seqid and ``ValueError`` if the region
        runs past the contig end (pyfaidx would silently truncate).
        """
        if region.seqid not in self._fa:
            raise KeyError(f"seqid {region.seqid!r} not in {self.path}")
        clen = len(self._fa[region.seqid])
        if region.end > clen:
            raise ValueError(
                f"region {region.seqid}:{region.start + 1}-{region.end} runs past "
                f"contig end ({clen} bp)"
            )
        seq = self._fa[region.seqid][region.start:region.end].seq
        if region.strand == "-":
            seq = reverse_complement(seq)
        return seq

    def fetch_whole(self, seqid: str) -> str:
        """Whole sequence of *seqid* (case preserved)."""
        return self._fa[seqid][:].seq
