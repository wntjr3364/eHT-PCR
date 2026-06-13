"""SpecificityEngine ABC — must not leak bwa-specific assumptions.

Two responsibilities, split so a future thermodynamic/BLAST backend fits:
- ``prepare(ref)``        : build/locate whatever index this engine needs
                            (idempotent, file-locked, cached); returns a handle.
- ``find_amplicons(...)`` : given primer pairs + a prepared reference + engine
                            params, return all amplicons. The engine OWNS hit
                            parsing, F/R pairing and FR/RF orientation. The
                            pipeline applies the uniqueness policy on the result.

Engine-specific tuning is passed via a typed params object (see params.py), NOT
as loose scalars like ``max_mismatch`` (a bwa concept; a thermo engine takes a
ΔG/Tm threshold, BLAST takes e-value/word-size).
"""
from __future__ import annotations

from abc import ABC, abstractmethod

from ..records import Amplicon, PrimerPair


class PreparedReference:
    """Opaque, engine-specific handle to a prepared (indexed) reference."""

    def __init__(self, fasta: str, **meta) -> None:
        self.fasta = fasta
        self.meta = meta


class SpecificityEngine(ABC):
    name: str = "base"

    @abstractmethod
    def prepare(self, fasta: str, *, index_dir: str | None = None) -> PreparedReference:
        ...

    @abstractmethod
    def find_amplicons(
        self,
        pairs: list[PrimerPair],
        ref: PreparedReference,
        params,
    ) -> list[Amplicon]:
        """Map *pairs* against *ref* and return all amplicons.

        All tuning (incl. the max product size — a bwa/PCR concept, not a generic
        one) lives in the engine-specific ``params`` object, so the signature
        stays backend-agnostic (a thermo/BLAST backend takes different params).
        """
        ...
