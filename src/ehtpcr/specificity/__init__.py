"""Specificity engines: map primers to a reference and report amplicons.

Pluggable: ``base.SpecificityEngine`` is the ABC, ``bwa_aln.BwaAlnEngine`` is the
default (pure Python over real bwa). Future backends (BLAST, thermodynamic)
implement the same ABC and are selected by ``specificity.engine`` in the config.
"""
from __future__ import annotations


def get_engine(name: str):
    """Return a SpecificityEngine instance by config name."""
    if name == "bwa_aln":
        from .bwa_aln import BwaAlnEngine
        return BwaAlnEngine()
    raise ValueError(f"unknown specificity engine: {name!r}")
