"""Sequence utilities: IUPAC reverse-complement, GC content, base validation.

Single source of truth. The legacy code had **two** divergent reverse-complement
implementations: ``Sequence.py`` had a full IUPAC table with the bug ``M -> "L"``
(M's complement is K), and ``Primer.py`` had an ACGTUN-only table that raised
``KeyError`` on any ambiguity code. This module fixes both.

Rules:
- Reverse-complement **preserves case** (soft-masking survives revcomp).
- GC content is computed **case-insensitively** over canonical bases A/C/G/T.
- Ambiguous / IUPAC degenerate bases are rejected from primer candidates by
  default (callers check :func:`is_canonical`).
"""
from __future__ import annotations

# Full IUPAC complement (uppercase). U complements to A.
_IUPAC_UPPER = {
    "A": "T", "T": "A", "U": "A", "C": "G", "G": "C",
    "W": "W", "S": "S", "R": "Y", "Y": "R", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}
# Case-preserving table (both cases).
_COMPLEMENT = {}
for _k, _v in _IUPAC_UPPER.items():
    _COMPLEMENT[_k] = _v
    _COMPLEMENT[_k.lower()] = _v.lower()

_TRANSLATION = str.maketrans(_COMPLEMENT)

_CANONICAL = frozenset("ACGTacgt")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of *seq*, preserving case.

    Raises ``ValueError`` on any character outside the IUPAC alphabet (gaps,
    whitespace, ``*`` etc. are rejected, not silently passed). Note: this is a
    DNA complement, so RNA ``U`` -> ``A``; a U-containing input is therefore NOT
    round-trip stable (``rc(rc("U")) == "T"``).
    """
    bad = set(seq) - _COMPLEMENT.keys()
    if bad:
        raise ValueError(f"non-IUPAC base(s) in sequence: {sorted(bad)!r}")
    return seq.translate(_TRANSLATION)[::-1]


def is_canonical(seq: str) -> bool:
    """True iff *seq* contains only A/C/G/T (case-insensitive), non-empty."""
    return len(seq) > 0 and all(b in _CANONICAL for b in seq)


def gc_content(seq: str) -> float:
    """GC fraction over canonical bases A/C/G/T (case-insensitive).

    For all-canonical sequences this equals the legacy ``(G+C)/len``. The result
    is **undefined/misleading for sequences containing non-canonical bases**
    (e.g. ``gc_content("GCNN") == 1.0`` because N is excluded from the
    denominator), so callers must restrict primers to canonical input via
    :func:`is_canonical` (``design.allow_degenerate: false``). Returns 0.0 if
    there are no canonical bases (which fails the GC-min filter, as intended).
    """
    s = seq.upper()
    g = s.count("G")
    c = s.count("C")
    a = s.count("A")
    t = s.count("T")
    canonical = a + t + g + c
    if canonical == 0:
        return 0.0
    return (g + c) / canonical
