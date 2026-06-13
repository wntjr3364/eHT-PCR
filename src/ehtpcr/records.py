"""Core domain dataclasses for eHT-PCR.

All genomic/template coordinates are **0-based, half-open** (``seq[start:end]``;
``length == end - start``). See docs/legacy_bwa_semantics.md.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

Strand = Literal["+", "-"]
StrandClass = Literal["FR", "RF"]
HitSource = Literal["primary", "XA"]


@dataclass(frozen=True)
class Region:
    """A 0-based, half-open interval on a named sequence."""

    seqid: str
    start: int  # 0-based inclusive
    end: int    # 0-based exclusive
    strand: Strand = "+"

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"Region.start must be >= 0, got {self.start}")
        if self.end < self.start:
            raise ValueError(f"Region.end ({self.end}) < start ({self.start})")
        if self.strand not in ("+", "-"):
            raise ValueError(f"Region.strand must be '+'/'-', got {self.strand!r}")

    @property
    def length(self) -> int:
        return self.end - self.start

    def overlap_bp(self, other: "Region") -> int:
        """Overlap length in bp with *other* on the same seqid (0 if none/diff seqid)."""
        if self.seqid != other.seqid:
            return 0
        lo = max(self.start, other.start)
        hi = min(self.end, other.end)
        return max(0, hi - lo)

    def overlaps(self, other: "Region", min_bp: int = 1) -> bool:
        return self.overlap_bp(other) >= min_bp


@dataclass(frozen=True)
class Primer:
    """A primer candidate, located on the design template (0-based half-open)."""

    seq: str       # as it anneals (forward primer = template sense; reverse = revcomp)
    start: int     # 0-based on the design template
    end: int       # 0-based half-open on the design template
    tm: float
    gc: float

    def __post_init__(self) -> None:
        if not self.seq:
            raise ValueError("Primer.seq must be non-empty")
        if self.start < 0:
            raise ValueError(f"Primer.start must be >= 0, got {self.start}")
        if self.end < self.start:
            raise ValueError(f"Primer.end ({self.end}) < start ({self.start})")

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class PrimerPair:
    """A forward+reverse primer pair designed for one target."""

    pair_id: str
    target_id: str
    forward: Primer
    reverse: Primer        # reverse.seq is the reverse-complement (legacy ``r_seq``)
    product_size: int      # on the design template
    tm_diff: float
    gc_diff: float


@dataclass(frozen=True)
class Amplicon:
    """A predicted PCR product on the specificity reference (0-based half-open)."""

    pair_id: str
    seqid: str
    start: int
    end: int
    strand_class: StrandClass   # FR or RF
    product_size: int
    f_mismatches: int
    r_mismatches: int
    f_hit_source: HitSource
    r_hit_source: HitSource

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"Amplicon.start must be >= 0, got {self.start}")
        if self.end < self.start:
            raise ValueError(f"Amplicon.end ({self.end}) < start ({self.start})")
        if self.strand_class not in ("FR", "RF"):
            raise ValueError(f"Amplicon.strand_class must be 'FR'/'RF', got {self.strand_class!r}")
        if self.f_hit_source not in ("primary", "XA") or self.r_hit_source not in ("primary", "XA"):
            raise ValueError("Amplicon hit_source must be 'primary' or 'XA'")
        if self.product_size < 0 or self.f_mismatches < 0 or self.r_mismatches < 0:
            raise ValueError("Amplicon product_size/mismatches must be >= 0")

    def as_region(self) -> Region:
        return Region(self.seqid, self.start, self.end, "+")
