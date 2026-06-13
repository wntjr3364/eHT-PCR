"""Target resolution: turn the user's target spec into design templates.

One ``resolve_targets(config, design_ref) -> list[ResolvedTarget]`` entry point
dispatches on ``target.mode`` to a strategy module:

- ``from_name``   : substring match against FASTA headers (legacy ``--contain``)
- ``from_region`` : explicit coordinates ``seqid:start-end[:strand]``
- ``from_gff``    : a locus id resolved via a GFF (mode ``locus``)
- ``from_fasta``  : locate a query sequence by alignment (mode ``target_fasta``)

A ResolvedTarget carries a stable ``target_id``, the design ``Region`` (where it
came from in the reference, when known), and the extracted template sequence.
"""
from __future__ import annotations

from dataclasses import dataclass

from ..records import Region


@dataclass(frozen=True)
class ResolvedTarget:
    # target_id is the STABLE key for outputs and per-target overrides. Derivation:
    #   - name mode  : the matched FASTA header name (sanitized; whitespace->_).
    #                  Multiple headers match a substring -> one target each;
    #                  duplicates are de-duplicated deterministically.
    #   - region mode: "seqid:start1-end1:strand".
    #   - locus mode : the locus id.
    # It is NOT the --contain substring (a substring can match many headers).
    target_id: str
    sequence: str               # design template (case preserved)
    region: Region | None       # location in the reference, when known
    source_mode: str            # "name" | "region" | "locus" | "target_fasta" (internal only)
    mismatch_override: int | None = None   # per-target bwa -n override; None => engine default


def resolve_targets(config, design_ref) -> list[ResolvedTarget]:
    """Dispatch on ``config.target.mode`` to the appropriate strategy."""
    mode = config.target.mode
    if mode == "name":
        from . import from_name
        return from_name.resolve(config, design_ref)
    if mode == "region":
        from . import from_region
        return from_region.resolve(config, design_ref)
    if mode == "locus":
        from . import from_gff
        return from_gff.resolve(config, design_ref)
    if mode == "target_fasta":
        from . import from_fasta
        return from_fasta.resolve(config, design_ref)
    raise ValueError(f"unknown target mode {mode!r}")
