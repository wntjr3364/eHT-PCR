"""Target mode: name — substring match against FASTA headers (legacy --contain).

The user's real legacy workflow (``Sequence.py locus --contain``): pull each whole
FASTA entry whose header contains the given substring and use it as a design
template. target_id = the matched header (sanitized), de-duplicated deterministically.
"""
from __future__ import annotations

from . import ResolvedTarget
from ..records import Region


def _sanitize(name: str) -> str:
    return "_".join(name.split())


def _override_for(target_id: str, substring: str, overrides: dict[str, int]) -> int | None:
    # exact target_id wins; fall back to the --contain substring (legacy keyed -m
    # on the substring, one file per substring).
    if target_id in overrides:
        return overrides[target_id]
    if substring in overrides:
        return overrides[substring]
    return None


def resolve(config, design_ref) -> list[ResolvedTarget]:
    overrides = config.specificity.params.bwa_aln.overrides
    seen: dict[str, ResolvedTarget] = {}
    for substring in config.target.name:
        for header in design_ref.find_names(substring):   # sorted, deterministic
            target_id = _sanitize(header)
            if target_id in seen:
                continue
            seq = design_ref.fetch_whole(header)
            seen[target_id] = ResolvedTarget(
                target_id=target_id,
                sequence=seq,
                region=Region(header, 0, len(seq), "+"),
                source_mode="name",
                mismatch_override=_override_for(target_id, substring, overrides),
            )
    return list(seen.values())
