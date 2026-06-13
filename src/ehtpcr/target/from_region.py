"""Target mode: region — explicit coordinates ``seqid:start-end[:strand]``.

CLI/config coordinates are 1-based inclusive (matching legacy ``queryParser``) and
are converted to a 0-based half-open Region here. The seqid may itself contain
colons (some accessions do); the optional trailing ``:+``/``:-`` is stripped first
and the position is taken from the last remaining colon.
"""
from __future__ import annotations

from . import ResolvedTarget
from ..records import Region


def parse_region(text: str) -> Region:
    """Parse ``seqid:start-end[:strand]`` (1-based inclusive) -> 0-based Region."""
    strand = "+"
    body = text
    if body.endswith(":+") or body.endswith(":-"):
        strand = body[-1]
        body = body[:-2]
    if ":" not in body:
        raise ValueError(f"bad region {text!r}; expected seqid:start-end[:strand]")
    seqid, pos = body.rsplit(":", 1)
    if not seqid:
        raise ValueError(f"bad region {text!r}; empty seqid")
    try:
        s, e = pos.split("-")
        start1, end1 = int(s), int(e)
    except ValueError:
        raise ValueError(f"bad position in region {text!r}; expected start-end") from None
    if start1 < 1 or end1 < start1:
        raise ValueError(f"bad coordinates in region {text!r} (1-based, start<=end)")
    return Region(seqid, start1 - 1, end1, strand)


def _apply_flank(region: Region, flank, design_ref) -> Region:
    if flank.left == 0 and flank.right == 0:
        return region
    clen = design_ref.length(region.seqid)
    if region.strand == "+":
        s, e = region.start - flank.left, region.end + flank.right
    else:  # flank is relative to gene orientation (legacy queryParser)
        s, e = region.start - flank.right, region.end + flank.left
    return Region(region.seqid, max(0, s), min(clen, e), region.strand)


def resolve(config, design_ref) -> list[ResolvedTarget]:
    if not config.target.region:
        raise ValueError("target.mode == 'region' requires target.region")
    region = parse_region(config.target.region)
    # validate the REQUESTED region independently of flank, so an out-of-bounds
    # request fails loud whether or not a flank is set (a flank otherwise silently
    # clamps it). length() raises on an unknown seqid.
    clen = design_ref.length(region.seqid)
    if region.end > clen:
        raise ValueError(
            f"region {config.target.region} runs past contig end ({clen} bp)")
    # target_id is the requested region (stable override/output key); the design
    # template is fetched from the flank-expanded region, but `region` stays the
    # pre-flank biological target so on-target overlap is not satisfied by a product
    # that sits entirely inside the added flank.
    target_id = f"{region.seqid}:{region.start + 1}-{region.end}:{region.strand}"
    fetch_region = _apply_flank(region, config.target.flank, design_ref)
    seq = design_ref.fetch(fetch_region)
    override = config.specificity.params.bwa_aln.overrides.get(target_id)
    return [ResolvedTarget(
        target_id=target_id, sequence=seq, region=region,
        source_mode="region", mismatch_override=override,
    )]
