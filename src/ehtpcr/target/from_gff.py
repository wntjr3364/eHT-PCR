"""Target mode: locus — gene id resolved via GFF, extracted from the reference.

Needs ``reference.gff``. Looks up each locus id's feature coordinates (gffutils),
applies strand-aware flank, and extracts the template from the design reference.
"""
from __future__ import annotations

from . import ResolvedTarget
from ..io.gff import GffRef
from ..records import Region


def _apply_flank(region: Region, flank, design_ref) -> Region:
    if flank.left == 0 and flank.right == 0:
        return region
    clen = design_ref.length(region.seqid)
    if region.strand == "+":
        s, e = region.start - flank.left, region.end + flank.right
    else:
        s, e = region.start - flank.right, region.end + flank.left
    return Region(region.seqid, max(0, s), min(clen, e), region.strand)


def resolve(config, design_ref) -> list[ResolvedTarget]:
    if not config.reference.gff:
        raise ValueError("target.mode == 'locus' requires reference.gff")
    if not config.target.locus:
        raise ValueError("target.mode == 'locus' requires target.locus")
    gff = GffRef(config.reference.gff, db_dir=config.reference.gff_db_dir)
    overrides = config.specificity.params.bwa_aln.overrides
    out: list[ResolvedTarget] = []
    for locus in config.target.locus:
        region = gff.region_for(locus, feature=config.target.feature)
        clen = design_ref.length(region.seqid)   # raises on a GFF seqid missing from the FASTA
        if region.end > clen:
            raise ValueError(
                f"locus {locus!r} ({region.seqid}:{region.start + 1}-{region.end}) runs "
                f"past contig end ({clen} bp) — GFF/reference mismatch?")
        fetch_region = _apply_flank(region, config.target.flank, design_ref)
        seq = design_ref.fetch(fetch_region)
        # `region` stays the pre-flank gene span (the biological target) so on-target
        # overlap isn't satisfied by a product sitting entirely in the added flank.
        out.append(ResolvedTarget(
            target_id=locus, sequence=seq, region=region,
            source_mode="locus", mismatch_override=overrides.get(locus),
        ))
    return out
