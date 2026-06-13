"""Primer candidate enumeration and pairing (port of legacy Primer.py).

Ports the validated logic from the original legacy ``Primer.py``:
- enumerate sub-sequences of length [min_len, max_len] across the template
  (replicating the legacy ``__main__`` boundary: scan starts at index 0 and the
  last base is excluded, ``end < len``; the legacy ``boundaryList=None`` default
  path was a dead/crashing branch and is not a fidelity target);
- filter by GC and Tm (legacy Breslauer Tm by default);
- keep only primers whose sequence is **unique within the template**;
- pair forward/reverse with a 3' GC-clamp (last ``gc_clamp`` bases all G/C),
  product size in [min, max], and tm_diff / gc_diff limits.

Fixes carried in: the single correct reverse-complement (sequtils), reject
degenerate bases by default, case-insensitive GC/Tm/clamp, a stable ``pair_id``,
and a hard ``max_pairs`` cap to guard the O(n^2) pairing explosion.
"""
from __future__ import annotations

from .logging import get_logger
from .records import Primer, PrimerPair
from .sequtils import gc_content, is_canonical, reverse_complement
from .thermo import tm as tm_of

log = get_logger("design")


def enumerate_primers(template: str, design) -> list[Primer]:
    """Candidate primers on *template* (0-based half-open), Tm/GC filtered.

    The template is upper-cased first, matching the legacy ``Primer.py`` which
    read FASTA with ``.upper()``. This keeps the unique-within-template dedup key
    and the emitted oligo sequences identical to the legacy (and to what bwa is
    fed), i.e. ``softmask_policy: allow`` == legacy case-folding.
    (``reject`` / ``lowercase_penalty`` are not implemented — see config.)
    """
    template = template.upper()
    min_len, max_len = design.primer_len.min, design.primer_len.max
    gc_min, gc_max = design.gc.min, design.gc.max
    tm_min, tm_max = design.tm.min, design.tm.max
    model = design.tm.model
    n = len(template)

    counts: dict[str, int] = {}
    candidates: list[Primer] = []
    for start in range(0, n - min_len + 1):
        for end in range(start + min_len, start + max_len + 1):
            if end >= n:           # legacy boundary: last base excluded
                break
            seq = template[start:end]
            if seq in counts:
                counts[seq] += 1
                continue
            counts[seq] = 0
            if not design.allow_degenerate and not is_canonical(seq):
                continue
            gc = gc_content(seq)
            t = tm_of(seq, model)
            if gc_min <= gc <= gc_max and tm_min <= t <= tm_max:
                candidates.append(Primer(seq=seq, start=start, end=end, tm=t, gc=gc))
    # keep only primers unique within the template (legacy countDict == 0)
    return [p for p in candidates if counts[p.seq] == 0]


def _three_prime_has_at(oligo: str, gc_clamp: int) -> bool:
    if gc_clamp <= 0:
        return False
    tail = oligo[-gc_clamp:].upper()
    return ("A" in tail) or ("T" in tail)


def pair_primers(target_id: str, primers: list[Primer], design, *,
                 max_pairs: int, fail_on_overflow: bool) -> list[PrimerPair]:
    """Pair forward+reverse primers under product-size / tm_diff / gc_diff / clamp."""
    gc_clamp = design.gc_clamp
    p_min, p_max = design.product_size.min, design.product_size.max
    tmd_min, tmd_max = design.pair.tm_diff.min, design.pair.tm_diff.max
    gcd_min, gcd_max = design.pair.gc_diff.min, design.pair.gc_diff.max

    n = len(primers)
    pairs: list[PrimerPair] = []
    for fi in range(n):
        f = primers[fi]
        if _three_prime_has_at(f.seq, gc_clamp):       # forward 3' GC-clamp
            continue
        for ri in range(fi + 1, n):
            r = primers[ri]
            product = r.end - f.start                   # 0-based half-open span
            if product < p_min or product > p_max:
                continue
            tm_diff = abs(f.tm - r.tm)
            gc_diff = abs(f.gc - r.gc)
            if tm_diff < tmd_min or tm_diff > tmd_max:
                continue
            if gc_diff < gcd_min or gc_diff > gcd_max:
                continue
            r_oligo = reverse_complement(r.seq)
            if _three_prime_has_at(r_oligo, gc_clamp):  # reverse oligo 3' GC-clamp
                continue
            reverse = Primer(seq=r_oligo, start=r.start, end=r.end, tm=r.tm, gc=r.gc)
            # stable + collision-free: both primer footprints (a target can't have
            # two distinct primers sharing the same start AND end).
            pair_id = f"{target_id}:{f.start}-{f.end}:{r.start}-{r.end}"
            pairs.append(PrimerPair(
                pair_id=pair_id, target_id=target_id, forward=f, reverse=reverse,
                product_size=product, tm_diff=tm_diff, gc_diff=gc_diff,
            ))
            if len(pairs) > max_pairs:
                if fail_on_overflow:
                    raise RuntimeError(
                        f"target {target_id}: more than {max_pairs} primer pairs; "
                        "raise runtime.max_pairs_per_target or tighten design filters"
                    )
                log.warning("target %s: capping primer pairs at %d", target_id, max_pairs)
                return pairs[:max_pairs]
    return pairs
