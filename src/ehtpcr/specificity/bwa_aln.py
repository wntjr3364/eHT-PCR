"""Default specificity engine: pure Python over real bwa (replaces the jar).

Reproduces the legacy jar's thin layer (it does NOT reimplement the aligner):

    bwa index <ref>                                       # prepare(), cached
    bwa aln -t T -N -i 30 -d 30 -m 200000000 -l 100 -n M  # primer.fastq -> .sai
    bwa samse -n 300000000 <ref> .sai primer.fastq        # -> SAM

then parses SAM + XA:Z multi-hit tags into per-oligo hits, and pairs
forward/reverse hits into amplicons. Every bwa flag is a tested
invariant (docs/legacy_bwa_semantics.md): -N required, -l 100 disables seeding,
-n integer, samse -n huge (XA is DROPPED, not truncated, on overflow), XA pos is
1-based signed (sign = strand); primary-hit strand = SAM FLAG 0x10.
"""
from __future__ import annotations

import fcntl
import hashlib
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from ..cache import cache_root, ensure_writable
from ..io.fasta import FastaRef
from ..logging import get_logger
from ..records import Amplicon, PrimerPair, Region
from .base import PreparedReference, SpecificityEngine
from .params import BwaAlnParams

log = get_logger("bwa_aln")

_QUAL = "I" * 48                      # jar's fixed 48-char Phred-Q40 quality string
_BWA_IDX_EXT = ("amb", "ann", "bwt", "pac", "sa")


@dataclass(frozen=True)
class BwaHit:
    """A single primer-oligo alignment (0-based leftmost)."""
    seqid: str
    strand: str          # '+' or '-'
    start0: int          # 0-based leftmost position on the reference
    nm: int              # edit distance (mismatches)
    source: str          # 'primary' | 'XA'


def _bwa_version(bwa: str = "bwa") -> str:
    try:
        out = subprocess.run([bwa], capture_output=True, text=True).stderr
    except FileNotFoundError as e:
        raise RuntimeError("bwa not found on PATH") from e
    for line in out.splitlines():
        if line.startswith("Version:"):
            return line.split(":", 1)[1].strip()
    return "unknown"


def _cache_root(index_dir: str | None) -> Path:
    # An explicit index_dir (specificity.reference.index_dir) is used verbatim;
    # otherwise indexes go under the shared cache root (see ehtpcr.cache).
    return Path(index_dir) if index_dir else cache_root() / "index"


def _cache_key(fasta: str, bwa_version: str) -> str:
    st = os.stat(fasta)
    # nanosecond mtime so a same-size edit within one second still re-keys.
    raw = f"{os.path.abspath(fasta)}:{st.st_size}:{st.st_mtime_ns}:{bwa_version}"
    return hashlib.sha1(raw.encode()).hexdigest()[:16]


def _tag(fields: list[str], tag: str) -> str | None:
    prefix = tag + ":"
    for t in fields:
        if t.startswith(prefix):
            return t.split(":", 2)[2]
    return None


def _parse_sam(sam_text: str, oligos: list[str]) -> dict[str, list[BwaHit]]:
    """Parse samse SAM into per-oligo hit lists (primary + XA), deduplicated.

    Read names are the oligo index (jar convention), so hits are keyed back to the
    exact oligo regardless of the reverse-complement bwa may have applied.
    """
    hits: dict[str, list[BwaHit]] = {o: [] for o in oligos}
    for line in sam_text.splitlines():
        if not line or line[0] == "@":
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        flag = int(f[1])
        if flag & 0x4:                       # unmapped read
            continue
        try:
            oligo = oligos[int(f[0])]
        except (ValueError, IndexError):
            continue
        pos = int(f[3])
        if pos < 1:                          # defensive: a mapped record needs POS>=1
            continue
        seen: set[tuple] = set()
        out = hits[oligo]

        # primary hit: strand from FLAG 0x10, POS is 1-based
        strand = "-" if (flag & 0x10) else "+"
        start0 = pos - 1
        nm = _tag(f[11:], "NM")
        key = (f[2], strand, start0)
        seen.add(key)
        out.append(BwaHit(f[2], strand, start0, int(nm) if nm and nm.isdigit() else 0, "primary"))

        # XA:Z: alternative hits — "rname,(+/-)pos,CIGAR,NM;..."
        xa = _tag(f[11:], "XA")
        if xa:
            for rec in xa.split(";"):
                if not rec:
                    continue
                p = rec.split(",")
                if len(p) < 2:
                    continue
                try:
                    signed = int(p[1])
                except ValueError:
                    continue                 # malformed XA position: skip this hit,
                                             # not the whole parse (cf. read-name above)
                xstrand = "-" if signed < 0 else "+"
                xstart0 = abs(signed) - 1
                k = (p[0], xstrand, xstart0)
                if k in seen:
                    continue
                seen.add(k)
                xnm = int(p[3]) if len(p) >= 4 and p[3].isdigit() else 0
                out.append(BwaHit(p[0], xstrand, xstart0, xnm, "XA"))
    return hits


def _drop_n_rich(hits: dict[str, list[BwaHit]], ref: PreparedReference,
                 threshold: int) -> dict[str, list[BwaHit]]:
    """Drop hits whose reference window has >= ``threshold`` N's (legacy D6)."""
    fa = getattr(ref, "_fasta_ref", None)   # reuse one FastaRef across targets
    if fa is None:
        fa = FastaRef(ref.fasta)
        ref._fasta_ref = fa
    out: dict[str, list[BwaHit]] = {}
    for oligo, hit_list in hits.items():
        length = len(oligo)
        kept = []
        for h in hit_list:
            end = min(h.start0 + length, fa.length(h.seqid))
            window = fa.fetch(Region(h.seqid, h.start0, end, "+"))
            if window.upper().count("N") >= threshold:
                continue
            kept.append(h)
        out[oligo] = kept
    return out


class BwaAlnEngine(SpecificityEngine):
    name = "bwa_aln"

    def prepare(self, fasta: str, *, index_dir: str | None = None) -> PreparedReference:
        """Build/locate the bwa index in a per-user cache (file-locked).

        The index is keyed by ``hash(fasta) + bwa_version`` and never written next
        to a possibly read-only reference. The original FASTA is symlinked into the
        cache so the ``.amb/.ann/.bwt/.pac/.sa`` files land in the cache.
        """
        version = _bwa_version()
        cache = ensure_writable(_cache_root(index_dir) / _cache_key(fasta, version))
        db = cache / Path(fasta).name
        abs_fasta = os.path.abspath(fasta)
        with open(cache / ".lock", "w") as lockf:
            fcntl.flock(lockf, fcntl.LOCK_EX)
            try:
                # (re)create the symlink if missing, broken, or pointing elsewhere
                if os.path.lexists(db):
                    if not (db.is_symlink() and os.path.realpath(db) == os.path.realpath(abs_fasta)):
                        db.unlink()
                        os.symlink(abs_fasta, db)
                else:
                    os.symlink(abs_fasta, db)
                # `.done` is written only after a fully successful build, so a
                # killed/interrupted `bwa index` (leaving truncated sidecars that
                # still "exist") is rebuilt instead of silently reused (mirrors the
                # blast/miniprot index sentinels in target/from_fasta.py).
                done = cache / f"{db.name}.done"
                if not (done.exists()
                        and all((cache / f"{db.name}.{e}").exists() for e in _BWA_IDX_EXT)):
                    log.info("building bwa index for %s in %s", fasta, cache)
                    r = subprocess.run(["bwa", "index", str(db)], capture_output=True, text=True)
                    if r.returncode != 0:
                        raise RuntimeError(f"bwa index failed: {r.stderr[-500:]}")
                    done.touch()
            finally:
                fcntl.flock(lockf, fcntl.LOCK_UN)
        return PreparedReference(str(db), bwa_version=version, cache_dir=str(cache))

    def find_hits(self, oligos, ref: PreparedReference, params: BwaAlnParams, *,
                  threads: int = 1, tmp_root: str | None = None,
                  keep_temp: bool = False) -> dict[str, list[BwaHit]]:
        """Run bwa aln/samse on the unique oligo set; return per-oligo hits."""
        oligos = list(dict.fromkeys(oligos))      # dedup, preserve order
        for o in oligos:
            if len(o) > len(_QUAL):
                raise ValueError(f"oligo > {len(_QUAL)} bp not supported by the legacy fastq: {o!r}")
        tmp = Path(tempfile.mkdtemp(prefix="ehtpcr_bwa_", dir=tmp_root))
        try:
            fastq = tmp / "primer.fastq"
            with open(fastq, "w") as fh:
                for i, o in enumerate(oligos):
                    fh.write(f"@{i}\n{o}\n+\n{_QUAL[:len(o)]}\n")
            sai, sam = tmp / "primer.sai", tmp / "primer.sam"
            with open(sai, "wb") as so, open(tmp / "aln.log", "wb") as se:
                r = subprocess.run(
                    ["bwa", "aln", "-t", str(threads), "-N", "-i", "30", "-d", "30",
                     "-m", "200000000", "-l", "100", "-n", str(params.max_mismatch),
                     ref.fasta, str(fastq)],
                    stdout=so, stderr=se,
                )
            if r.returncode != 0:
                raise RuntimeError("bwa aln failed: " + (tmp / "aln.log").read_text()[-500:])
            with open(sam, "wb") as so, open(tmp / "samse.log", "wb") as se:
                r = subprocess.run(
                    ["bwa", "samse", "-n", "300000000", ref.fasta, str(sai), str(fastq)],
                    stdout=so, stderr=se,
                )
            if r.returncode != 0:
                raise RuntimeError("bwa samse failed: " + (tmp / "samse.log").read_text()[-500:])
            hits = _parse_sam(sam.read_text(), oligos)
        finally:
            if keep_temp:
                log.info("kept bwa temp files in %s", tmp)
            else:
                shutil.rmtree(tmp, ignore_errors=True)
        if params.drop_n_rich_windows:
            hits = _drop_n_rich(hits, ref, params.n_window_threshold)
        return hits

    def find_amplicons(self, pairs: list[PrimerPair], ref: PreparedReference,
                       params: BwaAlnParams, *, threads: int = 1,
                       tmp_root: str | None = None, keep_temp: bool = False) -> list[Amplicon]:
        """Map every oligo once, then pair F/R hits into FR/RF amplicons.

        Per docs/legacy_bwa_semantics.md: ``-`` hits are normalized to their
        rightmost base; ``productSize = revPos - fwdPos + 1`` (here in 0-based
        half-open: amplicon = [leftmost+0, rightmost+0)); kept iff
        ``min_size <= productSize <= max_size`` with ``min_size = len_f + len_r``
        (the pair's two primer lengths) and ``max_size = params.max_product_size``.
        Unlike the jar, amplicons are de-duplicated (divergence D4).
        """
        oligos: list[str] = []
        for p in pairs:
            oligos.append(p.forward.seq)
            oligos.append(p.reverse.seq)
        hit_map = self.find_hits(oligos, ref, params, threads=threads,
                                 tmp_root=tmp_root, keep_temp=keep_temp)
        return self.pair_hits(pairs, hit_map, params)

    def pair_hits(self, pairs: list[PrimerPair], hit_map: dict[str, list[BwaHit]],
                  params: BwaAlnParams) -> list[Amplicon]:
        """Pair forward/reverse hits from a ``hit_map`` into FR/RF amplicons.

        The pairing half of :meth:`find_amplicons` (= ``find_hits`` then this), kept
        separate so hit-finding and pairing can be reasoned about independently.
        """
        max_size = params.max_product_size
        cap = params.max_amplicons_per_pair
        amplicons: list[Amplicon] = []
        for p in pairs:
            f_hits = hit_map.get(p.forward.seq, [])
            r_hits = hit_map.get(p.reverse.seq, [])
            len_f, len_r = p.forward.length, p.reverse.length
            min_size = len_f + len_r
            f_plus = [h for h in f_hits if h.strand == "+"]
            f_minus = [h for h in f_hits if h.strand == "-"]
            r_plus = [h for h in r_hits if h.strand == "+"]
            r_minus = [h for h in r_hits if h.strand == "-"]
            seen: set[tuple] = set()
            n_for_pair = 0
            overflow = False

            # FR: forward '+' (leftmost) x reverse '-' (rightmost)
            for fh in f_plus:
                if overflow:
                    break
                for rh in r_minus:
                    if fh.seqid != rh.seqid:
                        continue
                    start0, end0 = fh.start0, rh.start0 + len_r
                    ps = end0 - start0
                    if ps < min_size or ps > max_size:
                        continue
                    # dedup by span across FR/RF: a single physical product (only
                    # possible when forward.seq == reverse.seq) is not double-counted.
                    key = (fh.seqid, start0, end0)
                    if key in seen:
                        continue
                    seen.add(key)
                    amplicons.append(Amplicon(
                        p.pair_id, fh.seqid, start0, end0, "FR", ps,
                        fh.nm, rh.nm, fh.source, rh.source))
                    n_for_pair += 1
                    if cap and n_for_pair >= cap:
                        overflow = True
                        break

            # RF: reverse '+' (leftmost) x forward '-' (rightmost)
            for rh in r_plus:
                if overflow:
                    break
                for fh in f_minus:
                    if fh.seqid != rh.seqid:
                        continue
                    start0, end0 = rh.start0, fh.start0 + len_f
                    ps = end0 - start0
                    if ps < min_size or ps > max_size:
                        continue
                    key = (rh.seqid, start0, end0)
                    if key in seen:
                        continue
                    seen.add(key)
                    amplicons.append(Amplicon(
                        p.pair_id, rh.seqid, start0, end0, "RF", ps,
                        fh.nm, rh.nm, fh.source, rh.source))
                    n_for_pair += 1
                    if cap and n_for_pair >= cap:
                        overflow = True
                        break

            if overflow:
                log.warning("pair %s hit the amplicon cap (%d); non-specific, capping",
                            p.pair_id, cap)
        return amplicons
