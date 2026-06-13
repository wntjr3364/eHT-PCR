"""Target mode: target_fasta — locate query sequences in the reference.

Given a FASTA of query sequences (genes of interest), align each to the design
reference, keep loci passing identity/coverage thresholds, and **fail loud if
more than one locus passes** (don't silently pick the wrong paralog/homeolog).

``target.locator`` selects the aligner; default ``auto`` sniffs the query FASTA
(nucleotide vs protein — see :func:`_query_is_protein`) and routes nucleotide to
``blastn`` and protein to ``tblastn``:
- ``bwa``      : DNA query, near-identical (no extra dep).
- ``blastn``   : DNA query, sensitive + paralog-safe; also handles a spliced
                 transcript→genome (exon HSPs grouped). Needs BLAST+.
- ``minimap2`` : spliced transcript→genome, intron-aware (mappy). Paralog-blind.
- ``tblastn``  : protein→genome, 6-frame translated; paralog-safe. Needs BLAST+.
- ``miniprot`` : protein→genome, splice/frameshift-aware. Paralog-blind.

Note ``auto`` only distinguishes nucleotide from protein (reliable by composition);
a *gene* vs a *transcript* query are both nucleotide and indistinguishable by
sequence — both route to ``blastn``, which locates either. The splice-precise
locators (``minimap2``/``miniprot``) are explicit opt-ins.
"""
from __future__ import annotations

import fcntl
import hashlib
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

from . import ResolvedTarget
from ..cache import cache_root, ensure_writable
from ..logging import get_logger
from ..records import Region
from ..specificity.bwa_aln import BwaAlnEngine, _tag

log = get_logger("target.from_fasta")
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

# a locus hit is a 6-tuple: (seqid, start0, end0, strand, identity, coverage)


def _read_fasta_lengths(path: str) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name, n = None, 0
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                lengths[name] = n
            name = line[1:].split()[0] if line[1:].strip() else ""
            if name in lengths:   # the prior record with this name is already flushed
                raise ValueError(
                    f"{path}: duplicate query name {name!r}. Query FASTA headers "
                    f"(first whitespace token) must be unique — otherwise the targets "
                    f"silently collapse to one (and coverage uses the wrong length)."
                )
            n = 0
        else:
            n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def _merge_overlapping(hits: list[tuple]) -> list[tuple]:
    """Merge alignments on the same seqid AND strand whose ref ranges overlap (one
    locus reported as several chains) so they aren't miscounted as distinct loci.
    Strand matters: two overlapping *antisense* loci are distinct and must stay
    separate, or the >1-locus fail-loud check is silently bypassed."""
    merged: list[tuple] = []
    for h in sorted(hits, key=lambda x: (x[0], x[3], x[1])):
        if merged and merged[-1][0] == h[0] and merged[-1][3] == h[3] and h[1] < merged[-1][2]:
            p = merged[-1]
            merged[-1] = (p[0], min(p[1], h[1]), max(p[2], h[2]), p[3],
                          max(p[4], h[4]), max(p[5], h[5]))
        else:
            merged.append(h)
    return merged


# --- bwa locator -----------------------------------------------------------

def _cigar_stats(cigar: str, nm: int) -> tuple[int, int, float]:
    """(query_aligned_bases, ref_span, identity) from a CIGAR + NM."""
    q_aln = ref_span = block = 0
    for num, op in _CIGAR.findall(cigar):
        n = int(num)
        if op in "MI=X":
            q_aln += n
        if op in "MDN=X":
            ref_span += n
        if op in "MID=X":
            block += n
    identity = (block - nm) / block if block else 0.0
    return q_aln, ref_span, identity


def _locate_bwa(ref_fasta, query_fasta, qlen, threads, min_id, min_cov) -> dict[str, list]:
    prep = BwaAlnEngine().prepare(ref_fasta)          # cached bwa index of design ref
    tmp = Path(tempfile.mkdtemp(prefix="ehtpcr_locate_"))
    try:
        sam = tmp / "q.sam"
        with open(sam, "wb") as out, open(tmp / "mem.log", "wb") as err:
            r = subprocess.run(
                ["bwa", "mem", "-a", "-t", str(threads), prep.fasta, query_fasta],
                stdout=out, stderr=err,
            )
        if r.returncode != 0:
            raise RuntimeError("bwa mem failed: " + (tmp / "mem.log").read_text()[-400:])

        raw: dict[str, dict] = {q: {} for q in qlen}

        def _add(q, seqid, start0, ref_span, strand, q_aln, identity):
            cov = q_aln / qlen[q] if qlen.get(q) else 0.0
            if identity >= min_id and cov >= min_cov:
                raw[q][(seqid, start0, strand)] = (seqid, start0, start0 + ref_span, strand, identity, cov)

        for line in open(sam):
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            q, flag = f[0], int(f[1])
            if flag & 0x4 or q not in qlen:
                continue
            nm = _tag(f[11:], "NM")
            q_aln, ref_span, identity = _cigar_stats(f[5], int(nm) if nm and nm.isdigit() else 0)
            _add(q, f[2], int(f[3]) - 1, ref_span, "-" if flag & 0x10 else "+", q_aln, identity)

            xa = _tag(f[11:], "XA")          # alternative hits: rname,(+/-)pos,CIGAR,NM;...
            if xa:
                for rec in xa.split(";"):
                    if not rec:
                        continue
                    p = rec.split(",")
                    if len(p) < 4:
                        continue
                    signed, xnm = int(p[1]), int(p[3]) if p[3].isdigit() else 0
                    xq, xspan, xident = _cigar_stats(p[2], xnm)
                    _add(q, p[0], abs(signed) - 1, xspan, "-" if signed < 0 else "+", xq, xident)

        return {q: _merge_overlapping(list(raw[q].values())) for q in qlen}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --- blastn locator --------------------------------------------------------

def _blast_db(ref_fasta: str) -> str:
    """Build (once, cached + file-locked) a nucleotide BLAST DB for *ref_fasta*."""
    st = os.stat(ref_fasta)
    key = hashlib.sha1(f"{os.path.abspath(ref_fasta)}:{st.st_size}:{st.st_mtime_ns}".encode()).hexdigest()[:16]
    d = str(ensure_writable(cache_root() / "blastdb" / key))
    db = os.path.join(d, "db")
    done = db + ".done"            # written only after makeblastdb fully succeeds
    lock = os.path.join(d, ".lock")
    with open(lock, "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        if not os.path.exists(done):   # absent => never built or a partial/interrupted build
            r = subprocess.run(
                ["makeblastdb", "-in", ref_fasta, "-dbtype", "nucl", "-out", db],
                capture_output=True, text=True,
            )
            if r.returncode != 0:
                raise RuntimeError(f"makeblastdb failed: {r.stderr[-400:]}")
            open(done, "w").close()
    return db


def _overlap_len(intervals: list[tuple], iv: tuple) -> int:
    """Total overlap of [a,b) with the union of *intervals*."""
    a, b = iv
    return sum(max(0, min(e, b) - max(s, a)) for s, e in _disjoint(intervals))


def _disjoint(intervals: list[tuple]) -> list[tuple]:
    out: list[tuple] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def _group_hsps(hsps: list[tuple], qlen: int, gap_bound: int | None = None) -> list[list[tuple]]:
    """Group HSPs (same query) into loci. Each hsp = (sid,s0,e0,strand,ident,q0,q1).

    HSPs join a locus only if they are FRAGMENTS of one alignment — same
    seqid+strand, bounded subject gap, AND their **query** interval barely
    overlaps the group's existing query coverage (a fragment covers a *different*
    part of the query). Two paralogs each align the *whole* query, so their query
    intervals overlap heavily and they are kept as separate loci — this is what
    preserves the fail-loud-on-ambiguity guarantee (a plain subject-gap rule
    silently merges adjacent paralogs).

    *gap_bound* limits the subject gap bridged within a locus (default one query
    length; tblastn passes qlen*3 since query is protein, subject nucleotide)."""
    bound = qlen if gap_bound is None else gap_bound
    groups: list[list[tuple]] = []
    qivs: list[list[tuple]] = []
    for h in sorted(hsps, key=lambda x: (x[0], x[3], x[1])):
        sid, s0, e0, strand, _ident, q0, q1 = h
        if groups:
            g = groups[-1]
            same = g[-1][0] == sid and g[-1][3] == strand
            gap = s0 - max(x[2] for x in g)
            # competing-vs-fragment test, symmetric: overlap relative to the
            # SMALLER of (this HSP, group's covered query). Dividing by this HSP
            # alone let a tiny spurious upstream HSP seed the locus and then absorb
            # the real full-length HSP (small overlap/long-HSP), over-extending the
            # region; min() makes "the seed is fully re-covered" read as competing.
            covered = _union_len(qivs[-1])
            q_overlap = _overlap_len(qivs[-1], (q0, q1)) / max(1, min(q1 - q0, covered))
            if same and gap <= bound and q_overlap < 0.5:   # a distinct query fragment
                g.append(h)
                qivs[-1].append((q0, q1))
                continue
        groups.append([h])
        qivs.append([(q0, q1)])
    return groups


def _union_len(intervals: list[tuple]) -> int:
    out = 0
    cur_s = cur_e = None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                out += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        out += cur_e - cur_s
    return out


def _locate_blastn(ref_fasta, query_fasta, qlen, threads, min_id, min_cov) -> dict[str, list]:
    if shutil.which("blastn") is None or shutil.which("makeblastdb") is None:
        raise RuntimeError(
            "target.locator == 'blastn' requires BLAST+ (blastn/makeblastdb) on PATH "
            "— install it (conda install -c bioconda blast) or use locator: bwa"
        )
    db = _blast_db(ref_fasta)
    fmt = "6 qseqid sseqid pident qstart qend sstart send sstrand"
    # -task blastn (word size 11), NOT the default megablast (word size 28): a
    # divergent paralog's short exon HSPs (e.g. 150bp at ~92%) are seeded by
    # blastn but missed by megablast -> they MUST be found or the fail-loud
    # ambiguity check silently misses a paralog. No -perc_identity prefilter
    # either: a locus can pass on its AGGREGATE identity even when one HSP is low;
    # the real threshold is applied below on the grouped locus.
    r = subprocess.run(
        ["blastn", "-task", "blastn", "-query", query_fasta, "-db", db, "-outfmt", fmt,
         "-num_threads", str(threads), "-max_target_seqs", "100000"],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"blastn failed: {r.stderr[-400:]}")

    hsps: dict[str, list] = defaultdict(list)
    for line in r.stdout.splitlines():
        if not line.strip():
            continue
        qid, sid, pid, qs, qe, ss, se, sstrand = line.split("\t")
        if qid not in qlen:
            continue
        strand = "+" if sstrand == "plus" else "-"
        s0, e0 = min(int(ss), int(se)) - 1, max(int(ss), int(se))
        q0, q1 = min(int(qs), int(qe)) - 1, max(int(qs), int(qe))
        hsps[qid].append((sid, s0, e0, strand, float(pid) / 100.0, q0, q1))

    return _group_to_loci(hsps, qlen, min_id, min_cov)


def _aggregate(group: list[tuple]) -> tuple[int, float]:
    """(unioned query bases, identity) for a locus, weighting identity over the
    SAME non-overlapping query basis as coverage (so an HSP that re-covers query
    already counted doesn't double-weight identity)."""
    covered: list[tuple] = []
    total = 0
    isum = 0.0
    for h in sorted(group, key=lambda x: x[5]):
        q0, q1, ident = h[5], h[6], h[4]
        new = (q1 - q0) - _overlap_len(covered, (q0, q1))
        if new > 0:
            total += new
            isum += new * ident
            covered.append((q0, q1))
    return total, (isum / total if total else 0.0)


def _group_to_loci(hsps, qlen, min_id, min_cov, gap_factor=1) -> dict[str, list]:
    """Group per-query HSPs into loci and keep those whose aggregate identity +
    query coverage (both on the unioned query basis) pass the thresholds."""
    loci_by_q: dict[str, list] = {}
    for q in qlen:
        loci = []
        for g in _group_hsps(hsps.get(q, []), qlen[q], gap_bound=qlen[q] * gap_factor):
            covered, ident = _aggregate(g)
            qcov = covered / qlen[q] if qlen[q] else 0.0
            sid = g[0][0]
            s0 = min(h[1] for h in g)
            e0 = max(h[2] for h in g)
            if ident >= min_id and qcov >= min_cov:
                loci.append((sid, s0, e0, g[0][3], ident, qcov))
        loci_by_q[q] = loci
    return loci_by_q


def _locate_tblastn(ref_fasta, query_fasta, qlen, threads, min_id, min_cov) -> dict[str, list]:
    """Protein query -> nucleotide reference, 6-frame translated (tblastn).

    identity/coverage are at the PROTEIN level; the located subject region is
    nucleotide (the gene's coding span). gap_bound is x3 (protein->nucleotide)."""
    if shutil.which("tblastn") is None or shutil.which("makeblastdb") is None:
        raise RuntimeError(
            "target.locator == 'tblastn' requires BLAST+ (tblastn/makeblastdb) on PATH "
            "— install it (conda install -c bioconda blast) or use a DNA locator"
        )
    db = _blast_db(ref_fasta)
    fmt = "6 qseqid sseqid pident qstart qend sstart send"
    r = subprocess.run(
        ["tblastn", "-query", query_fasta, "-db", db, "-outfmt", fmt,
         "-num_threads", str(threads), "-max_target_seqs", "100000"],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"tblastn failed: {r.stderr[-400:]}")
    hsps: dict[str, list] = defaultdict(list)
    for line in r.stdout.splitlines():
        if not line.strip():
            continue
        qid, sid, pid, qs, qe, ss, se = line.split("\t")
        if qid not in qlen:
            continue
        ssi, sei = int(ss), int(se)
        strand = "+" if ssi <= sei else "-"          # subject order encodes strand
        s0, e0 = min(ssi, sei) - 1, max(ssi, sei)    # nucleotide
        q0, q1 = min(int(qs), int(qe)) - 1, max(int(qs), int(qe))   # protein
        hsps[qid].append((sid, s0, e0, strand, float(pid) / 100.0, q0, q1))
    return _group_to_loci(hsps, qlen, min_id, min_cov, gap_factor=3)


def _locate_minimap2(ref_fasta, query_fasta, qlen, threads, min_id, min_cov) -> dict[str, list]:
    """Spliced transcript -> genome (minimap2 'splice' preset, intron-aware).

    Each spliced alignment spans one genomic locus (exons + introns). NOTE:
    minimap2 suppresses divergent secondary alignments, so the fail-loud-on-
    paralog guarantee is weaker here than blastn/tblastn — use this for the splice
    capability; prefer blastn/tblastn when paralog disambiguation is critical."""
    try:
        import mappy
    except ImportError as e:
        raise RuntimeError(
            "target.locator == 'minimap2' requires the 'mappy' package "
            "(pip install mappy) or use locator: blastn"
        ) from e
    log.warning("locator=minimap2 is intron-aware but does NOT reliably flag paralogs "
                "(minimap2 suppresses divergent secondaries); use blastn/tblastn when "
                "paralog disambiguation matters")
    aligner = mappy.Aligner(ref_fasta, preset="splice", best_n=25, n_threads=threads)
    if not aligner:
        raise RuntimeError(f"minimap2 (splice) index build failed for {ref_fasta}")
    loci_by_q: dict[str, list] = {q: [] for q in qlen}
    for name, seq, _q in mappy.fastx_read(query_fasta):
        if name not in qlen:
            continue
        hits = []
        for h in aligner.map(seq):
            ident = h.mlen / h.blen if h.blen else 0.0
            cov = (h.q_en - h.q_st) / qlen[name] if qlen[name] else 0.0
            if ident >= min_id and cov >= min_cov:
                hits.append((h.ctg, h.r_st, h.r_en, "+" if h.strand == 1 else "-", ident, cov))
        loci_by_q[name] = _merge_overlapping(hits)
    return loci_by_q


def _miniprot_index(ref_fasta: str, threads: int) -> str:
    """Build (once, cached + file-locked) a miniprot .mpi index for *ref_fasta*."""
    st = os.stat(ref_fasta)
    key = hashlib.sha1(f"{os.path.abspath(ref_fasta)}:{st.st_size}:{st.st_mtime_ns}".encode()).hexdigest()[:16]
    d = str(ensure_writable(cache_root() / "miniprot" / key))
    idx = os.path.join(d, "ref.mpi")
    done = idx + ".done"
    with open(os.path.join(d, ".lock"), "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        if not (os.path.exists(done) and os.path.exists(idx)):
            r = subprocess.run(["miniprot", "-t", str(threads), "-d", idx, ref_fasta],
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"miniprot index failed: {r.stderr[-400:]}")
            open(done, "w").close()
    return idx


def _locate_miniprot(ref_fasta, query_fasta, qlen, threads, min_id, min_cov) -> dict[str, list]:
    """Protein query -> genome, splice + frameshift aware (miniprot).

    miniprot models introns internally (one alignment per genomic locus). Like
    minimap2 it suppresses divergent secondaries, so it does NOT reliably flag
    paralogs — logs a warning; use tblastn when paralog disambiguation matters."""
    if shutil.which("miniprot") is None:
        raise RuntimeError(
            "target.locator == 'miniprot' requires miniprot on PATH "
            "(conda install -c bioconda miniprot) or use locator: tblastn"
        )
    log.warning("locator=miniprot is splice/frameshift-aware but does NOT reliably flag "
                "paralogs (like minimap2); use tblastn when paralog disambiguation matters")
    idx = _miniprot_index(ref_fasta, threads)
    r = subprocess.run(
        ["miniprot", "-t", str(threads), "--outs=0.5", "-N", "25", idx, query_fasta],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"miniprot failed: {r.stderr[-400:]}")
    loci_by_q: dict[str, list] = {q: [] for q in qlen}
    for line in r.stdout.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 11:
            continue
        qn = f[0]
        if qn not in qlen:
            continue
        qs, qe = int(f[2]), int(f[3])              # protein coords
        tname, ts, te = f[5], int(f[7]), int(f[8])  # 0-based half-open genomic
        nmatch, alen = int(f[9]), int(f[10])
        identity = nmatch / alen if alen else 0.0
        coverage = (qe - qs) / qlen[qn] if qlen.get(qn) else 0.0
        if identity >= min_id and coverage >= min_cov:
            loci_by_q[qn].append((tname, ts, te, f[4], identity, coverage))
    return {q: _merge_overlapping(loci_by_q[q]) for q in qlen}


_LOCATORS = {
    "bwa": _locate_bwa,
    "blastn": _locate_blastn,
    "tblastn": _locate_tblastn,
    "minimap2": _locate_minimap2,
    "miniprot": _locate_miniprot,
}


# --- public entry point ----------------------------------------------------

def _query_is_protein(path: str) -> bool:
    """Sniff a query FASTA as nucleotide vs protein by composition.

    Classifies on **protein-exclusive** residues (E F I L P Q Z X J O — letters that
    are not nucleotide codes), so IUPAC-degenerate DNA (R/Y/S/W/K/M/B/D/H/V) is NOT
    misread as protein. DNA/RNA has 0 % of these; a real protein has ~30 %. This
    only tells nucleotide from protein — a gene vs a transcript (both nucleotide)
    is not distinguishable and does not need to be (both use blastn)."""
    nucleotide = set("ACGTUNRYSWKMBDHV")     # IUPAC nucleotide alphabet
    total = protein_only = 0
    for line in open(path):
        if line.startswith(">"):
            continue
        for ch in line.strip().upper():
            if ch.isalpha():
                total += 1
                protein_only += ch not in nucleotide
    if total == 0:
        raise ValueError(f"query FASTA {path!r} has no sequence residues")
    return protein_only / total > 0.05


def resolve(config, design_ref) -> list[ResolvedTarget]:
    if not config.target.query_fasta:
        raise ValueError("target.mode == 'target_fasta' requires target.query_fasta")
    qlen = _read_fasta_lengths(config.target.query_fasta)
    if not qlen or not any(qlen.values()):
        raise ValueError(f"query FASTA {config.target.query_fasta!r} has no usable records")
    min_id, min_cov = config.target.min_identity, config.target.min_coverage
    overrides = config.specificity.params.bwa_aln.overrides
    locator = config.target.locator
    if locator == "auto":
        protein = _query_is_protein(config.target.query_fasta)
        locator = "tblastn" if protein else "blastn"
        log.info("locator=auto -> %s (query looks like %s)",
                 locator, "protein" if protein else "DNA")
    locate = _LOCATORS[locator]
    loci = locate(config.reference.fasta, config.target.query_fasta,
                  qlen, config.runtime.threads, min_id, min_cov)

    out_targets: list[ResolvedTarget] = []
    for q in qlen:
        passing = loci.get(q, [])
        if not passing:
            raise ValueError(
                f"query {q!r} not located in {config.reference.fasta} "
                f"(no hit with identity>={min_id}, coverage>={min_cov})"
            )
        if len(passing) > 1:
            locs = ", ".join(f"{h[0]}:{h[1] + 1}-{h[2]}({h[3]})" for h in passing)
            raise ValueError(
                f"query {q!r} matched {len(passing)} loci ({locs}) — ambiguous; "
                f"refine the query or raise min_identity/min_coverage"
            )
        seqid, start0, end0, strand, identity, cov = passing[0]
        region = Region(seqid, start0, end0, strand)
        log.info("located %s -> %s:%d-%d:%s (id=%.3f cov=%.3f, %s)",
                 q, seqid, start0 + 1, end0, strand, identity, cov, locator)
        out_targets.append(ResolvedTarget(
            target_id=q, sequence=design_ref.fetch(region), region=region,
            source_mode="target_fasta", mismatch_override=overrides.get(q),
        ))
    return out_targets
