"""Pipeline orchestration: resolve targets -> design -> (specificity) -> filter.

Pure orchestration. Holds the **uniqueness policy** (not the engine):

    n_valid_amplicons = number of FR/RF amplicons within max_product_size
    on_target         = an amplicon on the target's seqid overlaps it by >= 1 bp
    unique            = (n_valid_amplicons == 1) and on_target
    keep: all | unique | nomatch

With ``specificity.enabled == False`` the mapping step is skipped (design-only).
"""
from __future__ import annotations

import hashlib
import os
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone

from . import __version__, report
from .design import enumerate_primers, pair_primers
from .io.fasta import FastaRef
from .logging import get_logger
from .records import Region
from .report import SCHEMA_VERSION
from .target import resolve_targets

log = get_logger("pipeline")


def _on_target(amp, target) -> bool:
    if target.region is None:
        return False
    if amp.seqid != target.region.seqid:
        return False
    # name mode: the WHOLE named entry is the target, so seqid match is enough.
    # (The design template and the specificity-reference entry of the same name may
    #  differ in length/offset — e.g. transcriptome vs curated transcriptome — so a
    #  coordinate-window check would wrongly reject the real on-target amplicon.)
    if target.source_mode == "name":
        return True
    return amp.as_region().overlap_bp(target.region) >= 1


def _keep(mode: str, unique: bool, n_valid: int, on_target: bool) -> tuple[bool, str]:
    if mode == "all":
        return True, ""
    if mode == "unique":
        if unique:
            return True, ""
        if n_valid == 0:
            return False, "no_amplicon"
        if not on_target:
            return False, "off_target"
        return False, "non_unique"
    if mode == "nomatch":
        return (n_valid == 0), ("" if n_valid == 0 else "has_amplicon")
    return True, ""


def _candidate_row(target, pair, *, design_only, max_mismatch_used=None,
                   n_valid=None, unique=None, keep_decision="", failure_reason="") -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "target_id": target.target_id,
        "pair_id": pair.pair_id,
        "f_seq": pair.forward.seq,
        "r_seq": pair.reverse.seq,
        "f_start0": pair.forward.start,
        "f_end0": pair.forward.end,
        "r_start0": pair.reverse.start,
        "r_end0": pair.reverse.end,
        "product_size_design": pair.product_size,
        "f_tm": round(pair.forward.tm, 3),
        "r_tm": round(pair.reverse.tm, 3),
        "f_gc": round(pair.forward.gc, 4),
        "r_gc": round(pair.reverse.gc, 4),
        "max_mismatch_used": "" if design_only else max_mismatch_used,
        "n_valid_amplicons": "" if design_only else n_valid,
        "unique": "" if design_only else unique,
        "keep_decision": "design_only" if design_only else keep_decision,
        "failure_reason": failure_reason,
    }


def _amplicon_row(target, amp) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "target_id": target.target_id,
        "pair_id": amp.pair_id,
        "amplicon_id": f"{amp.pair_id}|{amp.seqid}:{amp.start}-{amp.end}:{amp.strand_class}",
        "seqid": amp.seqid,
        "strand_class": amp.strand_class,
        "start0": amp.start,
        "end0": amp.end,
        "product_size": amp.product_size,
        "f_mismatches": amp.f_mismatches,
        "r_mismatches": amp.r_mismatches,
        "f_hit_source": amp.f_hit_source,
        "r_hit_source": amp.r_hit_source,
        "on_target": _on_target(amp, target),
    }


_HASH_MAX_BYTES = 50_000_000   # full sha256 only for files up to this size


def _file_stamp(path) -> dict | None:
    """Reproducibility stamp for a reference file (size+mtime always; sha256 if small)."""
    if not path or not os.path.exists(path):
        return None
    st = os.stat(path)
    stamp = {"path": os.path.abspath(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    if st.st_size <= _HASH_MAX_BYTES:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        stamp["sha256"] = h.hexdigest()
    else:
        stamp["sha256"] = None   # too large to hash every run; size+mtime recorded
    return stamp


# Per-process cache of (engine, prepared reference): reused across all targets in
# the same process so the bwa index + FastaRef are not rebuilt per target.
_PREPARED_CACHE: dict = {}


def _get_prepared(config, spec_fasta, spec_index_dir):
    st = os.stat(spec_fasta)
    # include size+mtime so an in-place reference change re-keys (mirrors the bwa
    # index cache key) — avoids stale results when run() is reused in one process.
    key = (config.specificity.engine, os.path.abspath(spec_fasta), spec_index_dir,
           st.st_size, st.st_mtime_ns)
    if key not in _PREPARED_CACHE:
        from .specificity import get_engine
        eng = get_engine(config.specificity.engine)
        _PREPARED_CACHE[key] = (eng, eng.prepare(spec_fasta, index_dir=spec_index_dir))
    return _PREPARED_CACHE[key]


def _rows_for_target(config, target, pairs, amplicons, mm):
    """Apply the uniqueness policy to one target's amplicons -> (cand, amp, n_unique)."""
    by_pair: dict[str, list] = defaultdict(list)
    for a in amplicons:
        by_pair[a.pair_id].append(a)
    cand, amp, n_unique = [], [], 0
    for p in pairs:
        a_list = by_pair.get(p.pair_id, [])
        n_valid = len(a_list)
        on_t = any(_on_target(a, target) for a in a_list)
        unique = (n_valid == 1) and on_t
        n_unique += int(unique)
        kept, reason = _keep(config.specificity.keep, unique, n_valid, on_t)
        cand.append(_candidate_row(
            target, p, design_only=False, max_mismatch_used=mm, n_valid=n_valid,
            unique=unique, keep_decision=("kept" if kept else "filtered"), failure_reason=reason))
        amp.extend(_amplicon_row(target, a) for a in a_list)
    return cand, amp, n_unique


def _process_target(config, target, spec_fasta, spec_index_dir):
    """Worker (process-safe): design + optional specificity for ONE target.

    Returns ``(candidate_rows, amplicon_rows, n_unique)``. The engine/index are
    cached per process (the on-disk index is built once in :func:`run`).
    """
    design_only = not config.specificity.enabled
    primers = enumerate_primers(target.sequence, config.design)
    pairs = pair_primers(
        target.target_id, primers, config.design,
        max_pairs=config.runtime.max_pairs_per_target,
        fail_on_overflow=config.runtime.fail_on_overflow,
    )
    if design_only:
        return [_candidate_row(target, p, design_only=True) for p in pairs], [], 0
    if not pairs:
        return [], [], 0

    from .specificity.params import BwaAlnParams
    engine, prepared = _get_prepared(config, spec_fasta, spec_index_dir)
    bwa_cfg = config.specificity.params.bwa_aln
    mm = target.mismatch_override if target.mismatch_override is not None else bwa_cfg.max_mismatch
    params = BwaAlnParams(
        max_mismatch=mm, max_product_size=bwa_cfg.max_product_size,
        drop_n_rich_windows=bwa_cfg.drop_n_rich_windows,
        n_window_threshold=bwa_cfg.n_window_threshold,
        max_amplicons_per_pair=bwa_cfg.max_amplicons_per_pair,
    )
    amplicons = engine.find_amplicons(
        pairs, prepared, params, threads=config.runtime.threads,
        tmp_root=None, keep_temp=config.runtime.keep_temp,
    )
    return _rows_for_target(config, target, pairs, amplicons, mm)


def run(config) -> None:
    design_only = not config.specificity.enabled
    design_ref = FastaRef(config.reference.fasta)
    targets = resolve_targets(config, design_ref)
    if not targets:
        msg = ("no targets resolved (check target.mode / name / region) — a name "
               "substring that matches nothing is the usual cause")
        if config.runtime.fail_on_no_targets:
            raise ValueError(msg)
        log.warning(msg)

    out_dir = config.runtime.out_dir
    os.makedirs(out_dir, exist_ok=True)

    engine = prepared = None
    spec_fasta = None
    spec_names: set[str] = set()
    resolved_overrides: dict[str, int] = {}
    if not design_only and targets:    # no targets -> nothing to map; skip index build
        spec_fasta = config.specificity.reference.fasta or config.reference.fasta
        engine, prepared = _get_prepared(config, spec_fasta, config.specificity.reference.index_dir)
        # one FastaRef over the prepared (cache-symlinked) reference, reused for the
        # seqid set here and the N-window filter in the engine (.fai lands in the
        # writable cache even when the original db/ is read-only).
        spec_ref_fa = FastaRef(prepared.fasta)
        prepared._fasta_ref = spec_ref_fa
        spec_names = set(spec_ref_fa.names())
        log.info("specificity reference: %s (engine=%s)", spec_fasta, config.specificity.engine)

    # pre-pass (main process): seqid warnings + resolved per-target overrides
    if not design_only:
        for t in targets:
            if t.region and t.region.seqid not in spec_names:
                log.warning(
                    "target %s seqid %r is not in the specificity reference -> every "
                    "amplicon will be off_target (separate/renamed reference?)",
                    t.target_id, t.region.seqid,
                )
            if t.mismatch_override is not None:
                resolved_overrides[t.target_id] = t.mismatch_override

    spec_index_dir = config.specificity.reference.index_dir
    results: list = [None] * len(targets)     # by target index (ids may collide)

    def _collect(i, t, result):
        results[i] = result
        cand, amp, n_unique = result
        log.info("target %s: %d candidates, %d amplicons, %d unique",
                 t.target_id, len(cand), len(amp), n_unique)

    # Per target: design + (specificity) -> rows. Each target's oligos are aligned
    # in their own bwa call: that per-target grouping is the jar-parity contract
    # (batching oligos across targets changes bwa's per-read reporting), so speed
    # comes from running targets in parallel (runtime.jobs), not from merging them.
    jobs = config.runtime.jobs
    if jobs > 1 and len(targets) > 1:
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        log.info("processing %d targets with %d parallel jobs", len(targets), jobs)
        # 'spawn' avoids fork-in-a-multi-threaded-process deadlocks (bwa subprocesses).
        workers = min(jobs, len(targets))
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
            futures = [ex.submit(_process_target, config, t, spec_fasta, spec_index_dir) for t in targets]
            for i, (t, fut) in enumerate(zip(targets, futures)):
                _collect(i, t, fut.result())
    else:
        for i, t in enumerate(targets):
            _collect(i, t, _process_target(config, t, spec_fasta, spec_index_dir))

    # emit in target order regardless of dispatch order
    cand_rows: list[dict] = []
    amp_rows: list[dict] = []
    for res in results:
        cand, amp, _ = res if res is not None else ([], [], 0)
        cand_rows.extend(cand)
        amp_rows.extend(amp)

    cand_path = os.path.join(out_dir, "candidates.tsv")
    n = report.write_candidates(cand_path, cand_rows)
    log.info("wrote %d candidate rows -> %s", n, cand_path)
    if not design_only:
        amp_path = os.path.join(out_dir, "amplicons.tsv")
        m = report.write_amplicons(amp_path, amp_rows)
        log.info("wrote %d amplicon rows -> %s", m, amp_path)

    # reproducibility artifacts
    manifest = {
        "ehtpcr_version": __version__,
        "python_version": platform.python_version(),
        "bwa_version": (prepared.meta.get("bwa_version") if prepared else None),
        "command": sys.argv,
        "timestamp": datetime.now(timezone.utc).astimezone().isoformat(),
        "design_only": design_only,
        "n_targets": len(targets),
        "reference": _file_stamp(config.reference.fasta),
        "specificity_reference": (None if design_only else _file_stamp(spec_fasta)),
        "engine": (None if design_only else config.specificity.engine),
        "keep": (None if design_only else config.specificity.keep),
        "mismatch_overrides_applied": resolved_overrides,
    }
    report.write_manifest(os.path.join(out_dir, "manifest.json"), manifest)
    report.write_run_yaml(os.path.join(out_dir, "run.yaml"), config.model_dump(mode="json"))
    log.info("wrote manifest.json + run.yaml -> %s", out_dir)
