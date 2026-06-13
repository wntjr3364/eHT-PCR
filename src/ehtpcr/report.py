"""Output writers: candidates.tsv, amplicons.tsv, manifest.json, run.yaml.

Schemas are versioned (``schema_version``) and additive (new columns are
non-breaking). Coordinate columns carry the convention in their names
(``start0``/``end0``). See docs/legacy_bwa_semantics.md for the exact column lists.
"""
from __future__ import annotations

import csv
import json
from typing import Iterable, Mapping

SCHEMA_VERSION = 1

CANDIDATES_COLUMNS = [
    "schema_version", "target_id", "pair_id", "f_seq", "r_seq",
    "f_start0", "f_end0", "r_start0", "r_end0", "product_size_design",
    "f_tm", "r_tm", "f_gc", "r_gc",
    "max_mismatch_used",   # per-target bwa -n actually applied (self-describing run)
    "n_valid_amplicons", "unique", "keep_decision", "failure_reason",
]

AMPLICONS_COLUMNS = [
    "schema_version", "target_id", "pair_id", "amplicon_id", "seqid",
    "strand_class", "start0", "end0", "product_size",
    "f_mismatches", "r_mismatches", "f_hit_source", "r_hit_source", "on_target",
]


def _write_tsv(path: str, columns: list[str], rows: Iterable[Mapping]) -> int:
    n = 0
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            n += 1
    return n


def write_candidates(path: str, rows: Iterable[Mapping]) -> int:
    """Write candidates.tsv (one row per primer pair). Returns the row count."""
    return _write_tsv(path, CANDIDATES_COLUMNS, rows)


def write_amplicons(path: str, rows: Iterable[Mapping]) -> int:
    return _write_tsv(path, AMPLICONS_COLUMNS, rows)


def write_manifest(path: str, manifest: Mapping) -> None:
    """Write the run manifest (versions, reference stamps, config, argv) as JSON."""
    with open(path, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)
        fh.write("\n")


def write_run_yaml(path: str, config_dict: Mapping) -> None:
    """Write the merged effective config (for reproducibility) as YAML."""
    import yaml
    with open(path, "w") as fh:
        yaml.safe_dump(dict(config_dict), fh, sort_keys=False, default_flow_style=False)
