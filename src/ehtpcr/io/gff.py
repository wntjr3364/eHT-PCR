"""GFF feature lookup via gffutils (replaces the hand-rolled GFF.py).

The gffutils SQLite DB is built once into a writable per-user cache (keyed by the
GFF path+size+mtime), never next to a possibly read-only reference.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

import gffutils

from ..cache import cache_root, ensure_writable
from ..records import Region


def _gff_db_path(gff_path: str, db_dir: str | None) -> str:
    d = ensure_writable(Path(db_dir) if db_dir else cache_root() / "gff")
    st = os.stat(gff_path)
    key = hashlib.sha1(
        f"{os.path.abspath(gff_path)}:{st.st_size}:{st.st_mtime_ns}".encode()
    ).hexdigest()[:16]
    return str(d / f"{key}.db")


class GffRef:
    """gffutils-backed lookup of a feature's coordinates by id."""

    def __init__(self, gff_path: str, db_dir: str | None = None) -> None:
        self.gff_path = gff_path
        db_path = _gff_db_path(gff_path, db_dir)
        if os.path.exists(db_path):
            self._db = gffutils.FeatureDB(db_path)
        else:
            # build to a temp path then atomically rename (avoids partial DB races)
            tmp = f"{db_path}.{os.getpid()}.tmp"
            gffutils.create_db(
                gff_path, tmp, force=True, keep_order=True,
                merge_strategy="create_unique", sort_attribute_values=True,
                disable_infer_genes=True, disable_infer_transcripts=True,
            )
            os.replace(tmp, db_path)
            self._db = gffutils.FeatureDB(db_path)

    def region_for(self, locus: str, feature: str = "gene") -> Region:
        """0-based half-open Region of the feature whose ID is *locus*.

        Raises ``ValueError`` if the resolved feature's type is not *feature*.
        """
        try:
            f = self._db[locus]
        except gffutils.exceptions.FeatureNotFoundError as e:
            raise KeyError(f"locus {locus!r} not found in {self.gff_path}") from e
        # gffutils' merge_strategy="create_unique" renames a duplicate ID to
        # "<id>_1"; if that sibling exists with the same featuretype, the locus is
        # multi-mapped -> fail loud rather than silently returning the first copy
        # (mirrors target_fasta's fail-loud-on->1-locus guarantee).
        # A synthetic rename KEEPS the original ID attribute (== [locus]); a
        # genuinely distinct feature literally named "<locus>_1" has ID == [locus_1].
        # Only the former is a real duplicate — checking the ID attribute avoids
        # rejecting a valid, unique locus that merely sits next to a "<locus>_1".
        try:
            dup = self._db[f"{locus}_1"]
        except gffutils.exceptions.FeatureNotFoundError:
            dup = None
        if (dup is not None and dup.featuretype == f.featuretype
                and dup.attributes.get("ID") == [locus]):
            raise ValueError(
                f"locus {locus!r} matches multiple {f.featuretype!r} features in "
                f"{self.gff_path} (duplicate ID) — resolve the ambiguity in the GFF "
                f"(the design would otherwise pick one copy arbitrarily)."
            )
        if feature and f.featuretype != feature:
            raise ValueError(
                f"locus {locus!r} is a {f.featuretype!r} feature, not {feature!r} "
                f"(set target.feature to {f.featuretype!r})"
            )
        strand = f.strand if f.strand in ("+", "-") else "+"
        # gffutils features are 1-based inclusive
        return Region(f.seqid, f.start - 1, f.end, strand)
