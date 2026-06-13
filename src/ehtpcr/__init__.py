"""eHT-PCR: general-purpose in-silico PCR — design unique primer pairs.

Organism-agnostic. Given a target region in a reference, enumerate candidate
primer pairs and keep those whose PCR product is unique against a (possibly
separate, possibly curated) specificity reference.

Coordinate convention: all internal Region/Amplicon objects are 0-based,
half-open. See docs/legacy_bwa_semantics.md.
"""

__version__ = "0.1.0"
