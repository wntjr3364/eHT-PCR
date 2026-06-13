"""Engine-specific specificity parameters (pydantic; selected by engine name).

In the config, ``specificity.params`` is a discriminated union keyed by the
``engine`` field, so the config model and the engine share one type (no
dataclass<->pydantic shim). Keeping params engine-specific avoids leaking
``max_mismatch`` / ``max_product_size`` (bwa/PCR concepts) into the generic
engine interface.
"""
from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, field_validator


class BwaAlnParams(BaseModel):
    """Parameters for the bwa_aln engine.

    ``max_mismatch`` is bwa ``-n`` and MUST be an integer edit distance (a float
    would mean a fraction / error-rate). ``max_product_size`` is the legacy jar's
    ``-n`` (max PCR product span). The other bwa flags are fixed invariants (see
    docs/legacy_bwa_semantics.md) and are not user-tunable.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    engine: Literal["bwa_aln"] = "bwa_aln"
    max_mismatch: int = 5
    max_product_size: int = 3000
    # Legacy jar drops a hit whose reference window has >= threshold N's
    # (docs/legacy_bwa_semantics.md D6). Default ON = jar-faithful; turn off to
    # keep such hits (more conservative on N-rich draft genomes).
    drop_n_rich_windows: bool = True
    n_window_threshold: int = 5
    # Guard against amplicon explosion on repetitive references (O(F_hits x R_hits)
    # per pair). A pair exceeding this is non-specific by definition; stop
    # materializing once reached. 0 disables the cap.
    max_amplicons_per_pair: int = 100_000

    @field_validator("max_mismatch", mode="before")
    @classmethod
    def _mismatch_is_nonneg_int(cls, v):
        # before-mode: pydantic would coerce bool->int before an after-validator,
        # so a bool check must run on the raw value (true/false != 1/0 here).
        if isinstance(v, bool):
            raise ValueError("max_mismatch must be an int (bwa -n edit distance), not bool")
        if isinstance(v, int) and v < 0:
            raise ValueError("max_mismatch must be >= 0")
        return v

    @field_validator("max_product_size")
    @classmethod
    def _product_size_positive(cls, v: int) -> int:
        if v <= 0:
            raise ValueError("max_product_size must be > 0")
        return v

    @field_validator("max_amplicons_per_pair")
    @classmethod
    def _cap_not_one(cls, v: int) -> int:
        if v < 0:
            raise ValueError("max_amplicons_per_pair must be >= 0 (0 = no cap)")
        if v == 1:
            # a cap of 1 truncates a non-specific pair to one amplicon, which reads
            # downstream as unique. 0 (no cap) or >= 2 only.
            raise ValueError(
                "max_amplicons_per_pair must be 0 (no cap) or >= 2 — a cap of 1 "
                "cannot distinguish a unique pair from a capped non-specific one")
        return v


# Future: BlastParams(engine: Literal["blast"], ...), ThermoParams(...) join the
# discriminated union; config selects by the ``engine`` discriminator.
SpecificityParams = BwaAlnParams
