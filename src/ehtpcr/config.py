"""Configuration model (pydantic v2) + loader.

Precedence: pydantic defaults < YAML file < environment (``EHTPCR_*``) < explicit
CLI flags. Only **non-None** CLI values override (an unset flag must not clobber a
YAML value). See config/example.yaml for the schema.
"""
from __future__ import annotations

import os
from typing import Literal, Optional

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class _Base(BaseModel):
    model_config = ConfigDict(extra="forbid")


class MinMaxF(_Base):
    min: float
    max: float

    @model_validator(mode="after")
    def _ordered(self):
        if self.max < self.min:
            raise ValueError(f"max ({self.max}) < min ({self.min})")
        return self


class MinMaxI(_Base):
    min: int
    max: int

    @model_validator(mode="after")
    def _ordered(self):
        if self.max < self.min:
            raise ValueError(f"max ({self.max}) < min ({self.min})")
        return self


class Reference(_Base):
    fasta: str
    gff: Optional[str] = None
    gff_db_dir: Optional[str] = None


class Flank(_Base):
    left: int = 0
    right: int = 0

    @model_validator(mode="after")
    def _nonneg(self):
        if self.left < 0 or self.right < 0:
            raise ValueError("target.flank left/right must be >= 0")
        return self


class Target(_Base):
    mode: Literal["name", "region", "locus", "target_fasta"] = "name"
    name: list[str] = Field(default_factory=list)
    region: Optional[str] = None
    locus: list[str] = Field(default_factory=list)
    feature: str = "gene"                # GFF feature type for mode: locus
    query_fasta: Optional[str] = None    # mode: target_fasta
    # how to locate the query in the reference:
    #   auto     - DNA query (gene/transcript) -> blastn, protein query -> tblastn
    #   bwa      - DNA query, near-identical (no extra dep)
    #   blastn   - DNA query, sensitive + paralog-safe (BLAST+)
    #   minimap2 - spliced transcript -> genome (intron-aware; paralog-blind)
    #   tblastn  - protein query -> genome (6-frame translated, paralog-safe; BLAST+)
    #   miniprot - protein query -> genome (splice/frameshift-aware; paralog-blind)
    locator: Literal["auto", "bwa", "blastn", "minimap2", "tblastn", "miniprot"] = "auto"
    min_identity: float = 0.90           # target_fasta: min alignment identity
    min_coverage: float = 0.90           # target_fasta: min query coverage
    flank: Flank = Field(default_factory=Flank)

    @model_validator(mode="after")
    def _check(self):
        for name in ("min_identity", "min_coverage"):
            if not 0.0 <= getattr(self, name) <= 1.0:
                raise ValueError(f"target.{name} must be in [0, 1]")
        # An empty/whitespace substring matches every header (name mode) and an
        # empty locus id matches nothing useful -> reject rather than silently
        # designing against the whole reference (fail-loud).
        for field in ("name", "locus"):
            for s in getattr(self, field):
                if not s or not s.strip():
                    raise ValueError(
                        f"target.{field} entries must be non-empty "
                        f"(an empty substring would match the whole reference)")
        return self


class TmCfg(_Base):
    min: float
    max: float
    # legacy (default, verbatim Breslauer/Borer port) | NN tables breslauer /
    # sugimoto / santalucia / santalucia2004 | wallace | gc | primer3 (biopython
    # powers the NN/Wallace/GC models, primer3-py powers 'primer3' — both bundled).
    model: Literal["legacy", "breslauer", "sugimoto", "santalucia",
                   "santalucia2004", "wallace", "gc", "primer3"] = "legacy"

    @model_validator(mode="after")
    def _check(self):
        if self.max < self.min:
            raise ValueError(f"tm.max ({self.max}) < tm.min ({self.min})")
        return self


class PairCfg(_Base):
    tm_diff: MinMaxF = Field(default_factory=lambda: MinMaxF(min=0.0, max=1.0))
    gc_diff: MinMaxF = Field(default_factory=lambda: MinMaxF(min=0.0, max=1.0))


class Design(_Base):
    primer_len: MinMaxI = Field(default_factory=lambda: MinMaxI(min=18, max=24))
    gc: MinMaxF = Field(default_factory=lambda: MinMaxF(min=0.30, max=0.70))
    tm: TmCfg = Field(default_factory=lambda: TmCfg(min=57, max=62))
    gc_clamp: int = 2
    product_size: MinMaxI = Field(default_factory=lambda: MinMaxI(min=150, max=250))
    pair: PairCfg = Field(default_factory=PairCfg)
    allow_degenerate: bool = False
    softmask_policy: Literal["allow", "reject", "lowercase_penalty"] = "allow"

    @model_validator(mode="after")
    def _check(self):
        if self.gc_clamp < 0:
            raise ValueError("design.gc_clamp must be >= 0")
        if self.primer_len.min < 1:
            raise ValueError("design.primer_len.min must be >= 1")
        if self.product_size.min < 1:
            raise ValueError("design.product_size.min must be >= 1")
        for b in ("min", "max"):
            if not 0.0 <= getattr(self.gc, b) <= 1.0:
                raise ValueError(f"design.gc.{b} must be a fraction in [0, 1]")
        if self.softmask_policy != "allow":
            raise ValueError(f"design.softmask_policy {self.softmask_policy!r} is deferred; use 'allow'")
        if self.allow_degenerate:
            raise ValueError("design.allow_degenerate is not implemented yet (GC/Tm/clamp/bwa are not IUPAC-aware); use false")
        return self


class BwaAlnConfig(_Base):
    """Config-side bwa_aln params (includes the per-target override map)."""

    max_mismatch: int = 5
    max_product_size: int = 3000
    drop_n_rich_windows: bool = True   # jar drops hits whose ref window has >= N's (D6)
    n_window_threshold: int = 5
    max_amplicons_per_pair: int = 100_000   # cap to avoid amplicon explosion (0 = no cap)
    overrides: dict[str, int] = Field(default_factory=dict)

    @field_validator("overrides", mode="before")
    @classmethod
    def _none_to_empty(cls, v):
        # an empty YAML `overrides:` parses to None; treat as {}.
        return {} if v is None else v

    @field_validator("max_mismatch", mode="before")
    @classmethod
    def _mismatch_not_bool(cls, v):
        # pydantic coerces YAML true/false -> 1/0 BEFORE an after-validator runs, so
        # the bool check must be a before-validator. `max_mismatch: true` must error,
        # not silently become 1 (a core specificity parameter).
        if isinstance(v, bool):
            raise ValueError("max_mismatch must be an int (bwa -n edit distance), not bool")
        return v

    @model_validator(mode="after")
    def _check(self):
        if self.max_mismatch < 0:
            raise ValueError("max_mismatch must be >= 0")
        if self.max_product_size < 1:
            raise ValueError("max_product_size must be >= 1")
        if self.n_window_threshold < 1:
            raise ValueError("n_window_threshold must be >= 1")
        if self.max_amplicons_per_pair < 0:
            raise ValueError("max_amplicons_per_pair must be >= 0 (0 = no cap)")
        if self.max_amplicons_per_pair == 1:
            # a cap of 1 truncates a non-specific pair's amplicons to one, which the
            # pipeline then reads as n_valid==1 -> unique. 0 (no cap) or >=2 only.
            raise ValueError(
                "max_amplicons_per_pair must be 0 (no cap) or >= 2 — a cap of 1 "
                "cannot distinguish a truly unique pair from a capped non-specific one")
        for k, val in self.overrides.items():
            if val < 0:
                raise ValueError(f"override for {k!r} must be >= 0")
        return self


class EngineParams(_Base):
    bwa_aln: BwaAlnConfig = Field(default_factory=BwaAlnConfig)


class SpecRef(_Base):
    fasta: Optional[str] = None
    index_dir: Optional[str] = None


class Specificity(_Base):
    enabled: bool = True
    reference: SpecRef = Field(default_factory=SpecRef)
    engine: Literal["bwa_aln"] = "bwa_aln"
    keep: Literal["all", "unique", "nomatch"] = "unique"
    params: EngineParams = Field(default_factory=EngineParams)


class Runtime(_Base):
    threads: int = 8
    jobs: int = 1
    out_dir: str = "result/"
    log_level: str = "INFO"
    max_pairs_per_target: int = 1_000_000
    fail_on_overflow: bool = True
    fail_on_no_targets: bool = True   # error out (non-zero) when 0 targets resolve
    keep_temp: bool = False

    @model_validator(mode="after")
    def _positive(self):
        for name in ("threads", "jobs", "max_pairs_per_target"):
            if getattr(self, name) < 1:
                raise ValueError(f"runtime.{name} must be >= 1")
        return self


class Config(_Base):
    version: int = 1
    reference: Reference
    target: Target = Field(default_factory=Target)
    design: Design = Field(default_factory=Design)
    specificity: Specificity = Field(default_factory=Specificity)
    runtime: Runtime = Field(default_factory=Runtime)


# --- loading + precedence ---------------------------------------------------

def _deep_merge(base: dict, over: dict) -> dict:
    out = dict(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def _env_layer() -> dict:
    """Minimal EHTPCR_* environment overrides (between YAML and CLI)."""
    layer: dict = {}
    rt: dict = {}
    if (v := os.environ.get("EHTPCR_OUT_DIR")) is not None:
        rt["out_dir"] = v
    if (v := os.environ.get("EHTPCR_THREADS")) is not None:
        rt["threads"] = int(v)
    if (v := os.environ.get("EHTPCR_LOG_LEVEL")) is not None:
        rt["log_level"] = v
    if rt:
        layer["runtime"] = rt
    return layer


# Flat CLI override key -> nested config path. Every config knob a flag touches
# is listed here so `ehtpcr run` can set it without editing the YAML.
_CLI_PATHS: dict[str, tuple[str, ...]] = {
    "reference": ("reference", "fasta"),
    "gff": ("reference", "gff"),
    "spec_reference": ("specificity", "reference", "fasta"),
    "keep": ("specificity", "keep"),
    "engine": ("specificity", "engine"),
    "max_mismatch": ("specificity", "params", "bwa_aln", "max_mismatch"),
    "max_product_size": ("specificity", "params", "bwa_aln", "max_product_size"),
    "target_mode": ("target", "mode"),
    "name": ("target", "name"),
    "region": ("target", "region"),
    "locus": ("target", "locus"),
    "query_fasta": ("target", "query_fasta"),
    "locator": ("target", "locator"),
    "min_identity": ("target", "min_identity"),
    "min_coverage": ("target", "min_coverage"),
    "tm_min": ("design", "tm", "min"),
    "tm_max": ("design", "tm", "max"),
    "tm_model": ("design", "tm", "model"),
    "gc_min": ("design", "gc", "min"),
    "gc_max": ("design", "gc", "max"),
    "len_min": ("design", "primer_len", "min"),
    "len_max": ("design", "primer_len", "max"),
    "product_min": ("design", "product_size", "min"),
    "product_max": ("design", "product_size", "max"),
    "gc_clamp": ("design", "gc_clamp"),
    "threads": ("runtime", "threads"),
    "jobs": ("runtime", "jobs"),
    "out_dir": ("runtime", "out_dir"),
    "max_pairs": ("runtime", "max_pairs_per_target"),
    "fail_on_no_targets": ("runtime", "fail_on_no_targets"),
    "keep_temp": ("runtime", "keep_temp"),
    "log_level": ("runtime", "log_level"),
}


def _set_path(d: dict, path: tuple[str, ...], value) -> None:
    for k in path[:-1]:
        d = d.setdefault(k, {})
    d[path[-1]] = value


def _cli_layer(cli: dict) -> dict:
    """Map flat CLI override keys to the nested config shape (None already dropped)."""
    layer: dict = {}
    for key, value in cli.items():
        if key == "no_spec":
            _set_path(layer, ("specificity", "enabled"), not value)
        elif key in _CLI_PATHS:
            _set_path(layer, _CLI_PATHS[key], value)
    return layer


def load_config(path: Optional[str] = None, **cli) -> Config:
    data: dict = {}
    if path:
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
    # Drop unset flags: None AND [] (an empty list flag must not clobber a YAML
    # list). This mirrors the CLI's own filter so a direct/programmatic
    # load_config(..., name=[]) behaves identically to the CLI.
    cli = {k: v for k, v in cli.items() if v not in (None, [])}
    # Start from a COMPLETE default skeleton so partial overrides (e.g. only
    # --tm-min) merge onto fully-formed sub-objects and still validate.
    skeleton = Config(reference=Reference(fasta="")).model_dump(mode="python")
    merged = _deep_merge(skeleton, data)
    merged = _deep_merge(merged, _env_layer())
    merged = _deep_merge(merged, _cli_layer(cli))
    if not merged.get("reference", {}).get("fasta"):
        raise ValueError("reference.fasta is required (set it in the config or pass --reference)")
    return Config.model_validate(merged)
