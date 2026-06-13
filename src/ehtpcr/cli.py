"""Command-line interface for eHT-PCR (typer).

Subcommands:
- ``run``    : the full pipeline (extract -> design -> specificity -> filter).
- ``index``  : eagerly build/cache the specificity index for a reference.

Design-only is ``run --no-spec``. Every ``run`` flag overrides the config field;
precedence is defaults < config.yaml < environment < explicit CLI flags.
"""
from __future__ import annotations

from typing import Optional

import typer

from . import __version__
from .logging import configure

app = typer.Typer(
    add_completion=False,
    pretty_exceptions_enable=False,   # don't dump rich tracebacks for expected errors
    help="General-purpose in-silico PCR — unique primer design.",
)

# Exceptions that represent user/input errors (shown as a clean one-line message).
_USER_ERRORS: tuple[type[BaseException], ...] = (
    FileNotFoundError, NotImplementedError, RuntimeError, KeyError, ValueError, OSError,
)


@app.callback(invoke_without_command=True)
def _root(
    ctx: typer.Context,
    version: bool = typer.Option(False, "--version", help="Show version and exit."),
) -> None:
    if version:
        typer.echo(f"ehtpcr {__version__}")
        raise typer.Exit()
    if ctx.invoked_subcommand is None:
        typer.echo(ctx.get_help())
        raise typer.Exit()


@app.command()
def run(
    config: Optional[str] = typer.Option(None, "--config", "-c", help="Run config YAML (optional; every field is also a flag below)."),
    # --- reference / specificity ---
    reference: Optional[str] = typer.Option(None, "--reference", "-r", help="reference.fasta (design sequence space).", rich_help_panel="Reference & specificity"),
    gff: Optional[str] = typer.Option(None, "--gff", help="reference.gff (for target-mode locus).", rich_help_panel="Reference & specificity"),
    spec_reference: Optional[str] = typer.Option(None, "--spec-reference", help="Separate specificity reference FASTA.", rich_help_panel="Reference & specificity"),
    no_spec: Optional[bool] = typer.Option(None, "--no-spec/--spec", help="Design only (skip specificity) / force it on.", rich_help_panel="Reference & specificity"),
    keep: Optional[str] = typer.Option(None, "--keep", help="Which pairs to keep: unique | all | nomatch.", rich_help_panel="Reference & specificity"),
    max_mismatch: Optional[int] = typer.Option(None, "--max-mismatch", "-m", help="bwa max mismatches for specificity.", rich_help_panel="Reference & specificity"),
    max_product_size: Optional[int] = typer.Option(None, "--max-product-size", help="Max amplicon size scanned for off-targets.", rich_help_panel="Reference & specificity"),
    # --- target ---
    target_mode: Optional[str] = typer.Option(None, "--target-mode", help="name | region | locus | target_fasta (inferred from the target flag if omitted).", rich_help_panel="Target"),
    name: Optional[list[str]] = typer.Option(None, "--name", help="Target name substring(s) (mode name). Repeatable.", rich_help_panel="Target"),
    region: Optional[str] = typer.Option(None, "--region", help="SEQID:START-END[:STRAND], 1-based (mode region).", rich_help_panel="Target"),
    locus: Optional[list[str]] = typer.Option(None, "--locus", help="Gene id(s) (mode locus, needs --gff). Repeatable.", rich_help_panel="Target"),
    query_fasta: Optional[str] = typer.Option(None, "--query-fasta", help="Query FASTA to locate (mode target_fasta).", rich_help_panel="Target"),
    locator: Optional[str] = typer.Option(
        None, "--locator", rich_help_panel="Target",
        help=(
            "Locate target_fasta queries. auto (default): blastn for DNA, tblastn for protein.  "
            "DNA gene/transcript: bwa | blastn | minimap2 (spliced to genome).  "
            "protein: tblastn | miniprot (spliced to genome).  "
            "minimap2/miniprot are intron-aware but paralog-blind."
        ),
    ),
    min_identity: Optional[float] = typer.Option(None, "--min-identity", help="target_fasta locate identity threshold.", rich_help_panel="Target"),
    min_coverage: Optional[float] = typer.Option(None, "--min-coverage", help="target_fasta locate coverage threshold.", rich_help_panel="Target"),
    # --- design ---
    tm_min: Optional[float] = typer.Option(None, "--tm-min", help="Min primer Tm.", rich_help_panel="Design"),
    tm_max: Optional[float] = typer.Option(None, "--tm-max", help="Max primer Tm.", rich_help_panel="Design"),
    tm_model: Optional[str] = typer.Option(
        None, "--tm-model", rich_help_panel="Design",
        help="Tm model: legacy (default) | breslauer | sugimoto | santalucia | "
             "santalucia2004 | wallace | gc | primer3."),
    gc_min: Optional[float] = typer.Option(None, "--gc-min", help="Min primer GC fraction.", rich_help_panel="Design"),
    gc_max: Optional[float] = typer.Option(None, "--gc-max", help="Max primer GC fraction.", rich_help_panel="Design"),
    len_min: Optional[int] = typer.Option(None, "--len-min", help="Min primer length.", rich_help_panel="Design"),
    len_max: Optional[int] = typer.Option(None, "--len-max", help="Max primer length.", rich_help_panel="Design"),
    product_min: Optional[int] = typer.Option(None, "--product-min", help="Min product size.", rich_help_panel="Design"),
    product_max: Optional[int] = typer.Option(None, "--product-max", help="Max product size.", rich_help_panel="Design"),
    gc_clamp: Optional[int] = typer.Option(None, "--gc-clamp", help="3' GC-clamp bases.", rich_help_panel="Design"),
    # --- runtime ---
    out_dir: Optional[str] = typer.Option(None, "--out", "-o", help="Output directory.", rich_help_panel="Runtime"),
    threads: Optional[int] = typer.Option(None, "--threads", "-t", help="bwa/aligner threads.", rich_help_panel="Runtime"),
    jobs: Optional[int] = typer.Option(None, "--jobs", "-j", help="Parallel target workers.", rich_help_panel="Runtime"),
    max_pairs: Optional[int] = typer.Option(None, "--max-pairs", help="Max primer pairs per target.", rich_help_panel="Runtime"),
    fail_on_no_targets: Optional[bool] = typer.Option(None, "--fail-on-no-targets/--no-fail-on-no-targets", help="Exit non-zero if no targets resolve (default on; catches name typos). Off = warn + empty output.", rich_help_panel="Runtime"),
    keep_temp: Optional[bool] = typer.Option(None, "--keep-temp/--no-keep-temp", help="Keep aligner temp files.", rich_help_panel="Runtime"),
    cache_dir: Optional[str] = typer.Option(None, "--cache-dir", help="Writable root for all eHT-PCR caches: bwa/BLAST/miniprot indexes, FASTA .fai, GFF DB (sets EHTPCR_CACHE_DIR). Default: $XDG_CACHE_HOME or ~/.cache.", rich_help_panel="Runtime"),
    log_level: Optional[str] = typer.Option(None, "--log-level", help="DEBUG|INFO|WARNING|ERROR.", rich_help_panel="Runtime"),
) -> None:
    """Run the full pipeline (extract -> design -> specificity -> filter).

    Provide a --config YAML and/or set fields directly with the flags below.
    Precedence: defaults < config < env (EHTPCR_*) < explicit CLI flags. An unset
    flag never clobbers a config value. With no --config, reference + target flags
    are enough (e.g. `ehtpcr run -r genome.fa --name MyGene -o result/`).
    """
    configure(log_level or "INFO")
    if cache_dir:                       # one writable root for every index cache,
        import os                       # inherited by spawned target workers too
        os.environ["EHTPCR_CACHE_DIR"] = cache_dir
    # Infer the target mode from the target flag the user typed, so --region /
    # --locus / --query-fasta just work (like --name, which is the default mode).
    # An explicit --target-mode always wins.
    if target_mode is None:
        given = [f for f, v in (("--name", name), ("--region", region),
                                ("--locus", locus), ("--query-fasta", query_fasta)) if v]
        if len(given) > 1:
            typer.secho(
                f"error: conflicting target flags {given} — pass exactly one, or set "
                "--target-mode explicitly to choose.", fg=typer.colors.RED, err=True)
            raise typer.Exit(1)
        # Use truthiness (not `is not None`) so it matches the conflict check above:
        # an empty-string `--region ''` must not silently flip the mode while the
        # conflict guard (which uses `if v`) treats it as "not given".
        if query_fasta:
            target_mode = "target_fasta"
        elif region:
            target_mode = "region"
        elif locus:
            target_mode = "locus"
    from .config import load_config
    from .pipeline import run as run_pipeline
    try:
        overrides = {k: v for k, v in {
            "reference": reference, "gff": gff, "spec_reference": spec_reference,
            "no_spec": no_spec, "keep": keep, "max_mismatch": max_mismatch,
            "max_product_size": max_product_size, "target_mode": target_mode,
            "name": name, "region": region, "locus": locus, "query_fasta": query_fasta,
            "locator": locator, "min_identity": min_identity, "min_coverage": min_coverage,
            "tm_min": tm_min, "tm_max": tm_max, "tm_model": tm_model,
            "gc_min": gc_min, "gc_max": gc_max,
            "len_min": len_min, "len_max": len_max, "product_min": product_min,
            "product_max": product_max, "gc_clamp": gc_clamp, "out_dir": out_dir,
            "threads": threads, "jobs": jobs, "max_pairs": max_pairs,
            "fail_on_no_targets": fail_on_no_targets,
            "keep_temp": keep_temp, "log_level": log_level,
        }.items() if v not in (None, [])}
        cfg = load_config(config, **overrides)
        configure(cfg.runtime.log_level)   # apply the resolved level (yaml/env/CLI)
        run_pipeline(cfg)
    except _USER_ERRORS as e:
        typer.secho(f"error: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    except Exception as e:  # pydantic ValidationError etc. -> clean one-liner
        typer.secho(f"error: {type(e).__name__}: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)


@app.command()
def index(
    fasta: str = typer.Argument(..., help="Reference FASTA to index."),
    engine: str = typer.Option("bwa_aln", "--engine"),
    cache_dir: Optional[str] = typer.Option(None, "--cache-dir", help="Writable cache root for the index (sets EHTPCR_CACHE_DIR)."),
) -> None:
    """Eagerly build/locate the specificity index for a reference."""
    configure("INFO")
    if cache_dir:
        import os
        os.environ["EHTPCR_CACHE_DIR"] = cache_dir
    from .specificity import get_engine
    try:
        prepared = get_engine(engine).prepare(fasta)
    except _USER_ERRORS as e:
        typer.secho(f"error: {e}", fg=typer.colors.RED, err=True)
        raise typer.Exit(1)
    typer.echo(f"index ready: {prepared.meta.get('cache_dir', prepared.fasta)}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
