# eHT-PCR

**English** | [한국어](README.ko.md)

Design PCR primers that amplify **one** locus and nothing else.

eHT-PCR is a command-line tool for *in-silico* (electronic) high-throughput PCR.
You give it a **reference** and a **target**; it enumerates candidate primer pairs,
aligns every primer back against the reference, and keeps only the pairs whose PCR
product is **unique** — pairs that won't also amplify a paralog, a homeolog, or
some other off-target. It's organism-agnostic: the reference can be a genome, a
transcriptome, or any multi-FASTA.

## How it works

For each target, eHT-PCR (1) **extracts** the target sequence, (2) **enumerates**
candidate primer pairs under length / GC / Tm / GC-clamp / product-size filters,
(3) **maps** every primer to a specificity reference with `bwa` and pairs
forward/reverse hits into amplicons, and (4) **keeps** the pairs that yield a
single on-target amplicon. "Unique" means *one* amplicon — not just a primer that
looked good on its own.

## Install

The bundled conda environment brings eHT-PCR and the external tools it can use
(`bwa`, BLAST+, miniprot):

```bash
conda env create -f environment.yml
conda activate ehtpcr
ehtpcr --help
```

Or install just the package into an existing environment (Python ≥ 3.10):

```bash
pip install .
```

Every Tm model then works out of the box — the nearest-neighbor models and
`primer3` ship as Python dependencies. The alignment steps, however, call
external binaries you put on your `PATH`: **`bwa`** for the specificity check
(needed unless you run design-only with `--no-spec`), plus **BLAST+** and/or
**miniprot** *only* if you locate targets by sequence (`--query-fasta`).
`conda install -c bioconda bwa blast miniprot` installs all three.

## Quick start — a complete run

eHT-PCR ships a tiny reference, `examples/mini.fa`, with two genes:

```
>geneX example gene geneX
TCGCTGCTGTCGGACTCCTAGTTACGTGGCGTTGCTCCACAGGTAGCCTGCCGTCGTGGTCCGCAACACT
CGCACGCTGTTTCAGGGCGATCCTCCGGATAACACCACCTCCACAAACGAAGACAACCCTCTGGTTCTTT
...
>geneY example gene geneY
CTTCTTTAGGCGAGAGTACCCTATTTTTGGCCCTATGAGCGCCTTGATGGACTCGTTACTTGGGACCAAT
...
```

The bundled `config/example.yaml` designs primers for every header containing
`gene` (so both genes) and checks each pair for uniqueness against the same file:

```bash
ehtpcr run --config config/example.yaml --out result/
```

`result/candidates.tsv` has one row per candidate pair. Your usable primers are the
**`kept`** rows — each amplifies a single product (`n_valid_amplicons == 1`) on its
target (a few columns shown):

```
target_id  pair_id             f_seq                    r_seq                  product_size  n_valid_amplicons  unique  keep_decision
geneX      geneX:3-26:132-153  CTGCTGTCGGACTCCTAGTTACG  TCTTACGGACGGGAAAGAACC  150           1                  True    kept
geneX      geneX:3-26:132-154  CTGCTGTCGGACTCCTAGTTACG  GTCTTACGGACGGGAAAGAACC 151           1                  True    kept
...
```

`result/` also gets `amplicons.tsv` (where each pair binds), `manifest.json`, and
`run.yaml`. To do real work, point the reference at your own genome or
transcriptome and choose a target — read on.

## Usage

In the examples below, `reference.fa` is your own FASTA — a genome, a
transcriptome, anything. Every run has the same shape: a reference, a way to pick
the target, an output directory.

```bash
ehtpcr run -r reference.fa  <how to pick the target>  -o result/
```

You can put everything in a YAML config (`--config`) instead of flags; flags
override the config. `ehtpcr run --help` lists every flag, grouped by reference /
target / design / runtime.

### 1. Pick your target

Four ways to say *what* to design primers for — pick the one that matches what you
have. The mode follows the flag you use (`--name` is the default), so you rarely
need `--target-mode`. If your selection matches **nothing** (usually a typo), the
run errors out instead of quietly writing an empty result; pass
`--no-fail-on-no-targets` for batch pipelines where empty is acceptable.

**By name** — you know the gene/transcript identifier (`--name`, the default
mode). This pulls **every** FASTA entry whose header *contains* the string, so all
isoforms of a gene come along:

```bash
ehtpcr run -r reference.fa --name Glyma.08G110400 -o result/
```

If the FASTA has `>Glyma.08G110400.1 …` and `>Glyma.08G110400.2 …`, both become
targets. Repeat `--name` for several genes:

```bash
ehtpcr run -r reference.fa --name Glyma.08G110400 --name Glyma.05G001200 -o result/
```

**By coordinates** — you know where it is (`--region`). Coordinates are 1-based
and inclusive, like a genome browser; the `:+`/`:-` strand is optional:

```bash
ehtpcr run -r genome.fa --region "Chr08:110340-112000:+" -o result/
```

**By gene ID + GFF** — you have a gene model in a GFF3 (`--locus` + `--gff`).
eHT-PCR finds the feature with that `ID=` and uses its span:

```bash
ehtpcr run -r genome.fa --gff genes.gff3 --locus AT1G01010 -o result/
```

**By sequence** — you only have the target's sequence, not its name or
coordinates (`--query-fasta`). eHT-PCR locates it in the reference by alignment and
**stops with an error if it matches more than one place** — that ambiguity is
something you want to be told about:

```bash
ehtpcr run -r genome.fa --query-fasta my_gene.fa -o result/
```

The query can be DNA or protein; the aligner is chosen automatically — see
[Locating a query sequence](#3-locating-a-query-sequence-target_fasta).

### 2. Check specificity (and tune it)

By default, uniqueness is checked against the design reference. But you often want
to check against a **different** file — the original soybean workflow designs on a
transcriptome and checks against a *curated* genome with known duplicate genes
removed. Pass `--spec-reference`:

```bash
ehtpcr run -r Gmax.transcript.fa --name Glyma.08G110400 \
  --spec-reference Gmax_curated.fa --max-mismatch 1 -o result/
```

`--max-mismatch` is the main specificity knob — the number of mismatches `bwa`
tolerates when looking for off-target binding sites (the original workflow used
`1` for tight gene clusters and `5` elsewhere). To tune it **per target**, use a
config file:

```yaml
# myrun.yaml
reference:
  fasta: Gmax.transcript.fa
target:
  mode: name
  name: [Glyma.08G110400, Glyma.05G001200]
specificity:
  reference:
    fasta: Gmax_curated.fa     # separate / curated reference (optional)
  keep: unique                 # unique | all | nomatch
  params:
    bwa_aln:
      max_mismatch: 5
      overrides:               # stricter mismatch limit for hard clusters
        Glyma.08G110400.1: 1
runtime:
  threads: 8
  jobs: 4                      # process targets in parallel
```

```bash
ehtpcr run --config myrun.yaml -o result/
```

`keep` chooses which pairs are written as `kept`: `unique` (the default) keeps only
pairs with a single on-target amplicon, `all` keeps every candidate with its verdict,
and `nomatch` keeps only pairs that produced no amplicon at all.

To skip specificity entirely and just design primers, add `--no-spec`.

### 3. Locating a query sequence (`target_fasta`)

When you give a query sequence (`--query-fasta` / `target.mode: target_fasta`),
`--locator` controls how it's aligned to the reference. The default is `auto`,
which sniffs the query and uses `blastn` for DNA, `tblastn` for protein (it picks
one locator for the whole file, so don't mix DNA and protein records — set
`--locator` explicitly if you must). Duplicate query names are rejected.

| `--locator` | for a query that is… | notes |
|-------------|----------------------|-------|
| `auto` | anything (default) | DNA → `blastn`, protein → `tblastn` |
| `bwa` | DNA, near-identical to the reference | no extra dependency |
| `blastn` | DNA, possibly divergent | sensitive + paralog-safe; also handles a transcript vs a genome |
| `minimap2` | a spliced transcript vs a genome | intron-aware; **not** reliable at flagging paralogs |
| `tblastn` | a protein vs a genome | 6-frame translation; paralog-safe |
| `miniprot` | a protein vs a genome | splice/frameshift-aware; **not** reliable at flagging paralogs |

```bash
ehtpcr run -r genome.fa --query-fasta my_protein.fa --locator tblastn -o result/
```

Use `blastn`/`tblastn` when telling paralogs apart matters; reach for
`minimap2`/`miniprot` only when you need a precise spliced gene structure.

### 4. Tune the design

Every candidate primer must pass these filters. The defaults suit most
genotyping / qPCR designs; override any of them on the CLI or in the `design:`
block of a config.

| Filter | Flags | Default | What it does |
|--------|-------|---------|--------------|
| Length | `--len-min` / `--len-max` | 18–24 nt | primer length window |
| GC content | `--gc-min` / `--gc-max` | 0.30–0.70 | fraction of G/C bases (too low → weak binding, too high → hard to denature) |
| Tm | `--tm-min` / `--tm-max` | 57–62 °C | melting temperature — see [Melting temperature](#melting-temperature) |
| Product size | `--product-min` / `--product-max` | 150–250 bp | amplicon length |
| 3′ GC-clamp | `--gc-clamp` | 2 | G/C bases required at the 3′ end to stabilize priming (0 = off) |

The forward and reverse primer of a pair must additionally be **well matched**:
by default their Tm may differ by at most 1 °C (`design.pair.tm_diff.max`), and
`design.pair.gc_diff.max` bounds their GC difference the same way. Tighten these
(in a config) for a uniform set of primers, or loosen them for more candidate pairs.

## Output

The output directory holds:

- **`candidates.tsv`** — one row per designed primer pair: both sequences, their
  coordinates, Tm/GC, how many valid amplicons it produced (`n_valid_amplicons`),
  whether it's `unique`, and the `keep_decision`. **Your hits are the rows with
  `keep_decision == kept`** (i.e. `unique == True`).
- **`amplicons.tsv`** — one row per amplicon found, with its location and the
  per-primer mismatch counts. This is the evidence behind the `unique` flag.
- **`manifest.json`** — eHT-PCR and `bwa` versions, reference checksums, the exact
  command line, a timestamp, and any per-target mismatch overrides.
- **`run.yaml`** — the merged effective config, enough to reproduce the run.

## Performance

Specificity (the `bwa` step) is the expensive part. Two knobs: `--threads`
(threads per `bwa` call) and `--jobs` (targets processed in parallel — the main
lever when you have many targets). The output is identical regardless of `--jobs`.

Indexes (`bwa`, BLAST+, miniprot) and the FASTA `.fai` / GFF databases are built
once and cached — keyed by each reference's path, size, and mtime — so repeated
runs against the same reference reuse them. Everything lives under one cache root:
`$XDG_CACHE_HOME/ehtpcr` (or `~/.cache/ehtpcr`). If that's read-only — common on
containers and shared HPC nodes — point the whole cache somewhere writable with
`--cache-dir /scratch/you/ehtpcr` (or `export EHTPCR_CACHE_DIR=...`); otherwise the
run stops with a clear "cache directory is not writable" error.

## Melting temperature

A nearest-neighbor Tm depends on salt and primer concentration; eHT-PCR fixes both
at **50 mM Na⁺ / 200 nM** for every model so the numbers are comparable. The
default `legacy` model is a verbatim port of the original tool's calculator
(Breslauer/Borer parameters), kept as the default so existing parameter sets
behave identically. For new designs, **`primer3`** (the de-facto-standard
SantaLucia 1998 calculation) or `santalucia2004` is more accurate; `breslauer`,
`sugimoto`, `santalucia`, `wallace`, and `gc` are also available via
`--tm-model` / `design.tm.model`. The legacy model reads a few °C hotter than the
SantaLucia models, so re-check your Tm window if you switch.

## Notes

**How specificity is judged.** Uniqueness is decided by *alignment*, not
thermodynamics: each primer is mapped with `bwa aln` allowing up to `max_mismatch`
mismatches counted uniformly. Real PCR is dominated by **3′-terminal** mismatches
(a primer that mismatches at its 3′ end usually won't extend), which a plain
mismatch count doesn't capture — so read the amplicon list as "places these
primers could plausibly bind," tune `max_mismatch` per target, and design against
a curated specificity reference. The engine is pluggable, so a thermodynamic
backend can be added later.

**On-target scope in `name` mode.** In `name` mode an amplicon counts as
*on-target* if it lands on a matched entry (whole-sequence, no sub-coordinate
check) — the right semantics for a transcriptome, where a transcript is the unit.
If you select a whole **chromosome** by name on a genome, a product elsewhere on
that chromosome would be treated as on-target; use `region`/`locus` for sub-contig
targets on a genome.

**Coordinates.** Internally everything is 0-based half-open, and the TSV columns
say so (`start0`, `end0`). Inputs that are conventionally 1-based (`--region` and
GFF features) are converted at the boundary, so you write coordinates the way you'd
read them in a browser.

## License

eHT-PCR's own code is MIT — see [LICENSE](LICENSE). Note that one core
dependency, [`primer3-py`](https://github.com/libnano/primer3-py) (used by the
`primer3` Tm model), is **GPLv2**; installing eHT-PCR pulls it in. If you need a
fully permissive dependency set, you can avoid the `primer3` model — every other
Tm model (the default `legacy`, the Biopython nearest-neighbor models, Wallace,
GC) relies only on permissively-licensed packages.
