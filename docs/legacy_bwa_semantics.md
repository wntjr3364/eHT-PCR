# Legacy bwa / SAM / XA semantics (the spec the Python engine must reproduce)

> Confirmed by decompiling the legacy `eHT-PCR_v8.0.1.jar` (`Linux_eHT_PCR_v08`,
> `eHT_PCR`, `PrimerPair`, `Primer`, `Hit`, `NoAscCompare`, `FastaIndex`) and
> cross-checking the legacy Python (`eHT-PCR.py`, `Primer.py`). Items marked
> **DIVERGENCE** are places the Python engine deliberately differs from the jar.

## Coordinate convention

- **Internal `Region`/`Amplicon`: 0-based, half-open** (`seq[start:end]`,
  `length == end - start`). Only convention in code + cross-process JSONL.
- Boundary conversions: GFF & `region` CLI (1-based inclusive) → `start0 = s-1,
  end0 = e`. SAM `POS` / `XA` pos (1-based) → convert at the parser. The jar
  prints **1-based** (`start1 = fwdHit.pos`, `end1 = revHit.pos`); the regression
  comparator normalizes **both** sides to 0-based half-open before set-comparison
  (`start0 = fwdHit.pos - 1`, `end0 = revHit.pos`).

## Flag-naming trap (read first)

The **jar CLI** uses `-i PrimerList -o outDir -d DB -n <maxProductSize> -m <maxMismatch>`.
Inside, the jar's `-m` (max mismatch) is passed to **`bwa aln -n`**, and the jar's
`-n` is the **max PCR product size** (the pairing ceiling). The `-m 200000000` in
the bwa command is **bwa's own** `-m` (queue-size cap), unrelated to the jar's
`-m`. In `ehtpcr`, `BwaAlnParams.max_mismatch` → `bwa aln -n`;
`BwaAlnParams.max_product_size` → the pairing ceiling.

## External commands (confirmed verbatim from `eHT_PCR.runBWA`)

```bash
# bwa index is NOT run by the jar — the legacy Python wrapper (eHT-PCR.py) ran it
# when .amb/.ann/.bwt/.pac/.sa are missing/empty, with the IS algorithm:
bwa index -a is <db>

# threads are HARDCODED to 32 in the jar (not derived from any flag):
bwa aln -t 32 -N -i 30 -d 30 -m 200000000 -l 100 -n <maxMismatch> <db> <out>/primer.fastq \
    1> <out>/primer.sai 2> <out>/primer.aln_log
bwa samse -n 300000000 <db> <out>/primer.sai <out>/primer.fastq \
    1> <out>/primer.sam 2> <out>/primer.samse_log
```
Run via `Runtime.exec(["/bin/bash","-c", <string>])`.
**DIVERGENCE:** `ehtpcr` makes `-t` = `runtime.threads` (jar hardcodes 32). Thread
count can change bwa-aln tie-break ordering → a real parity risk; the comparator
must be order-insensitive AND the oracle should be run at a pinned thread count.

## Invariants (each a tested contract)

1. **XA-drop, not truncate.** `bwa samse -n N`: if a read has > N hits the **whole
   `XA` tag is omitted**. `ehtpcr` keeps the same huge `-n 300000000`, so overflow
   is unreachable in practice; a missing XA is therefore "no other hits". (A
   hit-count ceiling that fails loud is NOT actively enforced — it relies on the
   huge `-n`; revisit if `-n` is ever lowered.)
2. **`-N` (find all hits) is required.** Without it bwa returns only the best hit.
3. **`-l 100` disables seeding** for 18–24 bp primers (seed > query). Reproduce.
4. **`-n` (→ `max_mismatch`) is an integer edit distance.** Float ⇒ a fraction.
5. **Strand decoding has two sources:**
   - **Primary hit**: strand from **SAM FLAG bit 0x10** (`flag & 16`). If set,
     strand `-` and the engine reverse-complements `SEQ` to recover the oligo.
   - **XA hit**: strand from the **sign of the XA position** (negative ⇒ `-`,
     `abs` ⇒ 1-based leftmost). XA tag = token starting `XA`; strip the 5-char
     `XA:Z:` prefix; split on `;`; each record split on `,` → `(rname, signedPos,
     CIGAR, NM)`; the jar uses only `rname` and `signedPos`.
6. **Strand-dependent hit-position normalization** (`Primer.addHit`): for `+`,
   `Hit.pos = pos` (1-based leftmost); for `-`, `Hit.pos = pos + primerLen - 1`
   (1-based rightmost). This is what makes the product-size formula work.
7. **bwa index cache key includes the bwa version** (`hash(fasta)+bwa_version`);
   record the bwa version string in `manifest.json`.
8. **Per-target/job tempdir** for `primer.fastq`/`.sai`/`.sam` (jobs > 1 collide).

## primer.fastq construction (`fastaToFastq` + `extractUniqPrimer`)

- Queries = the **deduplicated set of all oligos** (every forward `f_seq` and every
  reverse `r_seq`) collected in a `Hashtable` (so order is hash-order). bwa is run
  **once per unique oligo**; hits are then keyed back **by sequence**.
  → **Parity must key hits by oligo sequence, not by read index.**
- `r_seq` is **already reverse-complemented** upstream by `Primer.py`, so the query
  is the oligo as synthesized. `records.PrimerPair.reverse.seq` follows this.
- Record: `@<i>\n<SEQ>\n+\n<QUAL>` where `<i>` = array index, `<QUAL>` = the first
  `len(SEQ)` chars of a fixed 48-char `IIII…` (Phred Q40) string. **Oligos > 48 bp
  would crash the jar** (substring on a 48-char string); fine for ≤24-mers.

## Mismatch string (`getMisValue`)

- Reference window via `FastaIndex.getSequence(name, pos, pos+len-1)` (1-based).
- Per base: `'-'` (match) / `'#'` (mismatch). `+`: `primer[i]` vs `target[i]`. `-`:
  `rcNucl[primer[i]]` vs `target[len-1-i]`. `rcNucl` defines only A/T/G/C/N → any
  other base complements to `\0` (always mismatch).
- **N rule:** if the target window has **≥ 5 N's, `getMisValue` returns null → the
  hit is silently DROPPED.** Reproduce or mark **DIVERGENCE** with a test.
- `fmisvar`/`rmisvar` in the output are these `-`/`#` strings.

## Amplicon pairing + FR/RF (`PrimerPair.printParis`)

For each contig present in the forward primer's hit map, nested over
`i14 ∈ {F,R}` (the `+`-strand member) and `i15 ∈ {F,R}` (the `−`-strand member):
- `fwdHits = primer[i14].hits(name, "+")`, `revHits = primer[i15].hits(name, "-")`,
  each sorted ascending by `pos` (`NoAscCompare`).
- `direction = {F,R}[i14] + {F,R}[i15]` → the jar generates **all four**
  `FF/FR/RF/RR`. The **Python wrapper keeps only `FR`/`RF`**; `--nomatch` keeps only
  the `#N0` sentinel.
- For each `(fwdHit, revHit)`: **`productSize = revHit.pos − fwdHit.pos + 1`**
  (inclusive span using the strand-normalized positions: fwd = leftmost, rev =
  rightmost).
- **Keep iff `min_size ≤ productSize ≤ max_size`** (both inclusive), where
  **`min_size = fPrimer.len + rPrimer.len`** (the pair's two primer lengths,
  fixed) and `max_size` = the jar's `-n` (`max_product_size`).
- Output (`result_csgm`, 1-based): `name, direction, fwdHit.pos, revHit.pos,
  fmisvar, rmisvar` → `start0 = fwdHit.pos − 1`, `end0 = revHit.pos`,
  `strand_class = direction`.
- **No de-duplication in the jar.** The same amplicon can be emitted twice (e.g.
  via FR and RF, or via a primary hit and an XA hit at the same locus).
  **DIVERGENCE DECISION (Phase 4/5):** `ehtpcr` de-duplicates amplicons. Because
  the legacy "unique" notion is a *count* of kept rows, de-duplication can change
  the pair-level unique decision → the parity gate must compare the **retained
  primer-pair set** (not raw amplicon multiset) and the dedup behavior must be a
  documented divergence with its own test.
- `#N0`: a pair that produced **zero** kept rows across all contigs/directions
  gets a single `#N<flag>` sentinel (here `flag = 0`) written to the uniq writer.

## FastaIndex (`.idx` / `.seq`) — the jar's own coordinate source

- `.idx` lines: `name\tlength\toffset` (`[length, offset]` as ints).
- `.seq`: concatenated, **UPPERCASED**, newline-stripped sequence of all records in
  file order. → soft-masking is destroyed (**DIVERGENCE:** `ehtpcr`/pyfaidx
  preserves case).
- `getSequence(name, start, end)`: if `length < end` return null; seek
  **`offset + start − 1`** (proves `start` is **1-based**) and read a **fixed
  50 bytes** regardless of `end` (`end` only bounds-checks). Window `[start,
  start+49]` 1-based; oligos > 50 bp read garbage. Uses its **own** FastaIndex,
  not samtools — so jar coords and bwa coords are both 1-based and compose without
  conversion inside the jar.

## Enumerated divergences (Python deliberately != jar)

| # | Jar behavior | ehtpcr behavior | Why |
|---|--------------|-----------------|-----|
| D1 | revcomp `M → "L"` (`Sequence.py:95`) | `M → K` (correct IUPAC) | bug fix |
| D2 | `.upper()` everywhere (FastaIndex, getSequence) | case preserved (pyfaidx) | keep soft-mask info |
| D3 | XA overflow → silent under-count | rely on huge `-n` (unreachable); not enforced | overflow can't occur in practice |
| D4 | no amplicon de-duplication | de-duplicate amplicons by span | clearer; gate on pair set |
| D5 | `bwa aln -t 32` hardcoded | `-t = runtime.threads` | configurable (pin for oracle) |
| D6 | ≥5-N window → hit dropped | **replicated by default** (`drop_n_rich_windows: true`, threshold 5); optional off | jar-faithful default; opt-out for N-rich draft genomes |

## Known limitations (shared with the jar)

- **Gapped hits get an approximate reference span.** `specificity/bwa_aln.py` does
  not parse the SAM/XA CIGAR; it derives an amplicon's reference end from the primer
  length (`start0 + len`), which is exact only for an ungapped, full-length match.
  `bwa aln` defaults to `-o 1` (one gap open), so a rare indel-containing off-target
  hit gets a `product_size` (and span-dedup key) off by the net indel length — which
  can flip the size filter at the `max_product_size` boundary. The legacy jar makes
  the identical ungapped assumption (`getMisValue` is positional), so this is a
  shared latent limitation, not a port divergence. For 18–24 bp primers with a small
  integer `-n` and `-l 100` (seeding disabled), gapped hits are rare; if it matters,
  pass primers that don't tolerate indels or lower `-n`. (`target_fasta`'s bwa
  locator, by contrast, *does* read the CIGAR — `_cigar_stats` — so this is specific
  to the specificity engine.)
- **On-target window vs a differently-coordinated specificity reference.** For
  `locus`/`region`/`target_fasta` modes, on-target overlap is checked in the design
  reference's coordinates. If a separate `specificity.reference` shares seqid *names*
  but uses a *different* coordinate system (e.g. a genomic spec reference against a
  transcriptomic design reference), the window check is not meaningful — use a spec
  reference in the same coordinate system as the design reference (the curated
  same-assembly case, which is the supported workflow), or `name` mode (whose
  seqid-only on-target check is coordinate-agnostic).
