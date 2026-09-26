# ExpansionScout

Targeted repeat-expansion calling from Oxford Nanopore (ONT) BAM files,
with per-allele methylation.

ExpansionScout genotypes the disease-associated tandem repeat loci catalogued in
[STRchive](https://github.com/dashnowlab/STRchive). For each locus it sizes
every read by walking its repeat motifs, assigns reads to alleles, summarises
5mC over a window chosen for that locus' biology, and reports whether the
call is trustworthy given the read lengths in the sample.

It is a research tool, not a clinical assay.

## Get it

### Download the binary

Each [release](https://github.com/martinandclaude/expansionscout/releases) has
a statically linked binary for x86-64 Linux. It needs nothing installed --
no Python, no htslib, no conda environment -- and carries the locus
catalogue, so the file is the whole tool.

```bash
gh release download --repo martinandclaude/expansionscout   # the latest release
sha256sum -c expansionscout.sha256                            # expansionscout: OK
chmod +x expansionscout
./expansionscout --version
# expansionscout 0.1.0 (STRchive v2.26.1 compiled in)
```

Put it anywhere on your `PATH`, for instance
`install -m 755 expansionscout ~/.local/bin/`. The two files can also be
downloaded from the releases page in a browser. While this repository is
private, either way needs a GitHub account with access to it, and `gh` needs
to be logged in (`gh auth login`).

### Or build it

With [rustup](https://rustup.rs) and a C compiler installed:

```bash
git clone https://github.com/martinandclaude/expansionscout
cd expansionscout
./build-static.sh      # prints the binary's path and its SHA-256
```

This builds the same static x86-64 Linux binary as a release, and it is
reproducible: a checkout of a release tag gives a file with the SHA-256
published for that release, wherever the checkout is. The Rust compiler is
pinned in `rust-toolchain.toml` and `rustup` fetches it; the final link uses
the `rust-lld` and musl startup files that come with that toolchain rather
than whatever linker the machine has. The machine's `cc` links only the
dependencies' build scripts and macros, which run during the build and are
not part of the binary.

On another platform, `cargo install --locked --path .` builds for the
machine at hand and puts `expansionscout`, and `es` as a short alias, in
`~/.cargo/bin`. Only the static Linux build is tested and released.

## Run it

Requires a coordinate-sorted, indexed BAM and the indexed reference FASTA it
was aligned to. The default catalogue coordinates are GRCh38; for an
alignment to another assembly, pass `--build` (see below). Methylation output
additionally requires the BAM to carry MM/ML tags from modification-aware
basecalling.

```bash
expansionscout call --bam sample.bam --ref hg38.fa.gz --out results/sample
```

That calls the whole catalogue. For selected loci, by identifier, alias or
gene name:

```bash
expansionscout call --bam sample.bam --ref hg38.fa.gz --out results/sample --loci FMR1 C9ORF72 HTT
```

Everything is written next to the `--out` prefix: the two tables and the
HTML report (see [Output](#output)), and with `--vcf` a VCF.

Options for `call`:

| Flag | Meaning |
|---|---|
| `--sex XX\|XY` | ploidy on the sex chromosomes (default `XX`) |
| `--mosaic` | keep minor length modes instead of folding them into alleles |
| `--ignore-hp` | ignore `HP` haplotype tags and cluster on length alone |
| `--threads N` | call loci in parallel |
| `--sample NAME` | the sample name written to the outputs (default: the BAM's read-group `SM`, else its file name) |
| `--vcf` | also write `<prefix>.vcf` |
| `--no-report` | skip the HTML report |
| `--catalog PATH` | use another catalogue, in STRchive's JSON format (default: the release compiled in, see [Locus catalogue](#locus-catalogue)) |
| `--clinical` | report and classify under each locus's clinical size convention (see [What counts as the allele size](#what-counts-as-the-allele-size)) |
| `--build hg38\|t2t\|hg19` | which assembly's catalogue interval to use; `--ref` must be that assembly (default `hg38`) |
| `--anchor`, `--flank` | reference context required around the tract |
| `--min-mapq`, `--min-clip` | read filters |
| `--min-support`, `--min-frac` | minimum reads for an allele |
| `--meth-threshold` | modkit confidence threshold for a modification call |
| `--vcf-style NAME` | which VCF dialect to write (see below) |

Other subcommands:

```bash
expansionscout loci --regime expansion     # list the catalogue
expansionscout bed --loci FMR1 --out fmr1.bed
expansionscout qc results/sample.loci.tsv  # check a run against itself
expansionscout vcf-styles --checks         # the VCF dialects and their sources
expansionscout provenance --sources        # where every clinical band comes from
```

## Output

Three files per sample, and a fourth on request.

- `<prefix>.loci.tsv` — one row per locus: per-allele unit counts, spread,
  support, structure, interruptions, per-window methylation, size class,
  lower bound, detectability, evidence state, and `tagged_fraction`, the share
  of the clustered reads carrying an `HP` tag. Phasing can be thin without
  anything else saying so, and a per-haplotype split rests on it.
- `<prefix>.per_read.tsv` — one row per read: class, unit counts by label,
  structure string, per-CpG methylation.
- `<prefix>.report.html` — the sample as a worklist, and **the default
  output**. Self-contained: it opens from disk, needs no server and fetches
  nothing, so it can be attached to a case and will still render in five
  years. See below.
- `<prefix>.vcf` — written with `--vcf`, for tools that read one. One record
  per locus; per-haplotype values are `Number=.` FORMAT fields in genotype
  order, so two alleles of equal length but different methylation stay
  distinguishable.

### The report

```bash
expansionscout call --bam sample.bam --ref hg38.fa.gz --out results/sample
```

Every locus is a row in one table, with the same columns, ordered by what
needs looking at rather than by position — so rows can be compared by reading
down. Nothing is filtered away: a locus that could not be called and a locus
called normal on reads that could not have shown an expansion are both on the
page, in their own labelled sections. An invisible exclusion in diagnostic use
is worse than a visible one.

Expanding a row gives the per-read view, in the grammar the field already
reads — one read per row, longest first, each cell a repeat unit, interruptions
in a contrasting colour. That is the shape of a REViewer pileup and of TRGT's
waterfall, and there was no reason to invent a third.

What it adds to that grammar is what this tool measures and others do not:

- a **per-CpG methylation lane** under each read, sharing the unit axis;
- **partial reads drawn as open-ended bars**, because a clipped read says "at
  least this long" and drawing it as a point estimate would be a lie;
- the **spread that made the tool decline to classify** an allele, shown
  instead of the band the median happened to land in;
- the **clinical bands with the citation each one came from** — boundary and
  name separately, since they often come from different places, with the
  jurisdiction and a verbatim quote. Published ranges are never chained into a
  contiguous partition: at *AR* the benign range ends at 34 and the
  intermediate range begins at 36, so the page says that 35 belongs to no
  published band rather than inventing a boundary;
- the **checks from `expansionscout qc`**, at the top, before any call.

`--no-report` skips it, which is what a pipeline calling many samples wants.

### VCF dialects

A repeat genotype has no single correct spelling, so the dialect is data
([`data/vcf_styles.yaml`](data/vcf_styles.yaml)) rather than a branch in the writer,
and `--vcf-style` selects it. `expansionscout vcf-styles` lists them with the
source each was written against.

| Style | Writes |
|---|---|
| `native` (default) | VCF 4.5, explicit allele sequences, everything measured. The only style carrying methylation, interruption structure and the lower bound. |
| `tr` | VCF 4.5 `<CNV:TR>` symbolic alleles with `RN`, `RUS`, `RUL`, `RUC`, `RB` and `CIRUC`, as the specification defines them. Conformant rather than informative: allele sequences are not written. |
| `conservative` | VCF 4.2, sequence alleles, only fields a general reader will accept. |

Every style is written against a published specification and none is named
after a commercial platform: a style carrying a vendor's name makes a claim
only that vendor can settle.

Adding a dialect is a few lines of YAML in
[`data/vcf_styles.yaml`](data/vcf_styles.yaml) and a rebuild, including
renaming a field for a reader that insists on its own spelling. Field definitions stay in code, so a
style can select and rename but cannot invent a header line whose `Number` or
`Type` nobody checked.

The ExpansionHunter `<STRn>` convention is not implemented. REViewer, the
obvious consumer, needs the graph-realigned BAMlet ExpansionHunter emits
alongside its VCF, which this tool does not produce, so the VCF alone drives
nothing; and its `AD_SP` / `AD_FL` / `AD_IR` are per-allele read counts split
by read class, where this tool holds read class per locus.

### Checking a run against itself

```bash
expansionscout qc results/sample.loci.tsv
```

Six checks whose expected answer is known before the data arrives, so they
work on a clinical sample with no truth set. Each finding carries its basis,
because they are not equally well founded:

| Basis | Meaning |
|---|---|
| `mechanical` | a binary fact about the file — no threshold, no judgement |
| `biology` | follows from the sample's own karyotype or a cited standard |
| `consistency` | two of our own numbers that cannot both be right |
| `heuristic` | an unvalidated cut-point of ours, flagged so it is not quoted as a standard |

`skip` is not `pass`: a run with no modification tags has not passed the
methylation checks, and the report says so.

The methylation control is the clearest example. *FMR1* promoter methylation
should follow the sample's X complement — unmethylated on a single active X
with a normal-length allele, and split between one silenced and one active
allele on two X chromosomes. On the seven GIAB samples it passes in all seven
when phased; run without haplotags it warns in exactly the three samples with
two X chromosomes, because both alleles then report one pooled average that
describes neither. No truth set is involved in either direction.

## What is distinctive

A BAM-based caller can only measure an allele it can see. ExpansionScout is built
around being explicit about that boundary rather than silent at it.

- **Clipped reads still carry evidence.** A read the aligner soft-clipped at
  the repeat bounds the allele from below. A published benchmark of seven
  maintained ONT callers found that all of them missed a known C9orf72
  expansion because the expanded reads were clipped rather than aligned
  (Aliyev et al. 2026, preprint; see references). ExpansionScout classifies reads as
  spanning, split, partial or uninformative, and reports a lower bound with
  its supporting read count. The bound is the `min_reads`-th largest per-read
  count, not the largest: a per-read count is a noisy measurement inflated by
  insertion error in a repeat tract, so the maximum of several is biased
  upward and is not a bound at all. Measured against simulated alleles of
  known size, taking the maximum overshot the truth in 12.5 % of cells; the
  order statistic overshoots in none and gives up about two points of
  tightness.
- **A normal call is not automatically a negative result.** For every locus
  whose pathogenic range can exceed read length, the tool reports
  *detectability*: the fraction of reads long enough to have spanned a
  pathogenic allele. When that is low, `negative_reliable` is 0 and the
  output says so explicitly rather than returning a confident normal call.
- **Mosaicism is described, not hidden.** Each allele carries its median
  absolute deviation, central 90 % range, and the fraction of reads in an
  upward length tail. `--mosaic` additionally reports minor length modes, and
  each allele reports the fraction of its reads that are methylated.

Methylation windows are configured per locus in gene orientation, so a
promoter window lands on the correct flank for a gene on either strand.

## How a call is made

1. **Read classification.** Each alignment record is spanning, left or right
   partial, split (two records of one read that together span), or
   uninformative. Split records are paired by read name and the tract size is
   taken from read coordinates.
2. **Per-read sizing.** The tract segment is decomposed into units by a
   rotation-invariant walk that tolerates isolated base errors. Each unit is
   labelled canonical, pathogenic, benign, interruption, or unrecognised. A
   single-base indel costs at most one unit and never derails the walk.
   Motifs longer than 12 bp are rarely read without an error in some copy,
   so their units are also accepted at 80 % identity to the run's own
   reading frame: inside the run, and at its ends on a spanning read while
   the unit lies where the read's alignment places the tract and the run is
   no longer than the tract between the read's anchors -- the length its unit
   count is taken from. Not beyond, where a real tract is often flanked by
   degenerate copies of its own repeat. Where the flank continues the tract's
   period across a whole anchor, an aligner may put the tract's insertion or
   deletion beyond it, so that length is measured from the first flank base
   that breaks the period instead. Within those bounds a unit may also
   be matched by alignment, at most 20 % of it in edits, which carries the
   run across the insertions and deletions nanopore reads carry; a 61 or
   99 bp motif with no exact copy in the read is seeded the same way; and
   the tract window takes in an expansion that an aligner has moved into
   flank bases continuing the period. The bounded run only ever adds to the
   one exact matching gives.
3. **Allele assignment.** `HP` tags when present and well supported;
   otherwise a one-dimensional Gaussian mixture with the component count
   chosen by BIC. The mixture learns the spread from the data, which matters
   because sizing noise grows with tract length: a clean 17/18 heterozygote
   splits while a 500-unit homozygote with ±30 units of spread stays one
   allele.
4. **Allele sequence.** The sequence written for each allele -- the VCF ALT,
   and the structure and repeat unit reported with it -- is a consensus of
   the allele's reads rather than one read's copy, so that no single read's
   sequencing errors are written as the allele. Each read's run is cut into
   unit-sized pieces in the allele's reading frame, and each piece of the
   allele is decided by a vote, base by base where no two reads carry the
   same one. With fewer than three whole reads, or where the reads run on
   into degenerate copies of the motif that the representative read stops
   short of, the representative read is kept: the one nearest the allele's
   size with the least noise. On simulated nanopore reads over 77 real
   loci, reference alleles written as reference went from 692 of 2,256 to
   1,134, and variant alleles written exactly from 286 of 1,106 to 492.
   The rules are in [`src/consensus.rs`](src/consensus.rs).
5. **Methylation.** Per read from MM/ML tags, summarised over the tract and
   over gene-5' and gene-3' flank windows configured per locus.
6. **Classification.** Curated clinical size bands where they exist, catalogue
   thresholds otherwise, plus the evidence state from partial reads.

## Locus catalogue

The catalogue is [STRchive](https://github.com/dashnowlab/STRchive), vendored
under [`data/strchive/`](data/strchive/VERSION) pinned to a release
(currently v2.26.1) by commit and SHA-256, with its licence and citation. It supplies coordinates, gene strand, motifs in both
orientations, documented interruptions and clinical thresholds for 82 loci.

That release is compiled into the binary and is the default.
`--catalog PATH` runs against any other file in the same format -- a newer
STRchive release, or a laboratory's own edit of one -- without rebuilding
anything. Either way the run records which catalogue it used, by content
rather than by path: the report's provenance and the VCF header carry
`STRchive v2.26.1 (bundled), sha256:...` for the vendored release, wherever
the file was read from, and the path and SHA-256 of anything else.

One thing to know before editing motifs: above 12 bp, a unit is also matched
by identity to the run's reading frame -- its first exact unit, or where a
read holds none, the first catalogue rotation that fits -- and takes that
frame's class. So if a long motif is given more than one class of the same
length -- canonical and benign, say -- the classes are told apart only in
exact copies, and a noisy unit is counted with the frame's. The bundled
release has one class per long motif.

A small overlay, [`data/overlay.yaml`](data/overlay.yaml), adds only what the
catalogue does not carry: methylation windows and their clinical meaning,
which interruptions to score, aliases, and clinical size-band labels. It is
compiled in as well.

Adding a locus usually means adding nothing, because it is already in the
catalogue; one that is not can be added to a catalogue file given with
`--catalog`. An overlay entry is needed only to give a locus a methylation
window, interruption motifs, or finer size bands, and since the overlay is
compiled in, that means rebuilding.

### The interval belongs to one assembly

Every catalogue entry carries a separate interval per assembly, because the
tract itself is not the same in each. At `XYLT1` the repeat lies inside 241 bp
that GRCh38 does not contain and CHM13 does, so the annotated interval is
15 bp on GRCh38 and 93 bp on CHM13; at `HTT` the two assemblies simply carry
different alleles, 19 and 28 CAG. `--build` selects which interval is read,
and `--ref` must be a FASTA of that same assembly.

Nothing checks this for you at the level of a single locus, because a
mismatched pair does not error: the interval lands on unrelated sequence and
the locus quietly returns no spanning reads. `call` therefore tests the whole
selection at once and warns if an implausible share of intervals contain none
of their own motif in the reference given.

Clinical thresholds are a separate matter. They were derived once, on one
assembly, and the catalogue does not record which — its `ref_copies` field
tracks GRCh38. A threshold is only comparable with a measurement made on the
assembly it came from.

### Gene strand versus reference strand

Motifs used for sequence matching are on the **reference plus strand**;
motifs shown to a reader are on the **gene strand**. A BAM stores read
sequence plus-strand regardless of read orientation, so for a minus-strand
gene the reference carries the reverse complement. Getting this wrong is
silent: the locus finds no reads rather than raising an error. The catalogue
supplies both orientations, so nothing here is derived by hand.

The same applies to methylation windows, which is easy to get wrong. Windows
are declared **gene-5'** and **gene-3'** and mapped to reference sides by
strand. For C9orf72 and DMPK, both minus-strand genes, gene-5' is the
higher-coordinate flank. C9orf72 has a CpG island on each side of the repeat:
the gene-5' island at chr9:27,573,761-27,573,989 is the one reported as
hypermethylated on expanded alleles, while the lower-coordinate island is the
exon-1b promoter.

### What counts as the allele size

The annotated tract and the clinically quoted number are not always the same,
and the difference is a convention rather than an error.

An FMR1 allele of "30 CGG with 2 AGG interruptions" is thirty units in total:
the interruptions are counted. An HTT allele is quoted as its pure CAG count,
with the CAACAG excluded, so a tract annotated as 19 units is reported as 17.
Both numbers describe the same molecule and answer different questions.

By default `expansionscout` reports what it measured across the annotated tract.
`--clinical` reports and classifies under each locus's own convention, which
is declared per locus in the overlay. The measured count is kept either way,
in the `*_median_measured` columns, so the two are never confused and neither
is silently substituted.

This matters when comparing against other genotypers: several report the
pure repeat count by default, so a difference of a unit or two at an
interrupted locus is a difference of convention, not of accuracy.

### Clinical size bands

The catalogue's benign / intermediate / pathogenic split is coarser than
clinical vocabulary in places. It would report a 120-CGG FMR1 allele as
"intermediate", where a clinical reader expects "premutation".

At FMR1, HTT and DMPK the overlay draws sub-bands the catalogue does not:
gray zone and premutation, reduced versus full penetrance, and the
congenital range. At ATXN1 and C9orf72 it only renames the catalogue's
classes into the usual vocabulary, leaving the boundaries untouched.

These sub-band boundaries follow widely used clinical conventions, but they
are not themselves sourced from the catalogue. Check them against your own
laboratory's validated cut-points before relying on the labels.

## Scope and limits

Read this before trusting any output.

- **Alleles longer than the reads cannot be measured.** No method working on
  aligned sequence can size a tract no read crosses. What this tool adds is
  saying so per sample, through detectability and the partial-read lower
  bound, instead of returning a confident normal call.
- **A normal-range call is only negative evidence when detectability is
  high.** Reads carrying a long expansion often fail to align through it, so
  the result can look like a clean homozygous normal call. Check
  `negative_reliable` and `DETECT` before reading any normal result as a
  negative. Where they are low, the sample has not been excluded.
- **Somatic mosaicism is described, not quantified.** Spread, tail fractions
  and minor modes describe a distribution; none of them is a mosaic fraction.
- **Methylation depends entirely on modification-aware basecalling.** Without
  MM/ML tags every methylation column is `.`. Values are means of per-read
  modification probabilities, not calls.
- **Interpretation is locus-specific.** Composition loci such as RFC1 turn on
  which motif is expanded, not on length; at ATXN1 a long interrupted allele
  can be benign while a shorter pure one is not. Treat these as screens.
- **ADTKD-MUC1 is out of reach.** Its disease allele is one extra C in one
  unit of a VNTR kilobases long, which a length measurement cannot see; a
  result at that locus is the VNTR's length, not a finding for the disease.
  The catalogue also gives it a 61 bp motif, its pathogenic unit, where the
  repeat's own is 60 bp, so its counts run about 1.6 % low. That is read as
  published and recorded (`expansionscout provenance --discrepancies`).
- **Confirm every clinical threshold** against your own laboratory's
  validated cut-points. The size-band labels here follow common convention
  and are not themselves a validated standard.

This is a research tool. It is not a clinical assay, and no output is a
diagnostic claim.

## Reproducibility

A result depends on the binary and the inputs, and nothing else. The binary
has no runtime dependencies and carries the catalogue release, the clinical
overlay and the VCF dialects it was built with, and `--version` names the
catalogue. The build is reproducible: the same commit, `Cargo.lock` and
pinned toolchain, linker included, give the same SHA-256 in any checkout.
CI builds every commit twice, in two checkouts, and fails unless the two
match; a release is built a third time and must match again. So a result
can be tied to one file, and that file checked by rebuilding it. See Heng
Li, [The AI rewrite dilemma](https://lh3.github.io/2026/04/17/the-ai-rewrite-dilemma),
for why this matters.

`SOURCE_DATE_EPOCH`, when set, fixes the date printed on the report, so a
report can be reproduced byte for byte too.

The binary is a port of a Python reference implementation, which is where
the method is developed
([`martinandclaude/expansionscout-dev`](https://github.com/martinandclaude/expansionscout-dev)).
The two are required to agree byte for byte on every output: every table,
the VCF in each dialect, the report, `qc`, and the messages on standard
error. The reference's parity tests check this against the commit of this
repository it pins, on synthetic samples that vary read classes, strands,
clipping, phasing, sequencing noise, modification tags including malformed
ones, gzipped and plain references, BAI and CSI indexes, contig naming,
thread counts and catalogues, and on thousands of random cases built to
reach the edge cases of allele assignment and of the allele-sequence
consensus. Matching byte for byte meant pinning arithmetic that is easy to
leave implicit: the mixture model's `exp` and `log` come from the pure-Rust
`libm`, which gives the same bits on every CPU, where numpy's depend on
whether the processor has AVX-512.

## Status

The tests here are unit tests of the port, run with `cargo test`; the
engine's test suite and the parity tests live with the reference
implementation. CI (`.github/workflows/ci.yml`) runs `cargo fmt`, `clippy`
and the tests on the pinned toolchain, and the reproducible static build, on
every push and pull request, and publishes a release for a tag `vX.Y.Z`.

**It has not been benchmarked against orthogonal measurements.** Concordance
against PCR fragment analysis, repeat-primed PCR or Southern blot on samples
of known size has not been established. Treat the output accordingly.

## Licence

MIT, see [`LICENSE`](LICENSE).

That covers this project's own code. The bundled locus catalogue is separate
third-party content under CC BY 4.0 and carries its own attribution
requirement; see [`NOTICE.md`](NOTICE.md).

## Attribution

This tool bundles the STRchive locus catalogue by Hiatt et al., redistributed
under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). If you use
`expansionscout`, cite STRchive. Full provenance, the pinned commit and the citation
are in
[`NOTICE.md`](NOTICE.md),
[`data/strchive/VERSION`](data/strchive/VERSION) and
[`data/strchive/LICENSE.STRchive`](data/strchive/LICENSE.STRchive).

## References

- Hiatt L, Weisburd B, Dolzhenko E, Rubinetti V, Avvaru AK, VanNoy GE,
  Kurtas NE, Rehm HL, Quinlan AR, Dashnow H. STRchive: a dynamic resource
  detailing population-level and locus-specific insights at tandem repeat
  disease loci. *Genome Medicine* 2025;17(1):29.
  [doi:10.1186/s13073-025-01454-4](https://doi.org/10.1186/s13073-025-01454-4)
- Aliyev E, Avvaru A, De Coster W, et al. A comprehensive assessment of
  tandem repeat genotyping methods for Nanopore long-read genomes. *bioRxiv*
  2026. [doi:10.64898/2026.02.28.708646](https://doi.org/10.64898/2026.02.28.708646)
- Benarroch L, Pesovic J, Rossato M, et al. Toward the clinical application
  of long-read sequencing in repeat-expansion disorders. *Nature Genetics*
  2026.
- Russ J, Liu EY, Wu K, et al. Hypermethylation of repeat expanded C9orf72 is
  a clinical and molecular disease modifier. *Acta Neuropathologica*
  2015;129:39-52.
  [doi:10.1007/s00401-014-1365-0](https://doi.org/10.1007/s00401-014-1365-0)
- Udine E, et al. Targeted long-read sequencing to quantify methylation of
  the C9orf72 repeat expansion. *Molecular Neurodegeneration* 2024;19:99.
  doi:10.1186/s13024-024-00790-0
- Barbe L, Lanni S, Lopez-Castel A, et al. CpG methylation, a parent-of-origin
  effect for maternal-biased transmission of congenital myotonic dystrophy.
  *American Journal of Human Genetics* 2017;100:488-505.
  [doi:10.1016/j.ajhg.2017.01.033](https://doi.org/10.1016/j.ajhg.2017.01.033)
