"""Write the synthetic BAM/CRAM fixtures that tests/cram.rs reads.

    python tests/data/cram/make_fixtures.py ~/Documents/expansionscout-dev tests/data/cram

The reads are those of the development repository's tests/test_cram.py
(FMR1 moved onto a synthetic contig; phased, methylated, partial and split
reads), built with its tests/synth.py, so pysam is needed. The CRAMs are
written by htslib, as real ones are. Everything here is random sequence.
"""
import json
import re
import random
import sys
from pathlib import Path

import pysam

dev, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
sys.path[:0] = [str(dev / "tests"), str(dev)]
import synth  # noqa: E402
from expansionscout.catalog import DEFAULT_CATALOG  # noqa: E402

MOTIF = "CGG"
out.mkdir(parents=True, exist_ok=True)

ref = synth.make_reference(MOTIF)
# Soft-mask a stretch of left flank that the reads cover: a CRAM decoder
# copies matching bases from the reference, case and all.
masked = ref[:3500] + ref[3500:3900].lower() + ref[3900:]
synth.write_reference(out / "ref.fa", masked)
pysam.tabix_compress(str(out / "ref.fa"), str(out / "ref.fa.gz"), force=True)
pysam.faidx(str(out / "ref.fa.gz"))

rng = random.Random(3)
bb = synth.BamBuilder(ref, MOTIF)
for i in range(10):
    bb.spanning(f"n{i}", 30, hp=1, meth=0.1, reverse=(i % 2 == 0), noise=0.01, rng=rng)
for i in range(8):
    bb.spanning(f"p{i}", 90, hp=2, meth=0.9, reverse=(i % 3 == 0), noise=0.01, rng=rng)
for i in range(3):
    bb.left_partial(f"l{i}", 150)
    bb.split(f"s{i}", 120)
bb.write(out / "s.bam")
for p in out.glob("s.bam.unsorted.bam"):
    p.unlink()


def write_cram(path, fmt=()):
    with pysam.AlignmentFile(str(out / "s.bam")) as i, \
            pysam.AlignmentFile(str(path), "wc", template=i, reference_filename=str(out / "ref.fa"),
                                format_options=[o.encode() for o in fmt]) as o:
        for a in i:
            o.write(a)
    # htslib records the reference's absolute path as UR; keep it out of the fixture.
    with pysam.AlignmentFile(str(path), reference_filename=str(out / "ref.fa")) as f:
        text = re.sub(r"\tUR:[^\t\n]*", "", str(f.header))
    hdr = out / "header.sam"
    hdr.write_text(text)
    pysam.samtools.reheader("-i", str(hdr), str(path))
    hdr.unlink()
    pysam.index(str(path))


write_cram(out / "s.cram")
# Many containers of several slices each, so a region query has to pick
# containers from the index and read each once.
write_cram(out / "small_slices.cram", ["seqs_per_slice=4", "slices_per_container=2"])

# References the CRAM was not written against.
synth.write_reference(out / "other.fa", synth.make_reference(MOTIF, seed=99))
synth.write_reference(out / "short.fa", ref[:-100])
(out / "renamed.fa").write_text((out / "ref.fa").read_text().replace(f">{synth.CONTIG}", ">chrU"))
pysam.faidx(str(out / "renamed.fa"))

rec = next(r for r in json.load(open(DEFAULT_CATALOG)) if r["id"] == "FXS_FMR1")
rec = dict(rec, chrom=synth.CONTIG)
for b in ("hg38", "t2t", "hg19"):
    rec[f"start_{b}"], rec[f"stop_{b}"] = synth.TRACT_START, synth.TRACT_START + 3 * synth.REF_UNITS
(out / "catalog.json").write_text(json.dumps([rec]))
