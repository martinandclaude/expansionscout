//! CRAM input: the same reads as a BAM must give the same output, and a
//! reference the CRAM was not written against must be refused by name.
//!
//! The fixtures in tests/data/cram are synthetic (see make_fixtures.py there):
//! phased, methylated, partial and split reads at FMR1 moved onto a random
//! contig, written by htslib as a BAM and as CRAMs.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cram").join(name)
}

fn scratch(label: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("expansionscout-cram-{}-{label}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn call(flag: &str, aln: &str, reference: &str, out: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_expansionscout"))
        .arg("call")
        .arg(flag)
        .arg(data(aln))
        .arg("--ref")
        .arg(data(reference))
        .arg("--catalog")
        .arg(data("catalog.json"))
        .arg("--out")
        .arg(out)
        .args(["--vcf", "--no-report"])
        .output()
        .unwrap()
}

/// (loci.tsv, per_read.tsv, VCF without its ## header lines).
fn outputs(flag: &str, aln: &str, reference: &str, label: &str) -> (String, String, String) {
    let out = scratch(label).join("s");
    let run = call(flag, aln, reference, &out);
    assert!(run.status.success(), "{}", String::from_utf8_lossy(&run.stderr));
    let read = |ext: &str| std::fs::read_to_string(format!("{}.{ext}", out.display())).unwrap();
    let vcf = read("vcf")
        .lines()
        .filter(|l| !l.starts_with("##"))
        .collect::<Vec<_>>()
        .join("\n");
    (read("loci.tsv"), read("per_read.tsv"), vcf)
}

fn refused(aln: &str, reference: &str, label: &str) -> String {
    let run = call("--bam", aln, reference, &scratch(label).join("x"));
    assert!(!run.status.success());
    String::from_utf8_lossy(&run.stderr).into_owned()
}

#[test]
fn cram_gives_the_bam_output() {
    let from_bam = outputs("--bam", "s.bam", "ref.fa", "bam");
    assert_eq!(outputs("--cram", "s.cram", "ref.fa", "cram"), from_bam);
    // bgzipped reference; many containers of several slices each
    assert_eq!(outputs("--bam", "small_slices.cram", "ref.fa.gz", "small"), from_bam);

    // and the comparison was of something: phased, methylated, partial reads
    let mut loci = from_bam.0.lines();
    let (header, row) = (loci.next().unwrap(), loci.next().unwrap());
    let loc: Vec<(&str, &str)> = header.split('\t').zip(row.split('\t')).collect();
    assert!(loc.contains(&("method", "HP")));
    let reads = &from_bam.1;
    assert!(reads.contains("left_partial") && reads.contains("split"));
    let rows: Vec<Vec<&str>> = reads.lines().map(|l| l.split('\t').collect()).collect();
    let col = rows[0].iter().position(|c| *c == "frac_5mc").unwrap();
    assert!(
        rows[1..].iter().any(|r| !matches!(r[col], "" | ".")),
        "no read carried methylation"
    );
}

#[test]
fn cram_sample_comes_from_the_read_group() {
    let out = scratch("sample").join("s");
    assert!(call("--cram", "s.cram", "ref.fa", &out).status.success());
    let vcf = std::fs::read_to_string(format!("{}.vcf", out.display())).unwrap();
    assert!(vcf.lines().any(|l| l.starts_with("#CHROM") && l.ends_with("\tSYNTH")));
}

#[test]
fn a_cram_with_a_shorter_contig_in_the_reference_is_refused() {
    let err = refused("s.cram", "short.fa", "short");
    assert!(err.contains("written against a chrT of 8,000 bp"), "{err}");
    assert!(err.contains("has 7,900 bp"), "{err}");
}

#[test]
fn a_cram_whose_contig_is_not_in_the_reference_is_refused() {
    let err = refused("s.cram", "renamed.fa", "renamed");
    assert!(err.contains("contig chrT is not in"), "{err}");
}

#[test]
fn a_cram_decoded_against_other_sequence_says_so() {
    for aln in ["s.cram", "small_slices.cram"] {
        let err = refused(aln, "other.fa", "other");
        assert!(err.starts_with("FXS_FMR1: "), "{err}");
        assert!(err.contains("not the reference the CRAM was written against"), "{err}");
    }
}
