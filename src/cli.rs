//! Command line. Port of `expansionscout/cli.py`: the same subcommands, flags
//! and outputs. The one difference is where the default catalogue lives:
//! here it is compiled in, and `--catalog PATH` reads another at run time.

use std::io::Write;
use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand};

use crate::bam::{Bam, Fasta};
use crate::catalog::{bundled_release, squash, Catalog, DEFAULT_BUILD};
use crate::engine::{call_sample, loci_header, loci_row, Options, PER_READ_COLS};
use crate::methylation::{CONVENTIONS, DEFAULT_CONVENTION};
use crate::pyfmt::{fixed, percent0, repr, wrap, Json};
use crate::qc::{format_report, read_tsv, run_checks};
use crate::report::{build_report, report_date};
use crate::vcf::{get_style, load_styles, render_vcf};

const MISMATCH_MIN_COVER: f64 = 0.30;
const MISMATCH_MIN_SHARE: f64 = 0.50;
const MISMATCH_MIN_LOCI: usize = 8;

fn version() -> &'static str {
    let s = format!(
        "{} (STRchive {} compiled in)",
        env!("CARGO_PKG_VERSION"),
        bundled_release()
    );
    Box::leak(s.into_boxed_str())
}

#[derive(Parser)]
#[command(name = "expansionscout", version = version(), infer_long_args = true,
          about = "Targeted repeat-expansion calling from ONT BAMs, with methylation.\n\n  \
                   call   genotype one or more catalogue loci from a BAM\n  \
                   loci   list the catalogue (id, gene, coordinates, regime, thresholds)\n  \
                   bed    export a BED for selected loci (reference-strand motif in col 4)")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// genotype catalogue loci from a BAM
    Call(Box<CallArgs>),
    /// list catalogue loci
    Loci(LociArgs),
    /// export BED for loci
    Bed(BedArgs),
    /// where every clinical band label comes from
    Provenance(ProvArgs),
    /// check a run against what is knowable without a truth set
    Qc(QcArgs),
    /// list the VCF dialects and what each is for
    #[command(name = "vcf-styles")]
    VcfStyles,
}

const BUILD_HELP: &str = "reference build the catalogue intervals are read from; the --ref FASTA must be this assembly";
const CATALOG_HELP: &str =
    "locus catalogue JSON in STRchive's format (default: the STRchive release compiled into this binary)";

#[derive(Args)]
struct CallArgs {
    #[arg(long)]
    bam: PathBuf,
    /// indexed reference FASTA
    #[arg(long = "ref")]
    reference: PathBuf,
    /// output prefix
    #[arg(long)]
    out: String,
    /// locus ids, aliases or gene names (default: whole catalogue)
    #[arg(long, num_args = 0..)]
    loci: Option<Vec<String>>,
    #[arg(long, help = CATALOG_HELP)]
    catalog: Option<PathBuf>,
    #[arg(long, default_value = DEFAULT_BUILD, value_parser = ["hg38", "t2t", "hg19"], help = BUILD_HELP)]
    build: String,
    #[arg(long)]
    sample: Option<String>,
    #[arg(long, default_value = "XX", value_parser = ["XX", "XY"])]
    sex: String,
    /// keep minor length modes as reported minor alleles
    #[arg(long)]
    mosaic: bool,
    /// report and classify under each locus's clinical convention: at loci quoted as a pure
    /// repeat count, such as HTT and ATXN1, scored interruptions are excluded from the allele
    /// size. The measured tract count is kept in the *_median_measured columns.
    #[arg(long)]
    clinical: bool,
    /// do not use HP tags
    #[arg(long)]
    ignore_hp: bool,
    #[arg(long, default_value_t = 1)]
    threads: usize,
    #[arg(long, default_value_t = 50)]
    anchor: i64,
    #[arg(long, default_value_t = 100)]
    flank: i64,
    #[arg(long, default_value_t = 20)]
    min_mapq: i64,
    #[arg(long, default_value_t = 50)]
    min_clip: i64,
    /// reference context beyond the tract in which repeat units may be counted (default 12;
    /// also at least 2 motif lengths)
    #[arg(long, default_value_t = 12)]
    margin: i64,
    #[arg(long, default_value_t = 3)]
    min_support: i64,
    #[arg(long, default_value_t = 0.10)]
    min_frac: f64,
    /// modkit pass threshold on base-modification call confidence. Calls below it are
    /// discarded rather than averaged in, which is modkit's rule. Not estimated from the
    /// locus on purpose; pass modkit's own value from `modkit sample-probs` for exact agreement
    #[arg(long, default_value_t = 0.8)]
    meth_threshold: f64,
    /// how a call carrying both 5mC and 5hmC probabilities is resolved: 'modkit' (argmax over
    /// canonical, 5mC and 5hmC), 'combine' (5hmC counts as modified) or 'ignore_h' (5hmC mass
    /// split between canonical and 5mC)
    #[arg(long, default_value = DEFAULT_CONVENTION, value_parser = CONVENTIONS)]
    meth_convention: String,
    /// also write a VCF, for tools that read one. The report and the two tables are always written
    #[arg(long)]
    vcf: bool,
    /// skip the HTML report
    #[arg(long)]
    no_report: bool,
    /// VCF dialect, from data/vcf_styles.yaml: 'native' (default), 'tr' or 'conservative'.
    /// `expansionscout vcf-styles` lists them
    #[arg(long, value_name = "NAME")]
    vcf_style: Option<String>,
}

#[derive(Args)]
struct LociArgs {
    #[arg(long, help = CATALOG_HELP)]
    catalog: Option<PathBuf>,
    #[arg(long, default_value = DEFAULT_BUILD, value_parser = ["hg38", "t2t", "hg19"])]
    build: String,
    #[arg(long, value_parser = ["short", "expansion"])]
    regime: Option<String>,
}

#[derive(Args)]
struct BedArgs {
    #[arg(long, help = CATALOG_HELP)]
    catalog: Option<PathBuf>,
    #[arg(long, default_value = DEFAULT_BUILD, value_parser = ["hg38", "t2t", "hg19"])]
    build: String,
    #[arg(long, num_args = 0..)]
    loci: Option<Vec<String>>,
    /// bases to pad each side (default 0). Use a read length or more when subsetting
    /// alignments, so reads that align beside a locus rather than across it are included
    #[arg(long, default_value_t = 0)]
    flank: i64,
    /// merge intervals that overlap after padding
    #[arg(long)]
    merge: bool,
    #[arg(long)]
    out: Option<PathBuf>,
}

#[derive(Args)]
struct ProvArgs {
    #[arg(long, default_value = DEFAULT_BUILD, value_parser = ["hg38", "t2t", "hg19"])]
    build: String,
    /// list the cited standards and papers instead
    #[arg(long)]
    sources: bool,
    /// where STRchive disagrees with a guideline, the reference or itself. Recorded, never applied
    #[arg(long)]
    discrepancies: bool,
    /// include laboratory-local boundaries, which are off by default and marked [in-house]
    #[arg(long)]
    inhouse: bool,
}

#[derive(Args)]
struct QcArgs {
    /// a .loci.tsv written by `expansionscout call`
    loci: PathBuf,
    /// also write the findings as TSV
    #[arg(long)]
    tsv: Option<PathBuf>,
    /// exit non-zero on a warning too, not only on a failure
    #[arg(long)]
    strict: bool,
}

/// `str(Path(p))`: pathlib's normal form, without resolving anything.
fn py_path(p: &str) -> String {
    if p.is_empty() {
        return ".".into();
    }
    let root = if p.starts_with("//") && !p.starts_with("///") {
        "//"
    } else if p.starts_with('/') {
        "/"
    } else {
        ""
    };
    let parts: Vec<&str> = p.split('/').filter(|c| !c.is_empty() && *c != ".").collect();
    let body = parts.join("/");
    if root.is_empty() && body.is_empty() {
        ".".into()
    } else {
        format!("{root}{body}")
    }
}

fn eprint_exit(msg: impl std::fmt::Display) -> i32 {
    eprintln!("{msg}");
    1
}

pub fn main() -> i32 {
    let cli = Cli::parse();
    let r = match cli.cmd {
        Cmd::Call(a) => cmd_call(*a),
        Cmd::Loci(a) => cmd_loci(a),
        Cmd::Bed(a) => cmd_bed(a),
        Cmd::Provenance(a) => cmd_provenance(a),
        Cmd::Qc(a) => cmd_qc(a),
        Cmd::VcfStyles => cmd_vcf_styles(),
    };
    match r {
        Ok(code) => code,
        Err(e) => eprint_exit(e),
    }
}

fn check_reference_matches_build(reference: &Path, loci: &[&crate::catalog::Locus], build: &str) {
    let Ok(mut fa) = Fasta::open(reference) else { return };
    let mut covers = Vec::new();
    for loc in loci {
        if !fa.references.contains(&loc.chrom) {
            continue;
        }
        let Ok(seq) = fa.fetch(&loc.chrom, loc.start.max(0), loc.end) else {
            return;
        };
        let seq = seq.to_ascii_uppercase();
        if seq.is_empty() {
            continue;
        }
        let covered: usize = loc.motifset().find_regions(&seq, None).iter().map(|(s, e)| e - s).sum();
        covers.push(covered as f64 / seq.len() as f64);
    }
    if covers.len() < MISMATCH_MIN_LOCI {
        return;
    }
    let bad = covers.iter().filter(|&&c| c < MISMATCH_MIN_COVER).count();
    if bad as f64 / covers.len() as f64 > MISMATCH_MIN_SHARE {
        eprintln!(
            "WARNING: {bad} of {} --build {build} intervals are less than {} covered by their own motif in {}. \
                   That is what a build/reference mismatch looks like; check the FASTA is {build}.",
            covers.len(),
            percent0(MISMATCH_MIN_COVER),
            reference.display()
        );
    }
}

fn write_file(path: &str, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("{path}: {e}"))
}

fn cmd_call(a: CallArgs) -> Result<i32, String> {
    let cat = Catalog::load(a.catalog.as_deref(), &a.build, false)?;
    let loci = cat.select(a.loci.as_deref())?;
    if loci.is_empty() {
        return Ok(eprint_exit("no loci selected"));
    }
    check_reference_matches_build(&a.reference, &loci, &a.build);
    let bam = Bam::open(&a.bam).map_err(|e| format!("{}: {e}", a.bam.display()))?;
    let sample = match a.sample.clone().filter(|s| !s.is_empty()) {
        Some(s) => s,
        None => bam.sample().unwrap_or_else(|| {
            a.bam
                .file_name()
                .map(|n| n.to_string_lossy().replace(".bam", ""))
                .unwrap_or_default()
        }),
    };
    let opts = Options {
        anchor: a.anchor,
        flank: a.flank,
        min_mapq: a.min_mapq,
        min_clip: a.min_clip,
        min_support: a.min_support,
        min_frac: a.min_frac,
        mosaic: a.mosaic,
        margin: a.margin,
        use_hp: !a.ignore_hp,
        sex: a.sex.clone(),
        meth_threshold: a.meth_threshold,
        meth_convention: a.meth_convention.clone(),
        clinical: a.clinical,
    };
    let style = get_style(a.vcf_style.as_deref())?;
    let results = call_sample(&a.bam, &a.reference, &loci, &sample, &opts, a.threads)?;
    {
        // htslib's own messages, which the Python writes as it goes; here in
        // locus order whatever the thread count, so stderr is reproducible too
        let mut err = std::io::stderr().lock();
        for res in &results {
            for d in &res.diagnostics {
                let _ = writeln!(err, "{d}");
            }
        }
    }
    let out = py_path(&a.out);
    if let Some(parent) = Path::new(&out).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| format!("{}: {e}", parent.display()))?;
        }
    }
    let per_read_path = format!("{out}.per_read.tsv");
    let loci_path = format!("{out}.loci.tsv");
    let mut per_read = PER_READ_COLS.join("\t") + "\n";
    for res in &results {
        for row in &res.rows {
            per_read.push_str(&row.join("\t"));
            per_read.push('\n');
        }
    }
    write_file(&per_read_path, &per_read)?;
    let mut loci_tsv = loci_header().join("\t") + "\n";
    for res in &results {
        loci_tsv.push_str(&loci_row(res).join("\t"));
        loci_tsv.push('\n');
    }
    write_file(&loci_path, &loci_tsv)?;

    // Contig lengths from the BAM first, then the reference for any the BAM lacks.
    let mut contigs: Vec<(String, i64)> = bam
        .references
        .iter()
        .cloned()
        .zip(bam.lengths.iter().copied())
        .collect();
    if let Ok(fa) = Fasta::open(&a.reference) {
        for (name, len) in fa.references.iter().zip(&fa.lengths) {
            if !contigs.iter().any(|c| c.0 == *name) {
                contigs.push((name.clone(), *len));
            }
        }
    }
    let mut written = vec![per_read_path.clone(), loci_path.clone()];
    if a.vcf {
        let path = format!("{out}.vcf");
        write_file(
            &path,
            &render_vcf(&sample, &contigs, &results, &style, Some(&cat.label)),
        )?;
        written.push(path);
    }
    if !a.no_report {
        let loci_rows = read_tsv(&loci_tsv);
        let read_rows = read_tsv(&per_read);
        let findings = run_checks(&loci_rows);
        let meta: Vec<(&str, String)> = vec![
            ("sample", sample.clone()),
            ("build", a.build.clone()),
            ("reference", a.reference.to_string_lossy().into_owned()),
            ("alignment", a.bam.to_string_lossy().into_owned()),
            ("catalogue", cat.label.clone()),
            ("expansionscout", env!("CARGO_PKG_VERSION").to_string()),
            ("options", format!(
                "sex={} min_support={} min_frac={} mosaic={} use_hp={} clinical={} meth_threshold={} meth_convention={}",
                a.sex, a.min_support, repr(a.min_frac), py_bool(a.mosaic), py_bool(!a.ignore_hp),
                py_bool(a.clinical), repr(a.meth_threshold), a.meth_convention)),
        ];
        let path = format!("{out}.report.html");
        let page = build_report(&loci_rows, &read_rows, &findings, &sample, &cat, &meta, &report_date());
        write_file(&path, &page)?;
        written.push(path);
    }
    let mut err = std::io::stderr().lock();
    for res in &results {
        let h = &res.haplotypes;
        let cn = if h.is_empty() {
            ".".into()
        } else {
            h.iter().map(|x| fixed(x.median.f(), 0)).collect::<Vec<_>>().join("/")
        };
        let cls = if h.is_empty() {
            ".".into()
        } else {
            h.iter().map(|x| x.cls.clone()).collect::<Vec<_>>().join("/")
        };
        let lb = res.lb_units.map(|v| v.to_string()).unwrap_or_else(|| ".".into());
        let detect = res
            .detect
            .map(|d| repr(fixed(d, 2).parse().unwrap()))
            .unwrap_or_else(|| "None".into());
        let notes = if res.notes.is_empty() {
            String::new()
        } else {
            format!(" {}", res.notes.join("; "))
        };
        let _ = writeln!(
            err,
            "{:<22} span={:<4} split={:<3} partial={:<3} CN={:<10} {:<24} LB={:<5} detect={} {} [{}]{}",
            res.locus.id,
            res.n_spanning,
            res.n_split,
            res.n_partial,
            cn,
            cls,
            lb,
            detect,
            res.evidence,
            res.method,
            notes
        );
    }
    let _ = writeln!(err, "wrote {}", written.join(" "));
    Ok(0)
}

fn py_bool(b: bool) -> &'static str {
    if b {
        "True"
    } else {
        "False"
    }
}

fn cmd_loci(a: LociArgs) -> Result<i32, String> {
    let cat = Catalog::load(a.catalog.as_deref(), &a.build, false)?;
    let mut out = std::io::stdout().lock();
    let head = [
        "id",
        "gene",
        "chrom",
        "start",
        "end",
        "strand",
        "motif_gene",
        "motif_ref",
        "unit_len",
        "regime",
        "composition",
        "benign_max",
        "intermediate",
        "pathogenic_min",
        "pathogenic_max",
        "meth_relevant",
        "meth_up_bp",
        "meth_down_bp",
        "aliases",
    ];
    let _ = writeln!(out, "{}", head.join("\t"));
    use crate::pyfmt::opt_num_str;
    for loc in &cat.loci {
        if a.regime.as_ref().is_some_and(|r| *r != loc.regime) {
            continue;
        }
        let inter = match loc.intermediate_min {
            Some(m) => format!("{}-{}", m.py_str(), opt_num_str(loc.intermediate_max)),
            None => ".".into(),
        };
        let row = [
            loc.id.clone(),
            loc.gene.clone(),
            loc.chrom.clone(),
            loc.start.to_string(),
            loc.end.to_string(),
            loc.gene_strand.clone(),
            loc.gene_motif(),
            loc.ref_motif().to_string(),
            loc.unit_len.to_string(),
            loc.regime.clone(),
            (loc.composition() as i64).to_string(),
            opt_num_str(loc.benign_max),
            inter,
            opt_num_str(loc.pathogenic_min),
            opt_num_str(loc.pathogenic_max),
            (loc.meth_relevant as i64).to_string(),
            loc.meth_up_bp.to_string(),
            loc.meth_down_bp.to_string(),
            loc.aliases.join(","),
        ];
        let _ = writeln!(out, "{}", row.join("\t"));
    }
    if !cat.dropped.is_empty() {
        eprintln!(
            "dropped (no same-length reference motif, or no {} interval): {}",
            a.build,
            cat.dropped.join(", ")
        );
    }
    Ok(0)
}

fn cmd_bed(a: BedArgs) -> Result<i32, String> {
    let cat = Catalog::load(a.catalog.as_deref(), &a.build, false)?;
    let mut rows: Vec<(String, i64, i64, String, String)> = cat
        .select(a.loci.as_deref())?
        .iter()
        .map(|l| {
            (
                l.chrom.clone(),
                0.max(l.start - a.flank),
                l.end + a.flank,
                l.ref_motif().to_string(),
                l.id.clone(),
            )
        })
        .collect();
    rows.sort_by(|x, y| (&x.0, x.1).cmp(&(&y.0, y.1)));
    if a.merge {
        let mut merged: Vec<(String, i64, i64, String, String)> = Vec::new();
        for r in rows {
            match merged.last_mut() {
                Some(prev) if prev.0 == r.0 && r.1 <= prev.2 => {
                    prev.2 = prev.2.max(r.2);
                    prev.4 = format!("{},{}", prev.4, r.4);
                }
                _ => merged.push(r),
            }
        }
        rows = merged;
    }
    let text: String = rows
        .iter()
        .map(|(c, s, e, m, n)| format!("{c}\t{s}\t{e}\t{m}\t{n}\n"))
        .collect();
    match &a.out {
        Some(p) => std::fs::write(p, &text).map_err(|e| format!("{}: {e}", p.display()))?,
        None => {
            let _ = std::io::stdout().lock().write_all(text.as_bytes());
        }
    }
    eprintln!(
        "{} interval(s), flank {} bp{}",
        rows.len(),
        a.flank,
        if a.merge { ", overlaps merged" } else { "" }
    );
    Ok(0)
}

fn jstr(v: Option<&Json>) -> String {
    v.map(Json::py_str).unwrap_or_else(|| "None".into())
}

fn cmd_provenance(a: ProvArgs) -> Result<i32, String> {
    let cat = Catalog::load(None, &a.build, a.inhouse)?;
    let meta = &cat.provenance;
    let mut out = String::new();
    let mut p = |s: &str| {
        out.push_str(s);
        out.push('\n');
    };
    let source = |key: Option<&Json>| -> Json {
        key.and_then(|k| k.as_str())
            .and_then(|k| meta.sources.get(k))
            .cloned()
            .unwrap_or(Json::Obj(vec![]))
    };
    if a.discrepancies {
        p("Where STRchive disagrees with a guideline, with the reference, or with itself.");
        p("NONE of these are applied: the catalogue is read as published, so");
        p("results stay reproducible against it and the audit stays honest.\n");
        for d in &meta.discrepancies {
            let src = source(d.get("other_source"));
            p(&format!(
                "  {}  {}: catalogue {} vs {} ({})",
                jstr(d.get("locus")),
                jstr(d.get("field")),
                jstr(d.get("catalogue")),
                jstr(d.get("other")),
                jstr(d.get("other_source"))
            ));
            p(&format!(
                "    kind={}  applied={}",
                jstr(d.get("kind")),
                jstr(d.get("applied"))
            ));
            if let Some(j) = src.get("jurisdiction").filter(|j| j.truthy()) {
                p(&format!("    jurisdiction: {}", j.py_str()));
            }
            p(&format!(
                "    {}",
                squash(
                    &d.get("note")
                        .filter(|n| n.truthy())
                        .map(Json::py_str)
                        .unwrap_or_default()
                )
            ));
            p("");
        }
        print!("{out}");
        return Ok(0);
    }
    if a.sources {
        let mut items: Vec<&(String, Json)> = meta.sources.entries().iter().collect();
        let kind = |v: &Json| v.get("kind").map(Json::py_str).unwrap_or_default();
        items.sort_by(|x, y| (kind(&x.1), &x.0).cmp(&(kind(&y.1), &y.0)));
        for (key, s) in items {
            let oa = if s.get("open_access").is_some_and(Json::truthy) {
                "open access"
            } else {
                "paywalled"
            };
            let jur = s
                .get("jurisdiction")
                .filter(|j| j.truthy())
                .map(|j| format!("  |  {}", j.py_str()))
                .unwrap_or_default();
            p(&format!("[{key}]  {}{jur}  |  {oa}", jstr(s.get("kind"))));
            p(&format!(
                "   {}",
                squash(
                    &s.get("citation")
                        .filter(|c| c.truthy())
                        .map(Json::py_str)
                        .unwrap_or_default()
                )
            ));
            for f in ["pmid", "pmcid", "doi"] {
                if let Some(v) = s.get(f).filter(|v| v.truthy()) {
                    p(&format!("   {f}: {}", v.py_str()));
                }
            }
            p("");
        }
        print!("{out}");
        return Ok(0);
    }
    p(&format!(
        "overlay version {}{}",
        meta.version.py_str(),
        if meta.inhouse_enabled {
            "   IN-HOUSE BOUNDARIES ENABLED"
        } else {
            ""
        }
    ));
    p("Bands are emitted only where boundary AND label can both be cited.");
    p("Loci with no bands here fall back to the catalogue's own");
    p("benign / intermediate / pathogenic split, which is STRchive's.\n");
    let mut lids: Vec<&(String, Vec<Json>)> = meta.bands.iter().collect();
    lids.sort_by(|x, y| x.0.cmp(&y.0));
    for (lid, bands) in lids {
        let gene = cat
            .loci
            .iter()
            .find(|l| l.id == *lid)
            .map(|l| l.gene.clone())
            .unwrap_or_else(|| "?".into());
        if bands.is_empty() {
            p(&format!("{lid}  ({gene}): no bands -- catalogue split used"));
            continue;
        }
        p(&format!("{lid}  ({gene})"));
        for b in bands {
            let up = match b.get_some("upper") {
                None => "and above".to_string(),
                Some(u) => format!("<= {}", u.py_str()),
            };
            p(&format!("   {up:<12}{}", jstr(b.get("label"))));
            if b.get("_inhouse").is_some_and(Json::truthy) {
                p("     in-house  not supported by any standard we could find");
                if let Some(r) = b.get("rationale").filter(|r| r.truthy()) {
                    p(&format!("               {}", squash(&r.py_str())));
                }
                continue;
            }
            for (key, what) in [("boundary", "boundary"), ("naming", "label")] {
                let d = b.get(key).filter(|d| d.truthy()).cloned().unwrap_or(Json::Obj(vec![]));
                let src = source(d.get("source"));
                let mut bit = format!(
                    "     {what:<9}{}",
                    d.get("source").map(Json::py_str).unwrap_or_else(|| "?".into())
                );
                if let Some(sec) = d.get("section").filter(|s| s.truthy()) {
                    bit.push_str(&format!("  {}", sec.py_str()));
                }
                if let Some(j) = src.get("jurisdiction").filter(|j| j.truthy()) {
                    bit.push_str(&format!("  [{}]", j.py_str()));
                }
                p(&bit);
                if let Some(q) = d.get("quote").filter(|q| q.truthy()) {
                    p(&format!("                \"{}\"", squash(&q.py_str())));
                }
                if let Some(n) = d.get("note").filter(|n| n.truthy()) {
                    p(&format!("                {}", squash(&n.py_str())));
                }
            }
        }
        p("");
    }
    print!("{out}");
    Ok(0)
}

fn cmd_qc(a: QcArgs) -> Result<i32, String> {
    let text = std::fs::read_to_string(&a.loci).map_err(|e| format!("{}: {e}", a.loci.display()))?;
    let rows = read_tsv(&text);
    let sample = match rows.first() {
        Some(r) => crate::qc::get(r, "sample").map(String::from),
        None => a.loci.file_stem().map(|s| s.to_string_lossy().into_owned()),
    };
    let findings = run_checks(&rows);
    if let Some(tsv) = &a.tsv {
        let mut t = String::from("check\tstatus\tbasis\tmessage\n");
        for x in &findings {
            t.push_str(&format!("{}\t{}\t{}\t{}\n", x.check, x.status, x.basis, x.message));
        }
        std::fs::write(tsv, t).map_err(|e| format!("{}: {e}", tsv.display()))?;
        eprintln!("wrote {}", tsv.display());
    }
    println!("{}", format_report(&findings, sample.as_deref()));
    if findings.iter().any(|x| x.status == "fail") {
        return Ok(2);
    }
    if a.strict && findings.iter().any(|x| !x.ok()) {
        return Ok(1);
    }
    Ok(0)
}

fn cmd_vcf_styles() -> Result<i32, String> {
    let (styles, sources) = load_styles()?;
    let mut out = String::new();
    for st in &styles {
        let verified = st.spec.get("verified").filter(|v| v.truthy()).map(Json::py_str);
        let mark = match &verified {
            Some(v) => sources
                .get(v)
                .and_then(|s| s.get("citation"))
                .map(Json::py_str)
                .unwrap_or_else(|| v.clone()),
            None => "no external source claimed".into(),
        };
        out.push_str(&format!(
            "{}{}\n",
            st.name,
            if st.is_default() { "  [default]" } else { "" }
        ));
        out.push_str(&format!(
            "    VCF {}, {}, {} INFO / {} FORMAT\n",
            st.fileformat(),
            if st.alt() == "cnv_tr" {
                "symbolic <CNV:TR>"
            } else {
                "sequence alleles"
            },
            st.fields("info").len(),
            st.fields("format").len()
        ));
        out.push_str(&format!("    source: {mark}\n"));
        let desc = squash(
            &st.spec
                .get("description")
                .filter(|d| d.truthy())
                .map(Json::py_str)
                .unwrap_or_default(),
        );
        if !desc.is_empty() {
            for line in wrap(&desc, 74) {
                out.push_str(&format!("    {line}\n"));
            }
        }
        out.push('\n');
    }
    print!("{out}");
    Ok(0)
}
