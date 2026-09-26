//! VCF writer, one record per locus, in a dialect chosen at the command line.
//! Port of `expansionscout/vcf.py`; the dialects are data in
//! `data/vcf_styles.yaml`, compiled into the binary.

use std::fmt::Write as _;

use crate::engine::{fmt_opt_f, Hap, LocusResult, PyVal};
use crate::pyfmt::{fixed, Json};

pub const STYLES_YAML: &str = include_str!("../data/vcf_styles.yaml");
pub const DEFAULT_STYLE: &str = "native";

pub const INFO_DEFS: &[(&str, &str, &str, &str)] = &[
    ("END", "1", "Integer", "End of the reference tract (1-based, inclusive). Deprecated in VCF 4.5; kept because it is what makes a record readable without arithmetic"),
    ("SVLEN", "A", "Integer", "Length of the reference segment the repeat is defined over"),
    ("LOCUS", "1", "String", "STRchive locus id"),
    ("GENE", "1", "String", "Gene"),
    ("MOTIF", "1", "String", "Repeat motif on the gene strand"),
    ("REFMOTIF", "1", "String", "Repeat motif on the reference strand"),
    ("MOTIFLEN", "1", "Integer", "Motif length"),
    ("REGIME", "1", "String", "short or expansion (derived from pathogenic_max)"),
    ("REFCN", "1", "Float", "Reference copy number (STRchive)"),
    ("NSPAN", "1", "Integer", "Spanning reads"),
    ("NSPLIT", "1", "Integer", "Split-spanning reads (size from read coordinates)"),
    ("NPART", "1", "Integer", "Partial reads (one anchor aligned, clipped into the tract)"),
    ("LB", "1", "Integer", "Lower bound in units from partial/split reads, or missing"),
    ("LBSUP", "1", "Integer", "Reads supporting LB"),
    ("DETECT", "1", "Float", "Fraction of reads long enough to span a pathogenic_min allele"),
    ("EVIDENCE", "1", "String", "NONE, PARTIALS_PRESENT or EXPANSION_LB"),
    ("NEGREL", "1", "Integer", "1 if a normal-range call is reliable negative evidence"),
    ("METHOD", "1", "String", "Allele assignment method (HP, HP-haploid, GMM)"),
    ("MINOR", ".", "String", "Minor modes kept in mosaic mode as units:support"),
    ("RN", "A", "Integer", "Total number of repeat sequences in this allele"),
    ("RUS", ".", "String", "Repeat unit sequence of the corresponding repeat sequence"),
    ("RUL", ".", "Integer", "Repeat unit length of the corresponding repeat sequence"),
    ("RUC", ".", "Float", "Repeat unit count of corresponding repeat sequence"),
    ("RB", ".", "Integer", "Total number of bases in the corresponding repeat sequence"),
    ("CIRUC", ".", "Float", "Confidence interval around RUC"),
];

pub const FORMAT_DEFS: &[(&str, &str, &str, &str)] = &[
    ("GT", "1", "String", "Genotype"),
    ("CN", ".", "Float", "Per-haplotype median unit count"),
    (
        "CNR",
        ".",
        "String",
        "Per-haplotype central 90 percent range of unit counts (p5-p95)",
    ),
    ("SUP", ".", "Integer", "Per-haplotype supporting reads"),
    ("SIZE", ".", "Integer", "Per-haplotype median tract length in bp"),
    (
        "MAD",
        ".",
        "Float",
        "Per-haplotype median absolute deviation of unit counts",
    ),
    (
        "TAIL",
        ".",
        "Float",
        "Per-haplotype fraction of reads in an upward length tail (mosaic indicator)",
    ),
    (
        "NOISE",
        ".",
        "Float",
        "Per-haplotype mean fraction of unrecognised units",
    ),
    ("INT", ".", "Integer", "Per-haplotype interruption units"),
    (
        "INTPOS",
        ".",
        "String",
        "Per-haplotype consensus interruption positions (1-based units, / separated)",
    ),
    (
        "PM",
        ".",
        "Integer",
        "Per-haplotype pathogenic-motif units (composition loci)",
    ),
    ("CLS", ".", "String", "Per-haplotype size class"),
    ("MT", ".", "Float", "Per-haplotype mean 5mC over the tract"),
    ("MU", ".", "Float", "Per-haplotype mean 5mC over the gene-5prime window"),
    ("MD", ".", "Float", "Per-haplotype mean 5mC over the gene-3prime window"),
    (
        "MF",
        ".",
        "Float",
        "Per-haplotype fraction of reads with tract 5mC >= 0.5 (methylation mosaic indicator)",
    ),
];

const FORMAT_KEYS: &[(&str, &str, usize)] = &[
    ("CN", "median", 1),
    ("CNR", "cnr", 3),
    ("SUP", "support", 3),
    ("SIZE", "size_bp", 3),
    ("MAD", "mad", 2),
    ("TAIL", "tail_up_frac", 3),
    ("NOISE", "noise", 3),
    ("INT", "n_interruption", 3),
    ("INTPOS", "int_pos", 3),
    ("PM", "n_pathogenic", 3),
    ("CLS", "cls", 3),
    ("MT", "meth_tract", 3),
    ("MU", "meth_up", 3),
    ("MD", "meth_down", 3),
    ("MF", "frac_meth_reads_gt50", 3),
];

fn def<'a>(
    table: &'a [(&'a str, &'a str, &'a str, &'a str)],
    id: &str,
) -> Option<&'a (&'a str, &'a str, &'a str, &'a str)> {
    table.iter().find(|d| d.0 == id)
}

/// A resolved dialect: `like` followed, fields validated.
#[derive(Clone, Debug)]
pub struct Style {
    pub name: String,
    pub spec: Json,
}

impl Style {
    fn str_field(&self, key: &str) -> Option<String> {
        self.spec.get_some(key).map(|v| v.py_str())
    }
    pub fn fileformat(&self) -> String {
        self.str_field("fileformat").unwrap_or_else(|| "4.2".into())
    }
    pub fn alt(&self) -> String {
        self.str_field("alt").unwrap_or_else(|| "sequence".into())
    }
    pub fn is_default(&self) -> bool {
        self.spec.get("default").map(Json::truthy).unwrap_or(false)
    }
    /// [(emitted id, definition id)] for the style's INFO or FORMAT list.
    pub fn fields(&self, kind: &str) -> Vec<(String, String)> {
        self.spec
            .get(kind)
            .map(|l| {
                l.as_list()
                    .iter()
                    .map(|e| match e {
                        Json::Obj(_) => {
                            let id = e.get("id").map(Json::py_str).unwrap_or_default();
                            (e.get("as").map(Json::py_str).unwrap_or_else(|| id.clone()), id)
                        }
                        other => (other.py_str(), other.py_str()),
                    })
                    .collect()
            })
            .unwrap_or_default()
    }
}

/// Resolve `data/vcf_styles.yaml` into styles, in file order, and its sources.
pub fn load_styles() -> Result<(Vec<Style>, Json), String> {
    let docs = yaml_rust2::YamlLoader::load_from_str(STYLES_YAML).map_err(|e| format!("vcf_styles.yaml: {e}"))?;
    let doc = docs.first().map(Json::from_yaml).unwrap_or(Json::Null);
    let raw = doc.get("styles").cloned().unwrap_or(Json::Obj(vec![]));
    let mut out = Vec::new();
    for (name, spec) in raw.entries() {
        let mut merged: Vec<(String, Json)> = Vec::new();
        let put = |k: &str, v: &Json, merged: &mut Vec<(String, Json)>| match merged.iter_mut().find(|e| e.0 == k) {
            Some(e) => e.1 = v.clone(),
            None => merged.push((k.to_string(), v.clone())),
        };
        if let Some(parent) = spec.get("like").filter(|p| p.truthy()) {
            let pname = parent.py_str();
            let Some(p) = raw.get(&pname) else {
                return Err(format!("vcf style '{name}' aliases unknown style '{pname}'"));
            };
            for (k, v) in p.entries() {
                if k != "default" {
                    put(k, v, &mut merged);
                }
            }
        }
        for (k, v) in spec.entries() {
            if k != "like" {
                put(k, v, &mut merged);
            }
        }
        put("name", &Json::Str(name.clone()), &mut merged);
        out.push(Style {
            name: name.clone(),
            spec: Json::Obj(merged),
        });
    }
    for st in &out {
        for (kind, table) in [("info", INFO_DEFS), ("format", FORMAT_DEFS)] {
            for (_, id) in st.fields(kind) {
                if def(table, &id).is_none() {
                    return Err(format!(
                        "vcf style '{}' names unknown {} field '{id}'",
                        st.name,
                        kind.to_uppercase()
                    ));
                }
            }
        }
        if !st.fields("format").iter().any(|(_, id)| id == "GT") {
            return Err(format!("vcf style '{}' must emit GT", st.name));
        }
    }
    Ok((out, doc.get("sources").cloned().unwrap_or(Json::Obj(vec![]))))
}

pub fn get_style(name: Option<&str>) -> Result<Style, String> {
    let (styles, _) = load_styles()?;
    let name = match name {
        Some(n) => n.to_string(),
        None => styles
            .iter()
            .find(|s| s.is_default())
            .map(|s| s.name.clone())
            .unwrap_or(DEFAULT_STYLE.into()),
    };
    styles.iter().find(|s| s.name == name).cloned().ok_or_else(|| {
        let mut names: Vec<&str> = styles.iter().map(|s| s.name.as_str()).collect();
        names.sort();
        format!("unknown vcf style '{name}'; have {}", names.join(", "))
    })
}

/// `round(float(v), 3)`, as an allele key.
fn round3(v: f64) -> f64 {
    fixed(v, 3).parse().unwrap()
}

#[derive(PartialEq, Clone)]
enum Key {
    Num(f64),
    Seq(String),
}

struct Alleles<'r> {
    pad: String,
    ref_seq: String,
    alts: Vec<String>,
    alt_haps: Vec<&'r Hap>,
    gt: Vec<String>,
}

impl<'r> Alleles<'r> {
    fn new(res: &'r LocusResult, cnv_tr: bool) -> Alleles<'r> {
        let loc = res.locus;
        let pad = res
            .ref_pad_base
            .clone()
            .filter(|p| !p.is_empty())
            .unwrap_or_else(|| "N".into());
        let ref_seq = format!("{pad}{}", res.ref_tract_seq.clone().unwrap_or_default());
        let ref_key = if cnv_tr {
            loc.ref_copies.map(|c| Key::Num(round3(c.f())))
        } else {
            Some(Key::Seq(ref_seq.clone()))
        };
        let mut al = Alleles {
            pad: pad.clone(),
            ref_seq,
            alts: vec![],
            alt_haps: vec![],
            gt: vec![],
        };
        let mut keys: Vec<Key> = Vec::new();
        for h in &res.haplotypes {
            let key = if cnv_tr {
                Some(Key::Num(round3(h.median.f())))
            } else {
                h.seq.as_ref().map(|s| Key::Seq(format!("{pad}{s}")))
            };
            let Some(key) = key else {
                al.gt.push(".".into());
                continue;
            };
            if ref_key.as_ref() == Some(&key) {
                al.gt.push("0".into());
                continue;
            }
            if !keys.contains(&key) {
                al.alts.push(match &key {
                    Key::Num(_) => "<CNV:TR>".to_string(),
                    Key::Seq(s) => s.clone(),
                });
                al.alt_haps.push(h);
                keys.push(key.clone());
            }
            al.gt
                .push((keys.iter().position(|k| *k == key).unwrap() + 1).to_string());
        }
        al
    }
}

fn info_value(res: &LocusResult, al: &Alleles, id: &str) -> Option<String> {
    let loc = res.locus;
    let join = |f: &dyn Fn(&Hap) -> String| al.alt_haps.iter().map(|h| f(h)).collect::<Vec<_>>().join(",");
    let v = match id {
        "END" => loc.end.to_string(),
        "SVLEN" => {
            if al.alts.is_empty() {
                return None;
            }
            let tl = (al.ref_seq.chars().count() as i64 - 1).max(0);
            al.alts.iter().map(|_| tl.to_string()).collect::<Vec<_>>().join(",")
        }
        "LOCUS" => loc.id.clone(),
        "GENE" => loc.gene.clone(),
        "MOTIF" => loc.gene_motif(),
        "REFMOTIF" => loc.ref_motif().to_string(),
        "MOTIFLEN" => loc.unit_len.to_string(),
        "REGIME" => loc.regime.clone(),
        "REFCN" => PyVal::num(loc.ref_copies).fmt(1),
        "NSPAN" => res.n_spanning.to_string(),
        "NSPLIT" => res.n_split.to_string(),
        "NPART" => res.n_partial.to_string(),
        "LB" => PyVal::opt_i(res.lb_units).fmt(3),
        "LBSUP" => res.lb_support.to_string(),
        "DETECT" => fmt_opt_f(res.detect, 3),
        "EVIDENCE" => res.evidence.clone(),
        "NEGREL" => (res.negative_reliable as i64).to_string(),
        "METHOD" => res.method.clone(),
        "MINOR" => {
            if res.minor_modes.is_empty() {
                return None;
            }
            res.minor_modes
                .iter()
                .map(|m| format!("{}:{}", fixed(m.median.f(), 0), m.support))
                .collect::<Vec<_>>()
                .join(",")
        }
        _ if al.alts.is_empty() => return None,
        "RN" => al.alts.iter().map(|_| "1").collect::<Vec<_>>().join(","),
        "RUS" => join(&|h| h.motif.clone().unwrap_or_else(|| loc.ref_motif().to_string())),
        "RUL" => al
            .alts
            .iter()
            .map(|_| loc.unit_len.to_string())
            .collect::<Vec<_>>()
            .join(","),
        "RUC" => join(&|h| h.get("median").fmt(1)),
        "RB" => join(&|h| h.get("size_bp").fmt(0)),
        "CIRUC" => join(&|h| format!("{},{}", h.get("p5").fmt(1), h.get("p95").fmt(1))),
        _ => return None,
    };
    Some(v)
}

fn record(res: &LocusResult, style: &Style) -> String {
    let loc = res.locus;
    let cnv_tr = style.alt() == "cnv_tr";
    let al = Alleles::new(res, cnv_tr);
    let sep = if res.method.starts_with("HP") { "|" } else { "/" };
    let gt_str = if al.gt.is_empty() {
        ".".to_string()
    } else {
        al.gt.join(sep)
    };
    let mut info = Vec::new();
    for (emit, id) in style.fields("info") {
        match info_value(res, &al, &id) {
            Some(v) if v != "." => info.push(format!("{emit}={v}")),
            _ => {}
        }
    }
    let (mut keys, mut values) = (Vec::new(), Vec::new());
    for (emit, id) in style.fields("format") {
        keys.push(emit);
        if id == "GT" {
            values.push(gt_str.clone());
            continue;
        }
        let (_, key, nd) = FORMAT_KEYS.iter().find(|f| f.0 == id).unwrap();
        values.push(if res.haplotypes.is_empty() {
            ".".to_string()
        } else {
            res.haplotypes
                .iter()
                .map(|h| h.get(key).fmt(*nd))
                .collect::<Vec<_>>()
                .join(",")
        });
    }
    let reff = if cnv_tr { al.pad.clone() } else { al.ref_seq.clone() };
    [
        loc.chrom.clone(),
        loc.start.to_string(),
        loc.id.clone(),
        reff,
        if al.alts.is_empty() {
            ".".into()
        } else {
            al.alts.join(",")
        },
        ".".into(),
        if res.n_spanning + res.n_split as i64 > 0 {
            "PASS".into()
        } else {
            "NOREADS".into()
        },
        if info.is_empty() { ".".into() } else { info.join(";") },
        keys.join(":"),
        values.join(":"),
    ]
    .join("\t")
}

/// The VCF as a string. `contigs` in the order the lengths were collected;
/// only the contigs actually written are declared.
pub fn render_vcf(
    sample: &str,
    contigs: &[(String, i64)],
    results: &[LocusResult],
    style: &Style,
    catalogue: Option<&str>,
) -> String {
    let mut seen: Vec<&str> = Vec::new();
    for r in results {
        if !seen.contains(&r.locus.chrom.as_str()) {
            seen.push(&r.locus.chrom);
        }
    }
    let mut f = String::new();
    writeln!(f, "##fileformat=VCFv{}", style.fileformat()).unwrap();
    writeln!(f, "##source=expansionscout {}", env!("CARGO_PKG_VERSION")).unwrap();
    writeln!(f, "##expansionscoutVcfStyle={}", style.name).unwrap();
    if let Some(c) = catalogue {
        writeln!(f, "##expansionscoutCatalogue={c}").unwrap();
    }
    for name in seen {
        match contigs.iter().find(|c| c.0 == name) {
            Some((_, len)) => writeln!(f, "##contig=<ID={name},length={len}>").unwrap(),
            None => writeln!(f, "##contig=<ID={name}>").unwrap(),
        }
    }
    if style.alt() == "cnv_tr" {
        writeln!(f, "##ALT=<ID=CNV:TR,Description=\"Tandem repeat\">").unwrap();
    }
    writeln!(
        f,
        "##FILTER=<ID=NOREADS,Description=\"No spanning or split reads at this locus\">"
    )
    .unwrap();
    for (kind, table, head) in [("info", INFO_DEFS, "INFO"), ("format", FORMAT_DEFS, "FORMAT")] {
        for (emit, id) in style.fields(kind) {
            let d = def(table, &id).unwrap();
            writeln!(
                f,
                "##{head}=<ID={emit},Number={},Type={},Description=\"{}\">",
                d.1, d.2, d.3
            )
            .unwrap();
        }
    }
    writeln!(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}").unwrap();
    for r in results {
        f.push_str(&record(r, style));
        f.push('\n');
    }
    f
}
