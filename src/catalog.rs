//! Locus catalogue: STRchive records plus the local overlay.
//!
//! Port of `expansionscout/catalog.py`; see it for the reasoning. The
//! difference that matters here is where the default catalogue comes from:
//! the STRchive release vendored in the repository is compiled into the
//! binary, so a binary always carries the catalogue it was built with and
//! `--catalog` substitutes another one at run time.

use std::collections::HashMap;
use std::path::Path;

use sha2::{Digest, Sha256};

use crate::decompose::MotifSet;
use crate::pyfmt::{Json, Num};

/// The STRchive release this binary was built with.
pub const BUNDLED_CATALOG: &str = include_str!("../data/strchive/STRchive-loci.json");
pub const BUNDLED_VERSION: &str = include_str!("../data/strchive/VERSION");
pub const OVERLAY_YAML: &str = include_str!("../data/overlay.yaml");

pub const EXPANSION_REGIME_BP: f64 = 1000.0;
pub const BUILDS: [&str; 3] = ["hg38", "t2t", "hg19"];
pub const DEFAULT_BUILD: &str = "hg38";
pub const BAND_TOP: f64 = 1e9;

pub fn revcomp(seq: &str) -> String {
    seq.to_ascii_uppercase()
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            other => other as char,
        })
        .collect()
}

/// `key: value` from the VERSION file.
fn version_field(key: &str) -> Option<String> {
    BUNDLED_VERSION.lines().find_map(|l| {
        l.strip_prefix(key)
            .and_then(|r| r.strip_prefix(':'))
            .map(|v| v.trim().to_string())
    })
}

pub fn bundled_release() -> String {
    version_field("release").unwrap_or_else(|| "unknown".into())
}

pub fn sha256_hex(bytes: &[u8]) -> String {
    let d = Sha256::digest(bytes);
    d.iter().map(|b| format!("{:02x}", b)).collect()
}

/// How a run names the catalogue it used, in the report and the VCF header.
/// A file whose content is the bundled release is named as that release,
/// wherever it was read from.
pub fn catalogue_label(path: Option<&str>, bytes: &[u8]) -> String {
    let digest = sha256_hex(bytes);
    if Some(digest.as_str()) == version_field("sha256").as_deref() {
        format!("STRchive {} (bundled), sha256:{}", bundled_release(), digest)
    } else {
        format!("{}, sha256:{}", path.unwrap_or("?"), digest)
    }
}

#[derive(Clone, Debug)]
pub struct Band {
    pub upper: f64,
    pub label: String,
}

#[derive(Clone, Debug)]
pub struct Locus {
    pub id: String,
    pub gene: String,
    pub disease: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub gene_strand: String,
    pub unit_len: usize,
    pub ref_motifs: Vec<String>,
    pub pathogenic_motifs: Vec<String>,
    pub benign_motifs: Vec<String>,
    pub interruption_motifs: Vec<String>,
    pub benign_min: Option<Num>,
    pub benign_max: Option<Num>,
    pub intermediate_min: Option<Num>,
    pub intermediate_max: Option<Num>,
    pub pathogenic_min: Option<Num>,
    pub pathogenic_max: Option<Num>,
    pub ref_copies: Option<Num>,
    pub inheritance: String,
    pub mechanism: String,
    pub regime: String,
    pub meth_relevant: bool,
    pub meth_up_bp: i64,
    pub meth_down_bp: i64,
    pub meth_note: String,
    pub build: String,
    pub aliases: Vec<String>,
    pub caveats: Vec<String>,
    pub interruptions_in_size: bool,
    pub bands: Vec<Band>,
    motifset: MotifSet,
}

impl Locus {
    pub fn ref_motif(&self) -> &str {
        &self.ref_motifs[0]
    }

    pub fn gene_motif(&self) -> String {
        if self.gene_strand == "+" {
            self.ref_motif().to_string()
        } else {
            revcomp(self.ref_motif())
        }
    }

    pub fn to_gene_strand(&self, unit: &str) -> String {
        if self.gene_strand == "+" {
            unit.to_ascii_uppercase()
        } else {
            revcomp(unit)
        }
    }

    pub fn motifset(&self) -> &MotifSet {
        &self.motifset
    }

    pub fn composition(&self) -> bool {
        self.motifset.composition
    }

    /// Gene-5' flank window as reference [a, b), or None.
    pub fn up_window(&self) -> Option<(i64, i64)> {
        if self.meth_up_bp <= 0 {
            return None;
        }
        if self.gene_strand == "+" {
            Some((self.start - self.meth_up_bp, self.start))
        } else {
            Some((self.end, self.end + self.meth_up_bp))
        }
    }

    pub fn down_window(&self) -> Option<(i64, i64)> {
        if self.meth_down_bp <= 0 {
            return None;
        }
        if self.gene_strand == "+" {
            Some((self.end, self.end + self.meth_down_bp))
        } else {
            Some((self.start - self.meth_down_bp, self.start))
        }
    }

    pub fn ploidy(&self, sex: &str) -> usize {
        match self.chrom.as_str() {
            "chrX" | "X" => {
                if sex == "XY" {
                    1
                } else {
                    2
                }
            }
            "chrY" | "Y" => {
                if sex == "XY" {
                    1
                } else {
                    0
                }
            }
            _ => 2,
        }
    }
}

/// What the overlay adds to one locus, in the shape `build_locus` wants.
#[derive(Clone, Debug, Default)]
pub struct OverlayEntry {
    pub aliases: Vec<String>,
    pub interruptions_gene: Vec<String>,
    pub meth_relevant: bool,
    pub meth_up_bp: i64,
    pub meth_down_bp: i64,
    pub meth_note: String,
    pub caveats: Vec<String>,
    pub bands: Vec<Band>,
    pub interruptions_in_size: Option<bool>,
}

/// The overlay's provenance record, as `load_overlay` returns it in `meta`.
#[derive(Clone, Debug)]
pub struct Provenance {
    pub version: Json,
    pub sources: Json,
    /// Per locus, the band mappings as written in the YAML (with an
    /// `_inhouse` marker where an in-house band replaced them).
    pub bands: Vec<(String, Vec<Json>)>,
    pub discrepancies: Vec<Json>,
    pub inhouse_enabled: bool,
}

/// `" ".join(s.split())`.
pub fn squash(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}

fn str_list(v: Option<&Json>) -> Vec<String> {
    v.map(|l| {
        l.as_list()
            .iter()
            .filter_map(|x| x.as_str().map(String::from))
            .collect()
    })
    .unwrap_or_default()
}

/// Python's `int(x)` for the overlay's small integers.
fn py_int(v: Option<&Json>) -> i64 {
    match v {
        Some(Json::Int(i)) => *i,
        Some(Json::Float(f)) => f.trunc() as i64,
        Some(Json::Bool(b)) => *b as i64,
        Some(Json::Str(s)) => s.trim().parse().unwrap_or(0),
        _ => 0,
    }
}

pub fn load_overlay(inhouse: bool) -> Result<(HashMap<String, OverlayEntry>, Provenance), String> {
    let docs = yaml_rust2::YamlLoader::load_from_str(OVERLAY_YAML).map_err(|e| format!("overlay.yaml: {e}"))?;
    let doc = docs.first().map(Json::from_yaml).unwrap_or(Json::Null);
    let inhouse_doc = doc.get_some("inhouse").cloned().unwrap_or(Json::Obj(vec![]));
    let mut ov = HashMap::new();
    let mut prov_bands = Vec::new();
    for (lid, rec) in doc.get_some("loci").map(|l| l.entries()).unwrap_or(&[]) {
        let mut bands: Vec<Json> = rec
            .get("bands")
            .filter(|b| b.truthy())
            .map(|b| b.as_list().to_vec())
            .unwrap_or_default();
        if inhouse && inhouse_doc.has(lid) {
            let ih = inhouse_doc.get(lid).unwrap();
            let src = ih
                .get("bands")
                .filter(|b| b.truthy())
                .map(|b| b.as_list().to_vec())
                .unwrap_or(bands.clone());
            bands = src
                .into_iter()
                .map(|b| {
                    let mut items = b.entries().to_vec();
                    items.retain(|(k, _)| k != "_inhouse");
                    items.push(("_inhouse".into(), Json::Bool(true)));
                    Json::Obj(items)
                })
                .collect();
        }
        let meth = rec
            .get("methylation")
            .filter(|m| m.truthy())
            .cloned()
            .unwrap_or(Json::Obj(vec![]));
        let entry = OverlayEntry {
            aliases: str_list(rec.get("aliases")),
            interruptions_gene: str_list(rec.get("interruptions_gene")),
            meth_relevant: meth.get("relevant").map(Json::truthy).unwrap_or(false),
            meth_up_bp: py_int(meth.get("up_bp").filter(|v| v.truthy())),
            meth_down_bp: py_int(meth.get("down_bp").filter(|v| v.truthy())),
            meth_note: meth
                .get("note")
                .filter(|v| v.truthy())
                .and_then(|v| v.as_str())
                .unwrap_or("")
                .to_string(),
            caveats: rec
                .get("caveats")
                .map(|c| c.as_list().iter().filter_map(|x| x.as_str()).map(squash).collect())
                .unwrap_or_default(),
            bands: bands
                .iter()
                .map(|b| Band {
                    upper: b.get_some("upper").and_then(Json::as_f64).unwrap_or(BAND_TOP),
                    label: format!(
                        "{}{}",
                        b.get("label").and_then(Json::as_str).unwrap_or(""),
                        if b.get("_inhouse").map(Json::truthy).unwrap_or(false) {
                            " [in-house]"
                        } else {
                            ""
                        }
                    ),
                })
                .collect(),
            interruptions_in_size: rec.get("interruptions_in_size").map(Json::truthy),
        };
        ov.insert(lid.clone(), entry);
        prov_bands.push((lid.clone(), bands));
    }
    let meta = Provenance {
        version: doc.get("version").cloned().unwrap_or(Json::Null),
        sources: doc
            .get_some("sources")
            .cloned()
            .filter(Json::truthy)
            .unwrap_or(Json::Obj(vec![])),
        bands: prov_bands,
        discrepancies: doc
            .get_some("discrepancies")
            .map(|d| d.as_list().to_vec())
            .unwrap_or_default(),
        inhouse_enabled: inhouse,
    };
    Ok((ov, meta))
}

fn same_len(motifs: Option<&Json>, unit_len: usize) -> Vec<String> {
    motifs
        .map(|m| {
            m.as_list()
                .iter()
                .filter_map(|x| x.as_str())
                .filter(|s| s.chars().count() == unit_len)
                .map(|s| s.to_ascii_uppercase())
                .collect()
        })
        .unwrap_or_default()
}

fn regime(pathogenic_max: Option<Num>, unit_len: usize) -> &'static str {
    match pathogenic_max {
        None => "expansion",
        Some(p) => {
            if p.f() * unit_len as f64 > EXPANSION_REGIME_BP {
                "expansion"
            } else {
                "short"
            }
        }
    }
}

fn to_ref(motifs_gene: &[String], strand: &str) -> Vec<String> {
    motifs_gene
        .iter()
        .map(|m| {
            if strand == "+" {
                m.to_ascii_uppercase()
            } else {
                revcomp(m)
            }
        })
        .collect()
}

fn num(rec: &Json, key: &str) -> Option<Num> {
    rec.get(key).and_then(Num::from_json)
}

fn int_field(v: &Json) -> Option<i64> {
    match v {
        Json::Int(i) => Some(*i),
        Json::Float(f) => Some(f.trunc() as i64),
        _ => None,
    }
}

pub fn build_locus(rec: &Json, ov: Option<&OverlayEntry>, build: &str) -> Result<Option<Locus>, String> {
    let (start, stop) = match (
        rec.get_some(&format!("start_{build}")),
        rec.get_some(&format!("stop_{build}")),
    ) {
        (Some(a), Some(b)) => (a, b),
        _ => return Ok(None),
    };
    let strand = rec
        .get("gene_strand")
        .filter(|s| s.truthy())
        .and_then(Json::as_str)
        .unwrap_or("+")
        .to_string();
    let unit_len = rec.get("motif_len").and_then(int_field).ok_or("motif_len missing")? as usize;
    let mut ref_motifs = same_len(rec.get_some("reference_motif_reference_orientation"), unit_len);
    if ref_motifs.is_empty() {
        ref_motifs = same_len(rec.get_some("pathogenic_motif_reference_orientation"), unit_len);
    }
    if ref_motifs.is_empty() {
        return Ok(None);
    }
    let interruptions = match ov {
        Some(o) => to_ref(&o.interruptions_gene, &strand),
        None => same_len(rec.get_some("interruption_reference_orientation"), unit_len),
    };
    let (meth_relevant, meth_up_bp, meth_down_bp, meth_note) = match ov {
        Some(o) => (o.meth_relevant, o.meth_up_bp, o.meth_down_bp, o.meth_note.clone()),
        None => {
            let cpg = ref_motifs.iter().any(|m| m.repeat(2).contains("CG"));
            (
                cpg,
                0,
                0,
                if cpg {
                    "default: tract 5mC reported because the repeat unit contains CpG".to_string()
                } else {
                    String::new()
                },
            )
        }
    };
    let gene = rec
        .get("gene")
        .and_then(Json::as_str)
        .ok_or("gene missing")?
        .to_string();
    let mut aliases: Vec<String> = ov.map(|o| o.aliases.clone()).unwrap_or_default();
    aliases.push(gene.clone());
    let pathogenic_motifs = same_len(rec.get_some("pathogenic_motif_reference_orientation"), unit_len);
    let benign_motifs = same_len(rec.get_some("benign_motif_reference_orientation"), unit_len);
    let motifset = MotifSet::new(
        unit_len,
        &ref_motifs,
        &pathogenic_motifs,
        &benign_motifs,
        &interruptions,
    )?;
    Ok(Some(Locus {
        id: rec.get("id").and_then(Json::as_str).ok_or("id missing")?.to_string(),
        gene,
        disease: rec.get("disease").and_then(Json::as_str).unwrap_or("").to_string(),
        chrom: rec
            .get("chrom")
            .and_then(Json::as_str)
            .ok_or("chrom missing")?
            .to_string(),
        start: int_field(start).ok_or("bad start")?,
        end: int_field(stop).ok_or("bad stop")?,
        gene_strand: strand,
        unit_len,
        ref_motifs,
        pathogenic_motifs,
        benign_motifs,
        interruption_motifs: interruptions,
        benign_min: num(rec, "benign_min"),
        benign_max: num(rec, "benign_max"),
        intermediate_min: num(rec, "intermediate_min"),
        intermediate_max: num(rec, "intermediate_max"),
        pathogenic_min: num(rec, "pathogenic_min"),
        pathogenic_max: num(rec, "pathogenic_max"),
        ref_copies: num(rec, "ref_copies"),
        inheritance: rec
            .get("inheritance")
            .map(|l| {
                l.as_list()
                    .iter()
                    .map(|x| x.as_str().map(String::from).unwrap_or_else(|| x.py_str()))
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .unwrap_or_default(),
        mechanism: rec
            .get("mechanism")
            .filter(|m| m.truthy())
            .and_then(Json::as_str)
            .unwrap_or("")
            .to_string(),
        regime: regime(num(rec, "pathogenic_max"), unit_len).to_string(),
        meth_relevant,
        meth_up_bp,
        meth_down_bp,
        meth_note,
        build: build.to_string(),
        aliases,
        caveats: ov.map(|o| o.caveats.clone()).unwrap_or_default(),
        interruptions_in_size: ov.and_then(|o| o.interruptions_in_size).unwrap_or(true),
        bands: ov.map(|o| o.bands.clone()).unwrap_or_default(),
        motifset,
    }))
}

pub struct Catalog {
    pub build: String,
    pub loci: Vec<Locus>,
    alias: HashMap<String, usize>,
    pub dropped: Vec<String>,
    pub provenance: Provenance,
    /// How this catalogue is named in outputs; see `catalogue_label`.
    pub label: String,
}

impl Catalog {
    /// `path = None` reads the catalogue compiled into this binary.
    pub fn load(path: Option<&Path>, build: &str, inhouse: bool) -> Result<Catalog, String> {
        if !BUILDS.contains(&build) {
            return Err(format!(
                "unknown build '{build}'; expected one of {}",
                BUILDS.join(", ")
            ));
        }
        let owned;
        let bytes: &[u8] = match path {
            None => BUNDLED_CATALOG.as_bytes(),
            Some(p) => {
                owned = std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()))?;
                &owned
            }
        };
        let label = catalogue_label(path.map(|p| p.to_str().unwrap_or("?")), bytes);
        let value: serde_json::Value =
            serde_json::from_slice(bytes).map_err(|e| format!("catalogue is not valid JSON: {e}"))?;
        let records = Json::from_serde(&value);
        let (overlay, provenance) = load_overlay(inhouse)?;
        let mut cat = Catalog {
            build: build.to_string(),
            loci: Vec::new(),
            alias: HashMap::new(),
            dropped: Vec::new(),
            provenance,
            label,
        };
        for rec in records.as_list() {
            let id = rec.get("id").and_then(Json::as_str).unwrap_or("?").to_string();
            match build_locus(rec, overlay.get(&id), build).map_err(|e| format!("{id}: {e}"))? {
                None => cat.dropped.push(id),
                Some(loc) => {
                    // A later record with the same id replaces the earlier one,
                    // as a Python dict keyed by id does; aliases keep the first.
                    let idx = match cat.loci.iter().position(|l| l.id == loc.id) {
                        Some(i) => {
                            cat.loci[i] = loc;
                            i
                        }
                        None => {
                            cat.loci.push(loc);
                            cat.loci.len() - 1
                        }
                    };
                    let loc = &cat.loci[idx];
                    let keys: Vec<String> = std::iter::once(loc.id.clone())
                        .chain(loc.aliases.iter().cloned())
                        .collect();
                    for a in keys {
                        let id = cat.loci[idx].id.clone();
                        cat.alias
                            .entry(a.to_uppercase())
                            .or_insert_with(|| cat.loci.iter().position(|l| l.id == id).unwrap());
                    }
                }
            }
        }
        Ok(cat)
    }

    pub fn get(&self, key: &str) -> Result<&Locus, String> {
        self.alias
            .get(&key.to_uppercase())
            .map(|&i| &self.loci[i])
            .ok_or_else(|| format!("unknown locus '{key}'; try `expansionscout loci`"))
    }

    pub fn select(&self, keys: Option<&[String]>) -> Result<Vec<&Locus>, String> {
        match keys {
            None | Some([]) => Ok(self.loci.iter().collect()),
            Some(ks) => ks.iter().map(|k| self.get(k)).collect(),
        }
    }
}
