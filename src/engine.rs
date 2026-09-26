//! Per-locus orchestration: reads -> decomposition -> methylation ->
//! alleles -> classification. Port of `expansionscout/engine.py`.

use std::path::Path;

use crate::alleles::{cluster_counts, spread_stats};
use crate::bam::{Bam, Fasta};
use crate::catalog::Locus;
use crate::classify::{
    class_units, classify_units, detectability, evidence_state, lower_bound, negative_reliable, over_dispersed,
    size_units, SizeVal, EVIDENCE_NONE, NO_CALL,
};
use crate::consensus::{allele_consensus, shortfall, ConsensusRead};
use crate::decompose::{
    Decomposition, MotifSet, LABEL_BENIGN, LABEL_CANONICAL, LABEL_INTERRUPTION, LABEL_OTHER, LABEL_PATHOGENIC,
};
use crate::methylation::{read_methylation, ReadMeth, WindowMeth};
use crate::npstat::{mean, median, py_round};
use crate::pyfmt::{fixed, Num};
use crate::reads::{self as R, ReadObs};

pub const PER_READ_COLS: [&str; 34] = [
    "locus",
    "read_id",
    "read_class",
    "hp",
    "allele",
    "strand",
    "mapq",
    "read_length",
    "size_bp",
    "tract_len_bp",
    "n_total_units",
    "n_canonical",
    "n_interruption",
    "n_pathogenic",
    "n_benign",
    "n_other",
    "units_estimated",
    "lower_bound_units",
    "interruption_motifs_observed",
    "interruption_positions",
    "structure_string",
    "frac_5mc",
    "n_5mc_valid",
    "n_5mc_fail",
    "mean_5mc_prob",
    "n_5mc_sites",
    "cpg_5mc",
    "frac_5mc_up",
    "n_5mc_up_valid",
    "n_5mc_up_fail",
    "mean_5mc_up",
    "n_5mc_up",
    "mean_5mc_down",
    "n_5mc_down",
];

pub const HAP_KEYS: [&str; 23] = [
    "median",
    "median_measured",
    "p5",
    "p95",
    "mad",
    "support",
    "size_bp",
    "tail_up_frac",
    "noise",
    "n_interruption",
    "int_pos",
    "n_pathogenic",
    "cls",
    "dispersed",
    "structure",
    "meth_tract",
    "meth_tract_pooled",
    "meth_hydroxy",
    "meth_fail_frac",
    "meth_up",
    "meth_up_pooled",
    "meth_down",
    "frac_meth_reads_gt50",
];

const LOCI_FIXED: [&str; 20] = [
    "sample",
    "locus",
    "gene",
    "chrom",
    "start",
    "end",
    "gene_strand",
    "motif_gene",
    "motif_ref",
    "unit_len",
    "regime",
    "composition",
    "ref_copies",
    "ploidy",
    "method",
    "n_reads_total",
    "n_spanning",
    "n_split",
    "n_partial",
    "n_no_run",
];
const LOCI_TAIL: [&str; 9] = [
    "minor_modes",
    "lower_bound_units",
    "lower_bound_support",
    "detectability",
    "evidence",
    "negative_reliable",
    "tagged_fraction",
    "meth_convention",
    "notes",
];

#[derive(Clone, Debug)]
pub struct Options {
    pub anchor: i64,
    pub flank: i64,
    pub min_mapq: i64,
    pub min_clip: i64,
    pub margin: i64,
    pub min_support: i64,
    pub min_frac: f64,
    pub mosaic: bool,
    pub use_hp: bool,
    pub sex: String,
    pub meth_threshold: f64,
    pub meth_convention: String,
    pub clinical: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            anchor: 50,
            flank: 100,
            min_mapq: 20,
            min_clip: 50,
            margin: 12,
            min_support: 3,
            min_frac: 0.10,
            mosaic: false,
            use_hp: true,
            sex: "XX".into(),
            meth_threshold: 0.8,
            meth_convention: "modkit".into(),
            clinical: false,
        }
    }
}

/// A value as Python's dynamic typing holds it, for `_fmt`.
#[derive(Clone, Debug, PartialEq)]
pub enum PyVal {
    None,
    Int(i64),
    Float(f64),
    Str(String),
}

impl PyVal {
    pub fn opt_f(v: Option<f64>) -> PyVal {
        v.map(PyVal::Float).unwrap_or(PyVal::None)
    }
    pub fn opt_i(v: Option<i64>) -> PyVal {
        v.map(PyVal::Int).unwrap_or(PyVal::None)
    }
    pub fn opt_s(v: Option<&str>) -> PyVal {
        v.map(|s| PyVal::Str(s.to_string())).unwrap_or(PyVal::None)
    }
    pub fn size(v: SizeVal) -> PyVal {
        match v {
            SizeVal::Float(f) => PyVal::Float(f),
            SizeVal::Zero => PyVal::Int(0),
        }
    }
    pub fn num(v: Option<Num>) -> PyVal {
        match v {
            None => PyVal::None,
            Some(Num::Int(i)) => PyVal::Int(i),
            Some(Num::Float(f)) => PyVal::Float(f),
        }
    }
    /// `_fmt(v, nd)`.
    pub fn fmt(&self, nd: usize) -> String {
        match self {
            PyVal::None => ".".into(),
            PyVal::Float(f) => fixed(*f, nd),
            PyVal::Int(i) => i.to_string(),
            PyVal::Str(s) => s.clone(),
        }
    }
    pub fn is_none(&self) -> bool {
        matches!(self, PyVal::None)
    }
}

pub fn fmt_opt_f(v: Option<f64>, nd: usize) -> String {
    PyVal::opt_f(v).fmt(nd)
}

/// What `_summarise` returns for one allele.
#[derive(Clone, Debug)]
pub struct Hap {
    pub median: SizeVal,
    pub median_measured: f64,
    pub p5: f64,
    pub p95: f64,
    pub mad: f64,
    pub cnr: String,
    pub support: usize,
    pub size_bp: Option<i64>,
    pub tail_up_frac: f64,
    pub tail_down_frac: f64,
    pub noise: Option<f64>,
    pub n_interruption: Option<i64>,
    pub int_pos: Option<String>,
    pub n_pathogenic: Option<i64>,
    pub cls: String,
    pub dispersed: bool,
    pub structure: Option<String>,
    pub seq: Option<String>,
    pub motif: Option<String>,
    pub meth_tract: Option<f64>,
    pub meth_tract_pooled: Option<f64>,
    pub meth_up_pooled: Option<f64>,
    pub meth_hydroxy: Option<f64>,
    pub meth_fail_frac: Option<f64>,
    pub meth_up: Option<f64>,
    pub meth_down: Option<f64>,
    pub frac_meth_reads_gt50: Option<f64>,
    pub n_meth_reads: usize,
    pub hp: Option<i64>,
}

impl Hap {
    /// `h.get(key)` for the keys the tables and the VCF read.
    pub fn get(&self, key: &str) -> PyVal {
        match key {
            "median" => PyVal::size(self.median),
            "median_measured" => PyVal::Float(self.median_measured),
            "p5" => PyVal::Float(self.p5),
            "p95" => PyVal::Float(self.p95),
            "mad" => PyVal::Float(self.mad),
            "cnr" => PyVal::Str(self.cnr.clone()),
            "support" => PyVal::Int(self.support as i64),
            "size_bp" => PyVal::opt_i(self.size_bp),
            "tail_up_frac" => PyVal::Float(self.tail_up_frac),
            "tail_down_frac" => PyVal::Float(self.tail_down_frac),
            "noise" => PyVal::opt_f(self.noise),
            "n_interruption" => PyVal::opt_i(self.n_interruption),
            "int_pos" => PyVal::opt_s(self.int_pos.as_deref()),
            "n_pathogenic" => PyVal::opt_i(self.n_pathogenic),
            "cls" => PyVal::Str(self.cls.clone()),
            "dispersed" => PyVal::Int(self.dispersed as i64),
            "structure" => PyVal::opt_s(self.structure.as_deref()),
            "seq" => PyVal::opt_s(self.seq.as_deref()),
            "motif" => PyVal::opt_s(self.motif.as_deref()),
            "meth_tract" => PyVal::opt_f(self.meth_tract),
            "meth_tract_pooled" => PyVal::opt_f(self.meth_tract_pooled),
            "meth_up_pooled" => PyVal::opt_f(self.meth_up_pooled),
            "meth_hydroxy" => PyVal::opt_f(self.meth_hydroxy),
            "meth_fail_frac" => PyVal::opt_f(self.meth_fail_frac),
            "meth_up" => PyVal::opt_f(self.meth_up),
            "meth_down" => PyVal::opt_f(self.meth_down),
            "frac_meth_reads_gt50" => PyVal::opt_f(self.frac_meth_reads_gt50),
            "n_meth_reads" => PyVal::Int(self.n_meth_reads as i64),
            "hp" => PyVal::opt_i(self.hp),
            _ => PyVal::None,
        }
    }
}

pub struct LocusResult<'a> {
    pub locus: &'a Locus,
    pub sample: String,
    pub rows: Vec<Vec<String>>,
    pub haplotypes: Vec<Hap>,
    pub minor_modes: Vec<Hap>,
    pub method: String,
    pub n_reads_total: usize,
    pub n_spanning: i64,
    pub n_split: usize,
    pub n_partial: usize,
    pub n_no_run: usize,
    pub lb_units: Option<i64>,
    pub lb_support: usize,
    pub detect: Option<f64>,
    pub evidence: String,
    pub negative_reliable: bool,
    pub ref_tract_seq: Option<String>,
    pub ref_pad_base: Option<String>,
    pub ploidy: usize,
    pub meth_convention: String,
    pub tagged_fraction: Option<f64>,
    pub notes: Vec<String>,
    /// Messages htslib would have written to stderr while this locus was
    /// called (malformed modification tags), in the order it would have.
    pub diagnostics: Vec<String>,
}

fn empty_result<'a>(locus: &'a Locus, sample: &str, opts: &Options, note: String) -> LocusResult<'a> {
    LocusResult {
        locus,
        sample: sample.to_string(),
        rows: vec![],
        haplotypes: vec![],
        minor_modes: vec![],
        method: "none".into(),
        n_reads_total: 0,
        n_spanning: 0,
        n_split: 0,
        n_partial: 0,
        n_no_run: 0,
        lb_units: None,
        lb_support: 0,
        detect: None,
        evidence: EVIDENCE_NONE.into(),
        negative_reliable: false,
        ref_tract_seq: None,
        ref_pad_base: None,
        ploidy: locus.ploidy(&opts.sex),
        meth_convention: opts.meth_convention.clone(),
        tagged_fraction: None,
        notes: vec![note],
        diagnostics: vec![],
    }
}

fn col(name: &str) -> usize {
    PER_READ_COLS.iter().position(|c| *c == name).expect("known column")
}

struct Row(Vec<String>);

impl Row {
    fn set(&mut self, name: &str, v: impl Into<String>) {
        self.0[col(name)] = v.into();
    }
}

/// A spanning or split read that counts toward the alleles.
struct Call {
    idx: usize,
    units: i64,
    hp: i64,
    dec: Option<(Decomposition, Vec<u8>)>,
    meth: ReadMeth,
    seq: Option<String>,
    size_bp: i64,
    estimated: bool,
}

type Window = Option<(i64, i64)>;

fn meth_windows(o: &ReadObs, locus: &Locus) -> (Window, Window) {
    let rec = &o.record;
    let get = |win: Option<(i64, i64)>| {
        let (a, b) = win?;
        if rec.pos <= a && rec.reference_end().is_some_and(|e| e >= b) {
            o.pairs.window(a, b)
        } else {
            None
        }
    };
    (get(locus.up_window()), get(locus.down_window()))
}

fn crossed_units(locus: &Locus, o: &ReadObs) -> i64 {
    let Some(seg) = o.segment.as_ref().filter(|s| !s.is_empty()) else {
        return 0;
    };
    let upper = seg.to_ascii_uppercase();
    let ms = locus.motifset();
    match o.seg_tract {
        None => ms.count_units_loose(&upper) as i64,
        Some(w) => ms.count_units_crossed(&upper, Some(w), o.read_class != R::RIGHT_PARTIAL) as i64,
    }
}

fn fill_decomp(row: &mut Row, dec: &Decomposition, useq: &[u8], locus: &Locus, partial: bool) {
    row.set("n_canonical", dec.count(LABEL_CANONICAL).to_string());
    row.set("n_interruption", dec.count(LABEL_INTERRUPTION).to_string());
    row.set("n_pathogenic", dec.count(LABEL_PATHOGENIC).to_string());
    row.set("n_benign", dec.count(LABEL_BENIGN).to_string());
    row.set("n_other", dec.count(LABEL_OTHER).to_string());
    if !partial {
        row.set("n_total_units", dec.n_units().to_string());
    }
    let pos = dec.positions(LABEL_INTERRUPTION);
    if pos.is_empty() {
        row.set("interruption_positions", ".");
        row.set("interruption_motifs_observed", ".");
    } else {
        row.set(
            "interruption_positions",
            pos.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(","),
        );
        row.set(
            "interruption_motifs_observed",
            pos.iter()
                .map(|&p| locus.to_gene_strand(&String::from_utf8_lossy(dec.unit(useq, p - 1, locus.unit_len))))
                .collect::<Vec<_>>()
                .join(","),
        );
    }
    row.set(
        "structure_string",
        if dec.labels.is_empty() {
            ".".to_string()
        } else {
            dec.labels_str()
        },
    );
}

fn fill_meth(row: &mut Row, m: &ReadMeth) {
    row.set("frac_5mc", fmt_opt_f(m.tract.frac, 3));
    row.set("n_5mc_valid", m.tract.n_valid.to_string());
    row.set("n_5mc_fail", m.tract.n_fail.to_string());
    row.set("mean_5mc_prob", fmt_opt_f(m.tract.mean, 3));
    row.set("n_5mc_sites", m.tract.n.to_string());
    row.set(
        "cpg_5mc",
        if m.per_site.is_empty() {
            ".".to_string()
        } else {
            m.per_site
                .iter()
                .map(|(o, p)| format!("{o}:{p}"))
                .collect::<Vec<_>>()
                .join(",")
        },
    );
    row.set("frac_5mc_up", fmt_opt_f(m.up.frac, 3));
    row.set("n_5mc_up_valid", m.up.n_valid.to_string());
    row.set("n_5mc_up_fail", m.up.n_fail.to_string());
    row.set("mean_5mc_up", fmt_opt_f(m.up.mean, 3));
    row.set("n_5mc_up", m.up.n.to_string());
    row.set("mean_5mc_down", fmt_opt_f(m.down.mean, 3));
    row.set("n_5mc_down", m.down.n.to_string());
}

fn opt_mean(v: &[f64]) -> Option<f64> {
    if v.is_empty() {
        None
    } else {
        Some(mean(v))
    }
}

/// `int(np.median(values))`: the median, truncated toward zero.
fn int_median(v: &[f64]) -> i64 {
    median(v).trunc() as i64
}

/// The allele's repeat unit: the commonest unit; on a tie, an exact motif
/// rotation before a unit that is not, then the first seen (the Python's
/// `_modal_unit`).
fn modal_unit<'a>(units: impl Iterator<Item = &'a [u8]>, ms: &MotifSet) -> Option<&'a [u8]> {
    let mut counts: Vec<(&[u8], usize)> = Vec::new();
    for u in units {
        match counts.iter_mut().find(|e| e.0 == u) {
            Some(e) => e.1 += 1,
            None => counts.push((u, 1)),
        }
    }
    let key = |c: &(&[u8], usize)| (c.1, ms.match_unit(c.0).is_some());
    let mut best = *counts.first()?;
    for &c in &counts[1..] {
        if key(&c) > key(&best) {
            best = c;
        }
    }
    Some(best.0)
}

fn summarise(members: &[&Call], locus: &Locus, clinical: bool) -> Hap {
    let units: Vec<f64> = members.iter().map(|m| m.units as f64).collect();
    let st = spread_stats(&units).expect("an allele has reads");
    let med = st.median;
    let exact: Vec<&&Call> = members.iter().filter(|m| m.dec.is_some()).collect();
    fn dec(m: &Call) -> &Decomposition {
        &m.dec.as_ref().unwrap().0
    }
    // Nearest the median count; then not a fragment, by how far the
    // decomposition falls short of the read's own count beyond max(1, 10 %)
    // of it; then the least noise; then the decomposition nearest the
    // allele's typical one: see the Python's _summarise.
    let typical = if exact.is_empty() {
        0.0
    } else {
        median(&exact.iter().map(|m| dec(m).n_units() as f64).collect::<Vec<_>>())
    };
    let rep_key = |m: &Call| {
        let (u, n) = (m.units, dec(m).n_units());
        (
            (u as f64 - med).abs(),
            shortfall(u, n),
            dec(m).noise_frac(),
            (n as f64 - typical).abs(),
        )
    };
    let mut rep_i: Option<usize> = None;
    for (i, m) in exact.iter().enumerate() {
        match rep_i {
            None => rep_i = Some(i),
            Some(r) => {
                if rep_key(m).partial_cmp(&rep_key(exact[r])) == Some(std::cmp::Ordering::Less) {
                    rep_i = Some(i);
                }
            }
        }
    }
    let rep: Option<&&Call> = rep_i.map(|i| exact[i]);
    // The allele's sequence is a consensus of its reads' runs where one can
    // be formed, and its structure and repeat unit are that sequence's own
    // decomposition; else they are the representative read's. See the
    // Python's _summarise and consensus.rs.
    let ms = locus.motifset();
    let modal_of = |d: &Decomposition, useq: &[u8]| {
        if d.labels.is_empty() {
            return None;
        }
        let units = (0..d.labels.len()).map(|i| d.unit(useq, i, locus.unit_len));
        modal_unit(units, ms).map(|u| String::from_utf8_lossy(u).into_owned())
    };
    let (mut seq, mut structure, mut motif) = (None, None, None);
    if let Some(r) = rep {
        let (d, useq) = r.dec.as_ref().unwrap();
        seq = r.seq.clone();
        structure = Some(d.labels_str());
        motif = modal_of(d, useq);
        let reads: Vec<ConsensusRead> = exact
            .iter()
            .map(|m| {
                let (d, useq) = m.dec.as_ref().unwrap();
                ConsensusRead {
                    seq: m
                        .seq
                        .as_deref()
                        .expect("a read with a decomposition has its run")
                        .as_bytes(),
                    run_u: &useq[d.start..d.end],
                    dec: d,
                    units: m.units,
                }
            })
            .collect();
        if let Some(cons) = allele_consensus(&reads, ms, rep_i.unwrap()) {
            let cons_u = cons.to_ascii_uppercase();
            let n = cons_u.len();
            if let Some(cd) = ms.decompose(&cons_u, Some((0, n)), Some(n)) {
                if !cd.labels.is_empty() && 2 * cd.labels.len() >= d.labels.len() {
                    structure = Some(cd.labels_str());
                    motif = modal_of(&cd, &cons_u);
                }
            }
            seq = Some(String::from_utf8_lossy(&cons).into_owned());
        }
    }
    let rmed = py_round(med);
    let modal: Vec<&&&Call> = exact.iter().filter(|m| (m.units as f64 - rmed).abs() <= 0.0).collect();
    let mut pos_counts: Vec<(usize, usize)> = Vec::new();
    for m in &modal {
        for p in dec(m).positions(LABEL_INTERRUPTION) {
            match pos_counts.iter_mut().find(|e| e.0 == p) {
                Some(e) => e.1 += 1,
                None => pos_counts.push((p, 1)),
            }
        }
    }
    let need = 0.5 * 1.max(modal.len()) as f64;
    let mut int_pos: Vec<usize> = pos_counts.iter().filter(|e| e.1 as f64 >= need).map(|e| e.0).collect();
    int_pos.sort();
    let counts = |label: u8| -> Vec<f64> { exact.iter().map(|m| dec(m).count(label) as f64).collect() };
    let n_int = if exact.is_empty() {
        None
    } else {
        Some(int_median(&counts(LABEL_INTERRUPTION)))
    };
    let n_path = if exact.is_empty() {
        None
    } else {
        Some(int_median(&counts(LABEL_PATHOGENIC)))
    };
    let n_ben = if exact.is_empty() {
        0
    } else {
        int_median(&counts(LABEL_BENIGN))
    };
    let noise = if exact.is_empty() {
        None
    } else {
        Some(mean(&exact.iter().map(|m| dec(m).noise_frac()).collect::<Vec<_>>()))
    };
    let size_bp = if members.is_empty() {
        None
    } else {
        Some(int_median(
            &members.iter().map(|m| m.size_bp as f64).collect::<Vec<_>>(),
        ))
    };
    let tract: Vec<f64> = members.iter().filter_map(|m| m.meth.tract.frac).collect();
    let hyd: Vec<f64> = members
        .iter()
        .filter(|m| m.meth.tract.n_valid > 0)
        .map(|m| m.meth.tract.n_hydroxy as f64 / m.meth.tract.n_valid as f64)
        .collect();
    let tract_reads: Vec<f64> = members
        .iter()
        .filter(|m| m.meth.tract.frac.is_some() && m.meth.tract.n_valid >= 3)
        .map(|m| m.meth.tract.frac.unwrap())
        .collect();
    let up: Vec<f64> = members.iter().filter_map(|m| m.meth.up.frac).collect();
    let down: Vec<f64> = members.iter().filter_map(|m| m.meth.down.frac).collect();
    let pooled = |w: fn(&ReadMeth) -> &WindowMeth| -> Option<f64> {
        let (mut nm, mut nv) = (0usize, 0usize);
        for m in members {
            let win = w(&m.meth);
            if win.n_valid == 0 {
                continue;
            }
            nm += win.n_mod;
            nv += win.n_valid;
        }
        if nv > 0 {
            Some(nm as f64 / nv as f64)
        } else {
            None
        }
    };
    let fails: Vec<f64> = members.iter().filter_map(|m| m.meth.tract.fail_frac()).collect();
    let reported = size_units(locus, med, n_int.unwrap_or(0), clinical);
    let cls_units = class_units(locus, reported.f(), n_path.unwrap_or(0), n_ben);
    let dispersed = over_dispersed(Some(med), Some(st.p95), Some(st.tail_up_frac), members.len());
    let cls = if dispersed {
        NO_CALL.to_string()
    } else {
        classify_units(locus, Some(cls_units))
    };
    Hap {
        median: reported,
        median_measured: med,
        p5: st.p5,
        p95: st.p95,
        mad: st.mad,
        cnr: format!("{}-{}", fixed(st.p5, 0), fixed(st.p95, 0)),
        support: members.len(),
        size_bp,
        tail_up_frac: st.tail_up_frac,
        tail_down_frac: st.tail_down_frac,
        noise,
        n_interruption: n_int,
        int_pos: if int_pos.is_empty() {
            None
        } else {
            Some(int_pos.iter().map(|p| p.to_string()).collect::<Vec<_>>().join("/"))
        },
        n_pathogenic: n_path,
        cls,
        dispersed,
        structure,
        seq,
        motif,
        meth_tract: opt_mean(&tract),
        meth_tract_pooled: pooled(|m| &m.tract),
        meth_up_pooled: pooled(|m| &m.up),
        meth_hydroxy: opt_mean(&hyd),
        meth_fail_frac: opt_mean(&fails),
        meth_up: opt_mean(&up),
        meth_down: opt_mean(&down),
        frac_meth_reads_gt50: if tract_reads.is_empty() {
            None
        } else {
            Some(tract_reads.iter().filter(|&&v| v >= 0.5).count() as f64 / tract_reads.len() as f64)
        },
        n_meth_reads: tract_reads.len(),
        hp: None,
    }
}

pub fn call_locus<'a>(
    bam: &mut Bam,
    fasta: Option<&mut Fasta>,
    locus: &'a Locus,
    sample: &str,
    opts: &Options,
) -> Result<LocusResult<'a>, String> {
    let Some(contig) = R::resolve_contig(&bam.references, &locus.chrom) else {
        return Ok(empty_result(
            locus,
            sample,
            opts,
            format!("contig {} not in the BAM header", locus.chrom),
        ));
    };
    let ms = locus.motifset();
    let margin = opts.margin.max(2 * locus.unit_len as i64);
    let mut fasta = fasta;
    let fa_contig = fasta.as_deref().and_then(|fa| {
        [&locus.chrom, &contig]
            .into_iter()
            .find(|c| fa.references.contains(c))
            .cloned()
    });
    // How far the reference flanks continue the tract's period: see the
    // Python's call_locus.
    let mut period = (0, 0);
    if let (Some(fa), Some(c)) = (fasta.as_deref_mut(), &fa_contig) {
        let reach = opts.flank.max(margin);
        let lo = 0.max(locus.start - reach);
        let s = fa.fetch(c, lo, locus.end + reach).map_err(|e| e.to_string())?;
        period = R::period_flanks(
            &s.to_ascii_uppercase(),
            lo,
            locus.start,
            locus.end,
            locus.unit_len as i64,
            reach,
        );
    }
    let lr = R::collect_reads(
        bam,
        &contig,
        locus.start,
        locus.end,
        opts.anchor,
        opts.flank,
        opts.min_mapq,
        opts.min_clip,
        margin,
        period,
        locus.unit_len > ms.approx_above_len,
    )?;
    let mut rows: Vec<Row> = Vec::new();
    let mut calls: Vec<Call> = Vec::new();
    let mut partial_lb: Vec<i64> = Vec::new();
    let mut n_no_run = 0usize;
    let mut diag: Vec<String> = Vec::new();
    for o in &lr.obs {
        let mut row = Row(vec![".".to_string(); PER_READ_COLS.len()]);
        row.set("locus", locus.id.clone());
        row.set("read_id", o.name.clone());
        row.set("read_class", o.read_class);
        row.set("hp", o.hp.to_string());
        row.set("strand", o.strand);
        row.set("mapq", o.mapq.to_string());
        row.set("read_length", o.read_length.to_string());
        row.set("size_bp", PyVal::opt_i(o.size_bp_aln).fmt(3));
        row.set("units_estimated", "0");
        let useq = o.segment.as_ref().map(|s| s.to_ascii_uppercase());
        // A spanning read's alignment bounds the tract, by where it places it
        // and by its length between the anchors, and within those bounds a
        // long motif's edge units may be taken by identity; see
        // MotifSet::extend. Partial and split reads extend exactly.
        let (window, max_bp) = if o.read_class == R::SPANNING {
            (
                o.seg_tract.map(|(a, b)| (a as usize, b as usize)),
                o.size_bp_aln.map(|n| n as usize),
            )
        } else {
            (None, None)
        };
        let dec = useq
            .as_ref()
            .filter(|s| !s.is_empty())
            .and_then(|s| ms.decompose(s, window, max_bp));
        let meth;
        if o.read_class == R::SPANNING {
            let Some(dec) = dec else {
                n_no_run += 1;
                row.set("read_class", R::NO_RUN);
                rows.push(row);
                continue;
            };
            let useq = useq.unwrap();
            let run_q = (o.seg_q_left + dec.start as i64, o.seg_q_left + dec.end as i64);
            let (up, down) = meth_windows(o, locus);
            let m = read_methylation(
                &o.record,
                Some(run_q),
                up,
                down,
                opts.meth_threshold,
                &opts.meth_convention,
                &mut diag,
            );
            let size_bp_aln = o.size_bp_aln.unwrap_or(0);
            let units = if size_bp_aln != 0 {
                row.set("units_estimated", "1");
                1.max(py_round(size_bp_aln as f64 / locus.unit_len as f64) as i64)
            } else {
                dec.n_units() as i64
            };
            fill_decomp(&mut row, &dec, &useq, locus, false);
            row.set("n_total_units", units.to_string());
            row.set("tract_len_bp", dec.n_bp().to_string());
            let seg = o.segment.as_ref().unwrap();
            let seq = String::from_utf8_lossy(&seg[dec.start..dec.end]).into_owned();
            let size_bp = dec.n_bp() as i64;
            calls.push(Call {
                idx: rows.len(),
                units,
                hp: o.hp,
                meth: m.clone(),
                seq: Some(seq),
                size_bp,
                estimated: false,
                dec: Some((dec, useq)),
            });
            meth = Some(m);
        } else if o.read_class == R::SPLIT {
            let size_bp_aln = o.size_bp_aln.unwrap_or(0);
            let units = if size_bp_aln != 0 {
                Some(py_round(size_bp_aln as f64 / locus.unit_len as f64) as i64)
            } else {
                None
            };
            row.set("units_estimated", "1");
            let lb = crossed_units(locus, o);
            row.set("lower_bound_units", lb.to_string());
            if let Some(d) = &dec {
                fill_decomp(&mut row, d, useq.as_ref().unwrap(), locus, true);
            }
            row.set("n_total_units", PyVal::opt_i(units).fmt(3));
            let (up, down) = meth_windows(o, locus);
            let m = read_methylation(
                &o.record,
                None,
                up,
                down,
                opts.meth_threshold,
                &opts.meth_convention,
                &mut diag,
            );
            if let Some(u) = units.filter(|&u| u != 0) {
                calls.push(Call {
                    idx: rows.len(),
                    units: u,
                    hp: o.hp,
                    dec: None,
                    meth: m.clone(),
                    seq: None,
                    size_bp: size_bp_aln,
                    estimated: true,
                });
            }
            partial_lb.push(lb.max(units.unwrap_or(0)));
            meth = Some(m);
        } else {
            let lb = crossed_units(locus, o);
            row.set("lower_bound_units", lb.to_string());
            if let Some(d) = &dec {
                fill_decomp(&mut row, d, useq.as_ref().unwrap(), locus, true);
            }
            let (up, down) = meth_windows(o, locus);
            meth = Some(read_methylation(
                &o.record,
                None,
                up,
                down,
                opts.meth_threshold,
                &opts.meth_convention,
                &mut diag,
            ));
            partial_lb.push(lb);
        }
        if let Some(m) = &meth {
            fill_meth(&mut row, m);
        }
        rows.push(row);
    }

    let ploidy = locus.ploidy(&opts.sex);
    let mut notes: Vec<String> = Vec::new();
    let mut haps: Vec<Hap> = Vec::new();
    let mut minor: Vec<Hap> = Vec::new();
    let mut method = "none".to_string();
    let mut tagged_fraction = None;
    if ploidy == 0 {
        notes.push("locus skipped for this sex".into());
    } else if !calls.is_empty() {
        let exact: Vec<&Call> = calls.iter().filter(|c| !c.estimated).collect();
        let est: Vec<&Call> = calls.iter().filter(|c| c.estimated).collect();
        let pool: &Vec<&Call> = if exact.is_empty() { &est } else { &exact };
        if !pool.is_empty() {
            tagged_fraction = Some(pool.iter().filter(|c| c.hp == 1 || c.hp == 2).count() as f64 / pool.len() as f64);
        }
        let values: Vec<f64> = pool.iter().map(|c| c.units as f64).collect();
        let hps: Vec<i64> = pool.iter().map(|c| c.hp).collect();
        let asg = cluster_counts(
            &values,
            &hps,
            ploidy,
            opts.mosaic,
            opts.min_support.max(0) as usize,
            opts.min_frac,
            opts.use_hp,
        );
        method = asg.method.clone();
        notes.extend(asg.notes.iter().cloned());
        let mut members: Vec<Vec<&Call>> = asg
            .clusters
            .iter()
            .map(|c| c.members.iter().map(|&i| pool[i]).collect())
            .collect();
        for (c, &lab) in pool.iter().zip(&asg.labels) {
            rows[c.idx].set(
                "allele",
                if lab >= 0 {
                    (lab + 1).to_string()
                } else {
                    ".".to_string()
                },
            );
        }
        if !exact.is_empty() && !est.is_empty() {
            for c in &est {
                if asg.clusters.is_empty() {
                    break;
                }
                let d: Vec<f64> = asg
                    .clusters
                    .iter()
                    .map(|cl| (c.units as f64 - cl.median).abs())
                    .collect();
                let k = crate::npstat::argmin(&d);
                let m = asg.clusters[k].median;
                let lim = 0.1 * m;
                if (c.units as f64 - m).abs() <= if lim > 2.0 { lim } else { 2.0 } {
                    members[k].push(c);
                    rows[c.idx].set("allele", (k + 1).to_string());
                }
            }
        }
        for (i, cl) in asg.clusters.iter().enumerate() {
            let mut summ = summarise(&members[i], locus, opts.clinical);
            summ.hp = cl.hp;
            if cl.minor {
                minor.push(summ);
            } else {
                haps.push(summ);
            }
        }
        if haps.len() == 1 && ploidy == 2 {
            let h = haps[0].clone();
            if !h.dispersed {
                haps.push(h);
            } else {
                notes.push(format!(
                    "second allele not called: the reads assigned to allele 1 run to {} units \
                     against a median of {}, which is not one allele",
                    fixed(h.p95, 0),
                    fixed(h.median.f(), 0)
                ));
            }
        }
    } else {
        notes.push("no spanning or split reads".into());
    }

    if n_no_run > 0 && calls.is_empty() {
        notes.push(format!(
            "WARNING: {n_no_run} read(s) span the locus but no repeat run matched motif {} - check \
             the motif definition for this locus",
            locus.ref_motif()
        ));
    } else if n_no_run > 0 && n_no_run >= calls.len() {
        notes.push(format!(
            "WARNING: {n_no_run} spanning read(s) contained no recognisable repeat run, as many as \
             were used for the call"
        ));
    }

    let dispersed = haps.iter().any(|h| h.dispersed);
    for h in &haps {
        if h.dispersed {
            let band = classify_units(locus, Some(h.p95));
            let mut note = format!(
                "size not reported as one allele: reads run {} to {} units around a median of {}",
                fixed(h.p5, 0),
                fixed(h.p95, 0),
                fixed(h.median.f(), 0)
            );
            if band != NO_CALL {
                note.push_str(&format!(", and the upper end is in the {band} range"));
            }
            notes.push(note);
        }
    }

    let n_partial = lr
        .obs
        .iter()
        .filter(|o| o.read_class == R::LEFT_PARTIAL || o.read_class == R::RIGHT_PARTIAL)
        .count();
    let n_spanning = lr.obs.iter().filter(|o| o.read_class == R::SPANNING).count() as i64 - n_no_run as i64;
    let (lb_units, lb_support) = lower_bound(&partial_lb);
    let detect = detectability(locus, &lr.read_lengths, opts.flank, opts.anchor);
    let medians: Vec<f64> = haps.iter().map(|h| h.median.f()).collect();
    let evidence = evidence_state(locus, &medians, lb_units, lb_support, n_partial);
    let negrel = negative_reliable(locus, detect, evidence, n_spanning, opts.min_support, dispersed);

    let (mut ref_tract, mut pad) = (None, None);
    if let Some(fa) = fasta {
        if let Some(fa_contig) = fa_contig {
            let t = fa
                .fetch(&fa_contig, locus.start, locus.end)
                .map_err(|e| e.to_string())?;
            ref_tract = Some(String::from_utf8_lossy(&t).to_uppercase());
            if locus.start > 0 {
                let p = fa
                    .fetch(&fa_contig, locus.start - 1, locus.start)
                    .map_err(|e| e.to_string())?;
                pad = Some(String::from_utf8_lossy(&p).to_uppercase());
            }
        } else {
            notes.push(format!("contig {} not in the reference FASTA", locus.chrom));
        }
    }

    Ok(LocusResult {
        locus,
        sample: sample.to_string(),
        rows: rows.into_iter().map(|r| r.0).collect(),
        haplotypes: haps,
        minor_modes: minor,
        method,
        n_reads_total: lr.read_lengths.len(),
        n_spanning,
        n_split: lr.obs.iter().filter(|o| o.read_class == R::SPLIT).count(),
        n_partial,
        n_no_run,
        lb_units,
        lb_support,
        detect,
        evidence: evidence.to_string(),
        negative_reliable: negrel,
        ref_tract_seq: ref_tract,
        ref_pad_base: pad,
        ploidy,
        meth_convention: opts.meth_convention.clone(),
        tagged_fraction,
        notes,
        diagnostics: diag,
    })
}

// ----------------------------------------------------------------------------
// tables
// ----------------------------------------------------------------------------

pub fn loci_header() -> Vec<String> {
    let mut cols: Vec<String> = LOCI_FIXED.iter().map(|s| s.to_string()).collect();
    for a in ["a1", "a2"] {
        cols.extend(HAP_KEYS.iter().map(|k| format!("{a}_{k}")));
    }
    cols.extend(LOCI_TAIL.iter().map(|s| s.to_string()));
    cols
}

pub fn loci_row(res: &LocusResult) -> Vec<String> {
    let loc = res.locus;
    let mut vals: Vec<String> = vec![
        res.sample.clone(),
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
        PyVal::num(loc.ref_copies).fmt(1),
        res.ploidy.to_string(),
        res.method.clone(),
        res.n_reads_total.to_string(),
        res.n_spanning.to_string(),
        res.n_split.to_string(),
        res.n_partial.to_string(),
        res.n_no_run.to_string(),
    ];
    for i in 0..2 {
        let h = res.haplotypes.get(i);
        for k in HAP_KEYS {
            let v = h.map(|h| h.get(k)).unwrap_or(PyVal::None);
            vals.push(v.fmt(if matches!(k, "median" | "p5" | "p95") { 1 } else { 3 }));
        }
    }
    let minor: Vec<String> = res
        .minor_modes
        .iter()
        .map(|m| format!("{}:{}", fixed(m.median.f(), 0), m.support))
        .collect();
    vals.push(if minor.is_empty() { ".".into() } else { minor.join(";") });
    vals.push(PyVal::opt_i(res.lb_units).fmt(3));
    vals.push(res.lb_support.to_string());
    vals.push(fmt_opt_f(res.detect, 3));
    vals.push(res.evidence.clone());
    vals.push((res.negative_reliable as i64).to_string());
    vals.push(fmt_opt_f(res.tagged_fraction, 2));
    vals.push(res.meth_convention.clone());
    vals.push(if res.notes.is_empty() {
        ".".into()
    } else {
        res.notes.join("; ")
    });
    vals
}

// ----------------------------------------------------------------------------
// multi-locus driver
// ----------------------------------------------------------------------------

/// Call every locus, in order. With `threads > 1` loci are shared out over
/// worker threads, each with its own file handles; the order of the results
/// does not depend on which thread called which locus.
pub fn call_sample<'a>(
    bam_path: &Path,
    ref_path: &Path,
    loci: &[&'a Locus],
    sample: &str,
    opts: &Options,
    threads: usize,
) -> Result<Vec<LocusResult<'a>>, String> {
    let open = || -> Result<(Bam, Fasta), String> {
        let bam = Bam::open(bam_path).map_err(|e| format!("{}: {e}", bam_path.display()))?;
        let fasta = Fasta::open(ref_path).map_err(|e| format!("{}: {e}", ref_path.display()))?;
        Ok((bam, fasta))
    };
    if threads <= 1 || loci.len() <= 1 {
        let (mut bam, mut fasta) = open()?;
        return loci
            .iter()
            .map(|l| call_locus(&mut bam, Some(&mut fasta), l, sample, opts))
            .collect();
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let slots: std::sync::Mutex<Vec<Option<Result<LocusResult<'a>, String>>>> =
        std::sync::Mutex::new((0..loci.len()).map(|_| None).collect());
    std::thread::scope(|s| {
        for _ in 0..threads.min(loci.len()) {
            s.spawn(|| {
                let mut handles = match open() {
                    Ok(h) => Some(h),
                    Err(e) => {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                        if i < loci.len() {
                            slots.lock().unwrap()[i] = Some(Err(e));
                        }
                        None
                    }
                };
                loop {
                    let i = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                    if i >= loci.len() {
                        break;
                    }
                    let r = match handles.as_mut() {
                        Some((bam, fasta)) => call_locus(bam, Some(fasta), loci[i], sample, opts),
                        None => Err("could not open inputs".into()),
                    };
                    slots.lock().unwrap()[i] = Some(r);
                }
            });
        }
    });
    slots
        .into_inner()
        .unwrap()
        .into_iter()
        .map(|r| r.expect("every locus called"))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // tests/test_engine.py's test_modal_unit_prefers_an_exact_rotation_only_on_a_tie
    #[test]
    fn modal_unit_prefers_an_exact_rotation_only_on_a_tie() {
        const CEL: &str = "GGCCCCCCCCGTGCCGCCCACGGGTGACTCCGG";
        let ms = MotifSet::new(33, &[CEL.to_string()], &[], &[], &[]).unwrap();
        let modal = |units: &[&str]| {
            let got = modal_unit(units.iter().map(|u| u.as_bytes()), &ms);
            got.map(|u| String::from_utf8(u.to_vec()).unwrap())
        };
        let noisy = format!("{}A{}", &CEL[..5], &CEL[6..]);
        let rotated = format!("{}{}", &CEL[1..], &CEL[..1]);
        assert_eq!(modal(&[&noisy, &rotated, CEL]), Some(rotated.clone()));
        let variant = format!("{}T{}", &CEL[..10], &CEL[11..]);
        assert_eq!(modal(&[&variant, CEL, &variant, &variant]), Some(variant.clone()));
        let rot3 = format!("{}{}", &CEL[3..], &CEL[..3]);
        assert_eq!(modal(&[CEL, &rot3]), Some(CEL.to_string()));
        assert_eq!(modal(&[]), None);
    }
}
