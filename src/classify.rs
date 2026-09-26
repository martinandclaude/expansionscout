//! Size classification, detectability and expansion evidence.
//! Port of `expansionscout/classify.py`; the reasoning lives there.

use crate::catalog::Locus;
use crate::npstat::py_round;

pub const BENIGN: &str = "benign";
pub const INTERMEDIATE: &str = "intermediate";
pub const PATHOGENIC: &str = "pathogenic";
pub const UNCERTAIN: &str = "uncertain";
pub const NO_CALL: &str = "no_call";

pub const EVIDENCE_NONE: &str = "NONE";
pub const EVIDENCE_PARTIALS: &str = "PARTIALS_PRESENT";
pub const EVIDENCE_EXPANSION_LB: &str = "EXPANSION_LB";

pub fn classify_units(locus: &Locus, units: Option<f64>) -> String {
    let Some(units) = units else { return NO_CALL.into() };
    if !locus.bands.is_empty() {
        for b in &locus.bands {
            if units <= b.upper {
                return b.label.clone();
            }
        }
        return NO_CALL.into();
    }
    let f = |v: Option<crate::pyfmt::Num>| v.map(|n| n.f());
    let (pm, pmax) = (f(locus.pathogenic_min), f(locus.pathogenic_max));
    let (bmin, bmax) = (f(locus.benign_min), f(locus.benign_max));
    let (imin, imax) = (f(locus.intermediate_min), f(locus.intermediate_max));
    let in_benign = |n: f64| bmax.is_some_and(|b| n <= b) && bmin.is_none_or(|b| n >= b);

    // Contraction loci: the pathogenic range lies below the benign one.
    if let (Some(_), Some(pmax), Some(bmin)) = (pm, pmax, bmin) {
        if pmax < bmin {
            if units <= pmax {
                return PATHOGENIC.into();
            }
            if in_benign(units) {
                return BENIGN.into();
            }
            return UNCERTAIN.into();
        }
    }
    if in_benign(units) {
        return BENIGN.into();
    }
    if let (Some(a), Some(b)) = (imin, imax) {
        if a <= units && units <= b {
            return INTERMEDIATE.into();
        }
    }
    if pm.is_some_and(|p| units >= p) {
        return PATHOGENIC.into();
    }
    if bmax.is_some_and(|b| units <= b) {
        return BENIGN.into();
    }
    UNCERTAIN.into()
}

/// Allele size under the reporting convention in force. Returns `None` for
/// Python's integer 0, which `max(0, x)` yields when nothing is left; callers
/// print that as "0" rather than "0.0".
pub fn size_units(locus: &Locus, n_total: f64, n_interruption: i64, clinical: bool) -> SizeVal {
    if !clinical || locus.interruptions_in_size {
        return SizeVal::Float(n_total);
    }
    let v = n_total - n_interruption as f64;
    if v > 0.0 {
        SizeVal::Float(v)
    } else {
        SizeVal::Zero
    }
}

/// A reported size is a float, except that Python's `max(0, x)` returns the
/// integer 0 when x is not positive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SizeVal {
    Float(f64),
    Zero,
}

impl SizeVal {
    pub fn f(self) -> f64 {
        match self {
            SizeVal::Float(v) => v,
            SizeVal::Zero => 0.0,
        }
    }
}

/// Units that count toward the size class.
pub fn class_units(locus: &Locus, n_total: f64, n_pathogenic: i64, n_benign: i64) -> f64 {
    if locus.composition() {
        return n_pathogenic as f64;
    }
    let v = n_total - n_benign as f64;
    if v > 0.0 {
        v
    } else {
        0.0
    }
}

/// Are the reads behind one fitted component too spread out to be one allele?
pub fn over_dispersed(median: Option<f64>, p95: Option<f64>, tail_up_frac: Option<f64>, support: usize) -> bool {
    let (slack, floor, min_tail_reads) = (1.5, 3.0, 2.0);
    let (Some(median), Some(p95)) = (median, p95) else {
        return false;
    };
    let a = slack * median;
    let b = median + floor;
    if p95 <= if b > a { b } else { a } {
        return false;
    }
    if let Some(t) = tail_up_frac {
        if support > 0 {
            return py_round(t * support as f64) >= min_tail_reads;
        }
    }
    true
}

pub fn detectability(locus: &Locus, read_lengths: &[i64], flank: i64, anchor: i64) -> Option<f64> {
    if read_lengths.is_empty() {
        return None;
    }
    let pm = locus.pathogenic_min?;
    let needed = pm.f() * locus.unit_len as f64 + 2.0 * (flank + anchor) as f64;
    let n = read_lengths.iter().filter(|&&l| l as f64 >= needed).count();
    Some(n as f64 / read_lengths.len() as f64)
}

/// Corroborated lower bound: the second-largest per-read count.
pub fn lower_bound(partial_counts: &[i64]) -> (Option<i64>, usize) {
    let min_reads = 2;
    let mut vals: Vec<i64> = partial_counts.iter().copied().filter(|&v| v > 0).collect();
    vals.sort_by(|a, b| b.cmp(a));
    if vals.len() < min_reads {
        return (None, 0);
    }
    let v = vals[min_reads - 1];
    (Some(v), vals.iter().filter(|&&w| w >= v).count())
}

pub fn evidence_state(
    locus: &Locus,
    allele_medians: &[f64],
    lb_units: Option<i64>,
    lb_support: usize,
    n_partial: usize,
) -> &'static str {
    if let Some(lb) = lb_units {
        if lb_support >= 2 {
            let lbf = lb as f64;
            let explained = allele_medians.iter().any(|&m| m >= 0.9 * lbf);
            let mx = allele_medians
                .iter()
                .copied()
                .fold(None, |acc: Option<f64>, m| match acc {
                    Some(a) if a >= m => Some(a),
                    _ => Some(m),
                })
                .unwrap_or(0.0);
            if !explained && (locus.pathogenic_min.is_none_or(|p| lbf >= p.f()) || lbf > mx * 1.5) {
                return EVIDENCE_EXPANSION_LB;
            }
        }
    }
    if n_partial >= 2 {
        return EVIDENCE_PARTIALS;
    }
    EVIDENCE_NONE
}

pub fn negative_reliable(
    locus: &Locus,
    detect: Option<f64>,
    evidence: &str,
    n_spanning: i64,
    min_support: i64,
    dispersed: bool,
) -> bool {
    let min_detect = 0.8;
    if evidence == EVIDENCE_EXPANSION_LB || dispersed {
        return false;
    }
    if n_spanning < min_support {
        return false;
    }
    if locus.regime == "short" {
        return true;
    }
    detect.is_some_and(|d| d >= min_detect)
}
