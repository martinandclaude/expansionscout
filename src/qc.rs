//! Checks a run can make on itself, with no truth set.
//! Port of `expansionscout/qc.py`, which documents each check and its basis.

use crate::npstat::fsum;
use crate::pyfmt::{general, wrap};

pub const MAX_MODCALL_FAIL_FRAC: f64 = 0.50;
pub const MIN_TAGGED_FRACTION: f64 = 0.50;
pub const ZYGOSITY_MARGIN_UNITS: f64 = 1.0;
pub const FMR1_LOCUS: &str = "FXS_FMR1";
pub const FMR1_UNMETHYLATED_MAX: f64 = 0.20;

/// One row of a `.loci.tsv`, as `csv.DictReader` gives it.
pub type TsvRow = Vec<(String, String)>;

pub fn get<'a>(row: &'a TsvRow, key: &str) -> Option<&'a str> {
    row.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
}

/// Python's `float(str)`, or None.
pub fn py_float(s: &str) -> Option<f64> {
    let t = s.trim();
    if t.is_empty() {
        return None;
    }
    let bad_us = t.starts_with('_') || t.ends_with('_') || t.contains("__");
    let cleaned: String = t.chars().filter(|&c| c != '_').collect();
    if t.contains('_') && bad_us {
        return None;
    }
    match cleaned.to_ascii_lowercase().as_str() {
        "nan" | "+nan" | "-nan" => Some(f64::NAN),
        "inf" | "+inf" | "infinity" | "+infinity" => Some(f64::INFINITY),
        "-inf" | "-infinity" => Some(f64::NEG_INFINITY),
        other => {
            if other
                .chars()
                .all(|c| c.is_ascii_digit() || matches!(c, '.' | 'e' | '+' | '-'))
            {
                other.parse().ok()
            } else {
                None
            }
        }
    }
}

/// `_f(row, key)`.
pub fn f(row: &TsvRow, key: &str) -> Option<f64> {
    match get(row, key) {
        None | Some("") | Some(".") => None,
        Some(v) => py_float(v),
    }
}

/// `_i(row, key)`.
fn i(row: &TsvRow, key: &str) -> Option<i64> {
    f(row, key).map(|v| v.trunc() as i64)
}

fn s<'a>(row: &'a TsvRow, key: &str) -> &'a str {
    get(row, key).unwrap_or("")
}

#[derive(Clone, Debug)]
pub struct Finding {
    pub check: String,
    pub status: String,
    pub basis: String,
    pub message: String,
    pub detail: Vec<String>,
}

impl Finding {
    fn new(check: &str, status: &str, basis: &str, message: impl Into<String>, detail: Vec<String>) -> Finding {
        Finding {
            check: check.into(),
            status: status.into(),
            basis: basis.into(),
            message: message.into(),
            detail,
        }
    }
    pub fn ok(&self) -> bool {
        self.status == "pass" || self.status == "skip"
    }
}

fn pct1(v: f64) -> String {
    format!("{:.1}", 100.0 * v)
}

fn check_modifications_present(rows: &[TsvRow]) -> Finding {
    let n = rows
        .iter()
        .filter(|r| f(r, "a1_meth_tract").is_some() || f(r, "a1_meth_up").is_some())
        .count();
    if rows.is_empty() {
        return Finding::new(
            "modifications_present",
            "skip",
            "mechanical",
            "no loci in input",
            vec![],
        );
    }
    if n == 0 {
        return Finding::new(
            "modifications_present",
            "fail",
            "mechanical",
            format!("no 5mC calls at any of {} loci", rows.len()),
            vec![
                "the reads carry no MM/ML tags over these windows, or were basecalled without a modification model"
                    .into(),
                "every methylation column in this run is empty, and the genotype columns give no sign of it".into(),
            ],
        );
    }
    Finding::new(
        "modifications_present",
        "pass",
        "mechanical",
        format!("5mC calls present at {n} of {} loci", rows.len()),
        vec![],
    )
}

fn check_modcall_confidence(rows: &[TsvRow]) -> Finding {
    let vals: Vec<f64> = rows
        .iter()
        .flat_map(|r| [f(r, "a1_meth_fail_frac"), f(r, "a2_meth_fail_frac")])
        .flatten()
        .collect();
    if vals.is_empty() {
        return Finding::new(
            "modcall_confidence",
            "skip",
            "heuristic",
            "no modification calls to score",
            vec![],
        );
    }
    let worst = vals.iter().copied().fold(vals[0], |a, b| if b > a { b } else { a });
    let mean = fsum(vals.iter().copied()) / vals.len() as f64;
    let msg = format!(
        "mean {} % of 5mC calls discarded below the confidence threshold, worst locus {} %",
        pct1(mean),
        pct1(worst)
    );
    if mean > MAX_MODCALL_FAIL_FRAC {
        return Finding::new("modcall_confidence", "warn", "heuristic", msg, vec![
            format!("above our {:.0} % cut-point, which is unvalidated", 100.0 * MAX_MODCALL_FAIL_FRAC),
            "check that --meth-threshold matches the basecalling model; `modkit sample-probs` reports the value modkit would have chosen".into(),
        ]);
    }
    Finding::new("modcall_confidence", "pass", "heuristic", msg, vec![])
}

fn check_x_inactivation(rows: &[TsvRow]) -> Finding {
    let Some(r) = rows.iter().find(|r| get(r, "locus") == Some(FMR1_LOCUS)) else {
        return Finding::new(
            "x_inactivation",
            "skip",
            "biology",
            format!("{FMR1_LOCUS} not in this run"),
            vec![],
        );
    };
    let (a1, a2) = (f(r, "a1_meth_up"), f(r, "a2_meth_up"));
    let Some(a1) = a1 else {
        return Finding::new(
            "x_inactivation",
            "skip",
            "biology",
            format!("no promoter 5mC at {FMR1_LOCUS}"),
            vec![],
        );
    };
    let ploidy = i(r, "ploidy").filter(|&p| p != 0).unwrap_or(2);
    let cls1 = s(r, "a1_cls").to_lowercase();
    if ploidy == 1 {
        if cls1.contains("path") {
            return Finding::new(
                "x_inactivation",
                "skip",
                "biology",
                "single X carries an expanded allele; promoter methylation is expected and is not a control",
                vec![],
            );
        }
        let msg = format!("single X, normal-length allele, promoter 5mC {} %", pct1(a1));
        if a1 > FMR1_UNMETHYLATED_MAX {
            return Finding::new("x_inactivation", "warn", "biology", msg, vec![
                "an unexpanded allele on an active X should be unmethylated (ACMG fragile X standards, Maddalena et al. Genet Med 2001;3:200-205)".into(),
                format!("our cut-point for 'unmethylated' is {:.0} % and is a heuristic, so treat this as a prompt to look, not a result",
                        100.0 * FMR1_UNMETHYLATED_MAX),
            ]);
        }
        return Finding::new("x_inactivation", "pass", "biology", msg, vec![]);
    }
    let Some(a2) = a2 else {
        return Finding::new(
            "x_inactivation",
            "skip",
            "biology",
            "two X chromosomes but only one allele measured",
            vec![],
        );
    };
    let spread = (a1 - a2).abs();
    let msg = format!(
        "two X chromosomes, promoter 5mC {} % and {} % (difference {} points)",
        pct1(a1),
        pct1(a2),
        pct1(spread)
    );
    if a1 == a2 {
        return Finding::new("x_inactivation", "warn", "consistency", msg, vec![
            "both alleles report an identical value, which happens when reads could not be assigned to a haplotype".into(),
            "X-inactivation cannot be seen without phasing; the number shown is the average of an active and an inactive X and describes neither".into(),
        ]);
    }
    Finding::new("x_inactivation", "pass", "biology", msg, vec![])
}

fn check_per_allele_attribution(rows: &[TsvRow]) -> Finding {
    let tagged: Vec<f64> = rows.iter().filter_map(|r| f(r, "tagged_fraction")).collect();
    if tagged.is_empty() {
        return Finding::new(
            "per_allele_attribution",
            "skip",
            "consistency",
            "no tagged_fraction reported",
            vec![],
        );
    }
    let mean = fsum(tagged.iter().copied()) / tagged.len() as f64;
    let gmm = rows.iter().filter(|r| s(r, "method").starts_with("GMM")).count();
    let msg = format!(
        "mean tagged fraction {} %, {gmm} of {} loci called without tags",
        pct1(mean),
        rows.len()
    );
    if mean >= MIN_TAGGED_FRACTION && gmm == rows.len() && !rows.is_empty() {
        return Finding::new(
            "per_allele_attribution",
            "warn",
            "consistency",
            msg,
            vec![
                "tags are present on most reads but were used at no locus, which is what --ignore-hp looks like".into(),
                "per-allele methylation in this run is a pooled average reported twice, not two measurements".into(),
            ],
        );
    }
    if mean < MIN_TAGGED_FRACTION {
        return Finding::new(
            "per_allele_attribution",
            "warn",
            "heuristic",
            msg,
            vec![
                format!(
                    "below our {:.0} % cut-point, which is unvalidated",
                    100.0 * MIN_TAGGED_FRACTION
                ),
                "allele assignment falls back on read length alone, so alleles close in size may be pooled".into(),
            ],
        );
    }
    Finding::new("per_allele_attribution", "pass", "consistency", msg, vec![])
}

fn check_zygosity_margin(rows: &[TsvRow]) -> Finding {
    let (mut close, mut consequential): (Vec<String>, Vec<String>) = (vec![], vec![]);
    for r in rows {
        if i(r, "ploidy").filter(|&p| p != 0).unwrap_or(2) < 2 {
            continue;
        }
        let (Some(a1), a2) = (f(r, "a1_median"), f(r, "a2_median")) else {
            continue;
        };
        let (c1, c2) = (s(r, "a1_cls"), s(r, "a2_cls"));
        let locus = get(r, "locus").map(String::from).unwrap_or_else(|| "None".into());
        match a2 {
            Some(a2) if a1 != a2 && (a2 - a1).abs() <= ZYGOSITY_MARGIN_UNITS => {
                let entry = format!(
                    "{locus}: {}/{}, alleles {} unit apart",
                    general(a1),
                    general(a2),
                    general((a2 - a1).abs())
                );
                if c1 != c2 {
                    consequential.push(format!("{entry}, and they are called {c1} and {c2}"));
                } else {
                    close.push(entry);
                }
            }
            Some(a2) if a1 == a2 => {
                if let (Some(p5), Some(p95)) = (f(r, "a1_p5"), f(r, "a1_p95")) {
                    if p95 - p5 > 2.0 * ZYGOSITY_MARGIN_UNITS {
                        consequential.push(format!(
                            "{locus}: {}/{} called homozygous, but the reads span {}-{}, wide enough to hold a second allele",
                            general(a1), general(a2), general(p5), general(p95)));
                    }
                }
            }
            _ => {}
        }
    }
    if rows.is_empty() {
        return Finding::new("zygosity_margin", "skip", "heuristic", "no loci in input", vec![]);
    }
    let msg = format!("{} calls where the homozygote/heterozygote decision could change the answer, {} more that are close but land in the same class",
                      consequential.len(), close.len());
    if !consequential.is_empty() {
        let mut detail: Vec<String> = consequential.iter().take(20).cloned().collect();
        if consequential.len() > 20 {
            detail.push("...".into());
        }
        detail.push(format!("our margin is {} unit and is a heuristic; these calls may well be right, but a second assembly or one more read could move them",
                            general(ZYGOSITY_MARGIN_UNITS)));
        return Finding::new("zygosity_margin", "warn", "heuristic", msg, detail);
    }
    Finding::new("zygosity_margin", "pass", "heuristic", msg, vec![])
}

fn check_negative_evidence(rows: &[TsvRow]) -> Finding {
    let diploid: Vec<&TsvRow> = rows
        .iter()
        .filter(|r| !matches!(get(r, "negative_reliable"), None | Some("") | Some(".")))
        .collect();
    if diploid.is_empty() {
        return Finding::new(
            "negative_evidence",
            "skip",
            "consistency",
            "no negative_reliable column",
            vec![],
        );
    }
    let unreliable: Vec<String> = diploid
        .iter()
        .filter(|r| i(r, "negative_reliable") == Some(0))
        .map(|r| get(r, "locus").map(String::from).unwrap_or_else(|| "None".into()))
        .collect();
    let msg = format!(
        "{} of {} loci support a reliable negative",
        diploid.len() - unreliable.len(),
        diploid.len()
    );
    if !unreliable.is_empty() {
        let mut detail = vec!["a normal call at these loci is not evidence of a normal result; read length or depth could not have shown an expansion".to_string()];
        detail.extend(unreliable.iter().take(20).cloned());
        if unreliable.len() > 20 {
            detail.push("...".into());
        }
        return Finding::new("negative_evidence", "warn", "consistency", msg, detail);
    }
    Finding::new("negative_evidence", "pass", "consistency", msg, vec![])
}

pub fn run_checks(rows: &[TsvRow]) -> Vec<Finding> {
    vec![
        check_modifications_present(rows),
        check_modcall_confidence(rows),
        check_x_inactivation(rows),
        check_per_allele_attribution(rows),
        check_zygosity_margin(rows),
        check_negative_evidence(rows),
    ]
}

/// A tab-separated table with a header line, as `csv.DictReader` reads it.
pub fn read_tsv(text: &str) -> Vec<TsvRow> {
    let mut lines = text.split('\n').filter(|l| !l.is_empty());
    let Some(head) = lines.next() else { return vec![] };
    let head: Vec<&str> = head.trim_end_matches('\r').split('\t').collect();
    lines
        .map(|l| {
            let vals: Vec<&str> = l.trim_end_matches('\r').split('\t').collect();
            head.iter()
                .enumerate()
                .map(|(i, k)| (k.to_string(), vals.get(i).unwrap_or(&"").to_string()))
                .collect()
        })
        .collect()
}

/// Human-readable, one block per check.
pub fn format_report(findings: &[Finding], sample: Option<&str>) -> String {
    let width = 78;
    let mut out: Vec<String> = Vec::new();
    if let Some(s) = sample.filter(|s| !s.is_empty()) {
        out.push(format!("QC: {s}"));
        out.push(String::new());
    }
    for x in findings {
        let mark = match x.status.as_str() {
            "pass" => "ok  ",
            "warn" => "WARN",
            "fail" => "FAIL",
            "skip" => "--  ",
            _ => "?   ",
        };
        out.push(format!("{mark} {}  [{}]", x.check, x.basis));
        for line in wrap(&x.message, width - 6) {
            out.push(format!("      {line}"));
        }
        for d in &x.detail {
            for (k, line) in wrap(d, width - 10).iter().enumerate() {
                out.push(format!("        {}{line}", if k == 0 { "- " } else { "  " }));
            }
        }
        out.push(String::new());
    }
    let count = |st: &str| findings.iter().filter(|x| x.status == st).count();
    let (n_fail, n_warn, n_skip) = (count("fail"), count("warn"), count("skip"));
    out.push(format!(
        "{} checks: {n_fail} failed, {n_warn} warned, {n_skip} could not run",
        findings.len()
    ));
    if n_skip > 0 {
        out.push("A check that could not run has not passed.".into());
    }
    out.join("\n")
}
