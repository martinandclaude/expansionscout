//! Replays allele read sets written by `tests/parity/consensus_vectors.py` in the
//! development repository and compares them with the Python's answers:
//!
//!     python tests/parity/consensus_vectors.py /tmp/consensus.json 3000
//!     CONSENSUS_VECTORS=/tmp/consensus.json cargo test --release --test consensus_vectors -- --nocapture
//!
//! Skipped unless `CONSENSUS_VECTORS` names the file. Each read is decomposed
//! as the engine does it, the reads with a decomposition are the allele's
//! members, and the consensus and the consensus's own decomposition must be
//! the Python's exactly. A read's decomposition is compared too, so a
//! difference there is told apart from one in the consensus.

use expansionscout::consensus::{allele_consensus, ConsensusRead};
use expansionscout::decompose::{Decomposition, MotifSet};
use serde_json::Value;

fn strings(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap().to_string())
        .collect()
}

fn usize_of(v: &Value) -> usize {
    v.as_u64().unwrap() as usize
}

/// (start, end, labels) of a decomposition, as the Python writes it.
fn summary(d: &Option<Decomposition>) -> Value {
    match d {
        Some(d) => serde_json::json!([d.start, d.end, d.labels_str()]),
        None => Value::Null,
    }
}

/// Every part of the Rust's answer that differs from the Python's, and
/// whether a consensus was formed.
fn differences(c: &Value) -> (Vec<String>, bool) {
    let ms = MotifSet::new(
        usize_of(&c["unit_len"]),
        &strings(&c["canonical"]),
        &strings(&c["pathogenic"]),
        &strings(&c["benign"]),
        &strings(&c["interruption"]),
    )
    .unwrap();
    let e = &c["expect"];
    let mut out = Vec::new();
    let reads = c["reads"].as_array().unwrap();
    let segs: Vec<(&[u8], Vec<u8>)> = reads
        .iter()
        .map(|rd| {
            let s = rd["seg"].as_str().unwrap().as_bytes();
            (s, s.to_ascii_uppercase())
        })
        .collect();
    let decs: Vec<Option<Decomposition>> = reads
        .iter()
        .zip(&segs)
        .map(|(rd, (_, useq))| {
            let w = &rd["window"];
            ms.decompose(
                useq,
                Some((usize_of(&w[0]), usize_of(&w[1]))),
                Some(usize_of(&rd["max_bp"])),
            )
        })
        .collect();
    for (i, d) in decs.iter().enumerate() {
        let got = summary(d);
        if got != e["decs"][i] {
            out.push(format!("read {i} decomposition {got} vs {}", e["decs"][i]));
        }
    }
    let members: Vec<ConsensusRead> = decs
        .iter()
        .zip(reads.iter().zip(&segs))
        .filter_map(|(d, (rd, (seg, useq)))| {
            d.as_ref().map(|d| ConsensusRead {
                seq: &seg[d.start..d.end],
                run_u: &useq[d.start..d.end],
                dec: d,
                units: rd["units"].as_i64().unwrap(),
            })
        })
        .collect();
    let cons = allele_consensus(&members, &ms, usize_of(&c["rep_index"]));
    let got = cons.as_ref().map(|s| String::from_utf8_lossy(s).into_owned());
    if got.as_deref() != e["consensus"].as_str() {
        out.push(format!("consensus {got:?} vs {}", e["consensus"]));
    }
    if let Some(cons) = &cons {
        let cons_u = cons.to_ascii_uppercase();
        let n = cons_u.len();
        let labels = ms.decompose(&cons_u, Some((0, n)), Some(n)).map(|d| d.labels_str());
        if labels.as_deref() != e["cons_labels"].as_str() {
            out.push(format!("consensus structure {labels:?} vs {}", e["cons_labels"]));
        }
    }
    (out, cons.is_some())
}

#[test]
fn consensus_vectors_match_python() {
    let Ok(path) = std::env::var("CONSENSUS_VECTORS") else {
        eprintln!("CONSENSUS_VECTORS not set; skipping");
        return;
    };
    let cases: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut mismatches = 0;
    let mut formed = 0;
    for (i, c) in cases.iter().enumerate() {
        let (diffs, made) = differences(c);
        formed += usize::from(made);
        if !diffs.is_empty() {
            if mismatches < 5 {
                eprintln!("case {i} ({}): {}", c["locus"], diffs.join("; "));
            }
            mismatches += 1;
        }
    }
    eprintln!(
        "{} cases, {} mismatches; {} with a consensus",
        cases.len(),
        mismatches,
        formed
    );
    assert!(mismatches == 0, "{mismatches} cases differ from the Python");
}
