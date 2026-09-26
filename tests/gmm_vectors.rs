//! Replays allele-assignment cases written by `tests/parity/gmm_vectors.py` in the
//! development repository and compares them with the Python's answers:
//!
//!     python tests/parity/gmm_vectors.py /tmp/gmm.json 6000
//!     GMM_VECTORS=/tmp/gmm.json cargo test --release --test gmm_vectors -- --nocapture
//!
//! Skipped unless `GMM_VECTORS` names the file. The assignment itself --
//! method, notes, which read went to which allele, medians and spread -- must
//! match exactly. The raw mixture fits -- means, sds, weights, log-likelihood
//! and BIC -- must agree to within `FIT_TOLERANCE` relative, not exactly,
//! because they are where numpy's CPU-dependent `exp`/`log` and the Rust's
//! `libm` meet; the worst difference is printed.

use expansionscout::alleles::{cluster_counts, fit_gmm, spread_stats};
use serde_json::Value;

/// Relative agreement required of the raw fits: four orders of magnitude
/// above the largest difference seen, and far below anything that could move
/// an allele assignment.
const FIT_TOLERANCE: f64 = 1e-6;

fn f(v: &Value) -> f64 {
    v.as_f64().unwrap()
}

fn floats(v: &Value) -> Vec<f64> {
    v.as_array().unwrap().iter().map(f).collect()
}

fn ints(v: &Value) -> Vec<i64> {
    v.as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect()
}

fn strings(v: &Value) -> Vec<String> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_str().unwrap().to_string())
        .collect()
}

/// Every field of the Rust's assignment that differs from the Python's.
fn differences(c: &Value) -> Vec<String> {
    let values = floats(&c["values"]);
    let a = cluster_counts(
        &values,
        &ints(&c["hps"]),
        c["ploidy"].as_u64().unwrap() as usize,
        c["mosaic"].as_bool().unwrap(),
        c["min_support"].as_u64().unwrap() as usize,
        f(&c["min_frac"]),
        c["use_hp"].as_bool().unwrap(),
    );
    let e = &c["expect"];
    let mut out = Vec::new();
    let mut check = |what: String, same: bool| {
        if !same {
            out.push(what);
        }
    };
    check(
        format!("method {}", a.method),
        a.method == e["method"].as_str().unwrap(),
    );
    check(format!("notes {:?}", a.notes), a.notes == strings(&e["notes"]));
    check("labels".into(), a.labels == ints(&e["labels"]));
    let expected = e["clusters"].as_array().unwrap();
    check(
        format!("{} clusters", a.clusters.len()),
        a.clusters.len() == expected.len(),
    );
    for (k, (cl, ec)) in a.clusters.iter().zip(expected).enumerate() {
        let members: Vec<usize> = ints(&ec["members"]).iter().map(|&m| m as usize).collect();
        check(format!("cluster {k} members"), cl.members == members);
        check(
            format!("cluster {k} median {}", cl.median),
            cl.median == f(&ec["median"]),
        );
        check(format!("cluster {k} minor"), cl.minor == ec["minor"].as_bool().unwrap());
        check(format!("cluster {k} hp"), cl.hp == ec["hp"].as_i64());
        if let Some(s) = spread_stats(&cl.members.iter().map(|&m| values[m]).collect::<Vec<_>>()) {
            let es = &ec["spread"];
            for (name, got) in [
                ("median", s.median),
                ("mad", s.mad),
                ("p5", s.p5),
                ("p95", s.p95),
                ("tail_up_frac", s.tail_up_frac),
                ("tail_down_frac", s.tail_down_frac),
            ] {
                check(format!("cluster {k} {name} {got}"), got == f(&es[name]));
            }
        }
    }
    out
}

#[test]
fn gmm_vectors_match_python() {
    let Ok(path) = std::env::var("GMM_VECTORS") else {
        eprintln!("GMM_VECTORS not set; skipping");
        return;
    };
    let cases: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut mismatches = Vec::new();
    let mut worst_rel = 0.0f64;
    let mut bad_fits = 0;
    let mut n_fits = 0;
    let rel = |a: f64, b: f64| {
        if a == b {
            0.0
        } else {
            (a - b).abs() / a.abs().max(b.abs()).max(1e-300)
        }
    };
    for (i, c) in cases.iter().enumerate() {
        let diffs = differences(c);
        if !diffs.is_empty() {
            if mismatches.len() < 5 {
                eprintln!("case {i}: {}", diffs.join("; "));
            }
            mismatches.push(i);
        }
        let values = floats(&c["values"]);
        for fit in c.get("fits").and_then(Value::as_array).into_iter().flatten() {
            let g = fit_gmm(&values, fit["k"].as_u64().unwrap() as usize);
            n_fits += 1;
            let mut diffs = vec![rel(g.ll, f(&fit["ll"])), rel(g.bic, f(&fit["bic"]))];
            for (got, want) in [(&g.mu, &fit["mu"]), (&g.sd, &fit["sd"]), (&g.w, &fit["w"])] {
                let want = floats(want);
                assert_eq!(
                    got.len(),
                    want.len(),
                    "case {i}: a {}-component fit differs in size",
                    fit["k"]
                );
                diffs.extend(got.iter().zip(want).map(|(x, y)| rel(*x, y)));
            }
            for d in diffs {
                if d.is_nan() || d > FIT_TOLERANCE {
                    bad_fits += 1;
                }
                if d.is_nan() || d > worst_rel {
                    worst_rel = d;
                }
            }
        }
    }
    eprintln!(
        "{} cases, {} assignment mismatches; {} raw fits, worst relative difference {:e}",
        cases.len(),
        mismatches.len(),
        n_fits,
        worst_rel
    );
    assert!(
        mismatches.is_empty(),
        "assignment differs from the Python in {} cases",
        mismatches.len()
    );
    assert!(
        bad_fits == 0,
        "{bad_fits} raw fit values differ from the Python by more than {FIT_TOLERANCE:e} relative, or are NaN"
    );
}
