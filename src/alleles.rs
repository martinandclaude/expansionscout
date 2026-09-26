//! Allele assignment from per-read unit counts.
//!
//! Port of `expansionscout/alleles.py`, which explains the model. The EM fit
//! reproduces numpy's order of operations exactly, including which
//! reductions are pairwise and which are sequential, so that the only source
//! of difference from the Python is `exp` and `log` themselves (see
//! `npstat`).

use crate::npstat::{argmin, exp, ln, median, pairwise_sum, percentile, quantiles, var};

pub const QUANT_VAR: f64 = 1.0 / 12.0;
pub const REL_SD: f64 = 0.02;
pub const SEP_K: f64 = 2.0;
pub const SEP_MIN: f64 = 0.9;

#[derive(Clone, Debug)]
pub struct Cluster {
    pub index: i64,
    pub members: Vec<usize>,
    pub median: f64,
    pub hp: Option<i64>,
    pub minor: bool,
}

#[derive(Clone, Debug)]
pub struct Assignment {
    pub clusters: Vec<Cluster>,
    pub labels: Vec<i64>,
    pub method: String,
    pub notes: Vec<String>,
}

/// Python's `max(a, b)`: `a` unless `b` is strictly greater.
fn pmax(a: f64, b: f64) -> f64 {
    if b > a {
        b
    } else {
        a
    }
}

fn tol(m: f64) -> f64 {
    pmax(2.0, 0.10 * m)
}

fn var_floor(mu: f64) -> f64 {
    let s = REL_SD * pmax(mu, 1.0);
    pmax(QUANT_VAR, s * s)
}

fn resolve_tol(m: f64) -> f64 {
    pmax(SEP_MIN, SEP_K * var_floor(m).sqrt())
}

pub struct Gmm {
    pub mu: Vec<f64>,
    pub sd: Vec<f64>,
    pub w: Vec<f64>,
    pub ll: f64,
    pub bic: f64,
}

/// A column sum of an (n, k) C-ordered array, as numpy's `sum(axis=0)` does
/// it: pairwise when k == 1, because the axis is then contiguous, and a
/// row-by-row running sum otherwise.
fn colsum(m: &[Vec<f64>], j: usize, k: usize) -> f64 {
    if k == 1 {
        let col: Vec<f64> = m.iter().map(|row| row[j]).collect();
        pairwise_sum(&col)
    } else {
        let mut s = 0.0;
        for row in m {
            s += row[j];
        }
        s
    }
}

/// One-dimensional Gaussian mixture by EM.
pub fn fit_gmm(x: &[f64], k: usize) -> Gmm {
    let iters = 300;
    let tolerance = 1e-6;
    let n = x.len();
    let qs: Vec<f64> = (0..k).map(|i| (i as f64 + 0.5) / k as f64).collect();
    let mut mu = quantiles(x, &qs);
    let mut v = vec![pmax(var(x), QUANT_VAR); k];
    let mut w = vec![1.0 / k as f64; k];
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut ll_old = f64::NEG_INFINITY;
    let mut ll = f64::NEG_INFINITY;
    for _ in 0..iters {
        let a: Vec<f64> = (0..k).map(|j| ln(w[j]) - 0.5 * ln(two_pi * v[j])).collect();
        let lp: Vec<Vec<f64>> = x
            .iter()
            .map(|&xi| {
                (0..k)
                    .map(|j| {
                        let d = xi - mu[j];
                        a[j] - (d * d) / (2.0 * v[j])
                    })
                    .collect()
            })
            .collect();
        let mut lse = Vec::with_capacity(n);
        for row in &lp {
            let mx = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let mut s = 0.0;
            for &l in row {
                s += exp(l - mx);
            }
            lse.push(mx + ln(s));
        }
        let r: Vec<Vec<f64>> = lp
            .iter()
            .zip(&lse)
            .map(|(row, &l)| row.iter().map(|&p| exp(p - l)).collect())
            .collect();
        ll = pairwise_sum(&lse);
        let nk: Vec<f64> = (0..k).map(|j| colsum(&r, j, k) + 1e-9).collect();
        w = nk.iter().map(|&c| c / n as f64).collect();
        let rx: Vec<Vec<f64>> = r
            .iter()
            .zip(x)
            .map(|(row, &xi)| row.iter().map(|&p| p * xi).collect())
            .collect();
        mu = (0..k).map(|j| colsum(&rx, j, k) / nk[j]).collect();
        v = (0..k)
            .map(|j| {
                let terms: Vec<f64> = r
                    .iter()
                    .zip(x)
                    .map(|(row, &xi)| {
                        let d = xi - mu[j];
                        row[j] * (d * d)
                    })
                    .collect();
                pmax(pairwise_sum(&terms) / nk[j], var_floor(mu[j]))
            })
            .collect();
        if ll - ll_old < tolerance {
            break;
        }
        ll_old = ll;
    }
    let n_par = (3 * k - 1) as f64;
    let bic = -2.0 * ll + n_par * ln(n.max(2) as f64);
    Gmm {
        sd: v.iter().map(|x| x.sqrt()).collect(),
        mu,
        w,
        ll,
        bic,
    }
}

fn assign_nearest(x: &[f64], centres: &[f64]) -> Vec<usize> {
    x.iter()
        .map(|&xi| {
            let d: Vec<f64> = centres.iter().map(|&c| (xi - c).abs()).collect();
            argmin(&d)
        })
        .collect()
}

fn prune_by_support(x: &[f64], centres: Vec<f64>, min_n: usize) -> (Vec<f64>, Vec<usize>) {
    let mut centres = centres;
    let mut labels = assign_nearest(x, &centres);
    while centres.len() > 1 {
        let sizes: Vec<(usize, usize)> = (0..centres.len())
            .map(|k| (labels.iter().filter(|&&l| l == k).count(), k))
            .collect();
        let small = sizes.iter().filter(|s| s.0 < min_n).min();
        let Some(&(_, k)) = small else { break };
        centres.remove(k);
        labels = assign_nearest(x, &centres);
    }
    (centres, labels)
}

/// Fit k = 1..max_k, take the best BIC, then merge unresolvable components.
fn mixture_components(x: &[f64], max_k: usize, sep_k: f64) -> Vec<f64> {
    let mut best: Option<Gmm> = None;
    for k in 1..=max_k {
        if x.len() < 3 * k {
            break;
        }
        let cand = fit_gmm(x, k);
        if best.as_ref().is_none_or(|b| cand.bic < b.bic - 1e-9) {
            best = Some(cand);
        }
    }
    let Some(best) = best else { return vec![median(x)] };
    let mut order: Vec<usize> = (0..best.mu.len()).collect();
    order.sort_by(|&a, &b| best.mu[a].partial_cmp(&best.mu[b]).unwrap());
    let mu: Vec<f64> = order.iter().map(|&i| best.mu[i]).collect();
    let sd: Vec<f64> = order.iter().map(|&i| best.sd[i]).collect();
    let w: Vec<f64> = order.iter().map(|&i| best.w[i]).collect();
    let mut keep = vec![0usize];
    for j in 1..mu.len() {
        let i = *keep.last().unwrap();
        if (mu[j] - mu[i]).abs() < pmax(SEP_MIN, sep_k * pmax(sd[i], sd[j])) {
            if w[j] > w[i] {
                *keep.last_mut().unwrap() = j;
            }
        } else {
            keep.push(j);
        }
    }
    keep.iter().map(|&i| mu[i]).collect()
}

fn ceil_frac(min_frac: f64, n: usize) -> usize {
    (min_frac * n as f64).ceil() as usize
}

fn select(x: &[f64], idx: &[usize]) -> Vec<f64> {
    idx.iter().map(|&i| x[i]).collect()
}

fn hp_misses_an_allele(
    x: &[f64],
    hps: &[i64],
    good: &[i64],
    centres: &[(i64, f64)],
    min_support: usize,
    min_frac: f64,
    sep_k: f64,
) -> bool {
    let untagged: Vec<f64> = x
        .iter()
        .zip(hps)
        .filter(|(_, h)| !good.contains(h))
        .map(|(&v, _)| v)
        .collect();
    if untagged.is_empty() {
        return false;
    }
    let cs: Vec<f64> = good
        .iter()
        .map(|g| centres.iter().find(|c| c.0 == *g).unwrap().1)
        .collect();
    let n_orphan = untagged
        .iter()
        .filter(|&&u| {
            let a = (u - cs[0]).abs() - resolve_tol(cs[0]);
            let b = (u - cs[1]).abs() - resolve_tol(cs[1]);
            a.min(b) > 0.0
        })
        .count();
    let need = min_support.max(ceil_frac(min_frac, x.len()));
    if n_orphan < need {
        return false;
    }
    let (surviving, labels) = prune_by_support(x, mixture_components(x, 4, sep_k), need);
    let mut meds: Vec<f64> = (0..surviving.len())
        .filter_map(|k| {
            let m: Vec<usize> = (0..x.len()).filter(|&i| labels[i] == k).collect();
            if m.is_empty() {
                None
            } else {
                Some(median(&select(x, &m)))
            }
        })
        .collect();
    meds.sort_by(|a, b| a.partial_cmp(b).unwrap());
    meds.len() >= 2 && (meds[meds.len() - 1] - meds[0]) > resolve_tol(meds[meds.len() - 1])
}

fn members_of(labels: &[i64], k: i64) -> Vec<usize> {
    labels
        .iter()
        .enumerate()
        .filter(|(_, &l)| l == k)
        .map(|(i, _)| i)
        .collect()
}

/// Cluster per-read unit counts into alleles.
#[allow(clippy::too_many_arguments)]
pub fn cluster_counts(
    values: &[f64],
    hps: &[i64],
    ploidy: usize,
    mosaic: bool,
    min_support: usize,
    min_frac: f64,
    use_hp: bool,
) -> Assignment {
    let sep_k = SEP_K;
    let x = values;
    let n = x.len();
    if n == 0 {
        return Assignment {
            clusters: vec![],
            labels: vec![],
            method: "none".into(),
            notes: vec![],
        };
    }
    let mut notes = Vec::new();

    if use_hp {
        let mut tagged: Vec<i64> = hps.iter().copied().filter(|h| *h == 1 || *h == 2).collect();
        tagged.sort();
        tagged.dedup();
        let count = |h: i64| hps.iter().filter(|&&v| v == h).count();
        let good: Vec<i64> = tagged.iter().copied().filter(|&h| count(h) >= min_support).collect();
        if ploidy == 2 && good.len() == 2 {
            let centres: Vec<(i64, f64)> = good
                .iter()
                .map(|&h| {
                    let v: Vec<f64> = x.iter().zip(hps).filter(|(_, &g)| g == h).map(|(&v, _)| v).collect();
                    (h, median(&v))
                })
                .collect();
            let centre = |h: i64| centres.iter().find(|c| c.0 == h).unwrap().1;
            if hp_misses_an_allele(x, hps, &good, &centres, min_support, min_frac, sep_k) {
                notes.push(
                    "haplotype tags describe only one allele while the reads show two; \
                            used the mixture model instead"
                        .to_string(),
                );
            } else {
                let mut labels = vec![-1i64; n];
                for i in 0..n {
                    if let Some(k) = good.iter().position(|&g| g == hps[i]) {
                        labels[i] = k as i64;
                    } else {
                        // min() over `good` in order: the first nearest haplotype
                        let mut h = good[0];
                        for &g in &good[1..] {
                            if (x[i] - centre(g)).abs() < (x[i] - centre(h)).abs() {
                                h = g;
                            }
                        }
                        if (x[i] - centre(h)).abs() <= tol(centre(h)) {
                            labels[i] = good.iter().position(|&g| g == h).unwrap() as i64;
                        }
                    }
                }
                let clusters = good
                    .iter()
                    .enumerate()
                    .map(|(k, &h)| Cluster {
                        index: k as i64,
                        members: members_of(&labels, k as i64),
                        median: centre(h),
                        hp: Some(h),
                        minor: false,
                    })
                    .collect();
                return finalise(clusters, labels, "HP", notes);
            }
        }
        if ploidy == 1 && !good.is_empty() {
            if good.len() == 2 {
                notes.push("two HP labels on a haploid locus; merged".to_string());
            }
            let clusters = vec![Cluster {
                index: 0,
                members: (0..n).collect(),
                median: median(x),
                hp: None,
                minor: false,
            }];
            return finalise(clusters, vec![0; n], "HP-haploid", notes);
        }
        if !tagged.is_empty() && good.is_empty() {
            notes.push("HP tags present but below min_support; used the mixture model".to_string());
        }
    }

    if ploidy == 1 {
        let clusters = vec![Cluster {
            index: 0,
            members: (0..n).collect(),
            median: median(x),
            hp: None,
            minor: false,
        }];
        return finalise(clusters, vec![0; n], "GMM", notes);
    }

    let min_n = min_support.max(ceil_frac(min_frac, n));
    let (centres, labels) = prune_by_support(x, mixture_components(x, ploidy + 2, sep_k), min_n);
    let mut groups: Vec<(usize, Vec<usize>)> = (0..centres.len())
        .map(|k| (k, (0..n).filter(|&i| labels[i] == k).collect::<Vec<_>>()))
        .filter(|g| !g.1.is_empty())
        .collect();
    groups.sort_by(|a, b| {
        let ka = (-(a.1.len() as i64), median(&select(x, &a.1)));
        let kb = (-(b.1.len() as i64), median(&select(x, &b.1)));
        ka.0.cmp(&kb.0).then(ka.1.partial_cmp(&kb.1).unwrap())
    });
    let split = ploidy.min(groups.len());
    let (keep, extra) = groups.split_at(split);

    let mut clusters: Vec<Cluster> = Vec::new();
    let mut new_labels = vec![-1i64; n];
    for (k, (_, members)) in keep.iter().enumerate() {
        for &i in members {
            new_labels[i] = k as i64;
        }
        clusters.push(Cluster {
            index: k as i64,
            members: members.clone(),
            median: median(&select(x, members)),
            hp: None,
            minor: false,
        });
    }
    if !extra.is_empty() {
        if mosaic {
            for (off, (_, members)) in extra.iter().enumerate() {
                let k = (keep.len() + off) as i64;
                for &i in members {
                    new_labels[i] = k;
                }
                clusters.push(Cluster {
                    index: k,
                    members: members.clone(),
                    median: median(&select(x, members)),
                    hp: None,
                    minor: true,
                });
            }
        } else {
            let main: Vec<f64> = clusters.iter().map(|c| c.median).collect();
            for (_, members) in extra {
                for &i in members {
                    let d: Vec<f64> = main.iter().map(|&c| (x[i] - c).abs()).collect();
                    new_labels[i] = argmin(&d) as i64;
                }
            }
            notes.push(format!(
                "{} minor mode(s) folded into alleles (use --mosaic to keep)",
                extra.len()
            ));
            for c in clusters.iter_mut() {
                c.members = members_of(&new_labels, c.index);
                c.median = median(&select(x, &c.members));
            }
        }
    }
    finalise(clusters, new_labels, "GMM", notes)
}

/// Main alleles by size, minor modes after, per-read labels renumbered.
fn finalise(clusters: Vec<Cluster>, labels: Vec<i64>, method: &str, notes: Vec<String>) -> Assignment {
    let mut main: Vec<Cluster> = clusters.iter().filter(|c| !c.minor).cloned().collect();
    let mut minor: Vec<Cluster> = clusters.iter().filter(|c| c.minor).cloned().collect();
    main.sort_by(|a, b| a.median.partial_cmp(&b.median).unwrap());
    minor.sort_by(|a, b| a.median.partial_cmp(&b.median).unwrap());
    let mut ordered: Vec<Cluster> = main.into_iter().chain(minor).collect();
    let remap: Vec<(i64, i64)> = ordered.iter().enumerate().map(|(k, c)| (c.index, k as i64)).collect();
    let new_labels: Vec<i64> = labels
        .iter()
        .map(|&l| {
            if l >= 0 {
                remap.iter().rev().find(|r| r.0 == l).map(|r| r.1).unwrap_or(-1)
            } else {
                -1
            }
        })
        .collect();
    for (k, c) in ordered.iter_mut().enumerate() {
        c.index = k as i64;
        c.members = members_of(&new_labels, k as i64);
    }
    Assignment {
        clusters: ordered,
        labels: new_labels,
        method: method.to_string(),
        notes,
    }
}

pub struct Spread {
    pub median: f64,
    pub mad: f64,
    pub p5: f64,
    pub p95: f64,
    pub tail_up_frac: f64,
    pub tail_down_frac: f64,
}

/// Median, MAD, central 90 % range, and the tail fractions.
pub fn spread_stats(v: &[f64]) -> Option<Spread> {
    if v.is_empty() {
        return None;
    }
    let med = median(v);
    let t = tol(med);
    let (p5, p95) = if v.len() > 1 {
        (percentile(v, 5.0), percentile(v, 95.0))
    } else {
        (med, med)
    };
    let dev: Vec<f64> = v.iter().map(|&x| (x - med).abs()).collect();
    let n = v.len() as f64;
    Some(Spread {
        median: med,
        mad: median(&dev),
        p5,
        p95,
        tail_up_frac: v.iter().filter(|&&x| x > med + t).count() as f64 / n,
        tail_down_frac: v.iter().filter(|&&x| x < med - t).count() as f64 / n,
    })
}
