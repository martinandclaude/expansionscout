//! An allele's sequence as a consensus of its reads, by a vote over
//! repeat-unit tokens: a port of `expansionscout/consensus.py`, whose module
//! docstring states the rules. Every step is the
//! Python's, with the same integer scores and the same tie-breaks, so the
//! two give the same sequence.

use std::collections::{BTreeMap, BTreeSet};

use crate::decompose::{Decomposition, MotifSet, LABEL_OTHER};

/// Fewer reads than this form no consensus.
pub const MIN_READS: usize = 3;
/// The first and last columns are voted by the reads framed at that end
/// when at least this many are.
pub const EDGE_MIN: usize = 2;

/// One of the allele's reads with a decomposition.
pub struct ConsensusRead<'a> {
    /// The read's run, as the read has it: tokens are cut from it.
    pub seq: &'a [u8],
    /// The same run in upper case, which the decomposition was made of.
    pub run_u: &'a [u8],
    pub dec: &'a Decomposition,
    /// The read's own unit count.
    pub units: i64,
}

/// How far a read's decomposition falls short of the read's own count beyond
/// max(1, 10 %) of that count: the Python's `shortfall`.
pub fn shortfall(units: i64, n_units: usize) -> f64 {
    ((units - n_units as i64) as f64 - 1f64.max(0.1 * units as f64)).max(0.0)
}

/// round(g / l) with halves up, in integers.
fn ntok(g: usize, l: usize) -> usize {
    (2 * g + l) / (2 * l)
}

/// Segment offset of each unit's first base. `unit_at` repeats the offset of
/// a skipped span for each of its 'o' units, as the Python repeats the unit
/// string; the Python's `starts` puts them a unit apart, and so does this.
fn unit_starts(dec: &Decomposition, l: usize) -> Vec<usize> {
    let mut out: Vec<usize> = Vec::with_capacity(dec.unit_at.len());
    for k in 0..dec.unit_at.len() {
        if k > 0
            && dec.labels[k] == LABEL_OTHER
            && dec.labels[k - 1] == LABEL_OTHER
            && dec.unit_at[k] == dec.unit_at[k - 1]
        {
            out.push(out[k - 1] + l);
        } else {
            out.push(dec.unit_at[k]);
        }
    }
    out
}

/// (start, end, phase) of each unit of the read, relative to its run: the
/// Python's `read_units`.
fn read_units(ms: &MotifSet, r: &ConsensusRead) -> Vec<(usize, usize, Option<usize>)> {
    let l = ms.unit_len;
    let n = r.run_u.len();
    let rel: Vec<usize> = unit_starts(r.dec, l).iter().map(|&s| s - r.dec.start).collect();
    let wphase = r.dec.phase.as_ref().and_then(|(_, rot)| ms.unit_phase(rot));
    let mut out = Vec::with_capacity(rel.len());
    for (k, &s) in rel.iter().enumerate() {
        let e = if k + 1 < rel.len() { rel[k + 1] } else { n };
        if r.dec.labels[k] == LABEL_OTHER {
            out.push((s, e, None));
            continue;
        }
        let u = &r.run_u[s.min(n)..(s + l).min(n)];
        let p = if u.len() == l { ms.unit_phase(u) } else { None };
        match p {
            Some(p) => out.push((s, s + l, Some(p))),
            None => out.push((s, e, wphase)),
        }
    }
    out
}

/// The run cut into tokens at the unit boundaries its units imply in
/// `phase`: the Python's `cut_tokens`.
pub fn cut_tokens<'a>(seq: &'a [u8], units: &[(usize, usize, Option<usize>)], l: usize, phase: usize) -> Vec<&'a [u8]> {
    let n = seq.len();
    let mut cuts: BTreeSet<usize> = BTreeSet::new();
    for &(s, e, p) in units {
        match p {
            None => {}
            Some(p) if p == phase => {
                cuts.insert(s);
                cuts.insert(e);
            }
            Some(p) => {
                cuts.insert(s + (phase + l - p) % l);
            }
        }
    }
    let mut bounds = vec![0];
    bounds.extend(cuts.into_iter().filter(|&c| 0 < c && c < n));
    bounds.push(n);
    let mut pieces: Vec<(usize, usize)> = Vec::new();
    let mut pending: Option<usize> = None;
    for w in bounds.windows(2) {
        let (mut a, b) = (w[0], w[1]);
        if ntok(b - a, l) == 0 {
            if let Some(last) = pieces.last_mut() {
                last.1 = b;
            } else if pending.is_none() {
                pending = Some(a);
            }
            continue;
        }
        if let Some(p) = pending.take() {
            a = p;
        }
        pieces.push((a, b));
    }
    if pieces.is_empty() {
        return if n > 0 { vec![seq] } else { Vec::new() };
    }
    let mut toks = Vec::new();
    for (a, b) in pieces {
        let k = ntok(b - a, l).max(1);
        for j in 0..k - 1 {
            toks.push(&seq[a + j * l..a + (j + 1) * l]);
        }
        toks.push(&seq[a + (k - 1) * l..b]);
    }
    toks
}

/// Levenshtein distance.
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0; b.len() + 1];
    for (i, &ca) in a.iter().enumerate() {
        cur[0] = i + 1;
        for (j, &cb) in b.iter().enumerate() {
            cur[j + 1] = (prev[j] + usize::from(ca != cb)).min(prev[j + 1] + 1).min(cur[j] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

/// Global alignment of `s` to the centre `c`, unit costs: the base of `s`
/// aligned to each centre position (or `-`), and the bases inserted before
/// each (and after the last). The Python's `_align`, with its traceback.
fn align(c: &[u8], s: &[u8]) -> (Vec<u8>, Vec<Vec<u8>>) {
    let (n, m) = (c.len(), s.len());
    let w = m + 1;
    let mut d = vec![0usize; (n + 1) * w];
    for (j, x) in d.iter_mut().take(w).enumerate() {
        *x = j;
    }
    for i in 1..=n {
        d[i * w] = i;
        for j in 1..=m {
            let diag = d[(i - 1) * w + j - 1] + usize::from(c[i - 1] != s[j - 1]);
            d[i * w + j] = diag.min(d[(i - 1) * w + j] + 1).min(d[i * w + j - 1] + 1);
        }
    }
    let mut col = vec![b'-'; n];
    let mut ins: Vec<Vec<u8>> = vec![Vec::new(); n + 1];
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        if i > 0 && j > 0 && d[i * w + j] == d[(i - 1) * w + j - 1] + usize::from(c[i - 1] != s[j - 1]) {
            col[i - 1] = s[j - 1];
            i -= 1;
            j -= 1;
        } else if i > 0 && d[i * w + j] == d[(i - 1) * w + j] + 1 {
            i -= 1;
        } else {
            ins[i].push(s[j - 1]);
            j -= 1;
        }
    }
    for x in &mut ins {
        x.reverse();
    }
    (col, ins)
}

/// Distinct values in first-seen order, with their counts.
fn tally<T: PartialEq + Clone>(xs: impl IntoIterator<Item = T>) -> Vec<(T, usize)> {
    let mut out: Vec<(T, usize)> = Vec::new();
    for x in xs {
        match out.iter_mut().find(|e| e.0 == x) {
            Some(e) => e.1 += 1,
            None => out.push((x, 1)),
        }
    }
    out
}

/// Base-level consensus of a column whose tokens all differ: the Python's
/// `_star`, with the same centre, alignment and tie-breaks.
fn star(strings: &[&[u8]]) -> Vec<u8> {
    let k = strings.len();
    if k == 1 {
        return strings[0].to_vec();
    }
    let mut tot = vec![0usize; k];
    for a in 0..k {
        for b in a + 1..k {
            let d = edit_distance(strings[a], strings[b]);
            tot[a] += d;
            tot[b] += d;
        }
    }
    let ci = (0..k).min_by_key(|&i| (tot[i], i)).unwrap();
    let c = strings[ci];
    let mut cols: Vec<Vec<u8>> = vec![Vec::new(); c.len()];
    let mut slots: Vec<Vec<Vec<u8>>> = vec![Vec::new(); c.len() + 1];
    for s in strings {
        let (col, ins) = align(c, s);
        for (j, x) in col.into_iter().enumerate() {
            cols[j].push(x);
        }
        for (j, x) in ins.into_iter().enumerate() {
            slots[j].push(x);
        }
    }
    let rank = |x: u8| match x {
        b'A' => 0usize,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        b'-' => 4,
        _ => 5 + x as usize,
    };
    let mut out = Vec::new();
    for j in 0..=c.len() {
        let cnt = tally(slots[j].iter().cloned());
        let (best, _) = cnt
            .iter()
            .min_by_key(|(x, n)| (std::cmp::Reverse(*n), !x.is_empty(), x.len(), x.clone()))
            .unwrap();
        out.extend_from_slice(best);
        if j < c.len() {
            let cnt = tally(cols[j].iter().copied());
            let &(best, _) = cnt
                .iter()
                .min_by_key(|&&(x, n)| (std::cmp::Reverse(n), x != c[j], rank(x)))
                .unwrap();
            if best != b'-' {
                out.push(best);
            }
        }
    }
    out
}

/// Consensus run sequence of an allele, or None to keep the representative
/// read's. `reads`: the allele's members with a decomposition, in order;
/// `rep_index`: the representative read's index among them.
pub fn allele_consensus(reads: &[ConsensusRead], ms: &MotifSet, rep_index: usize) -> Option<Vec<u8>> {
    if reads.len() < MIN_READS {
        return None;
    }
    let keep: Vec<usize> = (0..reads.len())
        .filter(|&i| shortfall(reads[i].units, reads[i].dec.n_units()) == 0.0)
        .collect();
    if keep.len() < MIN_READS {
        return None;
    }
    let rep = keep.iter().position(|&i| i == rep_index)?;
    let reads: Vec<&ConsensusRead> = keep.iter().map(|&i| &reads[i]).collect();
    let l = ms.unit_len;
    let per: Vec<Vec<(usize, usize, Option<usize>)>> = reads.iter().map(|r| read_units(ms, r)).collect();
    let mut phases: BTreeMap<usize, usize> = BTreeMap::new();
    for units in &per {
        for &(_, _, p) in units {
            if let Some(p) = p {
                *phases.entry(p).or_insert(0) += 1;
            }
        }
    }
    // the commonest, the smaller on a tie
    let phase = *phases.iter().min_by_key(|&(&p, &n)| (std::cmp::Reverse(n), p))?.0;
    let mut toks: Vec<Vec<&[u8]>> = Vec::with_capacity(reads.len());
    let mut framed_start = Vec::with_capacity(reads.len());
    let mut framed_end = Vec::with_capacity(reads.len());
    for (r, units) in reads.iter().zip(&per) {
        toks.push(cut_tokens(r.seq, units, l, phase));
        framed_start.push(units.first().is_some_and(|u| u.2 == Some(phase) && u.0 == 0));
        framed_end.push(units.last().is_some_and(|u| u.2 == Some(phase) && u.1 == r.seq.len()));
    }
    let mut held: BTreeMap<usize, usize> = BTreeMap::new();
    let mut framed: BTreeMap<usize, usize> = BTreeMap::new();
    for (i, t) in toks.iter().enumerate() {
        *held.entry(t.len()).or_insert(0) += 1;
        if framed_start[i] && framed_end[i] {
            *framed.entry(t.len()).or_insert(0) += 1;
        }
    }
    let n_rep = toks[rep].len();
    let n = *held
        .keys()
        .min_by_key(|&&c| {
            (
                std::cmp::Reverse(framed.get(&c).copied().unwrap_or(0)),
                std::cmp::Reverse(held[&c]),
                c.abs_diff(n_rep),
                std::cmp::Reverse(c),
            )
        })
        .unwrap();
    let voters: Vec<usize> = (0..toks.len()).filter(|&i| toks[i].len() == n).collect();
    if n == 0 || voters.len() < MIN_READS {
        return None;
    }
    let rep_toks = if toks[rep].len() == n { Some(&toks[rep]) } else { None };
    let in_phase = |s: &[u8]| s.len() == l && ms.unit_phase(&s.to_ascii_uppercase()) == Some(phase);
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(n);
    for c in 0..n {
        let mut vs = voters.clone();
        if c == 0 || c == n - 1 {
            let sel: Vec<usize> = voters
                .iter()
                .copied()
                .filter(|&i| (c != 0 || framed_start[i]) && (c != n - 1 || framed_end[i]))
                .collect();
            if sel.len() >= EDGE_MIN {
                vs = sel;
            }
        }
        let col: Vec<&[u8]> = vs.iter().map(|&i| toks[i][c]).collect();
        // distinct tokens in first-seen order, so the index is `first`
        let cnt = tally(col.iter().copied());
        let rep_tok = rep_toks.map(|t| t[c]);
        let (best, best_n) = cnt
            .iter()
            .enumerate()
            .min_by_key(|&(first, &(s, k))| (std::cmp::Reverse(k), !in_phase(s), rep_tok != Some(s), first))
            .map(|(_, &x)| x)
            .unwrap();
        out.push(if best_n >= 2 { best.to_vec() } else { star(&col) });
    }
    let impure = |ts: &mut dyn Iterator<Item = &[u8]>| {
        ts.filter(|t| t.len() != l || ms.match_unit(&t.to_ascii_uppercase()).is_none())
            .count()
    };
    if n > n_rep && impure(&mut out.iter().map(|t| t.as_slice())) > impure(&mut toks[rep].iter().copied()) {
        return None;
    }
    Some(out.concat())
}

#[cfg(test)]
mod tests {
    use super::*;

    // The same sequences and expectations as tests/test_consensus.py.
    const LEFT: &str = "TTACT";
    const RIGHT: &str = "TTGACAT";
    const CEL: &str = "GGCCCCCCCCGTGCCGCCCACGGGTGACTCCGG";
    const FLANK_L: &str = "GCATCTTAGAGAGCCTAACGTGCCCTTTAAAGTCGCTTGCTTGTAAACAAGGCTGAGATCTCTTGT";
    const FLANK_R: &str = "CCGGACGTAACTCCCTATCAGGAATCGACCCCTACCCATAAGATGGGCTTCAAGTAATTGATTAGG";

    fn motifs(l: usize, m: &str) -> MotifSet {
        MotifSet::new(l, &[m.to_string()], &[], &[], &[]).unwrap()
    }

    fn sub(unit: &str, p: usize) -> String {
        let mut u: Vec<u8> = unit.bytes().collect();
        u[p] = match u[p] {
            b'A' => b'C',
            b'C' => b'A',
            b'G' => b'T',
            _ => b'G',
        };
        String::from_utf8(u).unwrap()
    }

    fn edited(tract: &str, kind: &str, pos: usize) -> String {
        let mut t: Vec<u8> = tract.bytes().collect();
        match kind {
            "sub" => t[pos] = b'T',
            "del" => {
                t.remove(pos);
            }
            _ => t.insert(pos, b'T'),
        }
        String::from_utf8(t).unwrap()
    }

    /// Each read's decomposition within the tract window, as the engine
    /// makes it, and the consensus of them with the first as representative.
    fn consensus(ms: &MotifSet, segments: &[String], window: (usize, usize), units: i64) -> Option<String> {
        let decs: Vec<Decomposition> = segments
            .iter()
            .map(|s| {
                ms.decompose(s.as_bytes(), Some(window), Some(window.1 - window.0))
                    .unwrap()
            })
            .collect();
        let reads: Vec<ConsensusRead> = segments
            .iter()
            .zip(&decs)
            .map(|(s, d)| ConsensusRead {
                seq: &s.as_bytes()[d.start..d.end],
                run_u: &s.as_bytes()[d.start..d.end],
                dec: d,
                units,
            })
            .collect();
        allele_consensus(&reads, ms, 0).map(|c| String::from_utf8(c).unwrap())
    }

    #[test]
    fn tokens_are_cut_in_the_allele_phase() {
        let ms = motifs(3, "CAG");
        let seg = format!("{LEFT}{}CA{}{RIGHT}", "CAG".repeat(3), "CAG".repeat(4));
        let d = ms.decompose(seg.as_bytes(), Some((5, 28)), Some(23)).unwrap();
        let run = &seg.as_bytes()[d.start..d.end];
        assert_eq!((d.start, d.end, run), (5, 26, &b"CAGCAGCAGCACAGCAGCAGC"[..]));
        let r = ConsensusRead {
            seq: run,
            run_u: run,
            dec: &d,
            units: 8,
        };
        let units = read_units(&ms, &r);
        let shown: Vec<(&[u8], Option<usize>)> = units.iter().map(|&(a, b, p)| (&run[a..b], p)).collect();
        assert_eq!(
            shown,
            vec![
                (&b"CAG"[..], Some(0)),
                (b"CAG", Some(0)),
                (b"CAG", Some(0)),
                (b"CAC", None),
                (b"AGC", Some(1)),
                (b"AGC", Some(1)),
                (b"AGC", Some(1)),
            ]
        );
        let toks: Vec<&[u8]> = cut_tokens(run, &units, 3, 0);
        let want: Vec<&[u8]> = vec![b"CAG", b"CAG", b"CAG", b"CAC", b"AG", b"CAG", b"CAGC"];
        assert_eq!(toks, want);
    }

    #[test]
    fn a_column_vote_restores_the_tract() {
        let ms = motifs(3, "CAG");
        let truth = "CAG".repeat(10);
        let edits = [("sub", 7), ("del", 10), ("sub", 17), ("ins", 13), ("sub", 19)];
        let reads: Vec<String> = edits
            .iter()
            .map(|&(k, p)| format!("{LEFT}{}{RIGHT}", edited(&truth, k, p)))
            .collect();
        assert_eq!(consensus(&ms, &reads, (5, 35), 10), Some(truth));
        assert_eq!(consensus(&ms, &reads[..2], (5, 35), 10), None);
    }

    #[test]
    fn a_long_motif_column_is_voted_base_by_base() {
        let ms = motifs(33, CEL);
        let reads: Vec<String> = (0..5)
            .map(|r| {
                let rest: String = (1..4).map(|k| sub(CEL, (7 * r + 5 * k + 3) % 33)).collect();
                format!("{FLANK_L}{CEL}{rest}{FLANK_R}")
            })
            .collect();
        assert_eq!(consensus(&ms, &reads, (66, 198), 4), Some(CEL.repeat(4)));
    }

    #[test]
    fn shortfall_matches_the_representative_rule() {
        assert_eq!(shortfall(30, 30), 0.0);
        assert_eq!(shortfall(30, 27), 0.0);
        assert_eq!(shortfall(30, 26), 1.0);
        assert_eq!(shortfall(5, 4), 0.0);
        assert_eq!(shortfall(5, 3), 1.0);
    }
}
