//! Motif-aware decomposition of a read segment into repeat units.
//!
//! Port of `expansionscout/decompose.py`, which documents the conventions and
//! the reasoning behind every threshold. Sequences here are uppercase byte
//! strings in reference orientation.
//!
//! The run finder is a regex over every rotation of every motif. Python's `re`
//! and the `regex` crate agree on it: all alternatives have the same length,
//! so leftmost-first and backtracking semantics give the same spans.

use std::collections::{HashMap, HashSet};

use crate::npstat::py_round;

pub const LABEL_CANONICAL: u8 = b'C';
pub const LABEL_PATHOGENIC: u8 = b'P';
pub const LABEL_BENIGN: u8 = b'B';
pub const LABEL_INTERRUPTION: u8 = b'I';
pub const LABEL_OTHER: u8 = b'o';

fn iupac(c: u8) -> Option<&'static [u8]> {
    Some(match c {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"GC",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => return None,
    })
}

fn is_ambiguous(m: &[u8]) -> bool {
    m.iter().any(|c| !matches!(c, b'A' | b'C' | b'G' | b'T'))
}

/// Rotations in first-seen order. Python iterates a set here, but nothing
/// downstream depends on the order: see `MotifSet::new`.
fn rotations(m: &[u8]) -> Vec<Vec<u8>> {
    let mut out: Vec<Vec<u8>> = Vec::new();
    for i in 0..m.len() {
        let r: Vec<u8> = m[i..].iter().chain(m[..i].iter()).copied().collect();
        if !out.contains(&r) {
            out.push(r);
        }
    }
    out
}

/// `(edits, length)` of the best alignment of all of `unit` against a prefix
/// of `text`, or None if every one needs more than `max_edits` edits. Of
/// equally good lengths the one nearest `unit.len()` is taken, then the
/// shorter. See the Python's `align_unit`.
pub fn align_unit(unit: &[u8], text: &[u8], max_edits: usize) -> Option<(usize, usize)> {
    let l = unit.len();
    let band = max_edits;
    let n = text.len().min(l + band);
    if n + band < l {
        return None;
    }
    let cap = max_edits + 1;
    let mut prev: Vec<usize> = (0..=n).map(|b| b.min(cap)).collect();
    let mut cur = vec![cap; n + 1];
    for a in 1..=l {
        let lo = a.saturating_sub(band);
        let hi = (a + band).min(n);
        cur.iter_mut().for_each(|c| *c = cap);
        if lo == 0 {
            cur[0] = a.min(cap);
        }
        let mut best = if lo == 0 { cur[0] } else { cap };
        let ua = unit[a - 1];
        for b in lo.max(1)..=hi {
            let c = (prev[b - 1] + usize::from(text[b - 1] != ua))
                .min(prev[b] + 1)
                .min(cur[b - 1] + 1)
                .min(cap);
            cur[b] = c;
            best = best.min(c);
        }
        if best == cap {
            return None;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    let mut choice: Option<(usize, usize, usize)> = None;
    for (b, &e) in prev.iter().enumerate().skip(l.saturating_sub(band)) {
        if e <= max_edits {
            let key = (e, b.abs_diff(l), b);
            if choice.is_none_or(|c| key < c) {
                choice = Some(key);
            }
        }
    }
    choice.map(|(e, _, b)| (e, b))
}

#[derive(Clone, Debug)]
pub struct Decomposition {
    pub start: usize,
    pub end: usize,
    pub labels: Vec<u8>,
    /// Segment offset of each unit's first base, parallel to `labels`. A unit
    /// is `seq[u..u + unit_len]`, cut short at the end of `seq` as a Python
    /// slice is; an 'o' unit repeats the offset it was recorded at, as Python
    /// repeats the unit string, and a unit taken by alignment, which may span
    /// a few bases more or less, is the `unit_len` bases from its start.
    pub unit_at: Vec<usize>,
    /// (label, rotation): the first exact unit the walk met, or the motif
    /// rotation an identity seed was made with. A long motif's units are
    /// accepted by identity to it, in `walk` and in `extend`.
    pub phase: Option<(u8, Vec<u8>)>,
    /// True once the run has been lengthened within a spanning read's
    /// bounds: every unit it takes after that must meet them too.
    pub bridged: bool,
}

impl Decomposition {
    pub fn n_units(&self) -> usize {
        self.labels.len()
    }
    pub fn n_bp(&self) -> usize {
        self.end - self.start
    }
    pub fn count(&self, label: u8) -> usize {
        self.labels.iter().filter(|&&l| l == label).count()
    }
    pub fn noise_frac(&self) -> f64 {
        if self.labels.is_empty() {
            0.0
        } else {
            self.count(LABEL_OTHER) as f64 / self.labels.len() as f64
        }
    }
    /// 1-based unit positions carrying `label`.
    pub fn positions(&self, label: u8) -> Vec<usize> {
        self.labels
            .iter()
            .enumerate()
            .filter(|(_, &l)| l == label)
            .map(|(i, _)| i + 1)
            .collect()
    }
    pub fn unit<'a>(&self, seq: &'a [u8], i: usize, unit_len: usize) -> &'a [u8] {
        let at = self.unit_at[i];
        &seq[at..(at + unit_len).min(seq.len())]
    }
    pub fn labels_str(&self) -> String {
        String::from_utf8_lossy(&self.labels).into_owned()
    }
}

#[derive(Clone, Debug)]
struct Ambiguous {
    allowed: Vec<&'static [u8]>,
    label: u8,
    /// The motif this rotation's pattern came from.
    base: Vec<u8>,
}

#[derive(Clone, Debug)]
pub struct MotifSet {
    pub unit_len: usize,
    pub min_region_units: usize,
    pub max_interruption: Option<usize>,
    pub approx_above_len: usize,
    pub approx_identity: f64,
    labels: HashMap<Vec<u8>, u8>,
    /// Rotation -> the motif it is a rotation of, beside `labels`.
    bases: HashMap<Vec<u8>, Vec<u8>>,
    ambiguous: Vec<Ambiguous>,
    exact_re: regex::bytes::Regex,
    pub min_seed: usize,
    pub composition: bool,
    /// The most edits a unit aligned to the phase may carry.
    pub max_edits: usize,
    /// First `SEED_K` bases -> [(label, rotation)], for an identity seed.
    seed_index: HashMap<Vec<u8>, Vec<(u8, Vec<u8>)>>,
}

/// A run is seeded by identity only above this motif length.
pub const SEED_ABOVE_LEN: usize = 40;
/// Rotations are looked up for an identity seed by their first bases.
pub const SEED_K: usize = 12;

impl MotifSet {
    pub fn new(
        unit_len: usize,
        canonical: &[String],
        pathogenic: &[String],
        benign: &[String],
        interruption: &[String],
    ) -> Result<MotifSet, String> {
        let ordered: [(u8, &[String]); 4] = [
            (LABEL_CANONICAL, canonical),
            (LABEL_PATHOGENIC, pathogenic),
            (LABEL_BENIGN, benign),
            (LABEL_INTERRUPTION, interruption),
        ];
        let mut labels: HashMap<Vec<u8>, u8> = HashMap::new();
        let mut bases: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
        let mut ambiguous = Vec::new();
        let mut patterns: Vec<String> = Vec::new();
        let mut seen_patterns: Vec<Vec<&'static [u8]>> = Vec::new();
        let mut has_canonical = false;
        for (label, motifs) in ordered {
            for m in motifs {
                let m = m.to_ascii_uppercase().into_bytes();
                if m.len() != unit_len {
                    continue;
                }
                if label == LABEL_CANONICAL {
                    has_canonical = true;
                }
                if is_ambiguous(&m) {
                    // First label wins for a pattern, as in the Python: the
                    // polyalanine loci write GCN with a pathogenic NGC, which is
                    // the same motif rotated.
                    for r in rotations(&m) {
                        let allowed: Vec<&'static [u8]> = r
                            .iter()
                            .map(|&c| {
                                iupac(c).ok_or_else(|| {
                                    format!(
                                        "not a nucleotide or IUPAC code: '{}' in '{}'",
                                        c as char,
                                        String::from_utf8_lossy(&m)
                                    )
                                })
                            })
                            .collect::<Result<_, _>>()?;
                        if seen_patterns.contains(&allowed) {
                            continue;
                        }
                        let pat: String = allowed
                            .iter()
                            .map(|a| {
                                if a.len() == 1 {
                                    (a[0] as char).to_string()
                                } else {
                                    format!("[{}]", String::from_utf8_lossy(a))
                                }
                            })
                            .collect();
                        seen_patterns.push(allowed.clone());
                        ambiguous.push(Ambiguous {
                            allowed,
                            label,
                            base: m.clone(),
                        });
                        patterns.push(pat);
                    }
                } else {
                    for r in rotations(&m) {
                        if let std::collections::hash_map::Entry::Vacant(e) = labels.entry(r.clone()) {
                            e.insert(label);
                            bases.insert(r.clone(), m.clone());
                            patterns.push(regex::escape(&String::from_utf8_lossy(&r)));
                        }
                    }
                }
            }
        }
        if !has_canonical {
            return Err("MotifSet needs at least one canonical motif of unit_len".into());
        }
        let min_seed = if unit_len <= 12 { 2 } else { 1 };
        let exact_re = regex::bytes::RegexBuilder::new(&format!("(?:{}){{{},}}", patterns.join("|"), min_seed))
            .unicode(false)
            .size_limit(1 << 28)
            .dfa_size_limit(1 << 28)
            .build()
            .map_err(|e| format!("motif regex: {e}"))?;
        let composition =
            labels.values().any(|&l| l == LABEL_PATHOGENIC) || ambiguous.iter().any(|a| a.label == LABEL_PATHOGENIC);
        let approx_identity = 0.80;
        let max_edits = (0..=unit_len)
            .filter(|&e| (unit_len - e) as f64 / unit_len as f64 >= approx_identity)
            .max()
            .unwrap_or(0);
        // Built in a fixed order (motif priority, then rotation offset), as
        // the Python builds it.
        let mut seed_index: HashMap<Vec<u8>, Vec<(u8, Vec<u8>)>> = HashMap::new();
        if unit_len > SEED_ABOVE_LEN {
            let mut seen: HashSet<Vec<u8>> = HashSet::new();
            for (label, motifs) in ordered {
                for m in motifs {
                    let m = m.to_ascii_uppercase().into_bytes();
                    if m.len() != unit_len || is_ambiguous(&m) {
                        continue;
                    }
                    for k in 0..m.len() {
                        let r: Vec<u8> = m[k..].iter().chain(m[..k].iter()).copied().collect();
                        if seen.insert(r.clone()) {
                            seed_index.entry(r[..SEED_K].to_vec()).or_default().push((label, r));
                        }
                    }
                }
            }
        }
        Ok(MotifSet {
            unit_len,
            min_region_units: 3,
            max_interruption: None,
            approx_above_len: 12,
            approx_identity,
            labels,
            bases,
            ambiguous,
            exact_re,
            min_seed,
            composition,
            max_edits,
            seed_index,
        })
    }

    /// Label for one unit, or None. Literal lookup first, then any ambiguous
    /// motif patterns in registration order.
    pub fn match_unit(&self, u: &[u8]) -> Option<u8> {
        self.match_base(u).map(|(l, _)| l)
    }

    /// Label and base motif for one unit, or None: the Python's
    /// `MotifSet.match`.
    pub fn match_base(&self, u: &[u8]) -> Option<(u8, &[u8])> {
        if let Some(&l) = self.labels.get(u) {
            return Some((l, &self.bases[u]));
        }
        for a in &self.ambiguous {
            if a.allowed.len() == u.len() && u.iter().zip(&a.allowed).all(|(c, set)| set.contains(c)) {
                return Some((a.label, &a.base));
            }
        }
        None
    }

    /// r such that the exact unit `u` is `base[r..] + base[..r]` of the motif
    /// it matches, the smallest if several, or None if `u` is not an exact
    /// unit: the Python's `consensus.unit_phase`.
    pub fn unit_phase(&self, u: &[u8]) -> Option<usize> {
        let (_, base) = self.match_base(u)?;
        let (n, amb) = (base.len(), is_ambiguous(base));
        if u.len() != n {
            return None;
        }
        (0..n).find(|&r| {
            (0..n).all(|j| {
                let c = base[(r + j) % n];
                if amb {
                    iupac(c).is_some_and(|set| set.contains(&u[j]))
                } else {
                    c == u[j]
                }
            })
        })
    }

    fn span(&self) -> usize {
        self.max_interruption.unwrap_or_else(|| 30.max(6 * self.unit_len))
    }

    /// Merged repeat regions in `seq`, in order, as (start, end).
    pub fn find_regions(&self, seq: &[u8], gap_tol: Option<usize>) -> Vec<(usize, usize)> {
        let gap_tol = gap_tol.unwrap_or(2 * self.unit_len + 2) as i64;
        let mut it = self.exact_re.find_iter(seq).map(|m| (m.start(), m.end()));
        let Some((mut cs, mut ce)) = it.next() else {
            return vec![];
        };
        let mut regions = Vec::new();
        for (a, b) in it {
            if a as i64 - ce as i64 <= gap_tol {
                ce = b;
            } else {
                regions.push((cs, ce));
                cs = a;
                ce = b;
            }
        }
        regions.push((cs, ce));
        regions
    }

    /// Envelope of the repeat tract in `seq`, or None.
    pub fn find_run(&self, seq: &[u8], gap_tol: Option<usize>) -> Option<(usize, usize)> {
        let regions = self.find_regions(seq, gap_tol);
        if regions.is_empty() {
            return None;
        }
        let need = if self.unit_len > self.approx_above_len {
            1
        } else {
            self.min_region_units
        };
        let solid: Vec<(usize, usize)> = regions
            .iter()
            .copied()
            .filter(|r| self.walk(seq, r.0, r.1).n_units() >= need)
            .collect();
        if solid.is_empty() {
            // max() returns the first of equal-length regions
            let mut best = regions[0];
            for &r in &regions[1..] {
                if r.1 - r.0 > best.1 - best.0 {
                    best = r;
                }
            }
            return Some(best);
        }
        let span = self.span() as i64;
        let (env_start, mut env_end) = solid[0];
        for &(a, b) in &solid[1..] {
            if a as i64 - env_end as i64 > span {
                break;
            }
            env_end = b;
        }
        Some((env_start, env_end))
    }

    fn identity(a: &[u8], b: &[u8]) -> f64 {
        if b.is_empty() {
            return 0.0;
        }
        let same = a.iter().zip(b).filter(|(x, y)| x == y).count();
        same as f64 / b.len() as f64
    }

    /// Walk units from `start` to `end`, tolerating isolated errors. For a
    /// motif longer than `approx_above_len`, a unit is also accepted by
    /// identity to the first exact unit of the walk, and the result carries
    /// that unit as its phase for `extend`.
    pub fn walk(&self, seq: &[u8], start: usize, end: usize) -> Decomposition {
        self.walk_bounded(seq, start, end, None, None)
    }

    /// `walk`, bounded as `extend` is by a spanning read's `window` and
    /// `max_bp`: a unit taken by identity must meet the bounds and may also
    /// be taken by alignment; once one has been, every later unit must meet
    /// them too. See the Python's `MotifSet.walk`.
    pub fn walk_bounded(
        &self,
        seq: &[u8],
        start: usize,
        end: usize,
        window: Option<(usize, usize)>,
        max_bp: Option<usize>,
    ) -> Decomposition {
        let l = self.unit_len;
        let approx = l > self.approx_above_len;
        let bounded = self.bounded(window, max_bp);
        let mut phase: Option<(u8, Vec<u8>)> = None;
        let mut bridged = false;
        let mut i = start;
        let mut labels = Vec::new();
        let mut unit_at = Vec::new();
        while i + l <= end {
            let u = &seq[i..i + l];
            let mut hit = self.match_unit(u);
            if let (Some(h), None) = (hit, &phase) {
                phase = Some((h, u.to_vec()));
            }
            let mut step = l;
            if hit.is_some() {
                if bridged && !self.may_approx(i, i + l, start, i + l, window, max_bp) {
                    break;
                }
            } else if approx {
                if let Some(n) = self.identity_unit(seq, i, end, phase.as_ref(), (start, i), window, max_bp, false) {
                    hit = phase.as_ref().map(|p| p.0);
                    step = n;
                    bridged = bounded;
                }
            }
            if let Some(h) = hit {
                labels.push(h);
                unit_at.push(i);
                i += step;
                continue;
            }
            let (ii, li, ei) = (i as i64, l as i64, end as i64);
            let mut nxt: Option<i64> = None;
            for j in [
                ii + li,
                ii + li - 1,
                ii + li + 1,
                ii + 1,
                ii + li - 2,
                ii + li + 2,
                ii + 2,
            ] {
                if j <= ii {
                    continue;
                }
                if j + li <= ei && self.match_unit(&seq[j as usize..(j + li) as usize]).is_some() {
                    nxt = Some(j);
                    break;
                }
            }
            if nxt.is_none() {
                let span = self.span() as i64;
                let hi = (ii + span + 1).min(ei - li + 1);
                let mut j = ii + 1;
                while j < hi {
                    let ju = j as usize;
                    if self.match_unit(&seq[ju..ju + l]).is_some()
                        && (j + 2 * li > ei || self.match_unit(&seq[ju + l..ju + 2 * l]).is_some())
                    {
                        nxt = Some(j);
                        break;
                    }
                    j += 1;
                }
            }
            let Some(nxt) = nxt else { break };
            let nxt_u = nxt as usize;
            if bridged && !self.may_approx(i, nxt_u + l, start, nxt_u + l, window, max_bp) {
                break;
            }
            let skipped = (nxt - ii) as usize;
            if skipped + 1 >= l {
                for _ in 0..1.max(skipped / l) {
                    labels.push(LABEL_OTHER);
                    unit_at.push(i);
                }
            }
            i = nxt_u;
        }
        Decomposition {
            start,
            end: i,
            labels,
            unit_at,
            phase,
            bridged,
        }
    }

    /// Does a spanning read's alignment bound the run? See `extend`.
    fn bounded(&self, window: Option<(usize, usize)>, max_bp: Option<usize>) -> bool {
        window.is_some() && max_bp.is_some_and(|m| m > 0) && self.unit_len > self.approx_above_len
    }

    /// May a bounded run take the unit at `[a, b)`, making the run
    /// `[new_start, new_end)`? See `extend`.
    fn may_approx(
        &self,
        a: usize,
        b: usize,
        new_start: usize,
        new_end: usize,
        window: Option<(usize, usize)>,
        max_bp: Option<usize>,
    ) -> bool {
        if !self.bounded(window, max_bp) {
            return false;
        }
        let ((w0, w1), max_bp) = (window.unwrap(), max_bp.unwrap());
        w0 <= a && b <= w1 && new_end - new_start <= max_bp
    }

    /// Length of the unit that starts at `at` (ends there if `backward`) and
    /// does not reach past `limit`, if it is taken by identity to the run's
    /// phase: position by position, or, within a spanning read's bounds, by
    /// alignment. `run` is the (start, end) it would lengthen. See the
    /// Python's `MotifSet._identity_unit`.
    #[allow(clippy::too_many_arguments)]
    fn identity_unit(
        &self,
        seq: &[u8],
        at: usize,
        limit: usize,
        phase: Option<&(u8, Vec<u8>)>,
        run: (usize, usize),
        window: Option<(usize, usize)>,
        max_bp: Option<usize>,
        backward: bool,
    ) -> Option<usize> {
        let (_, rot) = phase?;
        let l = self.unit_len;
        let bounded = self.bounded(window, max_bp);
        let fits = |n: usize| {
            let (a, b) = if backward { (at - n, at) } else { (at, at + n) };
            !bounded || self.may_approx(a, b, a.min(run.0), b.max(run.1), window, max_bp)
        };
        let whole = if backward { at >= limit + l } else { at + l <= limit };
        if whole {
            let u = if backward { &seq[at - l..at] } else { &seq[at..at + l] };
            if Self::identity(u, rot) >= self.approx_identity && fits(l) {
                return Some(l);
            }
        }
        if !bounded {
            return None;
        }
        let w = l + self.max_edits;
        let got = if backward {
            let lo = limit.max(at.saturating_sub(w));
            let text: Vec<u8> = seq[lo..at].iter().rev().copied().collect();
            let unit: Vec<u8> = rot.iter().rev().copied().collect();
            align_unit(&unit, &text, self.max_edits)
        } else {
            let hi = limit.min(at + w).max(at);
            align_unit(rot, &seq[at..hi], self.max_edits)
        };
        got.map(|g| g.1).filter(|&n| fits(n))
    }

    /// Grow a walked run outward while units still match.
    ///
    /// First the run takes exact motif rotations at both ends, without
    /// bounds. Then, for a motif longer than `approx_above_len` on a spanning
    /// read, it grows further by identity to the run's phase (position by
    /// position, or by alignment across indels), taking every unit, exact ones
    /// too, only within both bounds the read's alignment gives, because a
    /// real tract is often flanked by degenerate in-frame copies of its own
    /// repeat that identity would otherwise walk into: the unit must lie
    /// inside `window`, where the alignment places the tract (a unit
    /// straddling its edge is refused), and the run must stay within
    /// `max_bp`, the tract's length between the read's anchors, which unlike
    /// the window does not move with where the aligner places an indel. So a
    /// run this stage lengthened is never longer than `max_bp`, the length the
    /// read's unit count is taken from, and it only adds to the run exact
    /// extension gives. Without them (partial and split reads) extension is
    /// exact. A run already lengthened within the bounds skips the first
    /// stage.
    pub fn extend(
        &self,
        seq: &[u8],
        dec: Decomposition,
        window: Option<(usize, usize)>,
        max_bp: Option<usize>,
    ) -> Decomposition {
        let mut dec = dec;
        if !dec.bridged {
            dec = self.grow(seq, dec, None, None);
        }
        if self.bounded(window, max_bp) {
            dec = self.grow(seq, dec, window, max_bp);
        }
        dec
    }

    /// One stage of `extend`: exact units without bounds; or within them,
    /// exact units and those `identity_unit` takes, which may take an exact
    /// unit outside them at another length inside them when it is a
    /// different rotation from the phase.
    fn grow(
        &self,
        seq: &[u8],
        dec: Decomposition,
        window: Option<(usize, usize)>,
        max_bp: Option<usize>,
    ) -> Decomposition {
        let l = self.unit_len;
        let bounded = self.bounded(window, max_bp);
        let Decomposition {
            start: start0,
            end: end0,
            mut labels,
            mut unit_at,
            phase,
            bridged,
        } = dec;
        let (mut start, mut end) = (start0, end0);
        while end + l <= seq.len() {
            let mut hit = self.match_unit(&seq[end..end + l]);
            let mut step = l;
            if hit.is_some() && bounded && !self.may_approx(end, end + l, start, end + l, window, max_bp) {
                hit = None;
            }
            if hit.is_none() && bounded {
                if let Some(n) =
                    self.identity_unit(seq, end, seq.len(), phase.as_ref(), (start, end), window, max_bp, false)
                {
                    hit = phase.as_ref().map(|p| p.0);
                    step = n;
                }
            }
            let Some(h) = hit else { break };
            labels.push(h);
            unit_at.push(end);
            end += step;
        }
        let mut front_labels = Vec::new();
        let mut front_at = Vec::new();
        while start >= l {
            let mut hit = self.match_unit(&seq[start - l..start]);
            let mut step = l;
            if hit.is_some() && bounded && !self.may_approx(start - l, start, start - l, end, window, max_bp) {
                hit = None;
            }
            if hit.is_none() && bounded {
                if let Some(n) = self.identity_unit(seq, start, 0, phase.as_ref(), (start, end), window, max_bp, true) {
                    hit = phase.as_ref().map(|p| p.0);
                    step = n;
                }
            }
            let Some(h) = hit else { break };
            front_labels.push(h);
            front_at.push(start - step);
            start -= step;
        }
        front_labels.reverse();
        front_at.reverse();
        front_labels.extend(labels);
        front_at.extend(unit_at);
        let grew = (start, end) != (start0, end0);
        Decomposition {
            start,
            end,
            labels: front_labels,
            unit_at: front_at,
            phase,
            bridged: bridged || (bounded && grew),
        }
    }

    /// The first unit inside `window` that, with the unit after it, is taken
    /// by identity to one motif rotation, as a one-unit run already bounded;
    /// or None. See the Python's `MotifSet._identity_seed`.
    fn identity_seed(&self, seq: &[u8], window: (usize, usize), max_bp: usize) -> Option<Decomposition> {
        let l = self.unit_len;
        let (w0, w1) = (window.0, window.1.min(seq.len()));
        if 2 * l > max_bp || w1 < 2 * l {
            return None;
        }
        for j in w0..(w1 - 2 * l + 1).max(w0) {
            let Some(cands) = self.seed_index.get(&seq[j..j + SEED_K]) else {
                continue;
            };
            for (label, rot) in cands {
                let phase = (*label, rot.clone());
                let Some(n1) = self.identity_unit(seq, j, w1, Some(&phase), (j, j), Some(window), Some(max_bp), false)
                else {
                    continue;
                };
                if self
                    .identity_unit(
                        seq,
                        j + n1,
                        w1,
                        Some(&phase),
                        (j, j + n1),
                        Some(window),
                        Some(max_bp),
                        false,
                    )
                    .is_none()
                {
                    continue;
                }
                return Some(Decomposition {
                    start: j,
                    end: j + n1,
                    labels: vec![*label],
                    unit_at: vec![j],
                    phase: Some(phase),
                    bridged: true,
                });
            }
        }
        None
    }

    /// find_run, walk, then extend. `window` is the tract's `(start, end)`
    /// within `seq` as the read's alignment places it, and `max_bp` its
    /// length between the read's anchors, or None; see `extend`. With both,
    /// at a long motif, the seed is also walked within the bounds, and that
    /// run is kept when it contains the whole exact run, is no longer than
    /// `max_bp` and has no fewer units; where no exact copy seeds a run at
    /// all, a motif longer than `SEED_ABOVE_LEN` is seeded by identity inside
    /// the window. See the Python's `MotifSet.decompose`.
    pub fn decompose(
        &self,
        seq: &[u8],
        window: Option<(usize, usize)>,
        max_bp: Option<usize>,
    ) -> Option<Decomposition> {
        let bounded = self.bounded(window, max_bp);
        let Some(run) = self.find_run(seq, None) else {
            if bounded && self.unit_len > SEED_ABOVE_LEN {
                let seed = self.identity_seed(seq, window.unwrap(), max_bp.unwrap())?;
                return Some(self.extend(seq, seed, window, max_bp));
            }
            return None;
        };
        let mut at = self.seed_at_tract_start(seq, run, window);
        // Keep the seed where its exact run has more units, no more than the
        // tract holds and without leaving it: see the Python's decompose().
        if at != run {
            let exact = |(a, b): (usize, usize)| self.extend(seq, self.walk(seq, a, b), None, None);
            let own = exact(run);
            let (w0, w1) = window.unwrap();
            if own.end <= w1
                && own.n_units() as f64 <= py_round((w1 - w0) as f64 / self.unit_len as f64)
                && own.n_units() > exact(at).n_units()
            {
                at = run;
            }
        }
        Some(self.run_from_seed(seq, at, window, max_bp))
    }

    /// walk and extend from a find_run seed; see `decompose`.
    fn run_from_seed(
        &self,
        seq: &[u8],
        (a, b): (usize, usize),
        window: Option<(usize, usize)>,
        max_bp: Option<usize>,
    ) -> Decomposition {
        let bounded = self.bounded(window, max_bp);
        let exact = self.extend(seq, self.walk(seq, a, b), None, None);
        if !bounded {
            return exact;
        }
        let grown = self.extend(seq, exact.clone(), window, max_bp);
        let walked = self.extend(seq, self.walk_bounded(seq, a, b, window, max_bp), window, max_bp);
        if walked.start <= exact.start
            && walked.end >= exact.end
            && walked.n_bp() <= max_bp.unwrap()
            && walked.n_units() >= grown.n_units()
        {
            return walked;
        }
        grown
    }

    /// The seed, started where the read's alignment starts the tract when it
    /// begins less than a unit before that and a unit matches there, so an
    /// exact seed that starts in flank bases continuing the period does not
    /// shift the allele's sequence. See the Python's `_seed_at_tract_start`.
    fn seed_at_tract_start(&self, seq: &[u8], run: (usize, usize), window: Option<(usize, usize)>) -> (usize, usize) {
        let (r0, r1) = run;
        let Some((w0, _)) = window else {
            return run;
        };
        let l = self.unit_len;
        let n = seq.len();
        if r0 < w0 && w0 < r0 + l && self.match_unit(&seq[w0.min(n)..(w0 + l).min(n)]).is_some() {
            return (w0, r1.max(w0 + l));
        }
        run
    }

    /// Units the read demonstrably crossed, counted outward from the tract.
    pub fn count_units_crossed(&self, seq: &[u8], window: Option<(i64, i64)>, from_left: bool) -> usize {
        let min_units = 3;
        let n = seq.len() as i64;
        let (lo, hi) = window.unwrap_or((0, n));
        let lo = 0.max(lo.min(n));
        let hi = lo.max(hi.min(n));
        if hi - lo < (min_units * self.unit_len) as i64 {
            return 0;
        }
        let sub = &seq[lo as usize..hi as usize];
        let gap_tol = (6 * self.unit_len + 6) as i64;
        let mut regions = self.find_regions(sub, None);
        if regions.is_empty() {
            return 0;
        }
        if !from_left {
            regions.reverse();
        }
        let mut total = 0;
        let mut prev: Option<usize> = None;
        for (cs, ce) in regions {
            if let Some(p) = prev {
                let gap = if from_left {
                    cs as i64 - p as i64
                } else {
                    p as i64 - ce as i64
                };
                if gap > gap_tol {
                    break;
                }
            }
            total += self.walk(sub, cs, ce).n_units();
            prev = Some(if from_left { ce } else { cs });
        }
        if total >= min_units {
            total
        } else {
            0
        }
    }

    /// Longest repeat run in a soft-clipped tail, as a bound on the allele.
    pub fn count_units_loose(&self, seq: &[u8]) -> usize {
        let min_units = 3;
        let mut best = 0;
        for (cs, ce) in self.find_regions(seq, None) {
            let n = self.walk(seq, cs, ce).n_units();
            if n >= min_units {
                best = best.max(n);
            }
        }
        best
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The same sequences and expectations as tests/test_decompose.py, so the
    // two implementations pin the same numbers.
    const CEL: &str = "GGCCCCCCCCGTGCCGCCCACGGGTGACTCCGG";
    const FLANK_L: &str = "GCATCTTAGAGAGCCTAACGTGCCCTTTAAAGTCGCTTGCTTGTAAACAAGGCTGAGATCTCTTGT";
    const FLANK_R: &str = "CCGGACGTAACTCCCTATCAGGAATCGACCCCTACCCATAAGATGGGCTTCAAGTAATTGATTAGG";

    fn sub(unit: &str, positions: &[usize]) -> String {
        let mut u: Vec<u8> = unit.bytes().collect();
        for &p in positions {
            u[p] = match u[p] {
                b'A' => b'C',
                b'C' => b'A',
                b'G' => b'T',
                _ => b'G',
            };
        }
        String::from_utf8(u).unwrap()
    }

    fn cel() -> MotifSet {
        MotifSet::new(33, &[CEL.to_string()], &[], &[], &[]).unwrap()
    }

    fn span(ms: &MotifSet, seg: &str) -> (usize, usize, String) {
        let d = ms.decompose(seg.as_bytes(), None, None).unwrap();
        (d.start, d.end, d.labels_str())
    }

    /// Within a tract window and, unless given, a length equal to the
    /// window's, as a spanning read without boundary indels has.
    fn bounded(ms: &MotifSet, seg: &str, window: (usize, usize), max_bp: Option<usize>) -> (usize, usize, String) {
        let max_bp = max_bp.unwrap_or(window.1 - window.0);
        let d = ms.decompose(seg.as_bytes(), Some(window), Some(max_bp)).unwrap();
        (d.start, d.end, d.labels_str())
    }

    fn c(n: usize) -> String {
        "C".repeat(n)
    }

    #[test]
    fn walk_carries_the_phase_it_matched_by() {
        let ms = cel();
        let seg = format!(
            "{FLANK_L}{}{}{}{FLANK_R}",
            sub(CEL, &[32]),
            CEL.repeat(6),
            sub(CEL, &[5])
        );
        let (a, b) = ms.find_run(seg.as_bytes(), None).unwrap();
        let d = ms.walk(seg.as_bytes(), a, b);
        assert_eq!(d.phase, Some((LABEL_CANONICAL, seg.as_bytes()[99..132].to_vec())));
        assert_eq!(d.start, 99);
    }

    #[test]
    fn long_motif_edge_units_are_taken_by_identity_inside_the_tract() {
        let ms = cel();
        let seg = format!(
            "{FLANK_L}{}{}{}{FLANK_R}",
            sub(CEL, &[32]),
            CEL.repeat(6),
            sub(CEL, &[5])
        );
        assert_eq!(span(&ms, &seg), (99, 297, c(6)));
        assert_eq!(bounded(&ms, &seg, (66, 330), None), (66, 330, c(8)));
        // both bounds are needed
        let d = ms.decompose(seg.as_bytes(), Some((66, 330)), None).unwrap();
        assert_eq!((d.start, d.end), (99, 297));
        let d = ms.decompose(seg.as_bytes(), None, Some(264)).unwrap();
        assert_eq!((d.start, d.end), (99, 297));
        let fwd = format!("{FLANK_L}{}{}{FLANK_R}", CEL.repeat(6), sub(CEL, &[5]));
        assert_eq!(span(&ms, &fwd), (66, 264, c(6)));
        assert_eq!(bounded(&ms, &fwd, (66, 297), None), (66, 297, c(7)));
        let bwd = format!("{FLANK_L}{}{}{FLANK_R}", sub(CEL, &[32]), CEL.repeat(6));
        assert_eq!(bounded(&ms, &bwd, (66, 297), None), (66, 297, c(7)));
    }

    #[test]
    fn identity_extension_does_not_enter_degenerate_flank_copies() {
        let ms = cel();
        let copy = sub(CEL, &[3, 10, 17, 24, 31]);
        let seg = format!(
            "{FLANK_L}{}{}{}{FLANK_R}",
            copy.repeat(2),
            CEL.repeat(10),
            copy.repeat(2)
        );
        assert_eq!(span(&ms, &seg), (131, 461, c(10)));
        assert_eq!(bounded(&ms, &seg, (132, 462), None), (132, 462, c(10)));
        assert_eq!(bounded(&ms, &seg, (0, seg.len()), None), (65, 527, c(14)));
    }

    // tests/test_decompose.py's test_a_seed_that_starts_in_the_flank_starts_at_the_tract
    #[test]
    fn a_seed_that_starts_in_the_flank_starts_at_the_tract() {
        let ms = MotifSet::new(3, &["CAG".to_string()], &[], &[], &[]).unwrap();
        let seg = format!("ACGTTAG{}TTGACAT", "CAG".repeat(10));
        let run = |w: Option<(usize, usize)>| {
            let d = ms.decompose(seg.as_bytes(), w, None).unwrap();
            (d.start, d.end, d.labels_str())
        };
        assert_eq!(run(None), (5, 35, "C".repeat(10)));
        assert_eq!(run(Some((7, 37))), (7, 37, "C".repeat(10)));
        // a window starting a whole unit or more after the seed leaves it
        assert_eq!(run(Some((9, 37))), (5, 35, "C".repeat(10)));
        let long = cel();
        let seg = format!("{}G{}{FLANK_R}", &FLANK_L[..65], CEL.repeat(6));
        let d = long.decompose(seg.as_bytes(), Some((66, 264)), None).unwrap();
        assert_eq!((d.start, d.end), (66, 264));
        // an error in the tract's last base: the seed's own, fuller run is kept
        let at = |seg: String| {
            let d = ms.decompose(seg.as_bytes(), Some((7, 37)), None).unwrap();
            (d.start, d.end, d.labels_str())
        };
        assert_eq!(
            at(format!("ACGTTCG{}CATTTGACAT", "CAG".repeat(9))),
            (6, 36, "C".repeat(10))
        );
        assert_eq!(
            at(format!("ACGTTCG{}TTGACAT", "CAG".repeat(10))),
            (7, 37, "C".repeat(10))
        );
        // but not where its frame would take a unit from the flank after the
        // tract, nor make one of the flank before it and a base inserted in it
        assert_eq!(
            at(format!("ACGTTAG{}CATTGACAT", "CAG".repeat(10))),
            (7, 37, "C".repeat(10))
        );
        let seg = format!("ACGTTAG{}C{}TTGACAT", "CAG".repeat(5), "CAG".repeat(5));
        let d = ms.decompose(seg.as_bytes(), Some((7, 38)), None).unwrap();
        assert_eq!((d.start, d.n_units()), (7, 10));
    }

    #[test]
    fn an_identity_unit_straddling_the_window_edge_is_refused() {
        let ms = cel();
        let noisy = format!("{FLANK_L}{}{}{FLANK_R}", CEL.repeat(6), sub(CEL, &[5]));
        assert_eq!(bounded(&ms, &noisy, (66, 296), Some(1000)), (66, 264, c(6)));
        let clean = format!("{FLANK_L}{}{FLANK_R}", CEL.repeat(7));
        assert_eq!(bounded(&ms, &clean, (66, 296), None), (66, 297, c(7)));
    }

    #[test]
    fn the_length_between_anchors_caps_identity_extension() {
        let ms = cel();
        let seg = format!("{FLANK_L}{}{}A{FLANK_R}", sub(CEL, &[0]), CEL.repeat(6));
        assert_eq!(span(&ms, &seg), (67, 265, c(6)));
        assert_eq!(bounded(&ms, &seg, (67, 298), Some(199)), (67, 265, c(6)));
        assert_eq!(bounded(&ms, &seg, (67, 298), Some(231)), (67, 298, c(7)));
    }

    #[test]
    fn after_an_identity_unit_exact_units_obey_the_bounds_too() {
        let ms = cel();
        let seg = format!("{FLANK_L}{}{}{CEL}{FLANK_R}", CEL.repeat(6), sub(CEL, &[5]).repeat(7));
        assert_eq!(span(&ms, &seg), (66, 264, c(6)));
        assert_eq!(bounded(&ms, &seg, (66, 495), None), (66, 495, c(13)));
        assert_eq!(bounded(&ms, &seg, (66, 528), Some(429)), (66, 495, c(13)));
        assert_eq!(bounded(&ms, &seg, (66, 528), None), (66, 528, c(14)));
    }

    #[test]
    fn the_bounds_change_nothing_for_short_motifs() {
        let ms = MotifSet::new(3, &["CAG".to_string()], &[], &[], &[]).unwrap();
        let mut tract: Vec<u8> = "CAG".repeat(30).into_bytes();
        tract[1] = b'T';
        let n = tract.len();
        tract[n - 2] = b'T';
        let seg = format!("ACGTTAGC{}TTGACAGT", String::from_utf8(tract).unwrap());
        assert_eq!(bounded(&ms, &seg, (8, 98), None), span(&ms, &seg));
    }

    const MUC1: &str = "CCGGGGCCGAGGTGACACCGTGGGCTGGGGGGGGCGGTGGAGCCCGGGGCCGGCCTGGTGT";
    const PLIN4: &str = concat!(
        "TGGTGTCCACGCCGGTCTGGATGGTTCCTTTGGCCACATTCATGGCACCAGTCACCCCACTACAGACGGTGTCC",
        "TTGGTACCTGTTAGGACAGTCTTAC"
    );

    fn del(unit: &str, p: usize) -> String {
        format!("{}{}", &unit[..p], &unit[p + 1..])
    }

    fn ins(unit: &str, p: usize) -> String {
        format!("{}A{}", &unit[..p], &unit[p..])
    }

    #[test]
    fn align_unit_counts_an_indel_as_one_edit() {
        let ms = cel();
        assert_eq!(ms.max_edits, 6);
        let cut = format!("{}C", del(CEL, 16));
        assert!(MotifSet::identity(cut.as_bytes(), CEL.as_bytes()) < ms.approx_identity);
        let at = |t: String| align_unit(CEL.as_bytes(), t.as_bytes(), 6);
        assert_eq!(at(format!("{CEL}{FLANK_R}")), Some((0, 33)));
        assert_eq!(at(format!("{}{FLANK_R}", del(CEL, 16))), Some((1, 32)));
        assert_eq!(at(format!("{}{FLANK_R}", ins(CEL, 16))), Some((1, 34)));
        assert_eq!(at(FLANK_R.to_string()), None);
        let a20 = "A".repeat(20);
        assert_eq!(MotifSet::new(20, &[a20], &[], &[], &[]).unwrap().max_edits, 4);
    }

    #[test]
    fn an_indel_in_an_edge_unit_is_taken_by_alignment_inside_the_tract() {
        let ms = cel();
        let seg = format!("{FLANK_L}{}{}{FLANK_R}", CEL.repeat(6), del(CEL, 16));
        assert_eq!(span(&ms, &seg), (66, 264, c(6)));
        assert_eq!(bounded(&ms, &seg, (66, 296), None), (66, 296, c(7)));
        let seg = format!("{FLANK_L}{}{}{FLANK_R}", sub(&del(CEL, 5), &[31]), CEL.repeat(6));
        assert_eq!(span(&ms, &seg), (98, 296, c(6)));
        assert_eq!(bounded(&ms, &seg, (66, 296), None), (66, 296, c(7)));
        assert_eq!(bounded(&ms, &seg, (67, 296), Some(230)), (98, 296, c(6)));
    }

    #[test]
    fn the_bounded_walk_crosses_indels_the_exact_walk_turns_into_noise() {
        let ms = cel();
        let mid: String = (0..6)
            .map(|k| if k % 2 == 0 { del(CEL, 16) } else { ins(CEL, 20) })
            .collect();
        let seg = format!("{FLANK_L}{}{mid}{}{FLANK_R}", CEL.repeat(2), CEL.repeat(2));
        assert_eq!(span(&ms, &seg), (66, 382, "CCooooCCC".to_string()));
        assert_eq!(bounded(&ms, &seg, (66, 396), None), (66, 396, c(10)));
    }

    #[test]
    fn alignment_does_not_enter_degenerate_flank_copies() {
        let ms = cel();
        let copy = sub(CEL, &[3, 10, 17, 24, 31]);
        let seg = format!(
            "{FLANK_L}{}{}{}{}{}{FLANK_R}",
            copy.repeat(2),
            del(CEL, 16),
            CEL.repeat(8),
            ins(CEL, 20),
            copy.repeat(2)
        );
        let tract = (132, 132 + 32 + 8 * 33 + 34);
        assert_eq!(span(&ms, &seg), (148, 445, c(9)));
        assert_eq!(bounded(&ms, &seg, tract, None), (148, 445, c(9)));
        assert_eq!(bounded(&ms, &seg, (0, seg.len()), None), (82, 512, c(13)));
    }

    #[test]
    fn a_very_long_motif_with_no_exact_copy_is_seeded_by_identity() {
        for motif in [MUC1, PLIN4] {
            let (l, h) = (motif.len(), motif.len() / 2);
            let mut units: Vec<String> = (0..6).map(|_| sub(motif, &[20, 20 + h])).collect();
            units[1] = del(&units[1], 40);
            units[4] = del(&units[4], 40);
            let tract: String = units.concat();
            let seg = format!("{FLANK_L}{tract}{FLANK_R}");
            let window = (66, 66 + tract.len());
            let ms = MotifSet::new(l, &[motif.to_string()], &[], &[], &[]).unwrap();
            assert!(ms.find_run(seg.as_bytes(), None).is_none());
            assert!(ms.decompose(seg.as_bytes(), None, None).is_none());
            assert_eq!(bounded(&ms, &seg, window, None), (66, window.1, c(6)));
            assert!(ms.decompose(seg.as_bytes(), Some(window), Some(2 * l - 1)).is_none());
        }
        let seg = format!("{FLANK_L}{}{FLANK_R}", sub(CEL, &[3, 19]).repeat(6));
        let ms = cel();
        assert!(ms.find_run(seg.as_bytes(), None).is_none());
        assert!(ms.decompose(seg.as_bytes(), Some((66, 264)), Some(198)).is_none());
    }
}
