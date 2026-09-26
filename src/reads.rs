//! Read extraction and classification at one locus.
//! Port of `expansionscout/reads.py`; the read classes are described there.

use std::collections::HashMap;
use std::io;

use crate::bam::{Aln, Bam, CDEL, CDIFF, CEQUAL, CHARD_CLIP, CINS, CMATCH, CREF_SKIP, CSOFT_CLIP};

pub const SPANNING: &str = "spanning";
pub const SPLIT: &str = "split";
pub const LEFT_PARTIAL: &str = "left_partial";
pub const RIGHT_PARTIAL: &str = "right_partial";
pub const NO_RUN: &str = "no_run";

/// Matched (query, reference) positions of a record, sorted by reference.
#[derive(Clone, Debug)]
pub struct PairIndex {
    q: Vec<i64>,
    r: Vec<i64>,
}

impl PairIndex {
    pub fn new(read: &Aln) -> PairIndex {
        let (q, r) = read.aligned_pairs();
        PairIndex { q, r }
    }

    /// `np.searchsorted(r, rpos, side="left")`
    fn search(&self, rpos: i64) -> usize {
        self.r.partition_point(|&v| v < rpos)
    }

    pub fn q_at_or_after(&self, rpos: i64) -> Option<i64> {
        let i = self.search(rpos);
        self.q.get(i).copied()
    }

    /// Exclusive query end: one past the last query base aligned to ref < rpos.
    pub fn q_before(&self, rpos: i64) -> Option<i64> {
        let i = self.search(rpos);
        if i == 0 {
            None
        } else {
            Some(self.q[i - 1] + 1)
        }
    }

    /// Query [qL, qR) covering reference [ra, rb), or None.
    pub fn window(&self, ra: i64, rb: i64) -> Option<(i64, i64)> {
        let (ql, qr) = (self.q_at_or_after(ra)?, self.q_before(rb)?);
        if qr <= ql {
            None
        } else {
            Some((ql, qr))
        }
    }
}

#[derive(Clone, Debug)]
pub struct ReadObs {
    pub name: String,
    pub read_class: &'static str,
    pub strand: &'static str,
    pub mapq: i64,
    pub hp: i64,
    pub read_length: i64,
    pub record: std::rc::Rc<Aln>,
    pub pairs: std::rc::Rc<PairIndex>,
    pub segment: Option<Vec<u8>>,
    pub seg_q_left: i64,
    pub size_bp_aln: Option<i64>,
    pub orig_anchor: Option<i64>,
    pub seg_tract: Option<(i64, i64)>,
}

pub struct LocusReads {
    pub obs: Vec<ReadObs>,
    pub read_lengths: Vec<i64>,
}

struct Clips {
    left: i64,
    right: i64,
    left_soft: bool,
    right_soft: bool,
    hard_left: i64,
}

fn clips(read: &Aln) -> Clips {
    let c = &read.cigar;
    let first = c.first();
    let last = c.last();
    let is_clip = |op: u8| op == CSOFT_CLIP || op == CHARD_CLIP;
    Clips {
        left: first.filter(|x| is_clip(x.0)).map(|x| x.1 as i64).unwrap_or(0),
        right: last.filter(|x| is_clip(x.0)).map(|x| x.1 as i64).unwrap_or(0),
        left_soft: first.is_some_and(|x| x.0 == CSOFT_CLIP),
        right_soft: last.is_some_and(|x| x.0 == CSOFT_CLIP),
        hard_left: first.filter(|x| x.0 == CHARD_CLIP).map(|x| x.1 as i64).unwrap_or(0),
    }
}

/// Where the annotated tract begins or ends inside a partial's segment.
fn tract_window(
    read: &Aln,
    pairs: &PairIndex,
    segment: Option<&Vec<u8>>,
    seg_q_left: i64,
    tract_start: Option<i64>,
    tract_end: i64,
) -> Result<Option<(i64, i64)>, String> {
    let Some(segment) = segment else { return Ok(None) };
    let len = segment.len() as i64;
    if let Some(ts) = tract_start {
        let q = match pairs.q_at_or_after(ts) {
            Some(q) => q,
            None => {
                let Some(re) = read.reference_end() else {
                    return Ok(None);
                };
                read.query_alignment_end()? + 0.max(ts - re)
            }
        };
        return Ok(Some((0.max(len.min(q - seg_q_left)), len)));
    }
    let q = match pairs.q_before(tract_end) {
        Some(q) => q,
        // pysam's reference_start is -1, never None, for an unplaced record
        None => read.query_alignment_start()? - 0.max(read.pos - tract_end),
    };
    Ok(Some((0, 0.max(len.min(q - seg_q_left)))))
}

/// (ref_pos, q0, q1) for each insertion: query [q0, q1), placed before
/// reference base `ref_pos`.
fn insertions(read: &Aln) -> Vec<(i64, i64, i64)> {
    let mut out = Vec::new();
    let (mut r, mut q) = (read.pos, 0i64);
    for &(op, n) in &read.cigar {
        let n = n as i64;
        match op {
            CINS => {
                out.push((r, q, q + n));
                q += n;
            }
            CMATCH | CEQUAL | CDIFF => {
                r += n;
                q += n;
            }
            CDEL | CREF_SKIP => r += n,
            CSOFT_CLIP => q += n,
            _ => {}
        }
    }
    out
}

/// Where a spanning read's alignment puts the tract, within its segment: from
/// the first query base after the last base aligned to the left flank to the
/// first query base aligned to the right flank, so an expansion written as an
/// insertion against either boundary is inside, and widened to take in an
/// insertion placed among the flank bases that continue the tract's period
/// (`period`, from `period_flanks`). See the Python's `_span_window`.
fn span_window(
    pairs: &PairIndex,
    read: &Aln,
    tract_start: i64,
    tract_end: i64,
    seg_q_left: i64,
    seg_len: i64,
    period: (i64, i64),
) -> Option<(i64, i64)> {
    let (mut q0, mut q1) = (pairs.q_before(tract_start)?, pairs.q_at_or_after(tract_end)?);
    let (before, after) = period;
    if before != 0 || after != 0 {
        for (ref_pos, i0, i1) in insertions(read) {
            if tract_start - before <= ref_pos && ref_pos < tract_start {
                q0 = q0.min(i0);
            } else if tract_end < ref_pos && ref_pos <= tract_end + after {
                q1 = q1.max(i1);
            }
        }
    }
    let (a, b) = (0.max(q0 - seg_q_left), seg_len.min(q1 - seg_q_left));
    (a < b).then_some((a, b))
}

/// Flank bases beside the tract [start, end) that continue its period:
/// (left, right), each counted outward from the tract and at most `reach`.
/// `reference` is upper-case sequence beginning at `offset`; a flank base
/// continues the period when it equals the base one unit further in and is
/// A, C, G or T. See the Python's `period_flanks`.
pub fn period_flanks(reference: &[u8], offset: i64, start: i64, end: i64, unit_len: i64, reach: i64) -> (i64, i64) {
    let at = |x: i64| -> u8 {
        let i = x - offset;
        if i >= 0 && (i as usize) < reference.len() {
            reference[i as usize]
        } else {
            b'N'
        }
    };
    let continues = |x: i64, y: i64| matches!(at(x), b'A' | b'C' | b'G' | b'T') && at(x) == at(y);
    let mut left = 0;
    while left < reach && continues(start - 1 - left, start - 1 - left + unit_len) {
        left += 1;
    }
    let mut right = 0;
    while right < reach && continues(end + right, end + right - unit_len) {
        right += 1;
    }
    (left, right)
}

/// Reference [lo, hi) a spanning read's size is measured between: the
/// anchors, unless the whole anchor continues the tract's period, where an
/// aligner may place a tract's indel beyond the anchor. There the bound
/// moves out to the first flank base that breaks the period, but not past
/// `flank`. See the Python's `_size_bounds`.
fn size_bounds(start: i64, end: i64, anchor: i64, flank: i64, period: (i64, i64)) -> (i64, i64) {
    let (left, right) = period;
    let (mut lo, mut hi) = (start - anchor, end + anchor);
    if 0 < anchor && anchor <= left {
        lo = lo.min((start - flank).max(start - left - 1));
    }
    if 0 < anchor && anchor <= right {
        hi = hi.max((end + flank).min(end + right + 1));
    }
    (lo, hi)
}

/// `s[a:b]` for non-negative bounds: clamped, never out of range.
fn py_slice(s: &[u8], a: i64, b: i64) -> &[u8] {
    let n = s.len() as i64;
    let a = a.clamp(0, n) as usize;
    let b = b.clamp(0, n) as usize;
    if b <= a {
        &[]
    } else {
        &s[a..b]
    }
}

fn orig_coord(read: &Aln, q: i64, hard_left: i64, read_length: i64) -> i64 {
    let p = q + hard_left;
    if read.is_reverse() {
        read_length - 1 - p
    } else {
        p
    }
}

/// Map a catalogue contig name onto the BAM's own naming, or None.
pub fn resolve_contig(references: &[String], chrom: &str) -> Option<String> {
    let has = |c: &str| references.iter().any(|r| r == c);
    if has(chrom) {
        return Some(chrom.to_string());
    }
    let alt = match chrom.strip_prefix("chr") {
        Some(rest) => rest.to_string(),
        None => format!("chr{chrom}"),
    };
    if has(&alt) {
        return Some(alt);
    }
    for (a, b) in [("chrM", "MT"), ("MT", "chrM")] {
        if chrom == a && has(b) {
            return Some(b.to_string());
        }
    }
    None
}

/// `period` is (left, right) from `period_flanks`: flank bases beside the
/// tract that continue its period. They decide where a spanning read's size
/// is measured from (`size_bounds`) and, with `widen_window` (a long motif,
/// the only kind decomposed within the tract window), let its tract window
/// take in an insertion an aligner placed among the first `margin` of them
/// (`span_window`). (0, 0) changes neither.
#[allow(clippy::too_many_arguments)]
pub fn collect_reads(
    bam: &mut Bam,
    chrom: &str,
    start: i64,
    end: i64,
    anchor: i64,
    flank: i64,
    min_mapq: i64,
    min_clip: i64,
    margin: i64,
    period: (i64, i64),
    widen_window: bool,
) -> Result<LocusReads, String> {
    let (size_lo, size_hi) = size_bounds(start, end, anchor, flank, period);
    let win_period = if widen_window {
        (period.0.min(margin), period.1.min(margin))
    } else {
        (0, 0)
    };
    let mut obs = Vec::new();
    // partials by read name, in order of first appearance
    let mut partials: Vec<(String, Vec<ReadObs>)> = Vec::new();
    let mut partial_idx: HashMap<String, usize> = HashMap::new();
    let mut lengths: Vec<i64> = Vec::new();
    let mut seen_len: HashMap<String, ()> = HashMap::new();
    let mut seen_spanning: HashMap<String, ()> = HashMap::new();
    let records = bam
        .fetch(chrom, 0.max(start - flank - anchor), end + flank + anchor)
        .map_err(|e: io::Error| format!("{chrom}:{start}-{end}: {e}"))?;
    for read in records {
        if read.is_unmapped() || read.is_secondary() {
            continue;
        }
        if (read.mapq as i64) < min_mapq {
            continue;
        }
        let read_length = read.infer_read_length().unwrap_or(0);
        if seen_len.insert(read.name.clone(), ()).is_none() {
            lengths.push(read_length);
        }
        let ref_start = read.pos;
        let ref_end = read.reference_end();
        let covers_left = ref_start <= start - flank;
        let covers_right = ref_end.is_some_and(|e| e >= end + flank);
        let c = clips(&read);
        let strand = if read.is_reverse() { "-" } else { "+" };
        let name = read.name.clone();
        let mapq = read.mapq as i64;
        let hp = read.hp();
        let read = std::rc::Rc::new(read);
        if covers_left && covers_right {
            if seen_spanning.contains_key(&name) {
                continue;
            }
            let pairs = PairIndex::new(&read);
            let win = pairs.window(start - anchor, end + anchor);
            let (Some(win), Some(seq)) = (win, read.query_sequence()) else {
                continue;
            };
            // never None where `win` is not: the bounds are at least as wide
            let (ql, qr) = pairs.window(size_lo, size_hi).unwrap_or(win);
            let (dl, dr) = pairs.window(start - margin, end + margin).unwrap_or(win);
            let segment = py_slice(seq, dl, dr).to_vec();
            let seg_tract = span_window(&pairs, &read, start, end, dl, segment.len() as i64, win_period);
            seen_spanning.insert(name.clone(), ());
            obs.push(ReadObs {
                name,
                read_class: SPANNING,
                strand,
                mapq,
                hp,
                read_length,
                pairs: std::rc::Rc::new(pairs),
                segment: Some(segment),
                seg_q_left: dl,
                size_bp_aln: Some(0.max((qr - ql) - (start - size_lo) - (size_hi - end))),
                orig_anchor: None,
                seg_tract,
                record: read,
            });
        } else if covers_left && c.right >= min_clip {
            let pairs = PairIndex::new(&read);
            let Some(ql) = pairs.q_at_or_after(start - anchor) else {
                continue;
            };
            let seg = if c.right_soft {
                read.query_sequence().map(|s| py_slice(s, ql, s.len() as i64).to_vec())
            } else {
                None
            };
            let seg_tract = tract_window(&read, &pairs, seg.as_ref(), ql, Some(start), 0)?;
            let o = ReadObs {
                name: name.clone(),
                read_class: LEFT_PARTIAL,
                strand,
                mapq,
                hp,
                read_length,
                segment: seg,
                seg_q_left: ql,
                size_bp_aln: None,
                orig_anchor: Some(orig_coord(&read, ql, c.hard_left, read_length)),
                seg_tract,
                pairs: std::rc::Rc::new(pairs),
                record: read,
            };
            push_partial(&mut partials, &mut partial_idx, name, o);
        } else if covers_right && c.left >= min_clip {
            let pairs = PairIndex::new(&read);
            let Some(qr) = pairs.q_before(end + anchor) else {
                continue;
            };
            let seg = if c.left_soft {
                read.query_sequence().map(|s| py_slice(s, 0, qr).to_vec())
            } else {
                None
            };
            let seg_tract = tract_window(&read, &pairs, seg.as_ref(), 0, None, end)?;
            let o = ReadObs {
                name: name.clone(),
                read_class: RIGHT_PARTIAL,
                strand,
                mapq,
                hp,
                read_length,
                segment: seg,
                seg_q_left: 0,
                size_bp_aln: None,
                orig_anchor: Some(orig_coord(&read, qr, c.hard_left, read_length)),
                seg_tract,
                pairs: std::rc::Rc::new(pairs),
                record: read,
            };
            push_partial(&mut partials, &mut partial_idx, name, o);
        }
    }

    // pair partial records of the same read into split observations
    for (name, recs) in partials {
        if seen_spanning.contains_key(&name) {
            continue;
        }
        let lefts: Vec<&ReadObs> = recs.iter().filter(|r| r.read_class == LEFT_PARTIAL).collect();
        let rights: Vec<&ReadObs> = recs.iter().filter(|r| r.read_class == RIGHT_PARTIAL).collect();
        let mut paired = false;
        'outer: for l in &lefts {
            for r in &rights {
                if l.strand != r.strand {
                    continue;
                }
                let (la, ra) = (l.orig_anchor.unwrap(), r.orig_anchor.unwrap());
                let gap = if l.strand == "+" { ra - la } else { la - ra };
                if gap <= 0 {
                    continue;
                }
                let carrier = if l.segment.is_some() { l } else { r };
                obs.push(ReadObs {
                    name: name.clone(),
                    read_class: SPLIT,
                    strand: l.strand,
                    mapq: l.mapq.min(r.mapq),
                    hp: if l.hp != 0 { l.hp } else { r.hp },
                    read_length: l.read_length,
                    record: carrier.record.clone(),
                    pairs: carrier.pairs.clone(),
                    segment: carrier.segment.clone(),
                    seg_q_left: carrier.seg_q_left,
                    size_bp_aln: Some(0.max(gap - 2 * anchor)),
                    orig_anchor: None,
                    seg_tract: carrier.seg_tract,
                });
                paired = true;
                break 'outer;
            }
        }
        if !paired {
            // keep the best-supported single partial per read: max() keeps the first
            let mut best = &recs[0];
            for o in &recs[1..] {
                if (o.segment.is_some(), o.mapq) > (best.segment.is_some(), best.mapq) {
                    best = o;
                }
            }
            obs.push(best.clone());
        }
    }
    Ok(LocusReads {
        obs,
        read_lengths: lengths,
    })
}

fn push_partial(
    partials: &mut Vec<(String, Vec<ReadObs>)>,
    idx: &mut HashMap<String, usize>,
    name: String,
    o: ReadObs,
) {
    match idx.get(&name) {
        Some(&i) => partials[i].1.push(o),
        None => {
            idx.insert(name.clone(), partials.len());
            partials.push((name, vec![o]));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // tests/test_engine.py's test_size_bounds_move_only_where_the_whole_anchor_continues_the_period
    #[test]
    fn size_bounds_move_only_where_the_whole_anchor_continues_the_period() {
        let r = format!("{}TCGATC{}GAT{}", "C".repeat(14), "GATC".repeat(5), "G".repeat(17));
        let r = r.as_bytes();
        assert_eq!(period_flanks(r, 0, 20, 40, 4, 10), (6, 3));
        assert_eq!(period_flanks(&r[5..], 5, 20, 40, 4, 10), (6, 3));
        assert_eq!(period_flanks(r, 0, 20, 40, 4, 4), (4, 3));
        let mut n = r.to_vec();
        n[15] = b'N';
        assert_eq!(period_flanks(&n, 0, 20, 40, 4, 10), (4, 3));
        assert_eq!(period_flanks(&r[..42], 0, 20, 40, 4, 10), (6, 2));
        assert_eq!(size_bounds(20, 40, 3, 10, (6, 3)), (13, 44));
        assert_eq!(size_bounds(20, 40, 5, 10, (6, 3)), (13, 45));
        assert_eq!(size_bounds(20, 40, 3, 5, (6, 3)), (15, 44));
        assert_eq!(size_bounds(20, 40, 8, 10, (6, 3)), (12, 48));
        assert_eq!(size_bounds(20, 40, 3, 10, (0, 0)), (17, 43));
        assert_eq!(size_bounds(20, 40, 0, 10, (6, 3)), (20, 40));
    }
}
