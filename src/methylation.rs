//! Per-read 5mC from MM/ML tags, in query coordinates.
//!
//! Port of `expansionscout/methylation.py`, which explains the modkit
//! conventions this follows. The MM/ML parser is a port of htslib 1.24's
//! (`bam_parse_basemod2`, `bam_next_basemod`, `bam_mods_at_next_pos`) as
//! pysam's `modified_bases` drives it, so it returns what pysam returns:
//! nothing for a malformed or inconsistent tag, the sites found so far where
//! htslib stops part-way along a read, and htslib's messages about either.

use crate::bam::{Aln, Aux};
use crate::npstat::{fsum, py_round};

pub const DEFAULT_MOD_THRESHOLD: f64 = 0.8;
pub const CONVENTIONS: [&str; 3] = ["modkit", "combine", "ignore_h"];
pub const DEFAULT_CONVENTION: &str = "modkit";

/// One (canonical base, strand, code) key of pysam's `modified_bases`, with
/// its (position, ML byte or -1) list.
struct ModKey {
    base: u8,
    strand: u8,
    code: String,
    sites: Vec<(i64, i64)>,
}

const INT_MAX: i64 = i32::MAX as i64;
const MAX_BASE_MOD: usize = 256;
const SEQI_RC: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];
const SEQ_NT16_STR: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

fn nt16(c: u8) -> u8 {
    match c {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 15,
    }
}

/// htslib's `hts_str2uint(in, &end, 31, &failed)` on a NUL-terminated buffer.
fn str2uint(s: &[u8], at: usize) -> (i64, usize, bool) {
    let limit: u64 = (1u64 << 31) - 1;
    let mut v = at;
    if s[v] == b'+' {
        v += 1;
    }
    let mut n: u64 = 0;
    let mut fast = 31 * 1000 / 3322 + 1;
    loop {
        fast -= 1;
        if fast == 0 || !s[v].is_ascii_digit() {
            break;
        }
        n = n * 10 + (s[v] - b'0') as u64;
        v += 1;
    }
    let mut failed = false;
    if s[v].is_ascii_digit() && fast == 0 {
        let (d10, m10) = (limit / 10, limit - 10 * (limit / 10));
        while s[v].is_ascii_digit() {
            let d = (s[v] - b'0') as u64;
            if n < d10 || (n == d10 && d <= m10) {
                n = n * 10 + d;
                v += 1;
            } else {
                while s[v].is_ascii_digit() {
                    v += 1;
                }
                n = limit;
                failed = true;
                break;
            }
        }
    }
    (n as i64, v, failed)
}

/// `hts_base_mod_state`, with pointers held as indices.
struct ModState {
    typ: Vec<i64>,
    strand: Vec<u8>,
    canonical: Vec<u8>,
    count: Vec<i64>,
    mm: Vec<usize>,
    mm_end: Vec<usize>,
    ml: Vec<Option<i64>>,
    stride: Vec<i64>,
    seq_pos: usize,
}

/// `bam_parse_basemod2`, with pysam's flags (0). Err after logging, as htslib does.
fn parse_basemod(read: &Aln, mm: &[u8], ml: Option<&[u8]>, diag: &mut Vec<String>) -> Result<ModState, ()> {
    let q = &read.name;
    let rev = read.is_reverse();
    let mut err = |m: &str| diag.push(format!("[E::bam_parse_basemod2] {q}: {m}"));
    let mut st = ModState {
        typ: vec![],
        strand: vec![],
        canonical: vec![],
        count: vec![],
        mm: vec![],
        mm_end: vec![],
        ml: vec![],
        stride: vec![],
        seq_pos: 0,
    };
    let mut freq = [0i64; 16];
    if rev {
        for &c in &read.seq4 {
            freq[c as usize] += 1;
        }
        freq[15] = read.seq4.len() as i64;
    }
    let ml_end = ml.map(|m| m.len() as i64);
    let mut mlp: Option<i64> = ml.map(|_| 0);
    let mut cp = 0usize;
    while mm[cp] != 0 {
        while mm[cp] != 0 {
            let btype = mm[cp];
            cp += 1;
            if !matches!(btype, b'A' | b'C' | b'G' | b'T' | b'U' | b'N') {
                return Err(());
            }
            let btype = nt16(if btype == b'U' { b'T' } else { btype });
            if mm[cp] != b'+' && mm[cp] != b'-' {
                return Err(());
            }
            let strand = mm[cp];
            cp += 1;
            let mut ms = cp;
            let mut chebi = 0i64;
            if mm[cp].is_ascii_digit() {
                let (v, end, failed) = str2uint(mm, cp);
                if end == cp || failed {
                    err("malformed MM tag (invalid ChEBI code)");
                    return Err(());
                }
                chebi = v;
                cp = end;
                ms = cp - 1;
            } else {
                while mm[cp] != 0 && mm[cp].is_ascii_alphabetic() {
                    cp += 1;
                }
                if mm[cp] == 0 {
                    return Err(());
                }
            }
            let me = cp;
            if mm[cp] == b'.' || mm[cp] == b'?' {
                cp += 1;
            } else if mm[cp] != b',' && mm[cp] != b';' {
                return Err(());
            }
            let stride = (me - ms) as i64;
            let mut ndelta = 0i64;
            let delta: i64;
            let cp_end: usize;
            if rev {
                let mut total = 0i64;
                loop {
                    if mm[cp] == b',' {
                        cp += 1;
                    }
                    if mm[cp] == 0 || mm[cp] == b';' {
                        break;
                    }
                    let (d, end, failed) = str2uint(mm, cp);
                    if end == cp || failed {
                        err("Hit end of MM tag. Missing semicolon?");
                        return Err(());
                    }
                    cp = end;
                    total += d + 1;
                    ndelta += 1;
                }
                delta = freq[SEQI_RC[btype as usize] as usize] - total;
                cp_end = cp;
            } else if mm[cp] == b',' {
                let (d, end, failed) = str2uint(mm, cp + 1);
                if end == cp + 1 || failed {
                    err("Failed to parse integer from MM tag");
                    return Err(());
                }
                delta = d;
                cp_end = end;
            } else {
                delta = INT_MAX;
                cp_end = cp;
            }
            let mut n = 0i64;
            while ms < me {
                st.typ.push(if chebi != 0 { -chebi } else { mm[ms] as i64 });
                st.strand.push((strand == b'-') as u8);
                st.canonical.push(btype);
                st.stride.push(stride);
                if delta < 0 {
                    err("MM tag refers to bases beyond sequence length");
                    return Err(());
                }
                st.count.push(delta);
                if rev {
                    st.mm.push(me + 1);
                    st.mm_end.push(cp_end);
                    st.ml.push(mlp.map(|m| m + n + (ndelta - 1) * stride));
                } else {
                    st.mm.push(cp_end);
                    st.mm_end.push(0);
                    st.ml.push(mlp.map(|m| m + n));
                }
                if st.typ.len() >= MAX_BASE_MOD {
                    err("Too many base modification types");
                    return Err(());
                }
                ms += 1;
                n += 1;
            }
            if let Some(m) = mlp.as_mut() {
                if rev {
                    *m += ndelta * stride;
                } else {
                    while mm[cp] != 0 && mm[cp] != b';' {
                        if mm[cp] == b',' {
                            *m += stride;
                        }
                        cp += 1;
                    }
                }
                if *m > ml_end.unwrap() {
                    err("Insufficient number of entries in ML tag");
                    return Err(());
                }
            } else if rev {
                cp = cp_end;
            } else {
                while mm[cp] != 0 && mm[cp] != b';' {
                    cp += 1;
                }
            }
            if mm[cp] == 0 {
                err("Hit end of MM tag. Missing semicolon?");
                return Err(());
            }
            cp += 1;
        }
    }
    if let (Some(m), Some(e)) = (mlp, ml_end) {
        if m != e {
            err("Too many entries in ML tag");
            return Err(());
        }
    }
    Ok(st)
}

/// One hit reported by `bam_mods_at_next_pos`: (type, canonical, strand, qual).
type Hit = (i64, u8, u8, i64);

/// `bam_mods_at_next_pos`. Err after logging.
fn mods_at_next_pos(
    read: &Aln,
    st: &mut ModState,
    mm: &[u8],
    ml: Option<&[u8]>,
    out: &mut Vec<Hit>,
    diag: &mut Vec<String>,
) -> Result<(), ()> {
    let l_qseq = read.seq4.len();
    if st.seq_pos >= l_qseq {
        return Ok(());
    }
    let rev = read.is_reverse();
    let mut base = read.seq4[st.seq_pos];
    st.seq_pos += 1;
    if rev {
        base = SEQI_RC[base as usize];
    }
    let qual = |p: Option<i64>| -> i64 {
        match (p, ml) {
            (Some(i), Some(m)) => m[i as usize] as i64,
            _ => -1,
        }
    };
    let nmods = st.typ.len();
    let mut i = 0;
    while i < nmods {
        if st.canonical[i] != base && st.canonical[i] != 15 {
            i += 1;
            continue;
        }
        let old = st.count[i];
        st.count[i] -= 1;
        if old > 0 {
            i += 1;
            continue;
        }
        let mmptr = st.mm[i];
        out.push((st.typ[i], st.canonical[i], st.strand[i], qual(st.ml[i])));
        if let Some(p) = st.ml[i].as_mut() {
            *p += if rev { -st.stride[i] } else { st.stride[i] };
        }
        let mut failed = false;
        if rev {
            if st.mm_end[i] == 0 || st.mm_end[i] - 1 < st.mm[i] {
                diag.push("[E::bam_mods_at_next_pos] Assert failed while processing base modification states".into());
                return Err(());
            }
            let mut cp = st.mm_end[i] - 1;
            while cp != st.mm[i] {
                if mm[cp] == b',' {
                    break;
                }
                cp -= 1;
            }
            st.mm_end[i] = cp;
            if cp != st.mm[i] {
                let (v, _, f) = str2uint(mm, cp + 1);
                st.count[i] = v;
                failed = f;
            } else {
                st.count[i] = INT_MAX;
            }
        } else if mm[st.mm[i]] == b',' {
            let (v, end, f) = str2uint(mm, st.mm[i] + 1);
            st.count[i] = v;
            st.mm[i] = end;
            failed = f;
        } else {
            st.count[i] = INT_MAX;
        }
        if failed {
            diag.push(format!(
                "[E::bam_mods_at_next_pos] {}: Error parsing unsigned integer from MM tag",
                read.name
            ));
            return Err(());
        }
        let mut j = i + 1;
        while j < nmods && st.mm[j] == mmptr {
            out.push((st.typ[j], st.canonical[j], st.strand[j], qual(st.ml[j])));
            st.count[j] = st.count[i];
            st.mm[j] = st.mm[i];
            if let Some(p) = st.ml[j].as_mut() {
                *p += if rev { -st.stride[j] } else { st.stride[j] };
            }
            j += 1;
        }
        i = j;
    }
    Ok(())
}

/// `bam_next_basemod`: Ok(hits, position), Ok(empty) at the end, Err on error.
fn next_basemod(
    read: &Aln,
    st: &mut ModState,
    mm: &[u8],
    ml: Option<&[u8]>,
    diag: &mut Vec<String>,
) -> Result<(Vec<Hit>, usize), ()> {
    let rev = read.is_reverse();
    let l_qseq = read.seq4.len();
    let mut next = [0x7f7f_7f7fi64; 16];
    let mut freq = [0i64; 16];
    for i in 0..st.typ.len() {
        let mut base = st.canonical[i] as usize;
        if rev {
            base = SEQI_RC[base] as usize;
        }
        if next[base] > st.count[i] {
            next[base] = st.count[i];
        }
    }
    let mut i = st.seq_pos;
    while i < l_qseq {
        let bc = read.seq4[i] as usize;
        if next[bc] <= freq[bc] || next[15] <= freq[15] {
            break;
        }
        freq[bc] += 1;
        if bc != 15 {
            freq[15] += 1;
        }
        i += 1;
    }
    st.seq_pos = i;
    let pos = i;
    for k in 0..st.typ.len() {
        let c = st.canonical[k] as usize;
        st.count[k] -= freq[if rev { SEQI_RC[c] as usize } else { c }];
    }
    if l_qseq != 0 && st.seq_pos >= l_qseq {
        if !rev {
            for k in 0..st.typ.len() {
                if st.count[k] < 0x7f00_0000 || (mm[st.mm[k]] != 0 && mm[st.mm[k]] != b';') {
                    diag.push("[W::bam_next_basemod] MM tag refers to bases beyond sequence length".into());
                    return Err(());
                }
            }
        }
        return Ok((vec![], pos));
    }
    let mut hits = Vec::new();
    // An error here ends the iteration: htslib turns it into "no more", so
    // pysam keeps the positions before this one and drops this one's hits,
    // which were recorded before the error and are discarded with it.
    if mods_at_next_pos(read, st, mm, ml, &mut hits, diag).is_err() {
        return Ok((vec![], pos));
    }
    Ok((hits, pos))
}

/// pysam's `modified_bases`: None where it returns None, otherwise the keys
/// in first-emission order.
fn modified_bases(read: &Aln, diag: &mut Vec<String>) -> Option<Vec<ModKey>> {
    let mm_tag = read.tag(b"MM").or_else(|| read.tag(b"Mm"));
    let Some((mm_type, mm_val)) = mm_tag else {
        return Some(vec![]);
    };
    let Aux::Str(mm_str) = (if mm_type == b'Z' { mm_val } else { Aux::Array(0) }) else {
        diag.push(format!(
            "[E::bam_parse_basemod2] {}: MM tag is not of type Z",
            read.name
        ));
        return None;
    };
    let l_qseq = read.seq4.len() as i64;
    if let Some((_, v)) = read.tag(b"MN") {
        let mn = match v {
            Aux::Int(i) => i,
            _ => 0,
        };
        if mn != l_qseq && l_qseq != 0 {
            diag.push(format!(
                "[E::bam_parse_basemod2] {}: MM/MN data length is incompatible with SEQ length",
                read.name
            ));
            return None;
        }
    }
    let ml: Option<Vec<u8>> = match read.tag(b"ML").or_else(|| read.tag(b"Ml")) {
        None => None,
        Some((b'B', Aux::ArrayU8(v))) => Some(v.to_vec()),
        Some(_) => {
            diag.push(format!(
                "[E::bam_parse_basemod2] {}: ML tag is not of type B,C",
                read.name
            ));
            return None;
        }
    };
    let mut mm: Vec<u8> = mm_str.to_vec();
    mm.push(0);
    let mut st = parse_basemod(read, &mm, ml.as_deref(), diag).ok()?;
    let rev = read.is_reverse();
    let mut keys: Vec<ModKey> = Vec::new();
    loop {
        let (hits, pos) = next_basemod(read, &mut st, &mm, ml.as_deref(), diag).ok()?;
        if hits.is_empty() {
            break;
        }
        for (typ, canonical, strand, qual) in hits {
            let code = if typ > 0 {
                (typ as u8 as char).to_string()
            } else {
                (-typ).to_string()
            };
            let base = SEQ_NT16_STR[canonical as usize];
            let strand = if rev { 1 - strand } else { strand };
            match keys
                .iter_mut()
                .find(|k| k.base == base && k.strand == strand && k.code == code)
            {
                Some(k) => k.sites.push((pos as i64, qual)),
                None => keys.push(ModKey {
                    base,
                    strand,
                    code,
                    sites: vec![(pos as i64, qual)],
                }),
            }
        }
    }
    Some(keys)
}

/// (qpos, p_m, p_h) for every C with a modification call, by position.
pub fn iter_mods(read: &Aln, diag: &mut Vec<String>) -> Vec<(i64, f64, f64)> {
    let Some(keys) = modified_bases(read, diag) else {
        return vec![];
    };
    let mut m: Vec<(i64, f64)> = Vec::new();
    let mut h: Vec<(i64, f64)> = Vec::new();
    fn set(t: &mut Vec<(i64, f64)>, q: i64, v: f64) {
        match t.iter_mut().find(|e| e.0 == q) {
            Some(e) => e.1 = v,
            None => t.push((q, v)),
        }
    }
    for k in &keys {
        if k.base != b'C' {
            continue;
        }
        let target = match k.code.as_str() {
            "m" | "27551" => &mut m,
            "h" | "76792" => &mut h,
            _ => continue,
        };
        for &(q, prob) in &k.sites {
            set(target, q, (prob as f64 + 0.5) / 256.0);
        }
    }
    let mut qs: Vec<i64> = m.iter().chain(h.iter()).map(|e| e.0).collect();
    qs.sort();
    qs.dedup();
    let get = |t: &Vec<(i64, f64)>, q: i64| t.iter().find(|e| e.0 == q).map(|e| e.1).unwrap_or(0.0);
    qs.into_iter().map(|q| (q, get(&m, q), get(&h, q))).collect()
}

#[derive(Clone, Debug, Default)]
pub struct WindowMeth {
    pub mean: Option<f64>,
    pub n: usize,
    pub frac: Option<f64>,
    pub n_valid: usize,
    pub n_mod: usize,
    pub n_fail: usize,
    pub n_hydroxy: usize,
}

impl WindowMeth {
    pub fn fail_frac(&self) -> Option<f64> {
        if self.n > 0 {
            Some(self.n_fail as f64 / self.n as f64)
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct ReadMeth {
    pub tract: WindowMeth,
    pub up: WindowMeth,
    pub down: WindowMeth,
    pub per_site: Vec<(i64, i64)>,
}

/// One call -> 'C', 'm', 'h', or None for a call that fails the threshold.
fn classify(pm: f64, ph: f64, threshold: f64, convention: &str) -> Option<u8> {
    let pmax = |a: f64, b: f64| if b > a { b } else { a };
    if convention == "ignore_h" {
        let pm2 = pm + ph / 2.0;
        let pc = pmax(0.0, 1.0 - pm - ph) + ph / 2.0;
        let (best, label) = if pm2 >= pc { (pm2, b'm') } else { (pc, b'C') };
        return if best >= threshold { Some(label) } else { None };
    }
    let pc = pmax(0.0, 1.0 - pm - ph);
    // max() over (p, label) tuples: the highest p, and on a tie the label
    // that sorts last ('m' > 'h' > 'C').
    let mut best = (pc, b'C');
    for cand in [(pm, b'm'), (ph, b'h')] {
        if cand.0 > best.0 || (cand.0 == best.0 && cand.1 > best.1) {
            best = cand;
        }
    }
    if best.0 >= threshold {
        Some(best.1)
    } else {
        None
    }
}

fn window(calls: &[(i64, f64, f64)], q_left: i64, q_right: i64, threshold: f64, convention: &str) -> WindowMeth {
    let sel: Vec<(f64, f64)> = calls
        .iter()
        .filter(|c| q_left <= c.0 && c.0 < q_right)
        .map(|c| (c.1, c.2))
        .collect();
    if sel.is_empty() {
        return WindowMeth::default();
    }
    let kept: Vec<u8> = sel
        .iter()
        .filter_map(|&(pm, ph)| classify(pm, ph, threshold, convention))
        .collect();
    let n_mod = kept.iter().filter(|&&l| l == b'm').count();
    let n_hyd = kept.iter().filter(|&&l| l == b'h').count();
    let num = if convention == "combine" { n_mod + n_hyd } else { n_mod };
    WindowMeth {
        mean: Some(fsum(sel.iter().map(|s| s.0)) / sel.len() as f64),
        n: sel.len(),
        frac: if kept.is_empty() {
            None
        } else {
            Some(num as f64 / kept.len() as f64)
        },
        n_valid: kept.len(),
        n_mod,
        n_fail: sel.len() - kept.len(),
        n_hydroxy: n_hyd,
    }
}

/// Summarise 5mC over query windows.
pub fn read_methylation(
    read: &Aln,
    tract_q: Option<(i64, i64)>,
    up_q: Option<(i64, i64)>,
    down_q: Option<(i64, i64)>,
    threshold: f64,
    convention: &str,
    diag: &mut Vec<String>,
) -> ReadMeth {
    let calls = iter_mods(read, diag);
    let mut out = ReadMeth::default();
    if calls.is_empty() {
        return out;
    }
    if let Some((a, b)) = tract_q {
        out.tract = window(&calls, a, b, threshold, convention);
        let mut ps: Vec<(i64, i64)> = calls
            .iter()
            .filter(|c| a <= c.0 && c.0 < b)
            .map(|c| (c.0 - a, py_round(c.1 * 100.0) as i64))
            .collect();
        ps.sort();
        out.per_site = ps;
    }
    if let Some((a, b)) = up_q {
        out.up = window(&calls, a, b, threshold, convention);
    }
    if let Some((a, b)) = down_q {
        out.down = window(&calls, a, b, threshold, convention);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A read with SEQ `seq` (4-bit codes derived from it) and the given tags.
    fn read(seq: &str, reverse: bool, mm: &str, ml: Option<&[u8]>, mn: Option<i64>) -> Aln {
        let code = |b: u8| match b {
            b'A' => 1,
            b'C' => 2,
            b'M' => 3,
            b'G' => 4,
            b'T' => 8,
            _ => 15,
        };
        let mut aux = Vec::new();
        aux.extend_from_slice(b"MMZ");
        aux.extend_from_slice(mm.as_bytes());
        aux.push(0);
        if let Some(ml) = ml {
            aux.extend_from_slice(b"MLBC");
            aux.extend_from_slice(&(ml.len() as u32).to_le_bytes());
            aux.extend_from_slice(ml);
        }
        if let Some(mn) = mn {
            aux.extend_from_slice(b"MNi");
            aux.extend_from_slice(&(mn as i32).to_le_bytes());
        }
        Aln {
            name: "r".into(),
            flag: if reverse { 16 } else { 0 },
            pos: 0,
            mapq: 60,
            cigar: vec![(0, seq.len() as u32)],
            seq: seq.as_bytes().to_vec(),
            seq4: seq.bytes().map(code).collect(),
            aux,
        }
    }

    /// Keys written base-strand-code ("C0m"), each with its (pos, qual) list.
    type Mods = Vec<(String, Vec<(i64, i64)>)>;

    /// `modified_bases` as [(key, sites)], in key order.
    fn mods(r: &Aln) -> Option<Mods> {
        let mut diag = Vec::new();
        modified_bases(r, &mut diag).map(|keys| {
            keys.into_iter()
                .map(|k| (format!("{}{}{}", k.base as char, k.strand, k.code), k.sites))
                .collect()
        })
    }

    // Expected values are pysam 0.24.1's own output for the same records.
    const SEQ: &str = "ACGCGTCCA";

    #[test]
    fn forward_and_reverse_skip_counts() {
        let key = |s: &str, v: Vec<(i64, i64)>| Some(vec![(s.to_string(), v)]);
        assert_eq!(
            mods(&read(SEQ, false, "C+m,0,1;", Some(&[10, 200]), None)),
            key("C0m", vec![(1, 10), (6, 200)])
        );
        // on a reverse read the counts run over Gs from the end, and there are only two
        assert_eq!(mods(&read(SEQ, true, "C+m,0,1;", Some(&[10, 200]), None)), None);
        assert_eq!(
            mods(&read(SEQ, true, "C+m,0,0;", Some(&[10, 200]), None)),
            key("C1m", vec![(2, 200), (4, 10)])
        );
        assert_eq!(
            mods(&read(SEQ, false, "C+m?,0,1;", Some(&[10, 200]), None)),
            key("C0m", vec![(1, 10), (6, 200)])
        );
    }

    #[test]
    fn multiple_codes_and_groups() {
        assert_eq!(
            mods(&read(SEQ, false, "C+mh,0,1;", Some(&[1, 2, 3, 4]), None)),
            Some(vec![
                ("C0m".into(), vec![(1, 1), (6, 3)]),
                ("C0h".into(), vec![(1, 2), (6, 4)])
            ])
        );
        assert_eq!(
            mods(&read(SEQ, true, "C+mh,0,0;", Some(&[1, 2, 3, 4]), None)),
            Some(vec![
                ("C1m".into(), vec![(2, 3), (4, 1)]),
                ("C1h".into(), vec![(2, 4), (4, 2)])
            ])
        );
        // keys in order of first emission, not of the tag
        assert_eq!(
            mods(&read(SEQ, false, "C+h,1;C+m,0;", Some(&[5, 6]), None)),
            Some(vec![("C0m".into(), vec![(1, 6)]), ("C0h".into(), vec![(3, 5)])])
        );
        assert_eq!(
            mods(&read(SEQ, false, "C+27551,0;", Some(&[5]), None)),
            Some(vec![("C027551".into(), vec![(1, 5)])])
        );
        assert_eq!(
            mods(&read(SEQ, true, "N+n,0,2,1;", Some(&[1, 2, 3]), None)),
            Some(vec![("N1n".into(), vec![(3, 3), (5, 2), (8, 1)])])
        );
        // N canonical matches an IUPAC base in SEQ; C does not
        assert_eq!(
            mods(&read("ACNCGMC", false, "C+m,1,0;", Some(&[1, 2]), None)),
            Some(vec![("C0m".into(), vec![(3, 1), (6, 2)])])
        );
    }

    #[test]
    fn the_cases_pysam_refuses() {
        assert_eq!(mods(&read(SEQ, false, "C+m,0,5;", Some(&[1, 2]), None)), None); // overrun
        assert_eq!(mods(&read(SEQ, false, "C+m,0,1;", Some(&[10]), None)), None); // ML short
        assert_eq!(mods(&read(SEQ, false, "C+m,0,1;", Some(&[10, 200, 3]), None)), None); // ML long
        assert_eq!(mods(&read(SEQ, false, "C+m,0,1;", Some(&[10, 200]), Some(8))), None); // MN
        assert!(mods(&read(SEQ, false, "C+m,0,1;", Some(&[10, 200]), Some(9))).is_some());
        assert_eq!(mods(&read(SEQ, false, "C+m,0,x;", Some(&[1, 2]), None)), None);
        // no ML: positions with quality -1
        assert_eq!(
            mods(&read(SEQ, false, "C+m,0,1;", None, None)),
            Some(vec![("C0m".into(), vec![(1, -1), (6, -1)])])
        );
        // a delta that overflows later on stops the walk: pysam keeps what came before
        assert_eq!(
            mods(&read(SEQ, false, "C+m,0,0,99999999999;", Some(&[1, 2, 3]), None)),
            Some(vec![("C0m".into(), vec![(1, 1)])])
        );
    }
}
