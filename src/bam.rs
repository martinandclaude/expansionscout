//! Alignment and reference access, reproducing pysam's view of a record.
//!
//! noodles does the container work: BGZF, the BAM record framing, BAI/CSI
//! chunk lookup, FASTA indexes. Everything the caller can observe about a
//! record is decoded here from the raw fields, following htslib 1.24 and
//! pysam 0.24.1 (the versions the Python implementation runs on), because
//! small differences -- which records a region query returns, how a CIGAR
//! with a padding op advances the query, what a reverse read's MM tag points
//! at -- change the output.

use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::region::Interval;
use noodles_core::Position;
use noodles_csi::BinningIndex;
use noodles_fasta as fasta;
use noodles_sam as sam;

pub const CMATCH: u8 = 0;
pub const CINS: u8 = 1;
pub const CDEL: u8 = 2;
pub const CREF_SKIP: u8 = 3;
pub const CSOFT_CLIP: u8 = 4;
pub const CHARD_CLIP: u8 = 5;
pub const CPAD: u8 = 6;
pub const CEQUAL: u8 = 7;
pub const CDIFF: u8 = 8;

const FUNMAP: u16 = 0x4;
const FREVERSE: u16 = 0x10;
const FSECONDARY: u16 = 0x100;

/// One alignment record, decoded.
#[derive(Clone, Debug)]
pub struct Aln {
    pub name: String,
    pub flag: u16,
    pub pos: i64,
    pub mapq: u8,
    pub cigar: Vec<(u8, u32)>,
    /// Sequence as pysam's `query_sequence`; empty when SEQ is `*`.
    pub seq: Vec<u8>,
    /// The 4-bit codes, for the MM target test.
    pub seq4: Vec<u8>,
    pub aux: Vec<u8>,
}

/// An aux field value, as much of it as this tool reads.
#[derive(Clone, Debug)]
pub enum Aux<'a> {
    Int(i64),
    Float(f64),
    Char(u8),
    Str(&'a [u8]),
    ArrayU8(&'a [u8]),
    /// Any other array: the subtype and the element count.
    Array(u8),
}

fn aux_size(t: u8) -> Option<usize> {
    Some(match t {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        _ => return None,
    })
}

impl Aln {
    pub fn is_reverse(&self) -> bool {
        self.flag & FREVERSE != 0
    }
    pub fn is_unmapped(&self) -> bool {
        self.flag & FUNMAP != 0
    }
    pub fn is_secondary(&self) -> bool {
        self.flag & FSECONDARY != 0
    }

    /// `query_sequence`, or None for SEQ `*`.
    pub fn query_sequence(&self) -> Option<&[u8]> {
        if self.seq.is_empty() {
            None
        } else {
            Some(&self.seq)
        }
    }

    /// Reference length as htslib's `bam_cigar2rlen`.
    fn rlen(&self) -> i64 {
        self.cigar
            .iter()
            .filter(|(op, _)| matches!(*op, CMATCH | CDEL | CREF_SKIP | CEQUAL | CDIFF))
            .map(|&(_, l)| l as i64)
            .sum()
    }

    /// htslib's `bam_endpos`: never less than pos + 1.
    pub fn endpos(&self) -> i64 {
        let r = if self.is_unmapped() { 0 } else { self.rlen() };
        self.pos + r.max(1)
    }

    /// `reference_end`: None for an unmapped record or one with no CIGAR.
    pub fn reference_end(&self) -> Option<i64> {
        if self.is_unmapped() || self.cigar.is_empty() {
            None
        } else {
            Some(self.endpos())
        }
    }

    /// `infer_read_length()`: M, I, S, H, = and X; None without a CIGAR.
    pub fn infer_read_length(&self) -> Option<i64> {
        if self.cigar.is_empty() {
            return None;
        }
        Some(
            self.cigar
                .iter()
                .filter(|(op, _)| matches!(*op, CMATCH | CINS | CSOFT_CLIP | CHARD_CLIP | CEQUAL | CDIFF))
                .map(|&(_, l)| l as i64)
                .sum(),
        )
    }

    /// `query_alignment_start`.
    pub fn query_alignment_start(&self) -> Result<i64, String> {
        let l_qseq = self.seq.len() as i64;
        let mut off = 0i64;
        for &(op, len) in &self.cigar {
            match op {
                CHARD_CLIP => {
                    if off != 0 && off != l_qseq {
                        return Err("Invalid clipping in CIGAR string".into());
                    }
                }
                CSOFT_CLIP => off += len as i64,
                _ => break,
            }
        }
        Ok(off)
    }

    /// `query_alignment_end`. The backward walk stops at index 1, as pysam's does.
    pub fn query_alignment_end(&self) -> Result<i64, String> {
        let l_qseq = self.seq.len() as i64;
        let mut end = l_qseq;
        if end == 0 {
            for &(op, len) in &self.cigar {
                if matches!(op, CMATCH | CINS | CEQUAL | CDIFF) || (op == CSOFT_CLIP && end == 0) {
                    end += len as i64;
                }
            }
        } else {
            for k in (1..self.cigar.len()).rev() {
                let (op, len) = self.cigar[k];
                match op {
                    CHARD_CLIP => {
                        if end != l_qseq {
                            return Err("Invalid clipping in CIGAR string".into());
                        }
                    }
                    CSOFT_CLIP => end -= len as i64,
                    _ => break,
                }
            }
        }
        Ok(end)
    }

    /// `get_aligned_pairs(matches_only=True)`, as (query, reference) arrays.
    /// A padding op advances the query, as pysam's walk does.
    pub fn aligned_pairs(&self) -> (Vec<i64>, Vec<i64>) {
        let (mut q, mut r) = (Vec::new(), Vec::new());
        if self.cigar.is_empty() {
            return (q, r);
        }
        // pysam holds both counters as uint32, so pos = -1 wraps.
        let mut qpos: u32 = 0;
        let mut rpos: u32 = self.pos as u32;
        for &(op, len) in &self.cigar {
            match op {
                CMATCH | CEQUAL | CDIFF => {
                    for i in 0..len {
                        q.push(qpos.wrapping_add(i) as i64);
                        r.push(rpos.wrapping_add(i) as i64);
                    }
                    qpos = qpos.wrapping_add(len);
                    rpos = rpos.wrapping_add(len);
                }
                CINS | CSOFT_CLIP | CPAD => qpos = qpos.wrapping_add(len),
                CDEL | CREF_SKIP => rpos = rpos.wrapping_add(len),
                _ => {}
            }
        }
        (q, r)
    }

    /// Iterate aux fields as (tag, type, value).
    pub fn aux_fields(&self) -> Vec<([u8; 2], u8, Aux<'_>)> {
        let a = &self.aux;
        let mut out = Vec::new();
        let mut i = 0;
        while i + 3 <= a.len() {
            let tag = [a[i], a[i + 1]];
            let t = a[i + 2];
            i += 3;
            let rd = |b: &[u8], n: usize| -> Option<u64> {
                (b.len() >= n).then(|| b[..n].iter().rev().fold(0u64, |acc, &x| (acc << 8) | x as u64))
            };
            let v = match t {
                b'A' => {
                    let Some(&c) = a.get(i) else { break };
                    i += 1;
                    Aux::Char(c)
                }
                b'c' => {
                    let Some(&c) = a.get(i) else { break };
                    i += 1;
                    Aux::Int(c as i8 as i64)
                }
                b'C' => {
                    let Some(&c) = a.get(i) else { break };
                    i += 1;
                    Aux::Int(c as i64)
                }
                b's' => {
                    let Some(v) = rd(&a[i..], 2) else { break };
                    i += 2;
                    Aux::Int(v as u16 as i16 as i64)
                }
                b'S' => {
                    let Some(v) = rd(&a[i..], 2) else { break };
                    i += 2;
                    Aux::Int(v as i64)
                }
                b'i' => {
                    let Some(v) = rd(&a[i..], 4) else { break };
                    i += 4;
                    Aux::Int(v as u32 as i32 as i64)
                }
                b'I' => {
                    let Some(v) = rd(&a[i..], 4) else { break };
                    i += 4;
                    Aux::Int(v as i64)
                }
                b'f' => {
                    let Some(v) = rd(&a[i..], 4) else { break };
                    i += 4;
                    Aux::Float(f32::from_bits(v as u32) as f64)
                }
                b'd' => {
                    let Some(v) = rd(&a[i..], 8) else { break };
                    i += 8;
                    Aux::Float(f64::from_bits(v))
                }
                b'Z' | b'H' => {
                    let Some(z) = a[i..].iter().position(|&c| c == 0) else {
                        break;
                    };
                    let s = &a[i..i + z];
                    i += z + 1;
                    Aux::Str(s)
                }
                b'B' => {
                    let Some(&sub) = a.get(i) else { break };
                    let Some(n) = rd(&a[i + 1..], 4) else { break };
                    let Some(size) = aux_size(sub) else { break };
                    let start = i + 5;
                    let len = n as usize * size;
                    if start + len > a.len() {
                        break;
                    }
                    i = start + len;
                    if sub == b'C' {
                        Aux::ArrayU8(&a[start..start + len])
                    } else {
                        Aux::Array(sub)
                    }
                }
                _ => break,
            };
            out.push((tag, t, v));
        }
        out
    }

    pub fn tag(&self, name: &[u8; 2]) -> Option<(u8, Aux<'_>)> {
        self.aux_fields()
            .into_iter()
            .find(|(t, _, _)| t == name)
            .map(|(_, ty, v)| (ty, v))
    }

    /// The HP tag as `int(read.get_tag("HP"))`, 0 when absent or not an integer.
    pub fn hp(&self) -> i64 {
        match self.tag(b"HP") {
            Some((_, Aux::Int(v))) => v,
            Some((_, Aux::Float(f))) => {
                if f.is_finite() {
                    f.trunc() as i64
                } else {
                    0
                }
            }
            Some((_, Aux::Str(s))) => py_int_str(s).unwrap_or(0),
            Some((_, Aux::Char(c))) => py_int_str(&[c]).unwrap_or(0),
            _ => 0,
        }
    }
}

/// Python's `int(str)` for the simple cases: optional sign, digits,
/// surrounding whitespace, and underscores between digits.
fn py_int_str(s: &[u8]) -> Option<i64> {
    let t = std::str::from_utf8(s).ok()?.trim();
    let (neg, body) = match t.as_bytes().first() {
        Some(b'-') => (true, &t[1..]),
        Some(b'+') => (false, &t[1..]),
        _ => (false, t),
    };
    if body.is_empty() || body.starts_with('_') || body.ends_with('_') || body.contains("__") {
        return None;
    }
    let digits: String = body.chars().filter(|&c| c != '_').collect();
    if !digits.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let v: i64 = digits.parse().ok()?;
    Some(if neg { -v } else { v })
}

const SEQ_NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

fn decode(rec: &bam::Record) -> io::Result<Aln> {
    let flag = u16::from(rec.flags());
    let pos = match rec.alignment_start().transpose()? {
        Some(p) => usize::from(p) as i64 - 1,
        None => -1,
    };
    let mapq = rec.mapping_quality().map(u8::from).unwrap_or(255);
    let name = rec
        .name()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .unwrap_or_else(|| "*".into());
    let seq_len = rec.sequence().len();
    let packed = rec.sequence().as_bytes();
    let mut seq4 = Vec::with_capacity(seq_len);
    for i in 0..seq_len {
        let b = packed[i / 2];
        seq4.push(if i % 2 == 0 { b >> 4 } else { b & 0x0f });
    }
    let seq: Vec<u8> = seq4.iter().map(|&c| SEQ_NT16[c as usize]).collect();
    let raw = rec.cigar().as_bytes();
    let mut cigar: Vec<(u8, u32)> = raw
        .chunks_exact(4)
        .map(|c| {
            let v = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
            ((v & 0xf) as u8, v >> 4)
        })
        .collect();
    let mut aux = rec.data().as_bytes().to_vec();
    let tid_ok = rec.reference_sequence_id().transpose()?.is_some();
    let mut aln = Aln {
        name,
        flag,
        pos,
        mapq,
        cigar: Vec::new(),
        seq,
        seq4,
        aux: Vec::new(),
    };
    // htslib's bam_tag2cigar: a CIGAR too long for the record is stored in a
    // CG:B:I tag behind a placeholder that soft-clips the whole read.
    if !cigar.is_empty() && tid_ok && pos >= 0 && cigar[0].0 == CSOFT_CLIP && cigar[0].1 as usize == seq_len {
        aln.aux = aux.clone();
        let cg = aln
            .aux_fields()
            .into_iter()
            .find(|(t, _, _)| t == b"CG")
            .map(|(_, _, v)| v.clone());
        if let Some(Aux::Array(b'I' | b'i')) = cg {
            if let Some((start, n)) = find_cg(&aux) {
                let body = &aux[start + 8..start + 8 + n * 4];
                cigar = body
                    .chunks_exact(4)
                    .map(|c| {
                        let v = u32::from_le_bytes([c[0], c[1], c[2], c[3]]);
                        ((v & 0xf) as u8, v >> 4)
                    })
                    .collect();
                aux.drain(start..start + 8 + n * 4);
            }
        }
    }
    aln.cigar = cigar;
    aln.aux = aux;
    Ok(aln)
}

/// Byte offset and element count of a CG:B:I field in raw aux data.
fn find_cg(a: &[u8]) -> Option<(usize, usize)> {
    let mut i = 0;
    while i + 3 <= a.len() {
        let start = i;
        let t = a[i + 2];
        i += 3;
        match t {
            b'Z' | b'H' => i += a[i..].iter().position(|&c| c == 0)? + 1,
            b'B' => {
                let sub = *a.get(i)?;
                let n = u32::from_le_bytes(a.get(i + 1..i + 5)?.try_into().ok()?) as usize;
                if &a[start..start + 2] == b"CG" && (sub == b'I' || sub == b'i') {
                    return Some((start, n));
                }
                i += 5 + n * aux_size(sub)?;
            }
            other => i += aux_size(other)?,
        }
    }
    None
}

/// Open alignments: a coordinate-sorted BAM with its index.
pub struct Bam {
    reader: bam::io::Reader<bgzf::io::Reader<File>>,
    index: Box<dyn BinningIndex>,
    pub header: sam::Header,
    pub references: Vec<String>,
    pub lengths: Vec<i64>,
}

fn with_ext(p: &Path, ext: &str) -> PathBuf {
    let mut s = p.as_os_str().to_owned();
    s.push(ext);
    PathBuf::from(s)
}

fn read_index(path: &Path) -> io::Result<Box<dyn BinningIndex>> {
    let candidates = [
        with_ext(path, ".bai"),
        path.with_extension("bai"),
        with_ext(path, ".csi"),
        path.with_extension("csi"),
    ];
    for c in &candidates {
        if c.exists() {
            return if c.extension().is_some_and(|e| e == "csi") {
                Ok(Box::new(noodles_csi::fs::read(c)?))
            } else {
                Ok(Box::new(bam::bai::fs::read(c)?))
            };
        }
    }
    Err(io::Error::new(
        io::ErrorKind::NotFound,
        format!("no index for {} (looked for .bai and .csi)", path.display()),
    ))
}

impl Bam {
    pub fn open(path: &Path) -> io::Result<Bam> {
        let index = read_index(path)?;
        let mut reader = File::open(path).map(bam::io::Reader::new)?;
        let header = reader.read_header()?;
        let (references, lengths) = header
            .reference_sequences()
            .iter()
            .map(|(name, map)| {
                (
                    String::from_utf8_lossy(name.as_ref()).into_owned(),
                    usize::from(map.length()) as i64,
                )
            })
            .unzip();
        Ok(Bam {
            reader,
            index,
            header,
            references,
            lengths,
        })
    }

    /// The first read group's SM, as `_sample_from_bam` looks for it.
    pub fn sample(&self) -> Option<String> {
        use sam::header::record::value::map::read_group::tag::SAMPLE;
        self.header
            .read_groups()
            .values()
            .find_map(|rg| {
                rg.other_fields()
                    .get(&SAMPLE)
                    .map(|s| String::from_utf8_lossy(s.as_ref()).into_owned())
            })
            .filter(|s| !s.is_empty())
    }

    /// `fetch(contig, start, end)`: records with pos < end and endpos > start,
    /// in file order, as htslib's iterator returns them.
    pub fn fetch(&mut self, contig: &str, start: i64, end: i64) -> io::Result<Vec<Aln>> {
        let Some(tid) = self.references.iter().position(|r| r == contig) else {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "invalid contig"));
        };
        if start >= end {
            return Ok(vec![]);
        }
        // Chunks for a slightly wider interval, then htslib's own test.
        let lo = Position::new((start.max(0) + 1) as usize).unwrap();
        let hi = Position::new((end + 1) as usize).unwrap();
        let interval: Interval = (lo..=hi).into();
        let chunks = self.index.query(tid, interval)?;
        let query = noodles_csi::io::Query::new(self.reader.get_mut(), chunks);
        let mut r = bam::io::Reader::from(query);
        let mut rec = bam::Record::default();
        let mut out = Vec::new();
        while r.read_record(&mut rec)? != 0 {
            let rtid = rec.reference_sequence_id().transpose()?;
            if rtid != Some(tid) {
                if rtid.is_none() || rtid.unwrap() > tid {
                    break;
                }
                continue;
            }
            let a = decode(&rec)?;
            if a.pos >= end {
                break;
            }
            if a.endpos() > start {
                out.push(a);
            }
        }
        Ok(out)
    }
}

/// An indexed FASTA, plain or bgzip-compressed.
pub struct Fasta {
    reader: fasta::io::IndexedReader<fasta::io::BufReader<File>>,
    pub references: Vec<String>,
    pub lengths: Vec<i64>,
}

impl Fasta {
    pub fn open(path: &Path) -> io::Result<Fasta> {
        let reader = fasta::io::indexed_reader::Builder::default().build_from_path(path)?;
        let (references, lengths) = reader
            .index()
            .as_ref()
            .iter()
            .map(|r| {
                (
                    String::from_utf8_lossy(r.name().as_ref()).into_owned(),
                    r.length() as i64,
                )
            })
            .unzip();
        Ok(Fasta {
            reader,
            references,
            lengths,
        })
    }

    /// `fetch(contig, start, end)`, clamped as htslib clamps it; case preserved.
    pub fn fetch(&mut self, contig: &str, start: i64, end: i64) -> io::Result<Vec<u8>> {
        if start == end {
            return Ok(vec![]);
        }
        let Some(i) = self.references.iter().position(|r| r == contig) else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("sequence '{contig}' not present"),
            ));
        };
        let len = self.lengths[i];
        let start = start.max(0);
        if start >= len {
            return Ok(vec![]);
        }
        let end = end.min(len);
        if end <= start {
            return Ok(vec![]);
        }
        let region = noodles_core::Region::new(
            contig,
            Position::new(start as usize + 1).unwrap()..=Position::new(end as usize).unwrap(),
        );
        let rec = self.reader.query(&region)?;
        Ok(rec.sequence().as_ref().to_vec())
    }
}
