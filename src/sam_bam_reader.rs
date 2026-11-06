use anyhow::{Context, Error, bail};
use flate2::read::MultiGzDecoder;
use rust_htslib::bam::{Read, Reader as HtslibBamReader, ext::BamRecordExtensions, record::Aux};
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::cmd_cluster::taxid_from_id_str;

// ---------------------------------------------------------------------------
// Chunked name table
// ---------------------------------------------------------------------------

const CHUNK_SIZE: usize = 1_000_000;

/// Stores reference sequence names indexed by tid (zero-based @SQ order).
/// Split into fixed-size chunks to avoid one massive contiguous allocation.
pub struct NameTable {
    chunks: Vec<Vec<String>>,
}

impl NameTable {
    fn new() -> Self {
        Self { chunks: Vec::new() }
    }

    fn push(&mut self, name: String) {
        let total = self.chunks.iter().map(|c| c.len()).sum::<usize>();
        let chunk_idx = total / CHUNK_SIZE;
        if chunk_idx >= self.chunks.len() {
            self.chunks.push(Vec::with_capacity(CHUNK_SIZE));
        }
        self.chunks[chunk_idx].push(name);
    }

    pub fn get(&self, tid: usize) -> Option<&str> {
        let chunk_idx = tid / CHUNK_SIZE;
        let offset = tid % CHUNK_SIZE;
        self.chunks.get(chunk_idx)?.get(offset).map(String::as_str)
    }

    pub fn len(&self) -> usize {
        self.chunks.iter().map(|c| c.len()).sum()
    }
}

// ---------------------------------------------------------------------------
// Common record type
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct RawCigarOp {
    pub op: RawCigarKind,
    pub len: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RawCigarKind {
    Match,    // M — triggers error in filter_read_raw
    Equal,    // =
    Diff,     // X
    Ins,      // I
    SoftClip, // S
    Del,      // D
    RefSkip,  // N
    HardClip, // H
    Pad,      // P
}

#[derive(Default)]
pub struct RawRecord {
    pub qname: Vec<u8>,
    pub flags: u16,
    pub tid: i64, // -1 if unmapped
    pub pos: i64, // 0-based
    pub cigar: Vec<RawCigarOp>,
    pub seq_len: usize,
    pub ref_end: i64, // exclusive end on reference, computed from pos + cigar
    pub taxid: u32,
    pub rg_tag: String,
}

impl RawRecord {
    pub fn is_unmapped(&self) -> bool {
        self.flags & 0x4 != 0
    }
    pub fn is_reverse(&self) -> bool {
        self.flags & 0x10 != 0
    }
    pub fn is_secondary(&self) -> bool {
        self.flags & 0x100 != 0
    }
    pub fn is_supplementary(&self) -> bool {
        self.flags & 0x800 != 0
    }
}

// ---------------------------------------------------------------------------
// CIGAR parser (SAM text → Vec<RawCigarOp>)
// ---------------------------------------------------------------------------

fn parse_cigar(s: &[u8]) -> Result<(Vec<RawCigarOp>, i64), Error> {
    let mut ops = Vec::new();
    let mut num: u32 = 0;
    let mut ref_len: i64 = 0;

    for &b in s {
        if b.is_ascii_digit() {
            num = num * 10 + (b - b'0') as u32;
        } else {
            let kind = match b {
                b'M' => RawCigarKind::Match,
                b'=' => RawCigarKind::Equal,
                b'X' => RawCigarKind::Diff,
                b'I' => RawCigarKind::Ins,
                b'S' => RawCigarKind::SoftClip,
                b'D' => RawCigarKind::Del,
                b'N' => RawCigarKind::RefSkip,
                b'H' => RawCigarKind::HardClip,
                b'P' => RawCigarKind::Pad,
                other => bail!("Unknown CIGAR op: {}", other as char),
            };
            match kind {
                RawCigarKind::Equal
                | RawCigarKind::Diff
                | RawCigarKind::Match
                | RawCigarKind::Del
                | RawCigarKind::RefSkip => ref_len += num as i64,
                _ => {}
            }
            ops.push(RawCigarOp { op: kind, len: num });
            num = 0;
        }
    }

    Ok((ops, ref_len))
}

// ---------------------------------------------------------------------------
// SAM header scanner — pass 1
// ---------------------------------------------------------------------------

fn scan_sam_header<R: BufRead>(reader: &mut R, names: &mut NameTable) -> Result<(), Error> {
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        } // EOF
        if !line.starts_with('@') {
            break;
        }
        if line.starts_with("@SQ") {
            // Extract SN: field
            if let Some(sn) = line.split('\t').find(|f| f.starts_with("SN:")) {
                names.push(sn[3..].trim_end().to_string());
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// SAM record parser — pass 2
// ---------------------------------------------------------------------------

fn parse_sam_record_into(line: &str, out: &mut RawRecord) -> Result<(), Error> {
    let mut fields = line.splitn(12, '\t');

    macro_rules! next {
        ($name:literal) => {
            fields.next().with_context(|| format!("missing SAM field: {}", $name))?
        };
    }

    out.qname.clear();
    out.qname.extend_from_slice(next!("QNAME").as_bytes());
    out.flags = next!("FLAG").parse::<u16>().context("FLAG")?;
    let rname = next!("RNAME");
    out.pos = next!("POS").parse::<i64>().context("POS")? - 1;
    let _mapq = next!("MAPQ");
    let cigar_s = next!("CIGAR").as_bytes();
    let _rnext = next!("RNEXT");
    let _pnext = next!("PNEXT");
    let _tlen = next!("TLEN");
    let seq = next!("SEQ");
    let _qual = next!("QUAL");
    let tags = fields.next().unwrap_or("");

    out.tid = if rname == "*" { -1 } else { 0 };
    out.taxid = if rname == "*" {
        0
    } else {
        taxid_from_id_str(rname).with_context(|| format!("parsing taxid from RNAME: {rname}"))?
    };

    out.cigar.clear();
    let ref_consumed = if cigar_s == b"*" {
        0
    } else {
        let (ops, consumed) = parse_cigar(cigar_s)?;
        out.cigar = ops;
        consumed
    };

    out.seq_len = seq.len();
    out.ref_end = out.pos + ref_consumed;

    out.rg_tag.clear();
    if let Some(f) = tags.split('\t').find(|f| f.starts_with("RG:Z:")) {
        out.rg_tag.push_str(f[5..].trim_end());
    } else {
        out.rg_tag.push('0');
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Backends
// ---------------------------------------------------------------------------

pub struct BamBackend {
    reader: HtslibBamReader,
    names: NameTable,
    hts_rec: rust_htslib::bam::Record,
}

impl BamBackend {
    fn open(path: &Path) -> Result<Self, Error> {
        let mut reader = HtslibBamReader::from_path(path).context("open BAM")?;
        reader.set_threads(4).unwrap();

        let mut names = NameTable::new();
        for name in reader.header().target_names() {
            names.push(String::from_utf8_lossy(name).into_owned());
        }

        Ok(Self { reader, names, hts_rec: rust_htslib::bam::Record::new() })
    }
}

pub struct SamBackend {
    names: NameTable,
    sam_reader: SamRecordReader,
}

impl SamBackend {
    fn open(path: &Path) -> Result<Self, Error> {
        let mut names = NameTable::new();

        // Pass 1: scan header
        let f = std::fs::File::open(path).context("open SAM")?;
        let mut rdr = BufReader::new(MultiGzDecoder::new(f));
        scan_sam_header(&mut rdr, &mut names)?;

        eprintln!("SAM header scanned: {} reference sequences", names.len());

        // Pass 2: open streaming record reader
        let sam_reader = SamRecordReader::open(path)?;

        Ok(Self { names, sam_reader })
    }
}

// ---------------------------------------------------------------------------
// Unified reader
// ---------------------------------------------------------------------------

pub enum SamBamReader {
    Bam(BamBackend),
    Sam(SamBackend),
}

impl SamBamReader {
    pub fn open(path: &Path) -> Result<Self, Error> {
        let ext = path.file_name().and_then(|n| n.to_str()).unwrap_or("");

        if ext.ends_with(".bam") {
            Ok(Self::Bam(BamBackend::open(path)?))
        } else if ext.ends_with(".sam") || ext.ends_with(".sam.gz") {
            Ok(Self::Sam(SamBackend::open(path)?))
        } else {
            bail!("Unrecognised file extension for {:?}: expected .bam, .sam, or .sam.gz", path)
        }
    }

    pub fn tid_to_name(&self, tid: usize) -> Option<&str> {
        match self {
            Self::Bam(b) => b.names.get(tid),
            Self::Sam(s) => s.names.get(tid),
        }
    }
}

// ---------------------------------------------------------------------------
// BAM transcription helper
// ---------------------------------------------------------------------------

fn transcribe_bam_record_into(
    rec: &rust_htslib::bam::Record,
    names: &NameTable,
    out: &mut RawRecord,
) -> Result<(), Error> {
    use rust_htslib::bam::record::Cigar as HCigar;

    out.cigar.clear();
    for op in &rec.cigar().0 {
        let (kind, len) = match op {
            HCigar::Match(n) => (RawCigarKind::Match, *n),
            HCigar::Equal(n) => (RawCigarKind::Equal, *n),
            HCigar::Diff(n) => (RawCigarKind::Diff, *n),
            HCigar::Ins(n) => (RawCigarKind::Ins, *n),
            HCigar::SoftClip(n) => (RawCigarKind::SoftClip, *n),
            HCigar::Del(n) => (RawCigarKind::Del, *n),
            HCigar::RefSkip(n) => (RawCigarKind::RefSkip, *n),
            HCigar::HardClip(n) => (RawCigarKind::HardClip, *n),
            HCigar::Pad(n) => (RawCigarKind::Pad, *n),
        };
        out.cigar.push(RawCigarOp { op: kind, len });
    }

    out.taxid = if rec.tid() < 0 {
        0
    } else {
        let name = names.get(rec.tid() as usize).with_context(|| format!("tid {} not in name table", rec.tid()))?;
        taxid_from_id_str(name).with_context(|| format!("BAM tid {} name: {}", rec.tid(), name))?
    };

    out.qname.clear();
    out.qname.extend_from_slice(rec.qname());
    out.flags = rec.flags();
    out.tid = rec.tid() as i64;
    out.pos = rec.pos();
    out.seq_len = rec.seq_len();
    out.ref_end = rec.reference_end();
    out.rg_tag.clear();
    match rec.aux(b"RG") {
        Ok(Aux::String(s)) => out.rg_tag.push_str(s),
        _ => out.rg_tag.push('0'),
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// SAM record reader — pass 2, streaming, reuses caller-supplied RawRecord
// ---------------------------------------------------------------------------

pub struct SamRecordReader {
    inner: BufReader<MultiGzDecoder<std::fs::File>>,
    line: String,
}

impl SamRecordReader {
    fn open(path: &Path) -> Result<Self, Error> {
        let f = std::fs::File::open(path).context("open SAM for record pass")?;
        let mut inner = BufReader::new(MultiGzDecoder::new(f));

        // Skip header lines
        let mut line = String::new();
        loop {
            line.clear();
            let n = inner.read_line(&mut line)?;
            if n == 0 {
                break;
            }
            if !line.starts_with('@') {
                // First record line already in buffer — leave it there;
                // next read() call will parse it without re-reading.
                return Ok(Self { inner, line });
            }
        }
        Ok(Self { inner, line: String::new() })
    }

    /// Read the next record into `rec`, reusing its allocations.
    /// Returns `Some(Ok(()))` on success, `Some(Err(_))` on parse error, `None` at EOF.
    pub fn read(&mut self, rec: &mut RawRecord) -> Option<Result<(), Error>> {
        // If line is non-empty we already have content from the header-skip pass
        if self.line.is_empty() {
            self.line.clear();
            match self.inner.read_line(&mut self.line) {
                Err(e) => return Some(Err(e.into())),
                Ok(0) => return None,
                Ok(_) => {}
            }
        }

        let line = self.line.trim_end_matches('\n').trim_end_matches('\r');
        if line.is_empty() {
            self.line.clear();
            return None;
        }

        let result = parse_sam_record_into(line, rec);
        self.line.clear();
        Some(result)
    }
}

// ---------------------------------------------------------------------------
// SamBamReader public read() API
// ---------------------------------------------------------------------------

impl SamBamReader {
    /// Read the next record into `rec`, reusing its allocations.
    /// Returns `Some(Ok(()))` on success, `Some(Err(_))` on error, `None` at EOF.
    pub fn read(&mut self, rec: &mut RawRecord) -> Option<Result<(), Error>> {
        match self {
            Self::Bam(b) => {
                let result = b.reader.read(&mut b.hts_rec);
                let _ = result?; // None = EOF
                Some(transcribe_bam_record_into(&b.hts_rec, &b.names, rec))
            }
            Self::Sam(s) => s.sam_reader.read(rec),
        }
    }
}
