//! This script is designed to compile a multifasta database and accompanying files used by bamtax
//! from one or more directories of downloaded NCBI assemblies.
#![allow(clippy::unused_io_amount)]

use crate::assembly_dir_iterator::AssemblyDirIterator;
use crate::cmd_cluster::taxid_from_id_str;
use anyhow::{Context, Error};
use clap::Parser;
use seq_io::fasta::{Reader as FastaReader, Record};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use gzp::{
    deflate::Mgzip,
    par::compress::{Compression, ParCompress, ParCompressBuilder},
};

use flate2::read::GzDecoder;

const BUFWRITER_CAP: usize = 200 * 1024 * 1024;

#[derive(Parser)]
pub struct BuildDbArgs {
    #[arg(short = 'i', long, num_args = 1..)]
    pub inputs: Vec<String>,

    #[arg(short, long = "reheadered_fastas", num_args = 0..)]
    pub reheadered_fastas: Vec<String>,

    #[arg(short = 'o', long)]
    pub output_prefix: String,

    #[arg(short = 'n', long, default_value_t = 0.30)]
    pub max_frac_ambig: f32,

    #[arg(short = 't', long, num_args = 1..)]
    pub assembly_to_taxid_map: Vec<String>,

    #[arg(short = 'g', long)]
    pub gzip_fasta: bool,

    #[arg(short = 'l', long, default_value_t = 100)]
    pub min_len: usize,

    /// Approximate uncompressed FASTA size per chunk (for example, 500M)
    #[arg(long, value_parser = parse_chunksize)]
    pub chunksize: Option<u64>,
}

fn parse_chunksize(value: &str) -> Result<u64, String> {
    let Some((number, suffix)) = value.split_at_checked(value.len().saturating_sub(1)) else {
        return Err("chunk size must be a positive integer followed by K, M, or G".to_string());
    };
    let multiplier = match suffix {
        "K" => 1024_u64,
        "M" => 1024_u64.pow(2),
        "G" => 1024_u64.pow(3),
        _ => {
            return Err("chunk size must be a positive integer followed by K, M, or G".to_string());
        }
    };
    if number.is_empty() || !number.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err("chunk size must be a positive integer followed by K, M, or G".to_string());
    }
    let number = number
        .parse::<u64>()
        .map_err(|_| "chunk size is too large".to_string())?;
    if number == 0 {
        return Err("chunk size must be greater than zero".to_string());
    }
    number
        .checked_mul(multiplier)
        .ok_or_else(|| "chunk size is too large".to_string())
}

struct DatabaseWriter {
    output_prefix: String,
    gzip_fasta: bool,
    chunksize: Option<u64>,
    chunk_number: u64,
    chunk_bytes: u64,
    fasta_writer: Box<dyn Write>,
    header_writer: Box<dyn Write>,
}

impl DatabaseWriter {
    fn new(output_prefix: String, gzip_fasta: bool, chunksize: Option<u64>) -> Result<Self, Error> {
        let chunk_number = u64::from(chunksize.is_some());
        let fasta_writer =
            Self::create_fasta_writer(&output_prefix, gzip_fasta, chunksize.map(|_| chunk_number))?;
        let header_file = std::fs::File::create(format!("{output_prefix}_headers.txt"))?;
        header_file.try_lock()?;
        let header_writer = Box::new(std::io::BufWriter::with_capacity(
            BUFWRITER_CAP,
            header_file,
        ));
        Ok(Self {
            output_prefix,
            gzip_fasta,
            chunksize,
            chunk_number,
            chunk_bytes: 0,
            fasta_writer,
            header_writer,
        })
    }

    fn create_fasta_writer(
        output_prefix: &str,
        gzip_fasta: bool,
        chunk_number: Option<u64>,
    ) -> Result<Box<dyn Write>, Error> {
        let output_prefix = match chunk_number {
            Some(number) => format!("{output_prefix}.{number}"),
            None => output_prefix.to_string(),
        };
        if gzip_fasta {
            let output_file = std::fs::File::create(format!("{output_prefix}.fasta.gz"))?;
            output_file.try_lock()?;
            let writer: ParCompress<Mgzip, _> = ParCompressBuilder::new()
                .compression_level(Compression::new(4))
                .num_threads(num_cpus::get())?
                .from_writer(std::io::BufWriter::with_capacity(
                    BUFWRITER_CAP,
                    output_file,
                ));
            Ok(Box::new(writer))
        } else {
            let output_file = std::fs::File::create(format!("{output_prefix}.fasta"))?;
            output_file.try_lock()?;
            Ok(Box::new(std::io::BufWriter::with_capacity(
                BUFWRITER_CAP,
                output_file,
            )))
        }
    }

    fn write_record(&mut self, id: &str, seq: &[u8]) -> Result<(), Error> {
        let record_bytes = 1_u64
            .checked_add(id.len() as u64)
            .and_then(|size| size.checked_add(1))
            .and_then(|size| size.checked_add(seq.len() as u64))
            .and_then(|size| size.checked_add(1))
            .context("FASTA record is too large")?;

        if self.chunksize.is_some_and(|limit| {
            self.chunk_bytes > 0 && self.chunk_bytes.saturating_add(record_bytes) > limit
        }) {
            self.fasta_writer.flush()?;
            self.chunk_number += 1;
            self.fasta_writer = Self::create_fasta_writer(
                &self.output_prefix,
                self.gzip_fasta,
                Some(self.chunk_number),
            )?;
            self.chunk_bytes = 0;
        }

        self.fasta_writer.write_all(b">")?;
        self.fasta_writer.write_all(id.as_bytes())?;
        self.fasta_writer.write_all(b"\n")?;
        self.fasta_writer.write_all(seq)?;
        self.fasta_writer.write_all(b"\n")?;
        self.header_writer.write_all(id.as_bytes())?;
        self.header_writer.write_all(b"\n")?;
        self.chunk_bytes += record_bytes;
        Ok(())
    }

    fn flush(&mut self) -> Result<(), Error> {
        self.header_writer.flush()?;
        self.fasta_writer.flush()?;
        Ok(())
    }
}

#[inline(always)]
fn base_is_nonambig(base: u8) -> bool {
    let b = base.to_ascii_uppercase();
    b == b'A' || b == b'C' || b == b'G' || b == b'T'
}

fn seq_n_bases_ambig_and_total(seq: &[u8]) -> (usize, usize) {
    let mut ambig = 0;
    let mut total = 0;

    for b in seq.iter().copied() {
        if !base_is_nonambig(b) {
            ambig += 1;
        }
        total += 1;
    }

    (ambig, total)
}

fn construct_assembly_to_tid_db<P: AsRef<Path>>(
    paths: &[P],
) -> Result<HashMap<String, u32>, Error> {
    let mut line = String::new();
    let mut ret = HashMap::new();

    for p in paths {
        let f = std::fs::File::open(p)?;
        f.try_lock()?;
        let mut reader = BufReader::new(f);

        while reader.read_line(&mut line)? > 0 {
            let (assembly, taxid) = line
                .trim()
                .split_once("\t")
                .with_context(|| format!("Invalid assembly to tid line: {line}"))?;

            let taxid = taxid
                .parse::<u32>()
                .with_context(|| format!("Invalid taxid in line {taxid}"))?;

            if !ret.contains_key(assembly) {
                ret.insert(assembly.to_string(), taxid);
            }

            line.clear();
        }
    }

    Ok(ret)
}

#[inline(always)]
fn process_metadata_fasta(
    taxid: u32,
    fasta: &str,
    writer: &mut DatabaseWriter,
    seen_records: &mut HashSet<String>,
    min_len: usize,
    max_frac_ambig: f32,
    headers_already_changed: bool,
) -> Result<(), Error> {
    let reader: Box<dyn std::io::Read> = if fasta.ends_with(".gz") {
        Box::new(GzDecoder::new(std::fs::File::open(fasta)?))
    } else {
        Box::new(std::fs::File::open(fasta)?)
    };

    let mut reader = FastaReader::new(BufReader::new(reader));

    while let Some(rec) = reader.next() {
        let rec = rec?;
        let id = rec.id()?;

        if seen_records.contains(id) {
            continue;
        }

        let (ambig, total) = seq_n_bases_ambig_and_total(rec.seq());
        let frac_ambig = ambig as f32 / total as f32;

        if total < min_len || frac_ambig > max_frac_ambig {
            eprintln!(
                "skipping sequence {id}: failed filters. Length: {}. Fraction of sequence ambiguous: {}",
                total, frac_ambig
            );
            continue;
        }

        seen_records.insert(id.to_string());

        // let taxid = assembly_tid_map.get(&assembly).with_context(|| {
        //     format!("Assembly {assembly} not found in assembly to taxon id map!")
        // })?;

        let new_id = if !headers_already_changed {
            let (acc, desc) = match id.split_once(" ") {
                Some((acc, desc)) => (acc.to_string(), desc.replace(" ", "_")),
                None => (id.to_string(), "".to_string()),
            };
            format!("{}|taxid:{}|{}", acc, taxid, desc)
        } else {
            let _ = taxid_from_id_str(id).with_context(|| {
                format!("record in pre-reheadered fasta has invalid header: {id}")
            })?;
            id.to_string()
        };

        writer.write_record(&new_id, rec.seq())?;
    }

    Ok(())
}

pub fn build_db_main(args: BuildDbArgs) -> Result<(), Error> {
    let mut seen: HashSet<String> = HashSet::new();

    let assembly_tid_map = construct_assembly_to_tid_db(&args.assembly_to_taxid_map)?;

    let mut writer = DatabaseWriter::new(args.output_prefix, args.gzip_fasta, args.chunksize)?;

    for input in args.inputs {
        let mut iterator = AssemblyDirIterator::new(input)?;

        while let Some((assembly, fasta)) = iterator.next_item()? {
            let taxid = assembly_tid_map.get(&assembly).with_context(|| {
                format!("Assembly {assembly} not found in assembly to taxon id map!")
            })?;

            process_metadata_fasta(
                *taxid,
                fasta.as_str(),
                &mut writer,
                &mut seen,
                args.min_len,
                args.max_frac_ambig,
                false,
            )?;
        }
    }

    for file in &args.reheadered_fastas {
        process_metadata_fasta(
            0,
            file,
            &mut writer,
            &mut seen,
            args.min_len,
            args.max_frac_ambig,
            true,
        )?;
    }

    writer.flush()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{BuildDbArgs, DatabaseWriter, parse_chunksize};
    use clap::Parser;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn parses_strict_chunksizes() {
        assert_eq!(parse_chunksize("1K"), Ok(1024));
        assert_eq!(parse_chunksize("12M"), Ok(12 * 1024 * 1024));
        assert_eq!(parse_chunksize("2G"), Ok(2 * 1024 * 1024 * 1024));

        for invalid in ["", "0K", "1", "1k", "1KB", "1.5M", "+1M", " 1M"] {
            assert!(parse_chunksize(invalid).is_err(), "accepted {invalid:?}");
        }
    }

    #[test]
    fn clap_rejects_an_invalid_chunksize() {
        assert!(
            BuildDbArgs::try_parse_from(["build-db", "-o", "db", "--chunksize", "10MB"]).is_err()
        );
    }

    #[test]
    fn chunks_without_splitting_fasta_records() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!("bamtax-build-db-{unique}"));
        fs::create_dir(&root).unwrap();
        let prefix = root.join("db").to_string_lossy().into_owned();
        let mut writer = DatabaseWriter::new(prefix, false, Some(15)).unwrap();

        writer.write_record("one", b"AAAA").unwrap();
        writer.write_record("two", b"CCCC").unwrap();
        writer.write_record("oversized", b"GGGGGGGGGG").unwrap();
        writer.flush().unwrap();
        drop(writer);

        assert_eq!(
            fs::read_to_string(root.join("db.1.fasta")).unwrap(),
            ">one\nAAAA\n"
        );
        assert_eq!(
            fs::read_to_string(root.join("db.2.fasta")).unwrap(),
            ">two\nCCCC\n"
        );
        assert_eq!(
            fs::read_to_string(root.join("db.3.fasta")).unwrap(),
            ">oversized\nGGGGGGGGGG\n"
        );
        assert_eq!(
            fs::read_to_string(root.join("db_headers.txt")).unwrap(),
            "one\ntwo\noversized\n"
        );
        assert!(!root.join("db.1_headers.txt").exists());

        fs::remove_dir_all(root).unwrap();
    }
}
