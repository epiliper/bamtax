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
    header_writer: &mut Box<dyn std::io::Write>,
    fasta_writer: &mut Box<dyn std::io::Write>,
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

        fasta_writer.write(b">")?;
        fasta_writer.write(new_id.as_bytes())?;
        fasta_writer.write(b"\n")?;
        fasta_writer.write(rec.seq())?;
        fasta_writer.write(b"\n")?;

        header_writer.write(new_id.as_bytes())?;
        header_writer.write(b"\n")?;
    }

    Ok(())
}

pub fn build_db_main(args: BuildDbArgs) -> Result<(), Error> {
    let mut seen: HashSet<String> = HashSet::new();

    let assembly_tid_map = construct_assembly_to_tid_db(&args.assembly_to_taxid_map)?;

    let mut writer: Box<dyn std::io::Write> = if args.gzip_fasta {
        let output_file = std::fs::File::create(format!("{}.fasta.gz", args.output_prefix))?;
        output_file.try_lock()?;

        let writer: ParCompress<Mgzip, _> = ParCompressBuilder::new()
            .compression_level(Compression::new(4))
            .num_threads(num_cpus::get())?
            .from_writer(std::io::BufWriter::with_capacity(
                BUFWRITER_CAP,
                output_file,
            ));
        Box::new(writer)
    } else {
        let output_file = std::fs::File::create(format!("{}.fasta", args.output_prefix))?;
        output_file.try_lock()?;
        Box::new(std::io::BufWriter::with_capacity(
            BUFWRITER_CAP,
            output_file,
        ))
    };

    let mut header_writer: Box<dyn Write> = {
        let file = std::fs::File::create(format!("{}_headers.txt", args.output_prefix))?;
        file.try_lock()?;
        Box::new(std::io::BufWriter::with_capacity(BUFWRITER_CAP, file))
    };

    for input in args.inputs {
        let mut iterator = AssemblyDirIterator::new(input)?;

        while let Some((assembly, fasta)) = iterator.next_item()? {
            let taxid = assembly_tid_map.get(&assembly).with_context(|| {
                format!("Assembly {assembly} not found in assembly to taxon id map!")
            })?;

            process_metadata_fasta(
                *taxid,
                fasta.as_str(),
                &mut header_writer,
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
            &mut header_writer,
            &mut writer,
            &mut seen,
            args.min_len,
            args.max_frac_ambig,
            true,
        )?;
    }

    header_writer.flush()?;
    writer.flush()?;

    Ok(())
}
