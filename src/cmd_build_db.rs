//! This script is designed to compile a multifasta database and accompanying files used by bamtax
//! from one or more directories of downloaded NCBI assemblies.
#![allow(clippy::unused_io_amount)]

use crate::assembly_dir_iterator::AssemblyDirIterator;
use anyhow::{Context, Error};
use clap::Parser;
use seq_io::fasta::{Reader, Record};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use gzp::{
    deflate::Mgzip,
    par::compress::{ParCompress, ParCompressBuilder},
};

#[derive(Parser)]
pub struct BuildDbArgs {
    #[arg(short = 'i', long, value_delimiter = ',', num_args = 1..)]
    pub inputs: Vec<String>,

    #[arg(short = 'o', long)]
    pub output_prefix: String,

    #[arg(short = 'n', long, default_value_t = 0.30)]
    pub max_frac_ambig: f32,

    #[arg(short = 't', long, value_delimiter = ',', num_args = 1..)]
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
        if base_is_nonambig(b) {
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
        let mut reader = BufReader::new(std::fs::File::open(p)?);

        while reader.read_line(&mut line)? > 0 {
            let (assembly, taxid) = line
                .split_once("\t")
                .with_context(|| format!("Invalid assembly to tid line: {line}"))?;

            let taxid = taxid
                .parse::<u32>()
                .with_context(|| format!("Invalid taxid in line {line}"))?;

            if !ret.contains_key(assembly) {
                ret.insert(assembly.to_string(), taxid);
            }

            line.clear();
        }
    }

    Ok(ret)
}

pub fn build_db_main(args: BuildDbArgs) -> Result<(), Error> {
    let mut seen: HashSet<String> = HashSet::new();

    let assembly_tid_map = construct_assembly_to_tid_db(&args.assembly_to_taxid_map)?;

    let mut writer: Box<dyn std::io::Write> = if args.gzip_fasta {
        let output_file = std::fs::File::create(format!("{}.fasta.gz", args.output_prefix))?;
        let writer: ParCompress<Mgzip, _> =
            ParCompressBuilder::new().from_writer(std::io::BufWriter::new(output_file));
        Box::new(writer)
    } else {
        let output_file = std::fs::File::create(format!("{}.fasta", args.output_prefix))?;
        Box::new(std::io::BufWriter::new(output_file))
    };

    let mut header_writer = std::io::BufWriter::new(std::fs::File::create(format!(
        "{}_headers.txt",
        args.output_prefix
    ))?);

    for input in args.inputs {
        let mut iterator = AssemblyDirIterator::new(input)?;

        while let Some((assembly, fasta)) = iterator.next_item()? {
            let mut reader = Reader::new(BufReader::new(std::fs::File::open(fasta)?));
            while let Some(rec) = reader.next() {
                let rec = rec?;
                let id = rec.id()?;

                if seen.contains(id) {
                    continue;
                }

                let (ambig, total) = seq_n_bases_ambig_and_total(rec.seq());
                let frac_ambig = ambig as f32 / total as f32;

                if total < args.min_len || frac_ambig > args.max_frac_ambig {
                    eprintln!("skipping sequence {id}: failed filters");
                }

                seen.insert(id.to_string());

                let taxid = assembly_tid_map.get(&assembly).with_context(|| {
                    format!("Assembly {assembly} not found in assembly to taxon id map!")
                })?;

                let (acc, desc) = match id.split_once(" ") {
                    Some((acc, desc)) => (acc.to_string(), desc.replace(" ", "_")),
                    None => (id.to_string(), "".to_string()),
                };

                let new_id = format!("{}|taxid:{}|{}", acc, taxid, desc);

                writer.write(new_id.as_bytes())?;
                writer.write(b"\n")?;
                writer.write(rec.seq())?;
                writer.write(b"\n")?;

                header_writer.write(new_id.as_bytes())?;
                header_writer.write(b"\n")?;
            }
        }
    }

    header_writer.flush()?;
    writer.flush()?;

    Ok(())
}
