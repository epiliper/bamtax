#![allow(clippy::unused_io_amount)]
use crate::cmd_cluster::filter_read;
use anyhow::Error;
use clap::Parser;
use rust_htslib::bam::{IndexedReader, Read, Record};
use std::io::{BufWriter, Write};

#[derive(Parser)]
pub struct EmitNamesArgs {
    #[arg(short = 'i')]
    pub input: String,

    #[arg(short = 'o')]
    pub output: String,

    #[arg(short = 'f', default_value_t = 0.4)]
    pub min_frac_read_aligned: f32,

    #[arg(short = 'a', default_value_t = 0.4)]
    pub min_frac_read_matched: f32,
}

pub fn emit_names_main(args: EmitNamesArgs) -> Result<(), Error> {
    let mut reader = IndexedReader::from_path(&args.input).expect("create reader");
    reader.set_threads(4).expect("set reader threads");
    reader.fetch(".").expect("fetch everything");

    let mut tnames: Vec<Vec<u8>> = Vec::with_capacity(reader.header().target_count() as usize);
    for name in reader.header().target_names() {
        tnames.push(name.to_vec());
    }

    let mut rec = Record::new();
    let mut i = 0;
    let mut n_passed = 0;

    let mut writer = BufWriter::new(std::fs::File::create(args.output)?);

    while let Some(result) = reader.read(&mut rec) {
        result?;

        i += 1;
        if i % 1000 == 0 {
            eprintln!("processed {i} records from {}", &args.input);
        }

        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }

        if !filter_read(&rec, args.min_frac_read_aligned, args.min_frac_read_matched)? {
            continue;
        }

        n_passed += 1;

        writer.write(rec.qname())?;
        writer.write(b"\n")?;
    }

    writer.flush()?;

    eprintln!(
        "{n_passed} / {i} records in {} were determined to be legitimate maps",
        args.input
    );

    Ok(())
}
