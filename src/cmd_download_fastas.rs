#![allow(clippy::unused_io_amount)]

use clap::Parser;
use crate::cmd_download_accs::NCBIRequestTracker;
use std::io::{BufRead, BufWriter, BufReader, Read, Write};
use anyhow::Error;
use std::process::Command;
use std::sync::{Arc, Mutex};

#[derive(Parser)]
pub struct DownloadFastasArgs {
    // download input, should be TSV of TAXID, NUCLEOTIDE_ACCESSION
    #[arg(short = 'i', long)]
    input: String,

    #[arg(short = 'o', long, default_value_t = "-".to_string())]
    output_file: String,

    #[arg(short = 'n', long, default_value_t = 0.30)]
    pub max_frac_ambig: f32,

    #[arg(short, long = "max_seq_len", default_value_t = 500_000_000)]
    pub max_len: u64,

    #[arg(short = 'l', long, default_value_t = 100)]
    pub min_len: u64,
}

fn download_fasta_thread<W: std::io::Write + Send + 'static>(acc: String, tid: String, writer: Arc<Mutex<W>>) -> Result<(), Error> {

    std::thread::spawn(move || {
        let fetch = Command::new("efetch").args(["-db", "nucleotide", "-id", &acc, "-format", "fasta"]).output().unwrap();
        let mut lines = fetch.stdout.lines();

        // reheader
        let header = lines.next().unwrap().unwrap();
        assert!(header.starts_with(">"));

        let (first, second) = header.trim().split_once(" ").unwrap_or((&header, ""));
        let mut lock = writer.lock().unwrap();

        lock.write(first.as_bytes()).unwrap();
        lock.write(b"|").unwrap();
        lock.write(tid.as_bytes()).unwrap();
        lock.write(b"|").unwrap();
        lock.write(second.as_bytes()).unwrap();
        lock.write(b"\n").unwrap();
    });

    Ok(())
}

fn download_fasta<W: std::io::Write + Send>(acc: &str, tid: &str, writer: Arc<Mutex<W>>) -> Result<(), Error> {
    todo!();
}

pub fn cmd_download_fastas_main(args: DownloadFastasArgs) -> Result<(), Error> {
    Ok(())
}

