#![allow(clippy::unused_io_amount)]

use crate::cmd_build_db::{DBFilterArgs, DatabaseWriter, process_metadata_fasta};
use crate::cmd_download_accs::NCBIRequestTracker;
use anyhow::Error;
use clap::Parser;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

#[derive(Parser)]
pub struct DownloadFastasArgs {
    // download input, should be TSV of TAXID, NUCLEOTIDE_ACCESSION
    #[arg(short = 'i', long)]
    input: String,

    #[arg(short = 'o', long, default_value_t = "-".to_string())]
    output: String,

    #[arg(short = 'n', long, default_value_t = 0.30)]
    pub max_frac_ambig: f32,

    #[arg(short, long = "max_seq_len", default_value_t = 500_000_000)]
    pub max_len: u32,

    #[arg(short = 'l', long, default_value_t = 100)]
    pub min_len: u32,
}

fn download_fasta_thread(
    acc: String,
    tid: u32,
    writer: Arc<Mutex<DatabaseWriter>>,
    seen: Arc<Mutex<HashSet<String>>>,
    args: DBFilterArgs,
) -> JoinHandle<()> {
    std::thread::spawn(move || {
        let fetch = Command::new("efetch")
            .args(["-db", "nucleotide", "-id", &acc, "-format", "fasta"])
            .output()
            .unwrap();

        if !fetch.status.success() {
            eprintln!(
                "Downloading for acc {acc} failed: {}",
                std::str::from_utf8(&fetch.stderr).unwrap()
            );
        }

        let lock = &mut writer.lock().unwrap();
        let seen = &mut seen.lock().unwrap();

        process_metadata_fasta(
            tid,
            fetch.stdout.as_slice(),
            lock,
            seen,
            &args,
            false,
            false,
        )
        .expect("processing fasta");

        // // reheader
        // let header = lines.next().unwrap().unwrap();
        // assert!(header.starts_with(">"));

        // let (first, second) = header.trim().split_once(" ").unwrap_or((&header, ""));
        // let mut lock = writer.lock().unwrap();

        // lock.write(first.as_bytes()).unwrap();
        // lock.write(b"|").unwrap();
        // lock.write(tid.as_bytes()).unwrap();
        // lock.write(b"|").unwrap();
        // lock.write(second.as_bytes()).unwrap();
        // lock.write(b"\n").unwrap();
    })
}

pub fn download_fastas_main(args: DownloadFastasArgs) -> Result<(), Error> {
    let input_lines = BufReader::new(std::fs::File::open(&args.input)?)
        .lines()
        .count();

    let reader = BufReader::new(std::fs::File::open(&args.input)?);

    let output = Arc::new(Mutex::new(DatabaseWriter::new(
        args.output.clone(),
        args.output.ends_with(".gz"),
        None,
    )?));

    let mut tracker = NCBIRequestTracker::default();
    let seen: Arc<Mutex<HashSet<String>>> = Arc::new(Mutex::new(HashSet::new()));

    let args = DBFilterArgs {
        min_len: args.min_len,
        max_len: args.max_len,
        max_frac_ambig: args.max_frac_ambig,
    };

    const MAX_THREADS: usize = 100;
    let mut joins: Vec<JoinHandle<()>> = Vec::with_capacity(MAX_THREADS);

    for (i, line) in reader.lines().enumerate() {
        let l = line?;
        let (taxid, acc) = l.split_once("\t").unwrap();

        tracker.tick();

        let join = download_fasta_thread(
            acc.to_string(),
            taxid.parse::<u32>()?,
            Arc::clone(&output),
            Arc::clone(&seen),
            args.clone(),
        );

        joins.push(join);

        if joins.len() >= MAX_THREADS {
            joins.drain(..).for_each(|j| j.join().unwrap())
        }

        eprintln!("Processed {} of {input_lines}", i + 1);
    }

    Ok(())
}
