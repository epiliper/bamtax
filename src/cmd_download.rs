#![allow(clippy::unused_io_amount)]

use clap::Parser;
use std::io::{Read, Write, BufRead, BufReader};
use anyhow::Error;
use crate::{taxonomy::Taxonomy, cmd_report_species::open_file_reader};
use std::process::{Command, Stdio};
use std::collections::HashSet;

#[derive(Parser)]
pub struct DownloadArgs {
    /// File of taxon IDs, should be one per line.
    #[arg(short = 'i', long, default_value_t = "-".to_string())]
    input: String,

    #[arg(short = 'o', long, default_value_t = "-".to_string())]
    output: String,

    #[arg(short = 't', long)]
    taxonomy_dir: String,

    /// api key for NCBI. Required since downloading a metagenome-level DB from NCBI without one
    /// is slower than dying of old age.
    #[arg(short, long = "api-key")]
    api_key: String,

    #[arg(short, long = "species_blacklist")]
    species_blacklist: Option<String>,
}

fn nucleotide_query(taxid: u32) -> String {
    format!("txid{taxid}[Organism:exp] AND (\"complete genome\" OR \"partial genome\")")
}

pub fn download_genomes_for_taxid(taxid: u32) -> Result<Vec<u8>, Error> {
    let search = Command::new("esearch").stdout(Stdio::piped()).args(["-db", "nucleotide", "-query", &nucleotide_query(taxid)]).spawn()?; 
    let fetch = Command::new("efetch").stdin(search.stdout.unwrap()).args(["-format", "acc"]).output()?;
    if !fetch.status.success() {
        eprintln!("error fetching for {taxid}: {}", std::str::from_utf8(&fetch.stderr)?);
    }
    std::thread::sleep(std::time::Duration::from_secs(1));
    Ok(fetch.stdout)
}

pub fn download_main(args: DownloadArgs) -> Result<(), Error> {
    let input: Box<dyn Read> = if args.input == "-" {
        Box::new(std::io::stdin().lock())
    } else {
        open_file_reader(&args.input)?
    };

    let reader = BufReader::new(input);
    let mut taxo = Taxonomy::from_dir(args.taxonomy_dir)?;

    let blacklist: HashSet<u32> = if let Some(blacklist) = args.species_blacklist {
        let iter = BufReader::new(std::fs::File::open(blacklist)?).lines().map(|l| l.unwrap().trim().parse::<u32>().unwrap());
        HashSet::from_iter(iter)
    } else {
        HashSet::new()
    };

    // track refernce-level taxon ids to make sure we aren't doing repeated work
    let mut seen: HashSet<u32> = HashSet::new();

    for line in reader.lines() {
        let tid = line?.trim().parse::<u32>()?;

        if let Some(species) = taxo.species(tid) && !blacklist.contains(&species.tax_id) && !seen.contains(&tid) {
            seen.insert(tid);

            let bytes = download_genomes_for_taxid(tid)?;
            if bytes.is_empty() { continue; }
            std::io::stdout().write(&bytes)?;
        }
    }

    Ok(())
}
