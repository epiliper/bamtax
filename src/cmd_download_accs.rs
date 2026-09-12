#![allow(clippy::unused_io_amount)]

use crate::{cmd_report_species::open_file_reader, taxonomy::Taxonomy};
use anyhow::Error;
use clap::Parser;
use std::collections::HashSet;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

#[derive(Parser)]
pub struct DownloadAccessionsArgs {
    /// File of taxon IDs, should be one per line.
    #[arg(short = 'i', long)]
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

// I need some state to track requests so we can sleep() to avoid the dreaded 429
pub struct NCBIRequestTracker {
    nrequests: usize,
    time_last_clear: Instant,
}

impl NCBIRequestTracker {
    // rps when an API key is used.
    const NCBI_REQUESTS_PER_SECOND: usize = 10;
    const COOLDOWN: Duration = Duration::from_secs(1);

    #[inline(always)]
    pub fn tick(&mut self) {
        if self.time_last_clear >= Instant::now() - Self::COOLDOWN {
            if self.nrequests >= Self::NCBI_REQUESTS_PER_SECOND {
                std::thread::sleep(Self::COOLDOWN);
            }

            self.nrequests = 0;
            self.time_last_clear = Instant::now();
        }

        self.nrequests += 1;
    }

    #[inline(always)]
    pub fn default() -> Self {
        Self {
            nrequests: 0,
            time_last_clear: Instant::now(),
        }
    }
}

pub fn download_accs_for_taxid(
    taxid: u32,
    tracker: &mut NCBIRequestTracker,
    api_key: &str,
) -> Result<Vec<u8>, Error> {
    tracker.tick();
    let search = Command::new("esearch")
        .stdout(Stdio::piped())
        .args(["-db", "nucleotide", "-query", &nucleotide_query(taxid)])
        .env("NCBI_API_KEY", api_key)
        .spawn()?;

    tracker.tick();
    let fetch = Command::new("efetch")
        .stdin(search.stdout.unwrap())
        .args(["-format", "acc"])
        .env("NCBI_API_KEY", api_key)
        .output()?;
    if !fetch.status.success() {
        eprintln!(
            "error fetching for {taxid}: {}",
            std::str::from_utf8(&fetch.stderr)?
        );
    }

    Ok(fetch.stdout)
}

pub fn download_accs_main(args: DownloadAccessionsArgs) -> Result<(), Error> {
    let (totallines, input): (usize, Box<dyn Read>) = (
        BufReader::new(open_file_reader(&args.input)?)
            .lines()
            .count(),
        open_file_reader(&args.input)?,
    );

    let reader = BufReader::new(input);

    let mut writer: BufWriter<Box<dyn Write>> = {
        if args.output == "-" {
            BufWriter::new(Box::new(std::io::stdout().lock()))
        } else {
            BufWriter::new(Box::new(std::fs::File::create(&args.output)?))
        }
    };

    let mut taxo = Taxonomy::from_dir(args.taxonomy_dir)?;

    let blacklist: HashSet<u32> = if let Some(blacklist) = args.species_blacklist {
        let iter = BufReader::new(std::fs::File::open(blacklist)?)
            .lines()
            .map(|l| l.unwrap().trim().parse::<u32>().unwrap());
        HashSet::from_iter(iter)
    } else {
        HashSet::new()
    };

    // track reference-level taxon ids to make sure we aren't doing repeated work
    let mut seen: HashSet<u32> = HashSet::new();

    let mut tracker = NCBIRequestTracker {
        nrequests: 0,
        time_last_clear: Instant::now(),
    };

    for (i, line) in reader.lines().enumerate() {
        let tid = line?.trim().parse::<u32>()?;

        if let Some(species) = taxo.species(tid)
            && !blacklist.contains(&species.tax_id)
            && !seen.contains(&tid)
        {
            seen.insert(tid);
            let tidstr = tid.to_string();

            if let Ok(bytes) = download_accs_for_taxid(tid, &mut tracker, &args.api_key) {
                let output_lines = bytes.lines();

                for o in output_lines {
                    writer.write(tidstr.as_bytes())?;
                    writer.write(b"\t")?;
                    writer.write(o?.as_bytes())?;
                    writer.write(b"\n")?;
                }
            }
        }

        eprint!("Processed {} of {totallines} lines\r", i + 1);
    }

    writer.flush()?;
    Ok(())
}
