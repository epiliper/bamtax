
#![allow(clippy::unused_io_amount)]
use clap::Parser;
use anyhow::{Error, Context};
use std::io::{BufRead, BufReader, Read, Write, BufWriter};
use crate::taxonomy::Taxonomy;
use crate::cmd_cluster::taxid_from_id_str;
use flate2::read::MultiGzDecoder;
use std::collections::HashSet;

#[derive(Parser)]
pub struct ReportSpeciesArgs {
    // input fastas to report. Should be one file per line. Default is stdin.
    #[arg(short = 'i', long, default_value_t = "-".to_string())]
    pub input: String,

    #[arg(short = 't', long)]
    pub taxonomy_dir: String,

    #[arg(short = 'o', long, default_value_t = "-".to_string())]
    pub output: String,
}

pub fn open_file_reader(p: &str) -> Result<Box<dyn Read>, Error> {
        let inner = std::fs::File::open(p)?;

        let inner: Box<dyn Read> = if p.ends_with(".gz") {
            Box::new(MultiGzDecoder::new(inner))
        } else {
            Box::new(inner)
        };

        Ok(inner)
}

pub fn report_species_main(args: ReportSpeciesArgs) -> Result<(), Error> {

    let input: Box<dyn Read> = if &args.input == "-" {
        Box::new(std::io::stdin().lock())
    } else {
        open_file_reader(&args.input)?
    };

    let output: Box<dyn Write> = if &args.output == "-" {
        Box::new(std::io::stdout().lock())
    } else {
        Box::new(std::fs::File::create(&args.output)?)
    };

    let mut writer = BufWriter::new(output);

    let mut reader = BufReader::new(input);

    let mut taxo = Taxonomy::from_dir(&args.taxonomy_dir).context("create taxonomy")?;
    let mut seen: HashSet<u32> = HashSet::new();

    for f in reader.lines() {
        let f = f?;
        let mut inf = BufReader::new(open_file_reader(&f)?);
        eprintln!("Checking {f}...");

        for linebuf in inf.lines() {
            let l = linebuf?;
            if l.starts_with(">") { 
                let tid = taxid_from_id_str(&l[1..])?;
                if let Some(species) = taxo.species(tid) {

                    if seen.contains(&species.tax_id) { continue; }

                    writer.write(species.tax_id.to_string().as_bytes())?;
                    writer.write(b"\n")?;
                    seen.insert(species.tax_id);
                }
            }

        }
    }

    writer.flush()?;

    Ok(())

}
