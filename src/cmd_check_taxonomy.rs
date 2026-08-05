use crate::cmd_cluster::parse_delimiter;
use crate::taxonomy::Taxonomy;
use anyhow::Error;
use clap::Parser;
use serde::Serialize;
use std::collections::HashSet;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

#[derive(Parser)]
pub struct CheckTaxonomyArgs {
    /// taxon ID file to check for representatives in database, one per line
    #[arg(short = 'i')]
    pub query_headers: Option<String>,

    /// list of headers in database
    #[arg(short = 'r')]
    pub database_headers: String,

    #[arg(short = 'd', long, default_value = "\\t", value_parser = parse_delimiter)]
    pub delimiter: u8,

    /// output file. defaults to stdout
    #[arg(short = 'o')]
    pub output: Option<String>,

    #[arg(short = 't')]
    pub taxonomy_dir: String,
}

#[derive(Serialize)]
struct DatabaseTaxonQueryRow<'a> {
    tax_id: u32,
    reps: &'a [u32],
}

pub fn check_taxonomy_main(args: CheckTaxonomyArgs) -> Result<(), Error> {
    let outwriter: BufWriter<Box<dyn Write>> = if let Some(out) = args.output {
        BufWriter::new(Box::new(std::fs::File::create(&out)?))
    } else {
        BufWriter::new(Box::new(std::io::stdout()))
    };

    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(args.delimiter)
        .from_writer(outwriter);

    let mut query_reader: BufReader<Box<dyn Read>> = if let Some(input) = args.query_headers {
        BufReader::new(Box::new(std::fs::File::open(input)?))
    } else {
        BufReader::new(Box::new(std::io::stdin()))
    };

    let mut ref_reader = BufReader::new(std::fs::File::open(args.database_headers)?);

    let taxonomy = Taxonomy::from_dir(args.taxonomy_dir)?;

    let mut line = String::new();

    let mut db_contains: HashSet<u32> = HashSet::new();

    while ref_reader.read_line(&mut line)? > 0 {
        let tid = line.trim().parse::<u32>()?;
        line.clear();

        taxonomy.descendants(tid).iter().for_each(|t| {
            db_contains.insert(*t);
        });

        db_contains.insert(tid);
    }

    while query_reader.read_line(&mut line)? > 0 {
        let tid = line.trim().parse::<u32>()?;
        line.clear();

        let descendants = std::iter::chain(taxonomy.descendants(tid).iter(), std::iter::once(&tid))
            .filter(|d| db_contains.contains(d))
            .copied();

        csv_writer.serialize(DatabaseTaxonQueryRow {
            tax_id: tid,
            reps: descendants.collect::<Vec<u32>>().as_slice(),
        })?;
    }

    Ok(())
}
