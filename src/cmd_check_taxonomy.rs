use crate::cmd_cluster::parse_delimiter;
use crate::taxonomy::Taxonomy;
use anyhow::{Context, Error};
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
    query: u32,
    query_name: &'a str,
    representative: u32,
    representative_name: &'a str,
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

    let mut query_reader: BufReader<Box<dyn Read>> =
        if let Some(input) = args.query_headers.as_ref() {
            BufReader::new(Box::new(std::fs::File::open(input)?))
        } else {
            BufReader::new(Box::new(std::io::stdin()))
        };

    let mut ref_reader = BufReader::new(std::fs::File::open(&args.database_headers)?);

    let taxonomy = Taxonomy::from_dir(args.taxonomy_dir)?;

    let mut line = String::new();

    let mut db_contains: HashSet<u32> = HashSet::new();

    while ref_reader.read_line(&mut line)? > 0 {
        let tid = line
            .trim()
            .parse::<u32>()
            .with_context(|| format!("invalid taxon id in {}: {}", &args.database_headers, line))?;
        line.clear();

        let mut count = 0;
        taxonomy.descendants(tid).flatten().for_each(|t| {
            count += 1;
            db_contains.insert(*t);
        });

        // we queried a leaf node
        if taxonomy.get(tid).is_some() && count == 0 {
            db_contains.insert(tid);
        }
    }

    while query_reader.read_line(&mut line)? > 0 {
        let tid = line
            .trim()
            .parse::<u32>()
            .with_context(|| format!("invalid taxon id in query: {}", line))?;
        line.clear();

        let query_name = taxonomy
            .get(tid)
            .map(|d| d.name.as_str())
            .unwrap_or("not in taxonomy");

        let descendants =
            std::iter::chain(taxonomy.descendants(tid).flatten(), std::iter::once(&tid))
                .filter(|d| db_contains.contains(d))
                .copied();

        let mut count = 0;

        for d in descendants {
            count += 1;

            let representative_name = &taxonomy.get(d).unwrap().name;

            csv_writer.serialize(DatabaseTaxonQueryRow {
                query: tid,
                query_name,
                representative: d,
                representative_name,
            })?;
        }

        if count == 0 {
            csv_writer.serialize(DatabaseTaxonQueryRow {
                query: tid,
                query_name,
                representative: 0,
                representative_name: "MISSING",
            })?;
        }
    }

    Ok(())
}
