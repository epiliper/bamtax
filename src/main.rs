mod locus_tracker;
pub mod taxonomy;

use anyhow::{Context, Error};
use clap::Parser;
use rust_htslib::bam::{
    HeaderView, IndexedReader, Read, Record, ext::BamRecordExtensions, record::Cigar,
};
use taxonomy::{Taxon, Taxonomy};

use locus_tracker::LocusTracker;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Range {
    pub start: i64,
    pub end: i64,
}

#[derive(Parser)]
struct Args {
    #[arg(short = 'i')]
    pub input: String,

    #[arg(short = 'o')]
    pub output: String,

    #[arg(short = 'f', default_value_t = 0.5)]
    pub min_frac_read_aligned: f32,

    #[arg(short = 'a', default_value_t = 0.6)]
    pub min_frac_read_matched: f32,

    #[arg(short = 't')]
    pub taxonomy_dir: String,
}

fn filter_read(
    rec: &Record,
    min_frac_bases_aligned: f32,
    min_frac_bases_matched: f32,
) -> Result<bool, Error> {
    let mut bases_aligned: usize = 0;
    let mut bases_matched: usize = 0;

    let len = rec.seq_len();
    let mina = (min_frac_bases_aligned * len as f32) as usize;
    let minm = (min_frac_bases_matched * len as f32) as usize;

    for op in &rec.cigar().0 {
        match op {
            Cigar::Match(_) => {
                anyhow::bail!("Wrong SAM format! Need X/= instead of M. Use SAM format 1.4+")
            }

            Cigar::Equal(len) => {
                bases_aligned += *len as usize;
                bases_matched += *len as usize
            }

            Cigar::Diff(len) => bases_aligned += *len as usize,

            _ => (),
        }
    }

    Ok(bases_aligned >= mina && bases_matched >= minm)
}

fn record_get_taxid(header: &HeaderView, rec: &Record) -> Result<u32, Error> {
    let tname = header
        .target_names()
        .get(rec.tid() as usize)
        .with_context(|| format!("Invalid TID {} not in header", rec.tid()))
        .map(|bytes| std::str::from_utf8(bytes).unwrap())?;

    if let Some((_header, meta)) = tname.split_once("|taxid:") {
        let (num, _) = meta.split_once(" ").unwrap();
        num.parse::<u32>().map_err(|e| anyhow::anyhow!(e))
    } else {
        anyhow::bail!("no taxid pattern in read header!");
    }
}

fn main() {
    let args = Args::parse();

    let taxonomy = Taxonomy::from_dir(args.taxonomy_dir).expect("create taxonomy");

    let mut lt = LocusTracker::new();
    let mut reader = IndexedReader::from_path(args.input).expect("create reader");
    reader.set_threads(4).expect("set reader threads");

    let mut rec = Record::new();

    while let Some(result) = reader.read(&mut rec) {
        result.expect("read record");

        if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
            continue;
        }

        if filter_read(&rec, args.min_frac_read_aligned, args.min_frac_read_matched)
            .expect("filter record")
        {
            let taxid = record_get_taxid(reader.header(), &rec).expect("get read taxid");

            let species = taxonomy
                .species(taxid)
                .with_context(|| format!("failed to find species for taxid {}", taxid))
                .expect("failed to find species");

            lt.add(
                &species.name,
                Range {
                    start: rec.pos(),
                    end: rec.reference_end() - 1,
                },
            )
        }
    }

    lt.resolve();
}
