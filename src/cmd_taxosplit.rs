use crate::assembly_dir_iterator::{AssemblyDirIterator, file_name};
use crate::cmd_build_db::{
    DBFilterArgs, DatabaseWriter, construct_assembly_to_tid_db, seq_n_bases_ambig_and_total,
};
use crate::cmd_cluster::taxid_from_id_str;
use crate::taxonomy::{Rank, Taxonomy};
use anyhow::{Context, Error};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use seq_io::fasta::{Reader as FastaReader, Record};
use std::collections::{HashMap, HashSet};
use std::io::BufReader;
use std::path::Path;

#[derive(Parser)]
pub struct TaxoSplitArgs {
    #[arg(short = 'i', long)]
    pub inputs: Vec<String>,

    #[arg(short = 'o', long)]
    pub output_prefix: String,

    #[arg(short, long = "rank")]
    pub rank: Rank,

    #[arg(short = 't', long)]
    pub taxonomy_dir: String,

    #[arg(short = 'n', long, default_value_t = 0.30)]
    pub max_frac_ambig: f32,

    #[arg(short = 'l', long, default_value_t = 100)]
    pub min_len: u32,

    #[arg(short, long = "max_seq_len", default_value_t = 500_000_000)]
    pub max_len: u32,

    #[arg(short = 't', long, num_args = 1..)]
    pub assembly_to_taxid_map: Vec<String>,

    #[arg(short, long = "reheadered_fastas", num_args = 0..)]
    pub reheadered_fastas: Vec<String>,
}

fn taxosplit_fasta(
    taxid: u32,
    file: &str,
    seen_records: &mut HashSet<String>,
    filtargs: &DBFilterArgs,
    headers_already_changed: bool,
    split_map: &mut HashMap<u32, DatabaseWriter>,
    taxonomy: &mut Taxonomy,
    rank: Rank,
) -> Result<(), Error> {
    let inner: Box<dyn std::io::Read> = {
        let f = std::fs::File::open(file)?;

        if file.ends_with(".gz") {
            Box::new(MultiGzDecoder::new(f))
        } else {
            Box::new(f)
        }
    };

    let mut reader = FastaReader::new(BufReader::new(inner));

    while let Some(rec) = reader.next() {
        let rec = rec?;
        let id = rec.id()?;

        // already processed
        if seen_records.contains(id) {
            continue;
        }

        // check for quality sequence
        let (ambig, total) = seq_n_bases_ambig_and_total(rec.seq());
        let frac_ambig = ambig as f32 / total as f32;

        if total < filtargs.min_len || frac_ambig > filtargs.max_frac_ambig {
            eprintln!(
                "skipping sequence {id}: failed filters. Length: {}. Fraction of sequence ambiguous: {}",
                total, frac_ambig
            );
            continue;
        }

        // don't add anything too big.
        if total > filtargs.max_len {
            eprintln!(
                "skipping seuqence {id}: too long! {total} > {}",
                filtargs.max_len
            );
            continue;
        }

        let (tid, id) = if headers_already_changed {
            let tid = taxid_from_id_str(id)
                .with_context(|| format!("Record in file has invalid header: {id}"))?;
            let id = id.to_string();
            (tid, id)
        } else {
            let (acc, desc) = match id.split_once(" ") {
                Some((acc, desc)) => (acc.to_string(), desc.replace(" ", "_")),
                None => (id.to_string(), "".to_string()),
            };
            (taxid, format!("{}|taxid:{}|{}", acc, taxid, desc))
        };

        if let Some(taxon) = taxonomy.lookup(tid, rank) {
            let entry = split_map.entry(taxon.tax_id).or_insert_with(|| {
                DatabaseWriter::new(format!("{}_{}", taxon.name, taxon.tax_id), true, None)
                    .expect("Failed to create db writer")
            });

            entry.write_record(&id, rec.seq())?
        }

        seen_records.insert(id.to_string());
    }

    Ok(())
}

pub fn taxosplit_main(args: TaxoSplitArgs) -> Result<(), Error> {
    let mut taxonomy = Taxonomy::from_dir(&args.taxonomy_dir)?;
    let mut seen_records: HashSet<String> = HashSet::new();
    let mut splitmap: HashMap<u32, DatabaseWriter> = HashMap::new();

    let filtargs = DBFilterArgs {
        min_len: args.min_len,
        max_len: args.max_len,
        max_frac_ambig: args.max_frac_ambig,
    };

    let assembly_tid_map = construct_assembly_to_tid_db(&args.assembly_to_taxid_map)?;

    for file in &args.reheadered_fastas {
        if Path::new(file).is_file() {
            taxosplit_fasta(
                0,
                file,
                &mut seen_records,
                &filtargs,
                true,
                &mut splitmap,
                &mut taxonomy,
                args.rank,
            )?;
        } else {
            let entries = Path::new(file)
                .read_dir()
                .with_context(|| format!("error reading dir {}", file))?;

            for f in entries {
                let f = f?;
                if let Some(fname) = file_name(&f.path())? {
                    if !fname.contains(".fna")
                        && !fname.contains(".fa")
                        && !fname.contains(".fasta")
                    {
                        continue;
                    }

                    taxosplit_fasta(
                        0,
                        fname,
                        &mut seen_records,
                        &filtargs,
                        true,
                        &mut splitmap,
                        &mut taxonomy,
                        args.rank,
                    )?;
                }
            }
        }
    }

    for input in args.inputs {
        let mut iterator = AssemblyDirIterator::new(input)?;

        while let Some((assembly, fasta)) = iterator.next_item()? {
            let taxid = assembly_tid_map.get(&assembly).with_context(|| {
                format!("Assembly {assembly} not found in assembly to taxon id map!")
            })?;

            taxosplit_fasta(
                *taxid,
                &fasta,
                &mut seen_records,
                &filtargs,
                true,
                &mut splitmap,
                &mut taxonomy,
                args.rank,
            )?;
        }
    }

    for writer in splitmap.values_mut() {
        writer.flush()?;
    }

    Ok(())
}
