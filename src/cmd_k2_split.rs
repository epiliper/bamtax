use crate::k2_taxonomy::K2Taxonomy;
use anyhow::Error;
use clap::Parser;
use flate2::read::MultiGzDecoder;
use gzp::{
    deflate::Mgzip,
    par::compress::{Compression, ParCompress, ParCompressBuilder},
};
use seq_io::fastq::{Reader as FastqReader, Record, write_to as write_fastq};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
};

/////////////////////////////////////////////////////////////////////
pub enum FastqOut<'a> {
    Gzip(ParCompress<'a, Mgzip, BufWriter<File>>),
    Plaintext(BufWriter<File>),
}

impl<'a> std::io::Write for FastqOut<'a> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Gzip(w) => w.write(buf),
            Self::Plaintext(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Gzip(w) => w.flush(),
            Self::Plaintext(w) => w.flush(),
        }
    }
}

impl<'a> FastqOut<'a> {
    pub fn from_file(f: &str) -> Result<Self, Error> {
        let inner = BufWriter::new(File::create(f)?);

        if f.ends_with(".gz") {
            Ok(Self::Gzip(
                ParCompressBuilder::new()
                    .compression_level(Compression::new(4))
                    .num_threads(num_cpus::get())?
                    .from_writer(inner),
            ))
        } else {
            Ok(Self::Plaintext(inner))
        }
    }
}

pub enum FastqIn {
    Gzip(MultiGzDecoder<BufReader<File>>),
    Plaintext(BufReader<File>),
}

impl std::io::Read for FastqIn {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::Gzip(r) => r.read(buf),
            Self::Plaintext(r) => r.read(buf),
        }
    }
}

impl FastqIn {
    pub fn from_file(file: &str) -> Result<Self, Error> {
        let inner = BufReader::new(File::open(file)?);
        if file.ends_with(".gz") { Ok(Self::Gzip(MultiGzDecoder::new(inner))) } else { Ok(Self::Plaintext(inner)) }
    }
}

//////////////////////////////////////////////////////////////////

#[derive(Parser)]
pub struct K2SplitArgs {
    /// general k2 report
    #[arg(long)]
    pub report: String,

    /// Per-read k2 classifications and kmer lists
    #[arg(long)]
    pub classification: String,

    /// Optional output prefix to use for output report and fastqs.
    #[arg(short, long)]
    pub output_prefix: Option<String>,

    /// input fastq with reads mentioned in report and classification
    #[arg(short, long, num_args = 1..)]
    pub fastq: Vec<String>,

    /// whether or not to gzip output fastqs
    #[arg(short, long)]
    pub gzip: bool,
}

const HIGHEST_RANK_WE_PULL_READS_FROM: &str = "R2";
const SUB_HIGHEST_RANK_SPLIT_BY: &str = "G"; // genus

pub fn k2_split_main(args: K2SplitArgs) -> Result<(), Error> {
    let k2taxonomy = K2Taxonomy::build(
        BufReader::new(File::open(&args.classification)?),
        BufReader::new(File::open(&args.report)?),
    )?;

    let output_prefix = args
        .output_prefix
        .unwrap_or(args.report.split_once(".kreport2").expect("report should end in .kreport2").0.to_string());

    let mut reads_to_group: HashMap<String, u32> = HashMap::new();
    let mut group_to_out: HashMap<u32, FastqOut> = HashMap::new();

    let r2_iter = k2taxonomy.rank_nodes(HIGHEST_RANK_WE_PULL_READS_FROM);
    let mut report = BufWriter::new(File::create(format!("{output_prefix}_groups.tsv"))?);

    // 1. record which read names belong under which clusters, record genuses detected in report.
    for r2 in r2_iter {
        let fastq_name = format!("{output_prefix}_{}.fastq{}", r2.row.name, if args.gzip { ".gz" } else { "" });

        let writer = FastqOut::from_file(&fastq_name)?;
        group_to_out.insert(r2.row.taxid, writer);

        k2taxonomy.get_all_reads_under_node(r2.row.taxid).for_each(|r| {
            reads_to_group.insert(r.read_name.clone(), r2.row.taxid);
        });

        k2taxonomy.descendants_taxa(r2.row.taxid).filter(|l| l.row.rank == SUB_HIGHEST_RANK_SPLIT_BY).for_each(|g| {
            writeln!(report, "{fastq_name}\t{}\t{}", g.row.name, g.row.taxid).expect("writing to report");
        })
    }

    report.flush()?;

    // 2. write reads to their separate files.
    for fastq in args.fastq {
        let mut input = FastqReader::new(FastqIn::from_file(&fastq)?);
        while let Some(rec) = input.next() {
            let rec = rec?;
            let id = rec.id()?;

            if let Some(group) = reads_to_group.get(id) {
                let writer = group_to_out.get_mut(group).expect("getting writer");
                write_fastq(writer, id.as_bytes(), rec.seq(), rec.qual())?;
            }
        }
    }

    for writer in group_to_out.values_mut() {
        writer.flush()?;
    }

    Ok(())
}
