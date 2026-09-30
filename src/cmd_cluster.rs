#![allow(clippy::unused_io_amount)]
use crate::taxonomy::Taxonomy;
use anyhow::{Context, Error};
use clap::Parser;
use rust_htslib::bam::{
    Header, HeaderView, IndexedReader, Read, Reader, Record, Writer as BamWriter, ext::BamRecordExtensions,
    index::build as build_bam_index,
};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{self, Write};

use crate::filter_read::filter_read;
use crate::k2_taxonomy::K2Taxonomy;
use crate::locus_tracker::LocusTracker;

use std::collections::{HashSet, VecDeque};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Range {
    pub start: i64,
    pub end: i64,
}

#[allow(clippy::len_without_is_empty)]
impl Range {
    pub fn len(&self) -> usize {
        assert!(self.end > self.start);
        (self.end - self.start + 1) as usize
    }
}

#[derive(Parser)]
#[command(version, arg_required_else_help = true)]
pub struct ClusterArgs {
    #[arg(short = 'o', default_value = "-")]
    pub output: String,

    #[arg(short = 'd', long, default_value = "\\t", value_parser = parse_delimiter)]
    pub delimiter: u8,

    #[arg(short = 'f', default_value_t = 0.4)]
    pub min_frac_read_aligned: f32,

    #[arg(short = 'a', default_value_t = 0.4)]
    pub min_frac_read_matched: f32,

    #[arg(short = 't', long)]
    pub taxonomy_dir: String,

    #[arg(required = true, num_args = 1..)]
    pub inputs: Vec<String>,

    #[arg(short = 'm', default_value_t = 3)]
    pub min_loci_per_call: usize,
}

pub fn parse_delimiter(value: &str) -> Result<u8, String> {
    if value == "\\t" {
        return Ok(b'\t');
    }

    let bytes = value.as_bytes();
    if bytes.len() == 1 && bytes[0].is_ascii() {
        Ok(bytes[0])
    } else {
        Err("delimiter must be \\t or a single ASCII character".to_string())
    }
}

#[inline(always)]
pub fn taxid_from_id_str(id: &str) -> Result<u32, Error> {
    if let Some((_header, meta)) = id.split_once("|taxid:") {
        let digits = meta.bytes().take_while(|b| b.is_ascii_digit()).collect::<Vec<u8>>();

        std::str::from_utf8(&digits).expect("invalid taxid string").parse::<u32>().map_err(|e| anyhow::anyhow!(e))
    } else {
        anyhow::bail!("No taxid pattern in id {}", id)
    }
}

pub fn record_get_taxid(header: &HeaderView, rec: &Record) -> Result<u32, Error> {
    let tname = std::str::from_utf8(tid2name(header, rec.tid())?)?;

    taxid_from_id_str(tname).with_context(|| format!("Header: {tname}"))
}

pub fn is_broken_pipe(error: &Error) -> bool {
    error.chain().any(|cause| {
        cause.downcast_ref::<io::Error>().is_some_and(|error| error.kind() == io::ErrorKind::BrokenPipe)
            || cause.downcast_ref::<csv::Error>().is_some_and(|error| {
                matches!(
                    error.kind(),
                    csv::ErrorKind::Io(error) if error.kind() == io::ErrorKind::BrokenPipe
                )
            })
    })
}

pub fn tid2name(header: &HeaderView, tid: i32) -> Result<&[u8], Error> {
    Ok(header.tid2name(u32::try_from(tid)?))
}

pub fn cursory_header_equivalence_check(
    other: &HeaderView,
    headerfirst: &[u8],
    headerlast: &[u8],
    targetcount: u32,
) -> Result<bool, Error> {
    Ok(other.target_count() == targetcount
        && tid2name(other, 0)? == headerfirst
        && tid2name(other, i32::try_from(targetcount.saturating_sub(1))?)? == headerlast)
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct ReadHash {
    name_hash: u64,
    source_hash: u64,
}

impl ReadHash {
    pub fn from_read_and_source(rec: &Record, source: &str) -> Self {
        let mut hash = std::hash::DefaultHasher::new();
        rec.qname().hash(&mut hash);
        let name_hash = hash.finish();

        let mut hash = std::hash::DefaultHasher::new();
        source.hash(&mut hash);
        let source_hash = hash.finish();

        Self { name_hash, source_hash }
    }
}

pub struct ReadBucket {
    pub hash: ReadHash,
    pub rec: Record,
    pub mate: Option<Record>,
}

#[derive(Default)]
pub struct ReadTracker {
    q: VecDeque<ReadBucket>,
    held: HashSet<ReadHash>,
}

impl ReadTracker {
    pub fn intake(&mut self, record: &Record, src: &str) -> impl Iterator<Item = ReadBucket> {
        let hash = ReadHash::from_read_and_source(record, src);
        let nodes = if self.held.contains(&hash) {
            self.update_mate(record, hash);
            VecDeque::new()
        } else {
            std::mem::take(&mut self.q)
        };
        nodes.into_iter()
    }

    pub fn update_mate(&mut self, rec: &Record, hash: ReadHash) {
        for node in self.q.iter_mut() {
            if node.hash == hash {
                assert!(node.mate.is_none());
                node.mate = Some(rec.clone());
                return;
            }
        }

        panic!("Attempted to update mate of non-existent read")
    }

    pub fn flush(&mut self) -> impl Iterator<Item = ReadBucket> {
        std::mem::take(&mut self.q).into_iter()
    }
}

// #[inline(always)]
// fn process_read_bucket(
//     read_bucket: &mut ReadBucket,
//     taxonomy: &mut Taxonomy,
//     k2: &mut K2Taxonomy,
//     min_frac_aligned: f32,
//     min_frac_matched: f32,
// ) {
//     // let c1 = filter_read()
// }

pub fn cluster_main(args: ClusterArgs) -> Result<(), Error> {
    let mut taxonomy = Taxonomy::from_dir(&args.taxonomy_dir).context("create taxonomy")?;
    let to_stdout = args.output == "-";
    let output: Box<dyn Write> =
        if to_stdout { Box::new(io::stdout()) } else { Box::new(File::create(&args.output).context("create output")?) };

    let mut writer = csv::WriterBuilder::new().delimiter(args.delimiter).from_writer(output);

    // we expect headers across all input files to match. We just grab the first one.
    let header =
        Header::from_template(Reader::from_path(&args.inputs[0]).expect("create header sample reader").header());

    let headerview = HeaderView::from_header(&header);

    let headercount = headerview.target_count();
    let headerfirst = tid2name(&headerview, 0)?;
    let headerlast = tid2name(&headerview, i32::try_from(headercount)?.saturating_sub(1))?;

    for input in &args.inputs {
        let mut lt = LocusTracker::new();

        let mut reader = IndexedReader::from_path(input).expect("create reader");
        let mut seq_len = 0;
        reader.set_threads(4).context("set reader threads")?;
        reader.fetch(".").context("fetch everything")?;

        if !cursory_header_equivalence_check(reader.header(), headerfirst, headerlast, headercount)? {
            anyhow::bail!(
                "File {} has a different header from the first input file {}! All BAM headers should be the same.",
                input,
                &args.inputs[0],
            )
        }

        let basename = input.rsplit_once(".").unwrap_or((input, "")).0;

        let mut failed_read_writer = std::io::BufWriter::new(
            std::fs::File::create(format!("{basename}_failed_reads.txt")).context("create failed reads file")?,
        );

        let output_read_path = format!("{basename}_bamtax.bam");
        let mut read_writer = BamWriter::from_path(&output_read_path, &header, rust_htslib::bam::Format::Bam)?;

        let mut rec = Record::new();
        let mut i = 0;
        let mut n_passed = 0;

        while let Some(result) = reader.read(&mut rec) {
            result.context("read record")?;

            i += 1;
            if i % 1000 == 0 {
                eprintln!("processed {i} records from {input}");
            }

            if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
                continue;
            }

            let cover = filter_read(&rec).context("filter record")?;
            if cover.frac_aligned() >= args.min_frac_read_aligned && cover.frac_matched() >= args.min_frac_read_matched
            {
                let taxid = record_get_taxid(&headerview, &rec).expect("get read taxid");
                if let Some(species) = taxonomy.species(taxid) {
                    seq_len += rec.seq_len();
                    n_passed += 1;

                    lt.add(&species.name, rec.tid(), Range { start: rec.pos(), end: rec.reference_end() - 1 });

                    read_writer.write(&rec)?;
                } else {
                    eprintln!("Warning: failed to find species for taxon id {}. Skipping...", taxid);

                    failed_read_writer.write(rec.qname())?;
                    failed_read_writer.write(b"\t")?;
                    failed_read_writer.write(tid2name(&headerview, rec.tid())?)?;
                    failed_read_writer.write(b"\n")?;
                }
            }
        }

        failed_read_writer.flush().expect("writer flush");

        eprintln!("processed {i} records from {input}");

        let report = lt.resolve(args.min_loci_per_call, seq_len / n_passed);
        if let Err(error) = report.serialize(&mut writer, reader.header(), basename) {
            if to_stdout && is_broken_pipe(&error) {
                return Ok(());
            }

            return Err(error);
        }

        build_bam_index(&output_read_path, None, rust_htslib::bam::index::Type::Bai, num_cpus::get() as u32)?;
    }

    Ok(())
}
