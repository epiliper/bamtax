#![allow(clippy::unused_io_amount)]
use crate::taxonomy::Taxonomy;
use anyhow::{Context, Error};
use clap::Parser;
use rust_htslib::bam::{
    HeaderView, IndexedReader, Read, Reader, Record, ext::BamRecordExtensions, record::Cigar,
};
use std::fs::File;
use std::io::{self, Write};

use crate::cmd_build_db::base_is_nonambig;
use crate::locus_tracker::LocusTracker;

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
pub struct ClusterArgs {
    #[arg(short = 'o', default_value = "-")]
    pub output: String,

    #[arg(short = 'd', long, default_value = "\\t", value_parser = parse_delimiter)]
    pub delimiter: u8,

    #[arg(short = 'f', default_value_t = 0.4)]
    pub min_frac_read_aligned: f32,

    #[arg(short = 'a', default_value_t = 0.4)]
    pub min_frac_read_matched: f32,

    #[arg(short = 't')]
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
pub fn filter_read(
    rec: &Record,
    min_frac_bases_aligned: f32,
    min_frac_bases_matched: f32,
) -> Result<bool, Error> {
    let mut bases_aligned: usize = 0;
    let mut bases_matched: usize = 0;

    let len = rec.seq_len();
    let mina = (min_frac_bases_aligned * len as f32) as usize;
    let minm = (min_frac_bases_matched * len as f32) as usize;
    let seq = rec.seq();
    let mut seqi: usize = 0;

    for op in &rec.cigar().0 {
        match op {
            Cigar::Match(_) => {
                anyhow::bail!("Wrong SAM format! Need X/= instead of M. Use SAM format 1.4+")
            }

            Cigar::Equal(len) => {
                for i in seqi..seqi + (*len as usize) {
                    if base_is_nonambig(seq[i]) {
                        bases_matched += 1;
                    }

                    bases_aligned += 1;
                }
                seqi += *len as usize;
            }

            Cigar::Diff(len) => {
                bases_aligned += *len as usize;
                seqi += *len as usize;
            }

            Cigar::Ins(len) | Cigar::SoftClip(len) | Cigar::Pad(len) => {
                seqi += *len as usize;
            }

            _ => (),
        }
    }

    Ok(bases_aligned >= mina && bases_matched >= minm)
}

#[inline(always)]
pub fn taxid_from_id_str(id: &str) -> Result<u32, Error> {
    if let Some((_header, meta)) = id.split_once("|taxid:") {
        let digits = meta
            .bytes()
            .take_while(|b| b.is_ascii_digit())
            .collect::<Vec<u8>>();

        std::str::from_utf8(&digits)
            .expect("invalid taxid string")
            .parse::<u32>()
            .map_err(|e| anyhow::anyhow!(e))
    } else {
        anyhow::bail!("No taxid pattern in id {}", id)
    }
}

fn record_get_taxid(header: &HeaderView, rec: &Record) -> Result<u32, Error> {
    let tname = std::str::from_utf8(tid2name(header, rec.tid())?)?;

    taxid_from_id_str(tname).with_context(|| format!("Header: {tname}"))
}

fn is_broken_pipe(error: &Error) -> bool {
    error.chain().any(|cause| {
        cause
            .downcast_ref::<io::Error>()
            .is_some_and(|error| error.kind() == io::ErrorKind::BrokenPipe)
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

fn cursory_header_equivalence_check(
    other: &HeaderView,
    headerfirst: &[u8],
    headerlast: &[u8],
    targetcount: u32,
) -> Result<bool, Error> {
    Ok(other.target_count() == targetcount
        && tid2name(other, 0)? == headerfirst
        && tid2name(other, i32::try_from(targetcount.saturating_sub(1))?)? == headerlast)
}

pub fn cluster_main(args: ClusterArgs) -> Result<(), Error> {
    let mut taxonomy = Taxonomy::from_dir(&args.taxonomy_dir).context("create taxonomy")?;
    let to_stdout = args.output == "-";
    let output: Box<dyn Write> = if to_stdout {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(&args.output).context("create output")?)
    };

    let mut writer = csv::WriterBuilder::new()
        .delimiter(args.delimiter)
        .from_writer(output);

    // we expect headers across all input files to match. We just grab the first one.
    let header = Reader::from_path(&args.inputs[0])
        .expect("create header sample reader")
        .header()
        .clone();

    let headercount = header.target_count();
    let headerfirst = tid2name(&header, 0)?;
    let headerlast = tid2name(&header, i32::try_from(headercount)?.saturating_sub(1))?;

    for input in &args.inputs {
        let mut lt = LocusTracker::new();

        let mut reader = IndexedReader::from_path(input).expect("create reader");
        let mut seq_len = 0;
        reader.set_threads(4).context("set reader threads")?;
        reader.fetch(".").context("fetch everything")?;

        if !cursory_header_equivalence_check(reader.header(), headerfirst, headerlast, headercount)?
        {
            anyhow::bail!(
                "File {} has a different header from the first input file {}! All BAM headers should be the same.",
                input,
                &args.inputs[0],
            )
        }

        let basename = input.rsplit_once(".").unwrap_or((input, "")).0;

        let mut failed_read_writer = std::io::BufWriter::new(
            std::fs::File::create(format!("{basename}_failed_reads.txt"))
                .context("create failed reads file")?,
        );

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

            if filter_read(&rec, args.min_frac_read_aligned, args.min_frac_read_matched)
                .context("filter record")?
            {
                let taxid = record_get_taxid(&header, &rec).expect("get read taxid");
                if let Some(species) = taxonomy.species(taxid) {
                    seq_len += rec.seq_len();
                    n_passed += 1;

                    lt.add(
                        &species.name,
                        rec.tid(),
                        Range {
                            start: rec.pos(),
                            end: rec.reference_end() - 1,
                        },
                    )
                } else {
                    eprintln!(
                        "Warning: failed to find species for taxon id {}. Skipping...",
                        taxid
                    );

                    failed_read_writer.write(rec.qname())?;
                    failed_read_writer.write(b"\t")?;
                    failed_read_writer.write(tid2name(&header, rec.tid())?)?;
                    failed_read_writer.write(b"\n")?;
                }
            }
        }

        failed_read_writer.flush().expect("writer flush");

        eprintln!("processed {i} records from {input}");

        let report = lt.resolve(args.min_loci_per_call, seq_len / n_passed);
        if let Err(error) = report.serialize(&mut writer, reader.header(), input) {
            if to_stdout && is_broken_pipe(&error) {
                return Ok(());
            }

            return Err(error);
        }
    }

    Ok(())
}
