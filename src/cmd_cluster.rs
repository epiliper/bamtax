#![allow(clippy::unused_io_amount)]
use crate::taxonomy::Taxonomy;
use anyhow::{Context, Error, bail};
use clap::Parser;
use rust_htslib::bam::{
    Header, HeaderView, Read, Reader as BamReader, Record, Writer as BamWriter, ext::BamRecordExtensions,
    index::build as build_bam_index, record::Aux,
};
use std::fs::File;
use std::io::{self, BufReader, Write};

use crate::filter_read::filter_read;
use crate::k2_taxonomy::K2Taxonomy;
use crate::locus_tracker::LocusTracker;
use crate::sam_merge_buffer::{ReadBucket, SamMergeBuffer};

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
    pub bams: Vec<String>,

    #[arg(required = true, num_args = 1..)]
    pub k2_reports: Vec<String>,

    #[arg(required = true, num_args = 1..)]
    pub k2_classifications: Vec<String>,

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

const READ_SRC_TAG: &str = "RG";

pub fn get_read_src(record: &Record) -> &str {
    match record.aux(READ_SRC_TAG.as_bytes()) {
        Ok(Aux::String(str)) => str,
        Ok(_) => "0",
        _ => panic!("Read source tag {READ_SRC_TAG} overwritten!"),
    }
}

/// Given read pairs from a different alignment files (one file = one bucket), tie-break alignments
/// between mates in each bucket, then tie-break between buckets. Return the winning bucket and its hit
/// taxid if one was conclusively determined.
pub fn cluster_read_buckets(
    read_buckets: Vec<ReadBucket>,
    k2tax: &K2Taxonomy,
    headerview: &HeaderView,
    min_frac_aligned: f32,
    min_frac_matched: f32,
) -> Result<Option<(ReadBucket, u32)>, Error> {
    if read_buckets.is_empty() {
        return Ok(None);
    }

    let mut max_score = 0.0_f32;
    let mut winning_hit: u32 = 0;
    let mut winning_bucket: Option<ReadBucket> = None;

    for bucket in read_buckets {
        let r1_cov = filter_read(&bucket.rec)?;
        let r2_cov = bucket.mate.as_ref().map(filter_read).transpose()?;

        let checks =
            [Some((r1_cov, &bucket.rec, false)), r2_cov.map(|r2_cov| (r2_cov, bucket.mate.as_ref().unwrap(), true))];

        let (pair_hit, pair_score) = checks
            .iter()
            .flatten()
            .map(|(cov, read, r2)| {
                if cov.frac_aligned() >= min_frac_aligned && cov.frac_matched() >= min_frac_matched {
                    // get the taxid the mapper called
                    let map = record_get_taxid(headerview, read).expect("getting read taxonid");

                    // scale the mapper's call the number of kmers from that read that match its
                    // taxid.
                    let kmer_n = k2tax.get_kmer_calls(read.qname()).map(|row| row.get_taxid_kmer_count(map, *r2));
                    (map, cov.frac_matched() * kmer_n.unwrap_or(0) as f32)
                } else {
                    (0_u32, f32::MIN) // below alignment thresholds
                }
            })
            .max_by(|(_ahit, ascore), (_bhit, bscore)| {
                bscore.partial_cmp(ascore).expect("Invalid floating point score comparison")
            })
            .unwrap();

        // TODO: consider if epsilon is needed here for floating point comparisons?
        if pair_score > max_score {
            winning_bucket = Some(bucket);
            winning_hit = pair_hit;
            max_score = pair_score;
        }
    }

    Ok(winning_bucket.map(|bucket| (bucket, winning_hit)))
}

pub fn cluster_main(args: ClusterArgs) -> Result<(), Error> {
    if std::collections::HashSet::from([args.bams.len(), args.k2_reports.len(), args.k2_classifications.len()]).len()
        != 1
    {
        bail!("Mismatch between Number of BAMs, kraken2 reports, and kraken2 classifcations");
    }

    let taxonomy = Taxonomy::from_dir(&args.taxonomy_dir).context("create taxonomy")?;

    let to_stdout = args.output == "-";
    let output: Box<dyn Write> =
        if to_stdout { Box::new(io::stdout()) } else { Box::new(File::create(&args.output).context("create output")?) };
    let mut writer = csv::WriterBuilder::new().delimiter(args.delimiter).from_writer(output);

    for ((bam, k2_report), k2_classification) in
        std::iter::zip(&args.bams, &args.k2_reports).zip(&args.k2_classifications)
    {
        let mut merge = SamMergeBuffer::default();
        let basename = bam.rsplit_once(".").unwrap_or((bam, "")).0;
        let mut lt = LocusTracker::new();
        let k2taxon =
            K2Taxonomy::build(BufReader::new(File::open(k2_classification)?), BufReader::new(File::open(k2_report)?))?;

        let mut reader = BamReader::from_path(bam).expect("create reader");
        reader.set_threads(4).context("set reader threads")?;
        let header = reader.header();
        let mut seq_len = 0;

        let mut failed_read_writer = std::io::BufWriter::new(
            std::fs::File::create(format!("{basename}_failed_reads.txt")).context("create failed reads file")?,
        );
        let output_read_path = format!("{basename}_bamtax.bam");
        let mut read_writer =
            BamWriter::from_path(&output_read_path, &Header::from_template(header), rust_htslib::bam::Format::Bam)?;

        let mut rec = Record::new();
        let mut i = 0;
        let mut n_passed = 0;

        // copy header for now.
        let header = header.clone();

        /////////////// BEGIN CLOSURE
        let mut handle_buckets = |buckets: Vec<ReadBucket>| -> Result<(), Error> {
            if let Some((bucket, call_taxid)) = cluster_read_buckets(
                buckets,
                &k2taxon,
                &header,
                args.min_frac_read_aligned,
                args.min_frac_read_matched,
            )? {
                let reads = [Some(bucket.rec), bucket.mate];

                if let Some(species) = taxonomy.species(call_taxid) {
                    let hit_taxon = taxonomy.get(call_taxid).expect("Getting name for hit");

                    for rec in reads.iter().flatten() {
                        lt.add_and_hash(
                            &species.name,
                            &hit_taxon.name,
                            Range { start: rec.pos(), end: rec.reference_end() - 1 },
                        );

                        seq_len += rec.seq_len();
                        n_passed += 1;
                        read_writer.write(rec)?;
                    }
                } else {
                    eprintln!("Warning: failed to find species for taxon id {}. Skipping...", call_taxid);

                    for rec in reads.iter().flatten() {
                        failed_read_writer.write(rec.qname())?;
                        failed_read_writer.write(b"\t")?;
                        failed_read_writer.write(tid2name(&header, rec.tid())?)?;
                        failed_read_writer.write(b"\n")?;
                    }
                }
            }

            Ok(())
        };
        /////////////// END CLOSURE

        while let Some(result) = reader.read(&mut rec) {
            result.context("read record")?;

            if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
                continue;
            }

            if let Some(buckets) = merge.intake(&rec, get_read_src(&rec)) {
                handle_buckets(buckets)?;
            }

            i += 1;
            if i % 1000 == 0 {
                eprintln!("processed {i} records from {bam}");
            }
        }

        let remainder = merge.flush();
        if !remainder.is_empty() {
            i += remainder.len();
            handle_buckets(remainder)?;
        }

        failed_read_writer.flush().expect("writer flush");

        eprintln!("processed {i} records from {bam}");

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
