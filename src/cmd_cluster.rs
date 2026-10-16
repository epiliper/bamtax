#![allow(clippy::unused_io_amount)]
use crate::taxonomy::Taxonomy;
use anyhow::{Context, Error, bail};
use clap::Parser;
use rust_htslib::bam::{Header, HeaderView, Read, Reader as BamReader, Record, Writer as BamWriter, record::Aux};
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;

use crate::html_report::{HtmlSample, write_html_report};
use crate::k2_taxonomy::K2Taxonomy;
use crate::locus_tracker::LocusTracker;
use crate::sam_merge_buffer::{ClusterAlignment, ReadBucket, SamMergeBuffer};

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

const CLUSTER_HELP: &str = "
Generate species calls for one or more BAMs, optionally using output from kraken-style
classifiers to confirm read mapping specificity.

IMPORTANT: if specifying k2 files, the order of the k2 files must match the order of the input 
BAM files; these are not sorted internally.                                                    

Good:   --bams 1.bam 2.bam  --k2-reports 1.txt 2.txt --k2-classifications 1.txt 2.txt          
Bad:    --bams 1.bam 2.bam  --k2-reports 2.txt 1.txt --k2-classifications 1.txt 2.txt";

#[derive(Parser)]
#[command(version, arg_required_else_help = true, about = CLUSTER_HELP, long_about = CLUSTER_HELP, arg_required_else_help = true)]
pub struct ClusterArgs {
    /// output name of report
    #[arg(short = 'o', default_value = "-")]
    pub output: String,

    /// Write a combined, standalone HTML report
    #[arg(long, value_name = "PATH")]
    pub html_report: Option<String>,

    /// KronaTools ktImportTaxonomy executable
    #[arg(long, value_name = "PATH", requires = "html_report")]
    pub krona_command: Option<String>,

    /// Prepared KronaTools taxonomy database directory
    #[arg(long, value_name = "PATH", requires = "html_report")]
    pub krona_taxonomy: Option<String>,

    /// Delimiter to use in report
    #[arg(short = 'd', long, default_value = "\\t", value_parser = parse_delimiter)]
    pub delimiter: u8,

    /// Minimum fraction of considered reads that should be aligned, not matched
    #[arg(short = 'f', default_value_t = 0.4)]
    pub min_frac_read_aligned: f32,

    /// Minimum fraction of considered reads that should be matched to ref
    #[arg(short = 'a', default_value_t = 0.4)]
    pub min_frac_read_matched: f32,

    /// Directory with NCBI .dmp taxonomy files
    #[arg(short = 't', long)]
    pub taxonomy_dir: String,

    /// Input BAM files
    #[arg(short, long, required = true, num_args = 1..)]
    pub bams: Vec<String>,

    /// Input K2 report files. Using this option, k2-classifications are required
    #[arg(long, required = false, num_args = 1..)]
    pub k2_reports: Option<Vec<String>>,

    /// Input k2 per-read classifications. Should be with k2-reports
    #[arg(long, required = false, num_args = 1..)]
    pub k2_classifications: Option<Vec<String>>,

    /// Minimum # of loci per species to make a call
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

fn canonical_destination(path: &str) -> Result<std::path::PathBuf, Error> {
    let path = Path::new(path);
    let absolute = if path.is_absolute() { path.to_owned() } else { std::env::current_dir()?.join(path) };
    if absolute.exists() {
        return Ok(std::fs::canonicalize(absolute)?);
    }

    let parent = absolute.parent().context("output path has no parent directory")?;
    let name = absolute.file_name().context("output path has no file name")?;
    Ok(std::fs::canonicalize(parent)?.join(name))
}

#[inline(always)]
pub fn taxid_from_id_str(id: &str) -> Result<u32, Error> {
    let split = if let Some((_header, meta)) = id.split_once("|taxid:") {
        meta
    } else if let Some((_header, meta)) = id.split_once("kraken:taxid|") {
        meta
    } else {
        bail!("No taxid pattern in id {}", id);
    };

    let digits = split.bytes().take_while(|b| b.is_ascii_digit()).collect::<Vec<u8>>();
    std::str::from_utf8(digits.as_slice()).expect("invalid taxid string").parse::<u32>().map_err(|e| anyhow::anyhow!(e))
}

#[inline(always)]
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

#[inline(always)]
pub fn tid2name(header: &HeaderView, tid: i32) -> Result<&[u8], Error> {
    Ok(header.tid2name(u32::try_from(tid)?))
}

const READ_SRC_TAG: &str = "RG";

#[inline(always)]
pub fn get_read_src(record: &Record) -> &str {
    match record.aux(READ_SRC_TAG.as_bytes()) {
        Ok(Aux::String(str)) => str,
        _ => "0",
        // Ok(_) => "0",
        // _ => panic!("Read source tag {READ_SRC_TAG} overwritten!"),
    }
}

/// Given read pairs from a different alignment files (one file = one bucket), tie-break alignments
/// between mates in each bucket, then tie-break between buckets. Return the winning bucket and its hit
/// taxid if one was conclusively determined.
pub fn cluster_read_buckets(
    read_buckets: Vec<ReadBucket>,
    k2tax: Option<&K2Taxonomy>,
) -> Result<Option<(ReadBucket, u32)>, Error> {
    if read_buckets.is_empty() {
        return Ok(None);
    }

    let mut max_score = 0.0_f32;
    let mut winning_hit: u32 = 0;
    let mut winning_bucket: Option<ReadBucket> = None;

    #[allow(unused_variables)]
    for bucket in read_buckets {
        let checks = [Some((&bucket.rec, false)), bucket.mate.as_ref().map(|mate| (mate, true))];

        // tie-break between mates if they are chimeric
        let (pair_hit, pair_score) = checks
            .iter()
            .flatten()
            .map(|(read, r2)| {
                // scale the mapper's call the number of kmers from that read that match its
                // taxid.
                let kmer_n = if let Some(k2tax) = k2tax {
                    1
                    // k2tax.get_kmer_calls(read.qname()).map(|row| row.get_taxid_kmer_count(map, *r2)).unwrap_or(0)
                } else {
                    1 // if no k2 provided, then we just take read alignments at face value
                };

                (read.taxid, read.cover.frac_matched() * kmer_n as f32)
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
    let html_report = args.html_report.clone();
    let krona_command = args.krona_command.clone().unwrap_or_else(|| "ktImportTaxonomy".to_owned());
    let krona_taxonomy = args.krona_taxonomy.clone();
    if let Some(path) = html_report.as_deref()
        && args.output != "-"
        && canonical_destination(path)? == canonical_destination(&args.output)?
    {
        bail!("Table and HTML reports must use different output paths");
    }
    let n = args.bams.len();
    let main_iter: Vec<(String, Option<(String, String)>)> = match (args.k2_reports, args.k2_classifications) {
        (None, None) => std::iter::zip(args.bams, std::iter::repeat_n(None, n)).collect(),

        (None, Some(_)) | (Some(_), None) => bail!("Must provide both k2 reports and k2 classifications"),

        (Some(reports), Some(classifications)) => {
            if reports.len() != classifications.len() {
                bail!("# of k2 reports is not equal to # of k2 classifications!");
            }
            if reports.len() != n {
                bail!("# of k2 files is not equal to # BAMs!");
            }
            std::iter::zip(args.bams, classifications.into_iter().zip(reports).map(Some)).collect()
        }
    };

    let taxonomy = Taxonomy::from_dir(&args.taxonomy_dir).context("create taxonomy")?;

    let to_stdout = args.output == "-";
    let output: Box<dyn Write> =
        if to_stdout { Box::new(io::stdout()) } else { Box::new(File::create(&args.output).context("create output")?) };
    let mut writer = csv::WriterBuilder::new().delimiter(args.delimiter).from_writer(output);
    let mut html_samples = Vec::new();

    for (bam, k2) in main_iter {
        let mut merge = SamMergeBuffer::default();
        let basename = bam.rsplit_once(".").unwrap_or((&bam, "")).0;
        let mut lt = LocusTracker::new();

        let k2taxon = if let Some((c, r)) = k2 {
            Some(K2Taxonomy::build(BufReader::new(File::open(c)?), BufReader::new(File::open(r)?))?)
        } else {
            None
        };

        let mut reader = BamReader::from_path(&bam).expect("create reader");
        reader.set_threads(4).context("set reader threads")?;
        let header = reader.header();
        let mut seq_len = 0;

        let mut failed_read_writer = std::io::BufWriter::new(
            std::fs::File::create(format!("{basename}_failed_reads.txt")).context("create failed reads file")?,
        );
        // let output_read_path = format!("{basename}_bamtax.bam");
        // let mut read_writer =
        // BamWriter::from_path(&output_read_path, &Header::from_template(header), rust_htslib::bam::Format::Bam)?;

        let mut rec = Record::new();
        let mut i = 0;
        let mut n_passed = 0;

        // copy header for now.
        let header = header.clone();

        /////////////// BEGIN CLOSURE
        let mut handle_buckets = |buckets: Vec<ReadBucket>| -> Result<(), Error> {
            if let Some((bucket, call_taxid)) = cluster_read_buckets(buckets, k2taxon.as_ref())? {
                let reads = [Some(bucket.rec), bucket.mate];

                if let Some(species) = taxonomy.species(call_taxid) {
                    for rec in reads.iter().flatten() {
                        lt.add(species.tax_id, rec.tid, rec.map_interval);

                        seq_len += rec.seq_len;
                        n_passed += 1;
                        // read_writer.write(rec)?;
                    }
                } else {
                    eprintln!("Warning: failed to find species for taxon id {}. Skipping...", call_taxid);

                    // for rec in reads.iter().flatten() {
                    //     failed_read_writer.write(rec.qname())?;
                    //     failed_read_writer.write(b"\t")?;
                    //     failed_read_writer.write(tid2name(&header, rec.tid())?)?;
                    //     failed_read_writer.write(b"\n")?;
                    // }
                }
            }

            Ok(())
        };
        /////////////// END CLOSURE

        while let Some(result) = reader.read(&mut rec) {
            result.context("read record")?;
            i += 1;

            if rec.is_unmapped() || rec.is_secondary() || rec.is_supplementary() {
                continue;
            }

            let alignment =
                ClusterAlignment::from_read_and_source(&rec, record_get_taxid(&header, &rec)?, get_read_src(&rec))?;

            if alignment.cover.frac_aligned() < args.min_frac_read_aligned
                || alignment.cover.frac_matched() < args.min_frac_read_matched
            {
                continue;
            }

            if let Some(buckets) = merge.intake(alignment) {
                handle_buckets(buckets)?;
            }

            if i % 1000 == 0 {
                eprintln!("processed {i} records from {}", &bam);
            }
        }

        let remainder = merge.flush();
        if !remainder.is_empty() {
            i += remainder.len();
            handle_buckets(remainder)?;
        }

        failed_read_writer.flush().expect("writer flush");

        eprintln!("processed {i} records from {bam}");
        let krona_path = basename.to_string() + ".krona";
        let mut krona_writer = BufWriter::new(File::create(&krona_path)?);

        if n_passed > 0 {
            let report = lt.resolve(args.min_loci_per_call, seq_len / n_passed);

            if let Err(error) = report.serialize(&mut writer, &header, basename, Some(&taxonomy))
                && to_stdout
                && !is_broken_pipe(&error)
            {
                return Err(error);
            }

            report.serialize_krona(&mut krona_writer)?;

            if html_report.is_some() {
                html_samples.push(HtmlSample::from_report(basename, &report, &header));
            }
        } else if html_report.is_some() {
            html_samples.push(HtmlSample::empty(basename));
        }
    }

    if let Some(path) = html_report {
        write_html_report(path, &html_samples, &taxonomy, &krona_command, krona_taxonomy.as_deref().map(Path::new))?;
    }

    Ok(())
}
