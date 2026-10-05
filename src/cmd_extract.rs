#![allow(clippy::unused_io_amount)]

use crate::bam_utils::create_sub_bam;
use crate::locus_tracker::AlignmentReportRow;
use anyhow::{Error, bail};
use clap::Parser;
use rust_htslib::bam::{IndexedReader, Read as BamRead};
use rust_htslib::faidx::Reader as FaidxReader;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(version, arg_required_else_help = true)]
pub struct ExtractArgs {
    #[arg(long, conflicts_with = "row")]
    pub report: Option<String>,

    #[arg(long, conflicts_with = "report")]
    pub row: Option<String>,

    /// Provide a ref DB if you want to create ref fastas for each extracted BAM
    #[arg(long, long = "db")]
    pub db: Option<String>,

    /// if supplying rows, not report: directory of BAM files + report.
    #[arg(short, long, requires = "row")]
    pub data_dir: Option<String>,

    #[arg(short, long)]
    pub out_dir: Option<String>,
}

fn make_bam_name(dir: Option<&Path>, prefix: &str) -> String {
    let name = format!("{prefix}_bamtax.bam");
    let name = if let Some(dir) = dir { dir.join(Path::new(&name)).to_str().unwrap().to_string() } else { name };

    name.replace(" ", "_")
}

pub fn extract_main(args: ExtractArgs) -> Result<(), Error> {
    if args.report.is_none() && args.row.is_none() {
        bail!("Provide either an entire report file or one or more report rows");
    }

    let data_dir: PathBuf;

    let lines: Result<Vec<String>, std::io::Error> = if let Some(report) = &args.report {
        data_dir = PathBuf::from(report).parent().expect("failed to get parent dir of report file").to_path_buf();

        let input = File::open(report)?;
        BufReader::new(input).lines().skip(1).collect()
    } else {
        data_dir = PathBuf::from(args.data_dir.unwrap_or(".".to_string()));
        BufReader::new(args.row.as_ref().unwrap().as_bytes()).lines().collect()
    };

    let output_dir = data_dir.join(Path::new(&args.out_dir.unwrap_or(".".to_string())));
    if !output_dir.exists() {
        std::fs::create_dir_all(&output_dir)?;
    }

    let lines = lines?;
    let mut prev_reader: Option<(IndexedReader, String)> = None;

    for r in lines {
        // header
        if r.starts_with("source") {
            continue;
        };

        let mut rdr = csv::ReaderBuilder::new().has_headers(false).delimiter(b'\t').from_reader(r.trim().as_bytes());

        let alignrow: AlignmentReportRow = rdr.deserialize().next().unwrap()?;

        let reader = if let Some((bamreader, name)) = &mut prev_reader
            && *name == make_bam_name(Some(&data_dir), &alignrow.source)
        {
            bamreader
        } else {
            let bam_name = make_bam_name(Some(&data_dir), &alignrow.source);
            prev_reader = Some((IndexedReader::from_path(&bam_name)?, bam_name));
            &mut prev_reader.as_mut().unwrap().0
        };

        let references = alignrow.references.split(";");

        let mut tids = vec![];
        let mut names = vec![];

        let bam_name = make_bam_name(Some(&output_dir), &format!("{}_{}", &alignrow.target, &alignrow.source));

        let mut faidx: Option<(FaidxReader, BufWriter<File>)> = args.db.as_ref().map(|db| {
            let out = bam_name.replace(".bam", ".fasta");
            let reader = std::io::BufWriter::new(File::create(out).expect("create ref fasta"));
            let faidx = FaidxReader::from_path(db).expect("opening db file");
            (faidx, reader)
        });

        for r in references {
            let tid = reader.header().tid(r.as_bytes()).unwrap();
            let name = std::str::from_utf8(reader.header().tid2name(tid)).unwrap().to_string();
            tids.push(tid);

            if let Some((faidx_reader, writer)) = faidx.as_mut() {
                let len = faidx_reader.fetch_seq_len(&name);
                let seq = faidx_reader.fetch_seq(&name, 0, len as usize)?;
                writeln!(writer, ">{}", &name)?;
                writer.write(&seq)?;
                writer.write(b"\n")?;
            }
            names.push(name);
        }

        create_sub_bam(reader, &names, &tids, &bam_name)?;
    }

    Ok(())
}
