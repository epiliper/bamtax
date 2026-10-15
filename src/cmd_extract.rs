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
    pub row: Vec<String>,

    /// Read report rows from standard input
    #[arg(long, conflicts_with_all = ["report", "row"])]
    pub rows_stdin: bool,

    /// Provide a ref DB if you want to create ref fastas for each extracted BAM
    #[arg(long, long = "db")]
    pub db: Option<String>,

    /// Directory containing BAM files when supplying rows
    #[arg(short, long, conflicts_with = "report")]
    pub data_dir: Option<String>,

    #[arg(short, long)]
    pub out_dir: Option<String>,
}

fn make_bam_name(dir: Option<&Path>, prefix: &str) -> String {
    let name = format!("{prefix}_bamtax.bam");
    if let Some(dir) = dir { dir.join(Path::new(&name)).to_string_lossy().into_owned() } else { name }
}

fn make_output_bam_name(dir: &Path, target: &str, source: &str, references: &str) -> String {
    fn component(value: &str) -> String {
        let mut output = String::new();
        for character in value.chars() {
            let character =
                if character == '/' || character == '\\' || character.is_whitespace() { '_' } else { character };
            if output.len() + character.len_utf8() > 80 {
                break;
            }
            output.push(character);
        }
        output
    }

    let source_name = Path::new(source).file_name().and_then(|name| name.to_str()).unwrap_or(source);
    let mut hash = 0xcbf29ce484222325_u64;
    for byte in source.bytes().chain([0]).chain(target.bytes()).chain([0]).chain(references.bytes()) {
        hash = hash.wrapping_mul(0x100000001b3) ^ u64::from(byte);
    }
    let prefix = format!("{}_{}_{hash:016x}", component(target), component(source_name));
    make_bam_name(Some(dir), &prefix)
}

pub fn extract_main(args: ExtractArgs) -> Result<(), Error> {
    if args.report.is_none() && args.row.is_empty() && !args.rows_stdin {
        bail!("Provide either an entire report file or one or more report rows");
    }

    let data_dir: PathBuf;

    let lines = if let Some(report) = &args.report {
        data_dir = PathBuf::from(report).parent().expect("failed to get parent dir of report file").to_path_buf();

        let input = File::open(report)?;
        BufReader::new(input).lines().skip(1).collect::<Result<Vec<_>, _>>()?
    } else if args.rows_stdin {
        data_dir = PathBuf::from(args.data_dir.unwrap_or(".".to_string()));
        BufReader::new(std::io::stdin()).lines().collect::<Result<Vec<_>, _>>()?
    } else {
        data_dir = PathBuf::from(args.data_dir.unwrap_or(".".to_string()));
        args.row.iter().flat_map(|row| row.lines().map(str::to_owned)).collect()
    };

    let output_dir = data_dir.join(Path::new(&args.out_dir.unwrap_or(".".to_string())));
    if !output_dir.exists() {
        std::fs::create_dir_all(&output_dir)?;
    }

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

        let bam_name = make_output_bam_name(&output_dir, &alignrow.target, &alignrow.source, &alignrow.references);

        let mut faidx: Option<(FaidxReader, BufWriter<File>)> = args.db.as_ref().map(|db| {
            let out = PathBuf::from(&bam_name).with_extension("fasta");
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_repeated_report_rows() {
        let args = ExtractArgs::try_parse_from(["extract", "--row", "first", "--row", "second"]).unwrap();

        assert_eq!(args.row, ["first", "second"]);
        assert!(!args.rows_stdin);
    }

    #[test]
    fn accepts_rows_from_standard_input() {
        let args = ExtractArgs::try_parse_from(["extract", "--rows-stdin"]).unwrap();

        assert!(args.rows_stdin);
        assert!(args.row.is_empty());
    }

    #[test]
    fn preserves_input_paths_and_flattens_output_names() {
        assert_eq!(make_bam_name(Some(Path::new(".")), "reads/my sample"), "./reads/my sample_bamtax.bam");
        let output = make_output_bam_name(Path::new("extracted"), "Example species", "reads/my sample", "ref-a");
        assert!(output.starts_with("extracted/Example_species_my_sample_"));
        assert!(output.ends_with("_bamtax.bam"));
        assert_ne!(output, make_output_bam_name(Path::new("extracted"), "Example species", "reads/my_sample", "ref-a"));
        assert_ne!(output, make_output_bam_name(Path::new("extracted"), "Example species", "reads/my sample", "ref-b"));
    }
}
