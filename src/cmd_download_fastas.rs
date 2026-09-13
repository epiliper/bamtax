#![allow(clippy::unused_io_amount)]

use crate::cmd_build_db::{DBFilterArgs, DatabaseWriter, process_metadata_fasta};
use anyhow::Error;
use clap::Parser;
use seq_io::fasta::{Reader as FastaReader, Record};
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::process::Command;
use std::sync::Condvar;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{Receiver, Sender, channel};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;
use std::time::Duration;

const BATCH_SIZE: usize = 200;

pub struct ThreadOutput {
    fastabytes: Vec<u8>,
    taxid: u32,
}

pub struct ThreadSignal {
    n: Mutex<usize>,
    c: Condvar,
    jobs_done: AtomicUsize,
}

struct OutputThread {
    handle: JoinHandle<()>,
}

impl OutputThread {
    pub fn new(
        mut output_writer: DatabaseWriter,
        filtargs: DBFilterArgs,
    ) -> (Self, Sender<ThreadOutput>) {
        let (s, r): (Sender<ThreadOutput>, Receiver<ThreadOutput>) = channel();

        let j = std::thread::spawn(move || {
            let mut seen_records: HashSet<String> = HashSet::new();

            while let Ok(data) = r.recv() {
                process_metadata_fasta(
                    data.taxid,
                    data.fastabytes.as_slice(),
                    &mut output_writer,
                    &mut seen_records,
                    &filtargs,
                    false,
                    false,
                )
                .unwrap();
            }

            output_writer.flush().unwrap();
        });

        (Self { handle: j }, s)
    }
}

impl ThreadSignal {
    pub fn wait_while(&self) {
        let _l = self
            .c
            .wait_while(self.n.lock().unwrap(), |free| *free == 0)
            .unwrap();
    }

    pub fn mark_running(&self) {
        *self.n.lock().unwrap() -= 1;
        self.c.notify_one();
    }

    pub fn mark_done(&self) {
        *self.n.lock().unwrap() += 1;
        self.c.notify_one();
    }
}

struct Worker {
    handle: Option<JoinHandle<()>>,
    notify: Arc<ThreadSignal>,
}

impl Worker {
    fn is_finished(&mut self) -> bool {
        if let Some(ref handle) = self.handle {
            if handle.is_finished() {
                self.handle.take().unwrap().join().unwrap();
                return true;
            } else {
                return false;
            }
        }
        true
    }

    fn run(&mut self, batch: Vec<(u32, String)>, api_key: String, sender: Sender<ThreadOutput>) {
        assert!(self.handle.is_none());
        let notify = Arc::clone(&self.notify);
        self.handle = Some(std::thread::spawn(move || {
            notify.mark_running();
            download_fasta_thread(batch, api_key, sender);
            notify.mark_done();
            notify.jobs_done.fetch_add(1, Ordering::Relaxed);
        }));
    }
}

struct ThreadPool {
    workers: Vec<Worker>,
    notify: Arc<ThreadSignal>,
}

impl ThreadPool {
    pub fn new(n_threads: usize) -> Self {
        let notify = Arc::new(ThreadSignal {
            n: Mutex::new(n_threads),
            c: Condvar::new(),
            jobs_done: AtomicUsize::new(0),
        });

        let mut s = Self {
            notify,
            workers: Vec::with_capacity(n_threads),
        };

        (0..n_threads).for_each(|_| {
            s.workers.push(Worker {
                handle: None,
                notify: Arc::clone(&s.notify),
            })
        });

        s
    }

    pub fn get_available(&mut self) -> Option<&mut Worker> {
        self.notify.wait_while();
        eprint!(
            "Finished {} jobs so far...\r",
            self.notify.jobs_done.load(Ordering::Relaxed)
        );
        self.workers
            .iter_mut()
            .find_map(|w| w.is_finished().then_some(w))
    }

    // join all threads
    pub fn conclude(&mut self) {
        loop {
            if self.workers.iter_mut().all(|f| f.is_finished()) {
                break;
            }
        }
    }
}

#[derive(Parser)]
pub struct DownloadFastasArgs {
    // download input, should be TSV of TAXID, NUCLEOTIDE_ACCESSION
    #[arg(short = 'i', long)]
    input: String,

    #[arg(short = 'o', long, default_value_t = "-".to_string())]
    output_prefix: String,

    #[arg(short, long = "gzip")]
    gzip_output: bool,

    #[arg(short = 'n', long, default_value_t = 0.30)]
    pub max_frac_ambig: f32,

    #[arg(short, long = "max_seq_len", default_value_t = 500_000_000)]
    pub max_len: u32,

    #[arg(short = 'l', long, default_value_t = 100)]
    pub min_len: u32,

    #[arg(short, long)]
    pub api_key: String,
}

fn download_fasta_thread(batch: Vec<(u32, String)>, api_key: String, sender: Sender<ThreadOutput>) {
    let accessions = batch
        .iter()
        .map(|(_, accession)| accession.as_str())
        .collect::<Vec<_>>()
        .join(",");
    let fetch = Command::new("efetch")
        .args(["-db", "nucleotide", "-id", &accessions, "-format", "fasta"])
        .env("NCBI_API_KEY", api_key)
        .output()
        .unwrap();

    if !fetch.status.success() | fetch.stdout.is_empty() {
        eprintln!(
            "Downloading batch failed: {}",
            std::str::from_utf8(&fetch.stderr).unwrap()
        );
    }

    let mut reader = FastaReader::new(fetch.stdout.as_slice());
    while let Some(record) = reader.next() {
        let record = record.unwrap();
        let accession = record.id().unwrap();
        if let Some(taxid) = batch
            .iter()
            .find_map(|(taxid, requested)| (requested == accession).then_some(*taxid))
        {
            let mut fastabytes = Vec::new();
            record.write(&mut fastabytes).unwrap();
            sender.send(ThreadOutput { fastabytes, taxid }).unwrap();
        }
    }
}

pub fn download_fastas_main(args: DownloadFastasArgs) -> Result<(), Error> {
    let reader = BufReader::new(std::fs::File::open(&args.input)?);

    let output = DatabaseWriter::new(args.output_prefix.clone(), args.gzip_output, None)?;

    let filtargs = DBFilterArgs {
        min_len: args.min_len,
        max_len: args.max_len,
        max_frac_ambig: args.max_frac_ambig,
    };

    // spawn a thread to handle writes/validation
    let (output, sender) = OutputThread::new(output, filtargs);

    let mut input = reader.lines().peekable();
    let mut threadpool = ThreadPool::new(10);

    loop {
        if input.peek().is_none() {
            break;
        }

        if let Some(worker) = threadpool.get_available() {
            let mut batch = Vec::with_capacity(BATCH_SIZE);
            while batch.len() < BATCH_SIZE {
                let Some(line) = input.next() else {
                    break;
                };
                let line = line?;
                let (taxid, acc) = line.split_once("\t").unwrap();
                batch.push((taxid.parse::<u32>().unwrap(), acc.to_string()));
            }

            std::thread::sleep(Duration::from_millis(110));
            worker.run(batch, args.api_key.to_string(), sender.clone());
        }
    }

    threadpool.conclude();
    std::mem::drop(sender);
    output.handle.join().unwrap();

    Ok(())
}
