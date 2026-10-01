pub mod cmd_build_db;
pub mod cmd_check_taxonomy;
pub mod cmd_cluster;
pub mod cmd_download_accs;
pub mod cmd_download_fastas;
pub mod cmd_emit_read_names;
pub mod cmd_extract;
pub mod cmd_k2_split;
pub mod cmd_report_species;
pub mod cmd_taxosplit;

mod assembly_dir_iterator;
mod bam_utils;
mod filter_read;
mod k2_taxonomy;
mod k2_utils;
mod locus_tracker;
mod sam_merge_buffer;
mod taxonomy;

use clap::{Parser, Subcommand};
use cmd_build_db::{BuildDbArgs, build_db_main};
use cmd_check_taxonomy::{CheckTaxonomyArgs, check_taxonomy_main};
use cmd_cluster::{ClusterArgs, cluster_main};
use cmd_download_accs::{DownloadAccessionsArgs, download_accs_main};
use cmd_download_fastas::{DownloadFastasArgs, download_fastas_main};
use cmd_emit_read_names::{EmitNamesArgs, emit_names_main};
use cmd_extract::{ExtractArgs, extract_main};
use cmd_k2_split::{K2SplitArgs, k2_split_main};
use cmd_report_species::{ReportSpeciesArgs, report_species_main};
use cmd_taxosplit::{TaxoSplitArgs, taxosplit_main};

#[derive(Subcommand)]
enum Commands {
    Cluster(ClusterArgs),
    Names(EmitNamesArgs),
    TaxonCheck(CheckTaxonomyArgs),
    BuildDb(BuildDbArgs),
    ReportSpecies(ReportSpeciesArgs),
    DownloadAccs(DownloadAccessionsArgs),
    DownloadFastas(DownloadFastasArgs),
    TaxoSplit(TaxoSplitArgs),
    Extract(ExtractArgs),
    K2Split(K2SplitArgs),
}

#[derive(Parser)]
#[command(version, arg_required_else_help = true)]
pub struct Args {
    #[command(subcommand)]
    command: Commands,
}

pub fn main() {
    let args = Args::parse();
    let result = match args.command {
        Commands::Cluster(args) => cluster_main(args),
        Commands::Names(args) => emit_names_main(args),
        Commands::TaxonCheck(args) => check_taxonomy_main(args),
        Commands::BuildDb(args) => build_db_main(args),
        Commands::ReportSpecies(args) => report_species_main(args),
        Commands::DownloadAccs(args) => download_accs_main(args),
        Commands::DownloadFastas(args) => download_fastas_main(args),
        Commands::TaxoSplit(args) => taxosplit_main(args),
        Commands::Extract(args) => extract_main(args),
        Commands::K2Split(args) => k2_split_main(args),
    };

    if let Err(e) = result {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
