pub mod cmd_cluster;
pub mod cmd_emit_read_names;
mod locus_tracker;
pub mod taxonomy;

use clap::{Parser, Subcommand};
use cmd_cluster::{ClusterArgs, cluster_main};
use cmd_emit_read_names::{EmitNamesArgs, emit_names_main};

#[derive(Subcommand)]
enum Commands {
    Cluster(ClusterArgs),
    Names(EmitNamesArgs),
}

#[derive(Parser)]
#[command(version)]
pub struct Args {
    #[command(subcommand)]
    command: Commands,
}

pub fn main() {
    let args = Args::parse();
    let result = match args.command {
        Commands::Cluster(args) => cluster_main(args),
        Commands::Names(args) => emit_names_main(args),
    };

    result.unwrap();
}
