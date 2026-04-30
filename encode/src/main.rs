mod check;
mod gmm;
mod output;
mod pca;
mod prepare;
mod reference;
mod sample_meta;
mod vcf;

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "cloistr-encode", about = "Encode a cohort VCF into a query archive for control matching")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Extract counts, project to PCA space, fit GMM, and package output
    Prepare {
        /// Input VCF.gz (must have .tbi index)
        #[arg(long)]
        vcf: PathBuf,

        /// Reference PCA weights file (from ref_pack/)
        #[arg(long, default_value = "ref_pack/pca_weights.tsv.gz")]
        weights: PathBuf,

        /// Reference metadata JSON (age scaling, build info; from ref_pack/)
        #[arg(long, default_value = "ref_pack/manifest.json")]
        meta: PathBuf,

        /// Reference PCA eigenvalues file (from ref_pack/)
        #[arg(long, default_value = "ref_pack/db_pca.eigenval")]
        eigenval: PathBuf,

        /// Sample metadata TSV (columns: sample_id, sex, age)
        #[arg(long)]
        sample_meta: Option<PathBuf>,

        /// Output file path
        #[arg(long, default_value = "queries/query.enc.gz")]
        output: PathBuf,
    },

    /// Validate VCF coverage against reference SNPs
    Check {
        /// Input VCF.gz (must have .tbi index)
        #[arg(long)]
        vcf: PathBuf,

        /// Reference PCA weights file (from ref_pack/)
        #[arg(long, default_value = "ref_pack/pca_weights.tsv.gz")]
        weights: PathBuf,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::Prepare { vcf, weights, meta, eigenval, sample_meta, output } => {
            prepare::run(vcf, weights, meta, eigenval, sample_meta, output)
        }
        Command::Check { vcf, weights } => {
            check::run(vcf, weights)
        }
    }
}
