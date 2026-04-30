use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use serde::Deserialize;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// One SNP entry from pca_weights.tsv.gz
#[derive(Debug, Clone)]
pub struct SnpWeight {
    pub chrom: String,
    pub pos: u64,
    pub effect_allele: String,
    pub other_allele: String,
    pub effect_allele_freq: f64,
    pub weights: Vec<f64>,
}

/// Contents of manifest.json (from ref_pack/)
#[derive(Debug, Deserialize)]
pub struct GladMeta {
    pub reference_build: String,
    pub n_pcs: usize,
    pub n_samples: usize,
    pub age_mean: f64,
    pub age_sd: f64,
}

pub fn load_meta(path: &Path) -> Result<GladMeta> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let meta: GladMeta = serde_json::from_reader(BufReader::new(file))
        .with_context(|| format!("parsing {}", path.display()))?;
    Ok(meta)
}

pub fn load_weights(path: &Path) -> Result<Vec<SnpWeight>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let decoder = GzDecoder::new(file);
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(decoder);

    let headers = rdr.headers()?.clone();
    // Columns: chrom, pos, effect_allele, other_allele, effect_allele_freq, pc1..pcN
    let n_pcs = headers.len().saturating_sub(5);
    anyhow::ensure!(n_pcs > 0, "weights file has no PC columns");

    let mut snps = Vec::new();
    for (i, result) in rdr.records().enumerate() {
        let rec = result.with_context(|| format!("reading weights row {}", i + 1))?;
        let chrom = rec[0].to_string();
        let pos: u64 = rec[1]
            .parse()
            .with_context(|| format!("parsing pos at row {}", i + 1))?;
        let effect_allele = rec[2].to_string();
        let other_allele = rec[3].to_string();
        let effect_allele_freq: f64 = rec[4]
            .parse()
            .with_context(|| format!("parsing effect_allele_freq at row {}", i + 1))?;
        let mut weights = vec![0f64; n_pcs];
        for j in 0..n_pcs {
            weights[j] = rec[5 + j]
                .parse()
                .with_context(|| format!("parsing pc{} at row {}", j + 1, i + 1))?;
        }
        snps.push(SnpWeight {
            chrom,
            pos,
            effect_allele,
            other_allele,
            effect_allele_freq,
            weights,
        });
    }
    Ok(snps)
}

/// Load the first `n_pcs` eigenvalues from db_pca.eigenval (one per line, no header).
pub fn load_eigenvalues(path: &Path, n_pcs: usize) -> Result<Vec<f64>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let rdr = BufReader::new(file);
    let mut eigenvalues = Vec::with_capacity(n_pcs);
    for (i, line) in rdr.lines().take(n_pcs).enumerate() {
        let line = line.with_context(|| format!("reading eigenvalue line {}", i + 1))?;
        eigenvalues.push(
            line.trim()
                .parse::<f64>()
                .with_context(|| format!("parsing eigenvalue at line {}", i + 1))?,
        );
    }
    anyhow::ensure!(
        eigenvalues.len() == n_pcs,
        "eigenval file has {} values but weights file has {} PCs",
        eigenvalues.len(),
        n_pcs
    );
    Ok(eigenvalues)
}
