use crate::{gmm, output, pca, reference, sample_meta, vcf};
use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use std::path::PathBuf;

pub fn run(
    vcf_path: PathBuf,
    weights_path: PathBuf,
    meta_path: PathBuf,
    eigenval_path: PathBuf,
    sample_meta_path: Option<PathBuf>,
    output_path: PathBuf,
) -> Result<()> {
    let pb = ProgressBar::new_spinner();
    pb.set_style(ProgressStyle::default_spinner().template("{spinner} {msg}").unwrap());

    pb.set_message("Loading reference files...");
    let snps = reference::load_weights(&weights_path)?;
    let n_pcs = snps.first().map(|s| s.weights.len()).unwrap_or(0);
    let ref_meta = reference::load_meta(&meta_path)?;
    anyhow::ensure!(
        ref_meta.n_pcs == n_pcs,
        "manifest declares n_pcs={} but weights file has {} PC columns",
        ref_meta.n_pcs,
        n_pcs
    );
    let eigenvalues = reference::load_eigenvalues(&eigenval_path, n_pcs)?;
    pb.println(format!(
        "Reference: {} SNPs, {} PCs, build {}",
        snps.len(),
        n_pcs,
        ref_meta.reference_build
    ));

    let meta = sample_meta_path
        .as_deref()
        .map(sample_meta::load)
        .transpose()
        .context("loading sample metadata")?;

    pb.set_message("Extracting genotypes from VCF...");
    let vcf_data = vcf::extract(&vcf_path, &snps)?;

    let n_samples = vcf_data.samples.len();
    let n_found = vcf_data.snps.iter().filter(|s| s.is_some()).count();
    pb.println(format!(
        "Samples: {}  |  SNPs found: {}/{} ({:.1}%)",
        n_samples,
        n_found,
        snps.len(),
        n_found as f64 / snps.len() as f64 * 100.0
    ));

    anyhow::ensure!(n_samples >= 100, "fewer than 100 samples ({}); refusing to run", n_samples);

    // Validate sample metadata covers all VCF samples
    if let Some(ref m) = meta {
        let missing: Vec<&str> = vcf_data
            .samples
            .iter()
            .filter(|id| !m.map.contains_key(id.as_str()))
            .map(|s| s.as_str())
            .collect();
        anyhow::ensure!(
            missing.is_empty(),
            "{} VCF samples not found in metadata: {:?}...",
            missing.len(),
            &missing[..missing.len().min(5)]
        );
        pb.println(format!("Metadata mode: {:?}", m.mode));
    }

    pb.set_message("Projecting to PC space...");
    let pc_coords = pca::project(&vcf_data, &snps, &eigenvalues);

    pb.set_message("Fitting GMM...");
    let distributions = gmm::fit(&pc_coords, &vcf_data.samples, meta.as_ref(), &ref_meta)
        .context("fitting GMM")?;

    let per_sex_counts = meta.as_ref().and_then(|m| {
        matches!(
            m.mode,
            sample_meta::MetaMode::SexAndAge | sample_meta::MetaMode::SexOnly
        )
        .then(|| {
            let (female, male) = m.count_by_sex(&vcf_data.samples);
            output::PerSexCounts { female, male }
        })
    });

    pb.set_message("Writing output...");
    output::write(
        &output_path,
        &snps,
        &vcf_data,
        distributions,
        &ref_meta,
        per_sex_counts,
    )?;

    pb.finish_and_clear();
    println!("Output written to {}", output_path.display());

    Ok(())
}
