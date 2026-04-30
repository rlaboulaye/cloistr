use crate::reference;
use anyhow::Result;

pub fn run(vcf_path: std::path::PathBuf, weights_path: std::path::PathBuf) -> Result<()> {
    let snps = reference::load_weights(&weights_path)?;
    println!("Reference SNPs: {}", snps.len());

    let vcf_data = crate::vcf::extract(&vcf_path, &snps)?;
    let n_found = vcf_data.snps.iter().filter(|s| s.is_some()).count();
    let coverage = n_found as f64 / snps.len() as f64 * 100.0;

    println!("Samples in VCF: {}", vcf_data.samples.len());
    println!("SNPs found: {} / {} ({:.1}%)", n_found, snps.len(), coverage);

    if coverage < 90.0 {
        eprintln!("Warning: coverage below 90% — results may be unreliable");
    }

    Ok(())
}
