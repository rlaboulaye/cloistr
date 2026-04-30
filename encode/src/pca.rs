use crate::reference::SnpWeight;
use crate::vcf::VcfData;
use ndarray::Array2;

/// Project samples into PC space.
/// Returns an (n_samples × n_pcs) matrix where n_pcs = eigenvalues.len().
/// SNPs absent from the VCF contribute 0 to all projections.
pub fn project(vcf_data: &VcfData, snps: &[SnpWeight], eigenvalues: &[f64]) -> Array2<f64> {
    let n_samples = vcf_data.samples.len();
    let n_pcs = eigenvalues.len();
    let mut pc_coords = Array2::<f64>::zeros((n_samples, n_pcs));
    let n_snps = snps.len() as f64;

    for (snp_idx, snp_data_opt) in vcf_data.snps.iter().enumerate() {
        let snp_data = match snp_data_opt {
            Some(d) => d,
            None => continue,
        };
        let freq = snps[snp_idx].effect_allele_freq;
        let sqrt_2pq = (2.0 * freq * (1.0 - freq)).sqrt().max(1e-6);
        let weights = &snps[snp_idx].weights;
        for (sample_idx, &dosage) in snp_data.dosages.iter().enumerate() {
            if dosage == 255 {
                continue; // missing → mean-impute: centered = 0, no contribution
            }
            let standardized = (dosage as f64 - 2.0 * freq) / sqrt_2pq;
            for pc in 0..n_pcs {
                pc_coords[(sample_idx, pc)] += standardized * weights[pc];
            }
        }
    }

    // Scale by n_snps and eigenvalue
    for pc in 0..n_pcs {
        let scale = n_snps * eigenvalues[pc].sqrt().max(1e-6);
        for s in 0..n_samples {
            pc_coords[(s, pc)] /= scale;
        }
    }

    pc_coords
}
