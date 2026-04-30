use crate::reference::SnpWeight;
use anyhow::{Context, Result};
use noodles::{bgzf, tabix, vcf};
use noodles::vcf::variant::record::AlternateBases;
use noodles::vcf::variant::record::samples::Series;
use std::collections::HashMap;
use std::io::BufReader;
use std::path::Path;

pub struct SnpData {
    pub alt_count: u32,
    pub n_alleles: u32,
    pub dosages: Vec<u8>,
}

pub struct VcfData {
    pub samples: Vec<String>,
    pub snps: Vec<Option<SnpData>>,
}

pub fn extract(vcf_path: &Path, snps: &[SnpWeight]) -> Result<VcfData> {
    let index_path = vcf_path.with_extension("gz.tbi");
    let index = tabix::fs::read(&index_path)
        .or_else(|_| tabix::fs::read(vcf_path.with_extension("tbi")))
        .context("loading .tbi index — ensure it is present alongside the VCF")?;

    let file = std::fs::File::open(vcf_path)
        .with_context(|| format!("opening {}", vcf_path.display()))?;
    let mut reader = vcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(file)));
    let header = reader.read_header().context("reading VCF header")?;

    let samples: Vec<String> = header.sample_names().iter().map(|s| s.to_string()).collect();
    let n_samples = samples.len();

    // Group SNP indices by chromosome for one region query per chromosome
    let mut by_chrom: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, snp) in snps.iter().enumerate() {
        by_chrom.entry(snp.chrom.as_str()).or_default().push(i);
    }

    let mut results: Vec<Option<SnpData>> = (0..snps.len()).map(|_| None).collect();

    for (chrom, mut indices) in by_chrom {
        indices.sort_by_key(|&i| snps[i].pos);
        let min_pos = snps[indices[0]].pos;
        let max_pos = snps[indices[indices.len() - 1]].pos;

        let region = format!("{}:{}-{}", chrom, min_pos, max_pos)
            .parse()
            .with_context(|| format!("building region for {}", chrom))?;

        let pos_map: HashMap<u64, usize> = indices.iter().map(|&i| (snps[i].pos, i)).collect();

        // reader.query takes the index as a third argument; records() returns the iterator
        let query = reader
            .query(&header, &index, &region)
            .with_context(|| format!("querying {}", chrom))?;

        for result in query.records() {
            let record = result.context("reading VCF record")?;

            let pos = match record.variant_start() {
                Some(Ok(p)) => usize::from(p) as u64,
                _ => continue,
            };

            let snp_idx = match pos_map.get(&pos) {
                Some(&i) => i,
                None => continue,
            };
            let snp = &snps[snp_idx];

            let ref_allele = record.reference_bases().to_string();
            let alt_alleles: Vec<String> = record
                .alternate_bases()
                .iter()
                .map(|a| a.map(|s| s.to_string()).unwrap_or_default())
                .collect();

            if alt_alleles.len() != 1 {
                continue;
            }
            let alt_allele = &alt_alleles[0];

            let effect_is_alt =
                if &snp.effect_allele == alt_allele && &snp.other_allele == &ref_allele {
                    true
                } else if &snp.effect_allele == &ref_allele && &snp.other_allele == alt_allele {
                    false
                } else {
                    continue;
                };

            let mut alt_count = 0u32;
            let mut n_alleles = 0u32;
            let mut dosages = vec![0u8; n_samples];

            let samples_data = record.samples();
            let gt_series = match samples_data.select("GT") {
                Some(s) => s,
                None => continue,
            };

            for (sample_idx, gt_value) in gt_series.iter(&header).enumerate() {
                let gt = match gt_value {
                    Ok(Some(vcf::variant::record::samples::series::Value::Genotype(gt))) => gt,
                    _ => continue,
                };

                let mut sample_alt = 0u8;
                let mut sample_n = 0u8;
                for allele in gt.iter() {
                    let (allele_idx, _phasing) = allele.context("reading allele")?;
                    if let Some(idx) = allele_idx {
                        sample_n += 1;
                        let is_effect = if effect_is_alt { idx == 1 } else { idx == 0 };
                        if is_effect {
                            sample_alt += 1;
                        }
                    }
                }
                if sample_n > 0 {
                    alt_count += sample_alt as u32;
                    n_alleles += sample_n as u32;
                    dosages[sample_idx] = sample_alt;
                } else {
                    dosages[sample_idx] = 255; // missing → mean-impute in projection
                }
            }

            results[snp_idx] = Some(SnpData { alt_count, n_alleles, dosages });
        }
    }

    Ok(VcfData { samples, snps: results })
}
