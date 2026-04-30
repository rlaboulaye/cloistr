use crate::reference::GladMeta;
use crate::sample_meta::{MetaMode, SampleMeta, Sex};
use anyhow::{Context, Result};
use linfa::traits::Fit;
use linfa::DatasetBase;
use linfa_clustering::GaussianMixtureModel;
use ndarray::{Array2, ArrayView2, Axis};
use serde::Serialize;

// Query cohorts are expected to be relatively homogeneous; BIC handles the lower end naturally.
const MAX_COMPONENTS: usize = 5;
const GMM_MAX_ITER: u64 = 200;
const GMM_RUNS: usize = 5;

#[derive(Debug, Serialize)]
pub struct FittedGmm {
    pub n_components: usize,
    pub weights: Vec<f64>,
    pub means: Vec<Vec<f64>>,
    /// Full covariance matrices: n_components × n_dims × n_dims
    pub covariances: Vec<Vec<Vec<f64>>>,
}

#[derive(Debug, Serialize)]
pub struct Distributions {
    pub mode: String,
    pub n_dims: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub all: Option<FittedGmm>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub female: Option<FittedGmm>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub male: Option<FittedGmm>,
}

pub fn fit(
    pc_coords: &Array2<f64>,
    sample_order: &[String],
    meta: Option<&SampleMeta>,
    ref_meta: &GladMeta,
) -> Result<Distributions> {
    let mode = meta.map(|m| m.mode).unwrap_or(MetaMode::None);
    let use_sex = matches!(mode, MetaMode::SexAndAge | MetaMode::SexOnly);
    let use_age = matches!(mode, MetaMode::SexAndAge | MetaMode::AgeOnly);

    let mode_str = match mode {
        MetaMode::SexAndAge => "sex_and_age",
        MetaMode::SexOnly => "sex_only",
        MetaMode::AgeOnly => "age_only",
        MetaMode::None => "none",
    }
    .to_string();

    let features = build_features(pc_coords, sample_order, meta, ref_meta, use_age);
    let n_dims = features.ncols();

    if use_sex {
        let meta = meta.unwrap();
        let mut female_rows = Vec::new();
        let mut male_rows = Vec::new();
        for (i, id) in sample_order.iter().enumerate() {
            if meta.map.get(id).map(|s| s.sex == Sex::Female).unwrap_or(true) {
                female_rows.push(i);
            } else {
                male_rows.push(i);
            }
        }

        let female_features = select_rows(&features, &female_rows);
        let male_features = select_rows(&features, &male_rows);

        Ok(Distributions {
            mode: mode_str,
            n_dims,
            all: None,
            female: Some(fit_gmm(&female_features).context("fitting female GMM")?),
            male: Some(fit_gmm(&male_features).context("fitting male GMM")?),
        })
    } else {
        Ok(Distributions {
            mode: mode_str,
            n_dims,
            all: Some(fit_gmm(&features).context("fitting GMM")?),
            female: None,
            male: None,
        })
    }
}

fn build_features(
    pc_coords: &Array2<f64>,
    sample_order: &[String],
    meta: Option<&SampleMeta>,
    ref_meta: &GladMeta,
    use_age: bool,
) -> Array2<f64> {
    let n = sample_order.len();
    let n_pcs = pc_coords.ncols();
    let n_cols = if use_age { n_pcs + 1 } else { n_pcs };
    let mut features = Array2::<f64>::zeros((n, n_cols));

    for i in 0..n {
        for j in 0..n_pcs {
            features[(i, j)] = pc_coords[(i, j)];
        }
    }

    if use_age {
        let meta = meta.unwrap();
        let pc_sd = 1.0 / (ref_meta.n_samples as f64).sqrt();
        for (i, id) in sample_order.iter().enumerate() {
            if let Some(info) = meta.map.get(id) {
                let z_age = (info.age - ref_meta.age_mean) / ref_meta.age_sd;
                features[(i, n_pcs)] = z_age * pc_sd;
            }
        }
    }

    features
}

fn select_rows(matrix: &Array2<f64>, rows: &[usize]) -> Array2<f64> {
    let n = rows.len();
    let d = matrix.ncols();
    let mut result = Array2::<f64>::zeros((n, d));
    for (new_i, &old_i) in rows.iter().enumerate() {
        result.row_mut(new_i).assign(&matrix.row(old_i));
    }
    result
}

fn cholesky_log_det(mat: ArrayView2<f64>) -> Option<f64> {
    let d = mat.nrows();
    let mut l = Array2::<f64>::zeros((d, d));
    for i in 0..d {
        for j in 0..=i {
            let mut s = mat[(i, j)];
            for p in 0..j {
                s -= l[(i, p)] * l[(j, p)];
            }
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[(i, j)] = s.sqrt();
            } else {
                l[(i, j)] = s / l[(j, j)];
            }
        }
    }
    Some(2.0 * (0..d).map(|i| l[(i, i)].ln()).sum::<f64>())
}

fn gmm_log_likelihood(model: &GaussianMixtureModel<f64>, data: &Array2<f64>) -> Result<f64> {
    let weights = model.weights();
    let means = model.means();
    let precisions = model.precisions();
    let k = weights.len();
    let n = data.nrows();
    let d = data.ncols();
    let log_2pi = (2.0 * std::f64::consts::PI).ln();

    let mut log_norm = Vec::with_capacity(k);
    for c in 0..k {
        let prec = precisions.index_axis(Axis(0), c);
        let log_det_prec = cholesky_log_det(prec)
            .ok_or_else(|| anyhow::anyhow!("precision matrix not positive definite for component {}", c))?;
        log_norm.push(weights[c].ln() - 0.5 * (d as f64) * log_2pi + 0.5 * log_det_prec);
    }

    let mut total_ll = 0.0;
    for i in 0..n {
        let x = data.row(i);
        let mut log_probs = Vec::with_capacity(k);
        for c in 0..k {
            let prec = precisions.index_axis(Axis(0), c);
            let mu = means.row(c);
            let mut mahal = 0.0;
            for j in 0..d {
                let mut v = 0.0;
                for m in 0..d {
                    v += prec[(j, m)] * (x[m] - mu[m]);
                }
                mahal += (x[j] - mu[j]) * v;
            }
            log_probs.push(log_norm[c] - 0.5 * mahal);
        }
        let max_lp = log_probs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        total_ll += max_lp + log_probs.iter().map(|&lp| (lp - max_lp).exp()).sum::<f64>().ln();
    }

    Ok(total_ll)
}

fn gmm_n_params(k: usize, d: usize) -> usize {
    (k - 1) + k * d + k * d * (d + 1) / 2
}

fn bic_score(log_lik: f64, n_params: usize, n: usize) -> f64 {
    -2.0 * log_lik + (n_params as f64) * (n as f64).ln()
}

fn fit_gmm(data: &Array2<f64>) -> Result<FittedGmm> {
    let n = data.nrows();
    let d = data.ncols();
    anyhow::ensure!(n >= 10, "too few samples to fit GMM ({})", n);

    let dataset = DatasetBase::from(data.clone());
    let mut best_bic = f64::INFINITY;
    let mut best_model: Option<GaussianMixtureModel<f64>> = None;
    let mut best_k = 1usize;

    for k in 1..=MAX_COMPONENTS {
        let mut best_ll = f64::NEG_INFINITY;
        let mut best_model_for_k: Option<GaussianMixtureModel<f64>> = None;

        for _ in 0..GMM_RUNS {
            match GaussianMixtureModel::params(k)
                .max_n_iterations(GMM_MAX_ITER)
                .fit(&dataset)
            {
                Ok(model) => {
                    if let Ok(ll) = gmm_log_likelihood(&model, data) {
                        if ll > best_ll {
                            best_ll = ll;
                            best_model_for_k = Some(model);
                        }
                    }
                }
                Err(_) => {}
            }
        }

        if let Some(model) = best_model_for_k {
            let bic = bic_score(best_ll, gmm_n_params(k, d), n);
            if bic < best_bic {
                best_bic = bic;
                best_k = k;
                best_model = Some(model);
            }
        }
    }

    match best_model {
        Some(model) => serialize_gmm(&model, best_k, d),
        None => Err(anyhow::anyhow!("GMM fitting failed for all k in 1..={}", MAX_COMPONENTS)),
    }
}

fn serialize_gmm(
    gmm: &GaussianMixtureModel<f64>,
    k: usize,
    n_dims: usize,
) -> Result<FittedGmm> {
    let weights = gmm.weights().to_vec();
    let means_arr = gmm.means();
    let covs_arr = gmm.covariances();

    let mut means = Vec::with_capacity(k);
    let mut covariances = Vec::with_capacity(k);

    for i in 0..k {
        means.push(means_arr.row(i).to_vec());
        let cov = covs_arr.index_axis(Axis(0), i);
        let cov_rows: Vec<Vec<f64>> = (0..n_dims).map(|j| cov.row(j).to_vec()).collect();
        covariances.push(cov_rows);
    }

    Ok(FittedGmm { n_components: k, weights, means, covariances })
}
