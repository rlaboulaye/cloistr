use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Female = 0,
    Male = 1,
}

#[derive(Debug, Clone)]
pub struct SampleInfo {
    pub sex: Sex,
    pub age: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetaMode {
    SexAndAge,
    SexOnly,
    AgeOnly,
    None,
}

pub struct SampleMeta {
    pub map: HashMap<String, SampleInfo>,
    pub mode: MetaMode,
}

impl SampleMeta {
    /// Count `(female, male)` over the given sample IDs. Samples without sex
    /// info (or not in the metadata map) are silently skipped, so callers
    /// should only invoke this when `mode` includes sex.
    pub fn count_by_sex(&self, sample_ids: &[String]) -> (u32, u32) {
        let mut female = 0u32;
        let mut male = 0u32;
        for id in sample_ids {
            if let Some(info) = self.map.get(id) {
                match info.sex {
                    Sex::Female => female += 1,
                    Sex::Male => male += 1,
                }
            }
        }
        (female, male)
    }
}

pub fn load(path: &Path) -> Result<SampleMeta> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("opening {}", path.display()))?;

    let headers = rdr.headers()?.clone();
    let col = |name: &str| -> Option<usize> { headers.iter().position(|h| h == name) };

    let idx_id = col("sample_id").context("missing column: sample_id")?;
    let idx_sex = col("sex");
    let idx_age = col("age");

    let mut map: HashMap<String, SampleInfo> = HashMap::new();
    let mut any_missing_sex = false;
    let mut any_missing_age = false;

    for (i, result) in rdr.records().enumerate() {
        let rec = result.with_context(|| format!("reading metadata row {}", i + 1))?;
        let sample_id = rec[idx_id].to_string();

        let sex = match idx_sex {
            None => {
                any_missing_sex = true;
                Sex::Female // placeholder, mode will be downgraded
            }
            Some(j) => match rec[j].trim() {
                "0" => Sex::Female,
                "1" => Sex::Male,
                "" => {
                    any_missing_sex = true;
                    Sex::Female
                }
                other => bail!("invalid sex value '{}' for sample '{}' (expected 0 or 1)", other, sample_id),
            },
        };

        let age = match idx_age {
            None => {
                any_missing_age = true;
                0.0 // placeholder
            }
            Some(j) => match rec[j].trim() {
                "" => {
                    any_missing_age = true;
                    0.0
                }
                s => s
                    .parse::<f64>()
                    .with_context(|| format!("invalid age for sample '{}'", sample_id))?,
            },
        };

        map.insert(sample_id, SampleInfo { sex, age });
    }

    let mode = match (any_missing_sex || idx_sex.is_none(), any_missing_age || idx_age.is_none()) {
        (false, false) => MetaMode::SexAndAge,
        (false, true) => MetaMode::SexOnly,
        (true, false) => MetaMode::AgeOnly,
        (true, true) => MetaMode::None,
    };

    if any_missing_sex && idx_sex.is_some() {
        eprintln!("Warning: some samples have missing sex — using no-sex mode (single GMM)");
    }
    if any_missing_age && idx_age.is_some() {
        eprintln!("Warning: some samples have missing age — using no-age mode");
    }

    Ok(SampleMeta { map, mode })
}
