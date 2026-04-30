# sim

Simulation pipeline for generating synthetic genomic datasets to evaluate the cloistr matching algorithm against admixed Latin American populations.

This is an optional evaluation tool — it is not required to run `cloistr` or `cloistr-encode` with real data.

> **All commands in this document should be run from the repository root**, not from within `sim/`. This keeps default output paths (`raw/`, `pca/`, `db_pack/`, etc.) consistent with the rest of the pipeline.

---

## Purpose

Testing the matching approach requires a dataset where cases and controls have known, structured genetic ancestry, and a phenotype with realistic genetic architecture correlates with that ancestry. Real-world admixed Latin American cohorts are a natural target: they sit in a complex admixture space and are underrepresented in existing genomic resources.

The simulation provides a controlled ground truth to measure how much the algorithm reduces genomic inflation (λ) relative to unmatched control selection, before applying it to real data.

---

## Demographic Models

### Primary: MXB admixture model (m6-All-admixture.yml)

The default demographic model is taken from Medina-Muñoz et al. (2023), who inferred a 6-population admixture model from whole-genome sequences of Mexican and other Latin American populations. The model is stored in `m6-All-admixture.yml` in demes YAML format and is read at runtime via the `demes` library.

Reference repository (see README for full model description):
https://github.com/santiago1234/mxb-genomes

---

## Configuration

Simulation parameters are controlled by `config.yml`. Key sections:

```yaml
demographic_model: m6-All-admixture.yml   # path to demes YAML

populations:                              # per-population sample sizes
  YRI: { simulate: 1300, db_target: 1000 }
  ...

phenotypes:                               # list; one entry = one simulated phenotype
  - name: pheno1
    heritability: 0.25
    prevalence: 0.03
    n_causal_snps: 500
    beta_sex: 0.3
    beta_age: 0.02
```

Multiple phenotype entries simulate different genetic architectures in a single run, allowing comparison of matching difficulty across scenarios.

---

## Running

Install the simulation environment first:

```bash
uv pip install ".[sim]"   # from the repository root
```

Then run the simulation:

```bash
python sim/simulate.py   # writes to raw/: db.vcf.gz, query.vcf.gz, db_meta.parquet, query_meta.parquet, causal_variants.parquet
```

If the query cohort contains multiple populations and you want to restrict to one,
use `subset_query.py` before encoding:

```bash
python sim/subset_query.py --population MXL   # writes raw/mxl_query.vcf.gz + metadata
```

After simulation, continue with the standard database setup pipeline in `scripts/`:

```bash
bash scripts/plink_pca.sh          # reads raw/db.vcf.gz; writes to pca/
python scripts/build_db_pack.py    # reads from raw/ and pca/; writes to db_pack/
python scripts/build_ref_pack.py   # reads from db_pack/ and pca/; writes to ref_pack/
cloistr-encode prepare --vcf raw/query.vcf.gz --sample-meta raw/query_meta.parquet
                                   # writes to queries/query.enc.gz
cloistr run --query queries/query.enc.gz --db-pack db_pack/ --n-controls 500 \
            --out controls/controls.tsv.gz --summary controls/summary.json
```

---

## Output Files

| File | Contents |
|---|---|
| `raw/db.vcf.gz` (+ `.tbi`) | Genotypes for DB individuals |
| `raw/db_meta.parquet` | `sample_id`, `population`, `sex`, `age`, one `case_*` and `liability_*` column per phenotype |
| `raw/query.vcf.gz` (+ `.tbi`) | Genotypes for query cases |
| `raw/query_meta.parquet` | Same columns as `db_meta.parquet` |
| `raw/query_sample_meta.tsv` | `sample_id`, `sex`, `age` — direct input for `cloistr-encode --sample-meta` |
| `raw/causal_variants.parquet` | `phenotype`, `pos`, `ref`, `alt`, `effect` — true causal SNPs for GWAS benchmarking |

Both VCFs contain identical sites, which is required for projecting query samples onto a PCA space built from the DB.

---

## Phenotype Model

Case/control phenotypes are simulated via a liability threshold model:

```
L = G + β_sex · sex + β_age · (age − mean_age) + ε
```

- **G**: polygenic component (sum of dosage × effect over `n_causal_snps` random SNPs, scaled to `Var(G) = h²`)
- **Demographic effects**: sex and age with configurable β coefficients, drawn from the same distributions across all populations so genetic stratification is the sole driver of λ inflation
- **ε**: residual noise (floored to variance 0.05)

Cases are individuals whose liability exceeds the `(1 − prevalence)` quantile of the full population.

---

## References

- Medina-Muñoz et al., 2023. https://github.com/santiago1234/mxb-genomes
