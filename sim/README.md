# sim

Simulation pipeline for generating synthetic genomic datasets to evaluate the cloistr matching algorithm against admixed Latin American populations.

This is an optional evaluation tool — it is not required to run `cloistr` or `cloistr-encode` with real data.

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

### Background: Gravel / Browning model

The hand-rolled demographic model previously used in `simulate_db.py` is based on the three-population out-of-Africa model from **Gravel et al. (2011)**, adapted by **Browning et al. (2018)** to include an admixed American population. It uses YRI (African), IBS (European), and MXB (Native American proxy) as ancestral sources for four Latin American populations: PEL, MXL, CLM, and PUR.

This model remains documented in the code for reference but the `m6-All-admixture.yml` model should be preferred for new evaluations.

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
python sim/simulate_db.py   # writes to raw/: db.vcf.gz, query.vcf.gz, db_meta.parquet, query_meta.parquet
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
| `raw/db.vcf.gz` | Genotypes for DB individuals (bgzip compressed, tabix indexed) |
| `raw/db_meta.parquet` | `sample_id`, `population`, `sex`, `age`, one `case_*` column per phenotype, `liability_*` per phenotype |
| `raw/query.vcf.gz` | Genotypes for query cases |
| `raw/query_meta.parquet` | Same columns as `db_meta.parquet` |

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
- Browning et al., 2018. http://dx.doi.org/10.1371/journal.pgen.1007385
- Gravel et al., 2011. https://doi.org/10.1073/pnas.1019276108
