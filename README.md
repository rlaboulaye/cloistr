# CLOISTR
**C**ontrols via **L**atent-space **O**ptimization using **I**nferred **S**tatistical **TR**ansport

Select ancestry-, age-, and sex-matched genomic controls from a reference database for a query cohort, returning only aggregate allele and genotype counts — no individual-level data crosses either direction.

## Overview

CLOISTR has two binaries with distinct roles:

| Binary | Role | Who runs it |
|---|---|---|
| `cloistr` | Server-side matching — takes an encoded query, selects controls from a `db_pack`, returns aggregate counts | Database operator |
| `cloistr-encode` | Client-side encoding — takes a cohort VCF, projects to PC space, fits a GMM, and packages a query archive | Researcher |

The typical workflow is:

1. The database operator builds a `db_pack` from their VCF and runs `cloistr`.
2. The operator distributes a `ref_pack` + the `cloistr-encode` binary to researchers.
3. A researcher runs `cloistr-encode` on their cohort VCF to produce a `query.enc.gz` archive.
4. The researcher sends `query.enc.gz` to the operator (or a web service wrapping `cloistr`).
5. The operator runs `cloistr` and returns aggregate control counts — no individual records leave.

See [ALGORITHM.md](ALGORITHM.md) for a description of the optimal transport candidate selection and χ² greedy refinement.

## What you need

To run the database side you need:

- **`raw/db.vcf.gz`** — bgzipped, tabix-indexed VCF of all database samples (`.tbi` alongside).
- **`raw/db_meta.parquet`** — per-sample metadata with columns `sample_id`, `sex` (0/1), `age`, and optionally `population`.

If you do not have real data yet, `sim/` provides a simulation pipeline that writes these files automatically. See [sim/README.md](sim/README.md).

## Building

```bash
# Both binaries
cargo build --release
# Binaries at target/release/cloistr and target/release/cloistr-encode

# Or build individually
cargo build --release -p cloistr
cargo build --release -p cloistr-encode
```

Python scripts require a virtual environment with core dependencies:

```bash
uv pip install .
# For simulation only, additionally:
uv pip install ".[sim]"
```

## Database setup pipeline

### 0. Prepare raw inputs

Place your database VCF and metadata in `raw/`:

```
raw/db.vcf.gz        # bgzipped, tabix-indexed
raw/db.vcf.gz.tbi
raw/db_meta.parquet  # sample_id, sex, age, population
```

If using the simulation, `python sim/simulate_db.py` writes these files for you.

### 1. Run PCA

```bash
bash scripts/plink_pca.sh   # requires plink 1.9; edit NTHREADS as needed
```

Outputs: `pca/db.bim`, `pca/db.prune.in`, `pca/db_pca.eigenvec`, `pca/db_pca.eigenvec.var`, `pca/db_pca.eigenval`.

### 2. Build `db_pack`

```bash
python scripts/build_db_pack.py   # reads from raw/ and pca/; writes to db_pack/
```

### 3. Build `ref_pack` (for distribution to researchers)

```bash
python scripts/build_ref_pack.py   # reads from db_pack/ and pca/; writes to ref_pack/
```

Distribute `ref_pack/` and the `cloistr-encode` binary to researchers. See [encode/README.md](encode/README.md) for the researcher-facing workflow.

## Running `cloistr`

```bash
cloistr run \
  --query      queries/query.enc.gz \
  --db-pack    db_pack/ \
  --n-controls 500 \
  --out        controls/controls.tsv.gz \
  --summary    controls/summary.json
```

Run `cloistr run --help` for the full parameter list.

### Outputs

**`controls.tsv.gz`** — gzipped TSV of per-site aggregate control counts:

| Column | Description |
|---|---|
| `chrom`, `pos`, `ref`, `alt` | Site coordinates |
| `AC_ctrl` | Alt allele count across selected controls |
| `AN_ctrl` | Total allele count (excludes missing calls) |
| `n_00_ctrl`, `n_01_ctrl`, `n_11_ctrl`, `n_miss_ctrl` | Genotype counts |

**`summary.json`** — pipeline diagnostics: Sinkhorn convergence, initial/final genomic control λ, per-population counts, and a k-anonymized age histogram.

## Key tuning parameters

| Flag | Default | Notes |
|---|---|---|
| `--n-controls` | — | Number of controls to select (required) |
| `--pool-factor` | 4 | Candidate pool size as a multiple of `n-controls` |
| `--seed` | 42 | RNG seed for reproducibility |
| `--refine-tol` | 0.01 | Stop refinement when `\|log λ\| < tol` |
| `--sinkhorn-eps` | auto | OT regularization ε (0 → `median(cost)/50`) |
| `--sinkhorn-rho` | 0.1 | Marginal-KL penalty; lower = more mass may be discarded |
| `--exclude-population` | — | Comma-separated population labels to exclude |

## Privacy

The researcher submits only aggregates: per-site alt/total allele counts and a GMM fit over PCA (and optionally age) space. No individual genotypes, phenotype labels, or PCA coordinates leave the researcher's environment.

`cloistr` returns only aggregates: per-site genotype counts over the selected control set, and a k-anonymized demographic summary. The optional `--selected-out` flag writes individual-level records for internal operator QC only.
