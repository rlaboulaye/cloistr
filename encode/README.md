# cloistr-encode

Encode a cohort VCF into a query archive for control matching with `cloistr`.

## Overview

`cloistr-encode` prepares your cohort's genomic data for submission to a control matching service. Given an imputed VCF, it extracts allele counts at reference SNPs, projects your samples into a shared PC space, and summarizes the resulting ancestry distribution as a Gaussian mixture model. The output is a compact archive (`query.enc.gz`) that you send to the operator.

No individual-level genotype data leaves your machine. The archive contains only summary statistics (allele counts) and a parametric distribution over PC space.

## Prerequisites

- A genome-wide imputed VCF (GRCh38, bgzipped) with an accompanying `.tbi` tabix index.
- Reference files from the operator: a `ref_pack/` directory containing `pca_weights.tsv.gz` and `manifest.json`.
- *(Optional)* A sample metadata TSV with sex and/or age (see [Sample metadata](#sample-metadata)).

## Installation

Pre-built binaries are available from the database operator alongside the `ref_pack`.

To build from source, install [Rust](https://rustup.rs) and run from the repository root:

```bash
cargo build --release -p cloistr-encode
# binary at target/release/cloistr-encode
```

## Quick start

First, run `check` to confirm your VCF has adequate coverage of the reference SNPs:

```bash
cloistr-encode check \
  --vcf     your_cohort.vcf.gz \
  --weights ref_pack/pca_weights.tsv.gz
```

Then run `prepare` to generate the submission archive:

```bash
cloistr-encode prepare \
  --vcf         your_cohort.vcf.gz \
  --weights     ref_pack/pca_weights.tsv.gz \
  --meta        ref_pack/manifest.json \
  --sample-meta metadata.tsv
```

This writes `query.enc.gz` in the current directory. Send it to the operator.

## Command reference

### `check`

Validates that your VCF contains the reference SNPs and reports coverage.

```
cloistr-encode check --vcf <path> [--weights <path>]
```

### `prepare`

Extracts allele counts, projects samples to PC space, fits a GMM, and writes the query archive.

```
cloistr-encode prepare --vcf <path> [--weights <path>] [--meta <path>]
                       [--eigenval <path>] [--sample-meta <path>] [--output <path>]
```

| Flag | Default | Description |
|---|---|---|
| `--vcf` | *(required)* | Imputed VCF.gz with .tbi index |
| `--weights` | `pca_weights.tsv.gz` | Reference PCA weights (from ref_pack/) |
| `--meta` | `manifest.json` | Reference manifest (from ref_pack/) |
| `--eigenval` | `db_pca.eigenval` | Reference PCA eigenvalues (from ref_pack/) |
| `--sample-meta` | *(none)* | Sample metadata TSV (sex, age) |
| `--output` | `query.enc.gz` | Output archive path |

## Sample metadata

The `--sample-meta` file is a tab-separated file with the following columns:

```
sample_id	sex	age
SAMPLE_001	1	54.2
SAMPLE_002	0	61.8
```

- `sex`: `1` = male, `0` = female.
- `age`: continuous value in years.

All columns are optional, but if `sex` or `age` is present it must be present for every sample. Providing both enables the most precise ancestry + age matching. The minimum supported cohort size is 100 samples.

## Output

`query.enc.gz` is a gzip-compressed JSON archive containing:

- **Allele counts** at each reference SNP, used for genomic inflation (λ) calculation during matching.
- **Ancestry distribution** — a Gaussian mixture model in PC space (and optionally age space), used by the optimal transport matching algorithm.
