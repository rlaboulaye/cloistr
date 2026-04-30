#!/usr/bin/env python3
"""
build_ref_pack.py — build a `ref_pack/` for distribution to cloistr-encode users.

Reads all metadata from a completed db_pack/manifest.json. Also requires the
raw PLINK outputs that are not stored in db_pack.

Inputs:
  --db-pack <dir>          - completed db_pack/ directory
  --bim <file>             - PLINK BIM (db.bim)
  --eigenvec-var <file>    - PLINK eigenvec.var projection weights (db_pca.eigenvec.var)

Outputs (in --out-dir, default: ref_pack/):
  pca_weights.tsv.gz  - per-SNP PCA projection weights for cloistr-encode
  manifest.json            - copied from db_pack/ (build, n_pcs, age scaling, n_samples)
"""
import argparse
import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import polars as pl


def load_bim(path: Path) -> pl.DataFrame:
    return pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "snp_id", "cm", "pos", "a1", "a2"],
        schema_overrides={
            "chrom": pl.String, "snp_id": pl.Int64, "pos": pl.Int64,
            "a1": pl.String, "a2": pl.String,
        },
    )


def load_eigenvec_var(path: Path, n_pcs: int) -> pl.DataFrame:
    pc_cols = [f"pc{i+1}" for i in range(n_pcs)]
    # Format: CHROM SNP_ID EFFECT_ALLELE OTHER_ALLELE PC1..PCn (space-separated, no header)
    return pl.read_csv(
        path,
        separator=" ",
        has_header=False,
        new_columns=["chrom", "snp_id", "effect_allele", "other_allele"] + pc_cols,
        schema_overrides={
            "chrom": pl.String,
            "snp_id": pl.Int64,
            "effect_allele": pl.String,
            "other_allele": pl.String,
            **{c: pl.Float64 for c in pc_cols},
        },
    )


def compute_allele_freqs(db_pack_path: Path) -> pl.DataFrame:
    """Compute VCF-ALT allele frequency for each LD-independent site via dosage data."""
    sites = pl.read_parquet(db_pack_path / "sites.parquet").filter(pl.col("ld_indep"))
    print(f"  {len(sites):,} LD-independent sites")

    print("  Loading geno_ld_indep.parquet...")
    geno_np = pl.read_parquet(db_pack_path / "geno_ld_indep.parquet").to_numpy()
    if geno_np.shape[0] != len(sites):
        raise ValueError(
            f"geno_ld_indep rows ({geno_np.shape[0]}) != ld_indep sites ({len(sites)})"
        )
    print(f"  {geno_np.shape[0]:,} sites × {geno_np.shape[1]:,} samples")

    dosage_f = geno_np.astype(np.float32)
    dosage_f[geno_np == 255] = np.nan
    alt_freq = np.nanmean(dosage_f, axis=1) / 2.0

    return sites.select(["chrom", "pos", "alt"]).with_columns(
        pl.Series("alt_freq", alt_freq, dtype=pl.Float64),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db-pack", type=Path, default=Path("db_pack"),
                   help="Path to completed db_pack/ directory (default: db_pack/)")
    p.add_argument("--bim", type=Path, default=Path("pca/db.bim"),
                   help="PLINK BIM file (default: pca/db.bim)")
    p.add_argument("--eigenvec-var", type=Path, default=Path("pca/db_pca.eigenvec.var"),
                   help="PLINK eigenvec.var projection weights (default: pca/db_pca.eigenvec.var)")
    p.add_argument("--eigenval", type=Path, default=Path("pca/db_pca.eigenval"),
                   help="PLINK eigenval file (default: pca/db_pca.eigenval)")
    p.add_argument("--out-dir", type=Path, default=Path("ref_pack"),
                   help="Output directory (default: ref_pack/)")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading manifest from {args.db_pack / 'manifest.json'} ...")
    manifest = json.loads((args.db_pack / "manifest.json").read_text())
    n_pcs = manifest["n_pcs"]
    print(
        f"  build={manifest['reference_build']}, n_pcs={n_pcs}, "
        f"n_samples={manifest['n_samples']:,}, "
        f"age_mean={manifest['age_mean']:.3f}, age_sd={manifest['age_sd']:.3f}"
    )

    pc_cols = [f"pc{i+1}" for i in range(n_pcs)]

    print(f"\nLoading BIM from {args.bim} ...")
    bim = load_bim(args.bim)
    print(f"  {len(bim):,} SNPs")

    print(f"Loading eigenvec.var from {args.eigenvec_var} ...")
    evec = load_eigenvec_var(args.eigenvec_var, n_pcs)
    print(f"  {len(evec):,} SNPs, {n_pcs} PCs")

    # Join eigenvec.var onto BIM to resolve genomic positions
    weights = evec.join(bim.select(["snp_id", "pos"]), on="snp_id", how="left")
    n_unmatched = weights["pos"].null_count()
    if n_unmatched > 0:
        print(f"  WARNING: {n_unmatched:,} eigenvec.var SNPs not found in BIM — check snp_id format")
        print(f"  eigenvec.var snp_id sample: {evec['snp_id'].head(5).to_list()}")
        print(f"  BIM snp_id sample:          {bim['snp_id'].head(5).to_list()}")
    else:
        print(f"  All {len(weights):,} SNPs matched to BIM positions")

    # Verify allele consistency (effect/other should be {a1, a2} from BIM)
    check = evec.join(bim.select(["snp_id", "a1", "a2"]), on="snp_id", how="left")
    bad = check.filter(
        ~(
            ((pl.col("effect_allele") == pl.col("a1")) & (pl.col("other_allele") == pl.col("a2")))
            | ((pl.col("effect_allele") == pl.col("a2")) & (pl.col("other_allele") == pl.col("a1")))
        )
    )
    if len(bad) > 0:
        print(f"  WARNING: {len(bad):,} SNPs have allele mismatches between eigenvec.var and BIM")
        print(bad.head(5))
    else:
        print("  Allele check passed")

    weights = weights.select(["chrom", "pos", "effect_allele", "other_allele"] + pc_cols)

    # Compute per-site allele frequencies from geno_ld_indep.parquet
    print("\nComputing allele frequencies from db_pack...")
    freqs = compute_allele_freqs(args.db_pack)

    # Normalize chrom for joining (strip "chr" prefix from both sides)
    weights_j = weights.with_columns(pl.col("chrom").str.strip_prefix("chr").alias("_chrom_norm"))
    freqs_j = freqs.with_columns(pl.col("chrom").str.strip_prefix("chr").alias("_chrom_norm"))

    weights = weights_j.join(
        freqs_j.select(["_chrom_norm", "pos", "alt", "alt_freq"]),
        on=["_chrom_norm", "pos"],
        how="left",
    ).drop("_chrom_norm")

    # Flip frequency when the effect allele is REF (not VCF ALT)
    weights = weights.with_columns(
        pl.when(pl.col("effect_allele") == pl.col("alt"))
        .then(pl.col("alt_freq"))
        .otherwise(1.0 - pl.col("alt_freq"))
        .alias("effect_allele_freq")
    ).drop(["alt", "alt_freq"])

    n_missing_freq = weights["effect_allele_freq"].null_count()
    if n_missing_freq > 0:
        print(f"  WARNING: {n_missing_freq:,} SNPs had no frequency match — dropping them")
        weights = weights.filter(pl.col("effect_allele_freq").is_not_null())
    print(f"  Frequencies computed for {len(weights):,} SNPs")

    out_cols = ["chrom", "pos", "effect_allele", "other_allele", "effect_allele_freq"] + pc_cols
    weights_path = args.out_dir / "pca_weights.tsv.gz"
    with gzip.open(weights_path, "wb") as f:
        weights.select(out_cols).write_csv(f, separator="\t")
    print(f"\nWrote {weights_path}: {len(weights):,} SNPs x {n_pcs} PCs ({weights_path.stat().st_size / 1e6:.1f} MB)")

    shutil.copy(args.db_pack / "manifest.json", args.out_dir / "manifest.json")
    print(f"Copied manifest.json → {args.out_dir / 'manifest.json'}")

    shutil.copy(args.eigenval, args.out_dir / "db_pca.eigenval")
    print(f"Copied {args.eigenval} → {args.out_dir / 'db_pca.eigenval'}")


if __name__ == "__main__":
    main()
