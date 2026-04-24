#!/usr/bin/env python3
"""
Prepare GLAD reference files for distribution.

Inputs (run from project root):
  db.bim              - PLINK BIM (all db SNPs after MAF filter)
  db_pca.eigenvec.var - PCA projection weights (LD-pruned SNPs, no header)
  db_meta.parquet     - DB sample metadata (for age normalization params)
  query_meta.parquet  - Query sample metadata (to generate test TSV)

Outputs:
  glad_pca_weights.tsv.gz  - SNP weights for CLI distribution
  glad_meta.json           - Reference metadata (age scaling, build, etc.)
  query_meta.tsv           - Test TSV in the format the CLI expects from users
"""
import gzip
import json
from pathlib import Path
import polars as pl

N_PCS = 30
PC_COLS = [f"pc{i+1}" for i in range(N_PCS)]
ROOT = Path(__file__).parent


def load_bim(path: Path) -> pl.DataFrame:
    return pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "snp_id", "cm", "pos", "a1", "a2"],
        schema_overrides={"chrom": pl.String, "snp_id": pl.Int64, "pos": pl.Int64, "a1": pl.String, "a2": pl.String},
    )


def load_eigenvec_var(path: Path) -> pl.DataFrame:
    # Format: CHROM SNP_ID EFFECT_ALLELE OTHER_ALLELE PC1..PC30 (space-separated, no header)
    return pl.read_csv(
        path,
        separator=" ",
        has_header=False,
        new_columns=["chrom", "snp_id", "effect_allele", "other_allele"] + PC_COLS,
        schema_overrides={
            "chrom": pl.String,
            "snp_id": pl.Int64,
            "effect_allele": pl.String,
            "other_allele": pl.String,
            **{c: pl.Float64 for c in PC_COLS},
        },
    )


def add_chr_prefix(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.when(pl.col("chrom").str.starts_with("chr"))
        .then(pl.col("chrom"))
        .otherwise(pl.concat_str(pl.lit("chr"), pl.col("chrom")))
        .alias("chrom")
    )


def main():
    print("Loading BIM...")
    bim = load_bim(ROOT / "db.bim")
    print(f"  {len(bim):,} SNPs")

    print("Loading eigenvec.var...")
    evec = load_eigenvec_var(ROOT / "db_pca.eigenvec.var")
    print(f"  {len(evec):,} SNPs, {len(evec.columns) - 4} PCs")

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
        print(f"  Allele check passed")

    weights = weights.select(["chrom", "pos", "effect_allele", "other_allele"] + PC_COLS)

    out_path = ROOT / "glad_pca_weights.tsv.gz"
    with gzip.open(out_path, "wb") as f:
        weights.write_csv(f, separator="\t")
    print(f"\nWrote {out_path.name}: {len(weights):,} SNPs x {N_PCS} PCs ({out_path.stat().st_size / 1e6:.1f} MB)")

    # Age scaling parameters from the db
    print("\nComputing age scaling from db_meta...")
    db_meta = pl.read_parquet(ROOT / "db_meta.parquet")
    age_mean = float(db_meta["age"].mean())
    age_sd = float(db_meta["age"].std())
    print(f"  n={len(db_meta):,}  age_mean={age_mean:.3f}  age_sd={age_sd:.3f}")

    glad_meta = {
        "reference_build": "GRCh38",
        "n_pcs": N_PCS,
        "n_snps": len(weights),
        "age_mean": age_mean,
        "age_sd": age_sd,
    }
    meta_path = ROOT / "glad_meta.json"
    with open(meta_path, "w") as f:
        json.dump(glad_meta, f, indent=2)
    print(f"Wrote {meta_path.name}")

    # Test TSV from query_meta (columns the CLI expects from users)
    print("\nGenerating query_meta.tsv...")
    query_meta = pl.read_parquet(ROOT / "query_meta.parquet")
    tsv_path = ROOT / "query_meta.tsv"
    query_meta.select(["sample_id", "sex", "age"]).write_csv(tsv_path, separator="\t")
    print(f"Wrote {tsv_path.name}: {len(query_meta):,} samples")
    print(f"  sex distribution: {query_meta['sex'].value_counts().sort('sex').to_dict(as_series=False)}")


if __name__ == "__main__":
    main()
