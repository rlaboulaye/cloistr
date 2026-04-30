#!/usr/bin/env python3
"""Filter a simulated query cohort to a single population.

Reads {prefix}query_meta.parquet, extracts sample IDs for the chosen
population, then uses bcftools to produce a matching VCF subset. Writes:
  {prefix}{pop}_query.vcf.gz (+ .tbi)
  {prefix}{pop}_query_meta.parquet
  {prefix}{pop}_query_sample_meta.tsv

Usage:
    python sim/subset_query.py --population MXL [--data-dir raw] [--prefix toy_]
"""
import argparse
import subprocess
import tempfile
from pathlib import Path

import polars as pl


def select(population: str, data_dir: Path, prefix: str) -> None:
    pop_lower = population.lower()
    out_prefix = f"{prefix}{pop_lower}_"

    meta = pl.read_parquet(data_dir / f"{prefix}query_meta.parquet")
    subset = meta.filter(pl.col("population").str.to_lowercase() == pop_lower)
    if subset.is_empty():
        available = meta["population"].unique().sort().to_list()
        raise SystemExit(
            f"No samples found for population '{population}'. "
            f"Available: {available}"
        )
    print(f"Selected {len(subset)} samples from {population}")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp:
        tmp.write("\n".join(subset["sample_id"].to_list()))
        samples_file = Path(tmp.name)

    in_vcf = data_dir / f"{prefix}query.vcf.gz"
    out_vcf = data_dir / f"{out_prefix}query.vcf.gz"
    subprocess.run(
        [
            "bcftools", "view",
            "-S", str(samples_file),
            "--write-index=tbi",
            "-Oz", "-o", str(out_vcf),
            str(in_vcf),
        ],
        check=True,
    )
    samples_file.unlink()
    print(f"Wrote {out_vcf} (+ .tbi)")

    subset.write_parquet(data_dir / f"{out_prefix}query_meta.parquet")
    (
        subset
        .select(["sample_id", "sex", "age"])
        .write_csv(data_dir / f"{out_prefix}query_sample_meta.tsv", separator="\t")
    )
    print(f"Wrote {out_prefix}query_meta.parquet and {out_prefix}query_sample_meta.tsv")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--population", required=True, help="Population label to select (e.g. MXL)")
    p.add_argument("--data-dir", type=Path, default=Path("raw"), help="Directory containing query files (default: raw)")
    p.add_argument("--prefix", default="", help="Filename prefix used during simulation (default: \"\")")
    args = p.parse_args()
    select(args.population, args.data_dir, args.prefix)


if __name__ == "__main__":
    main()
