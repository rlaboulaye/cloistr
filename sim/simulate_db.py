#!/usr/bin/env python3
"""Simulate a synthetic genomic dataset for testing cloistr.

Reads all parameters from a YAML config. The demographic model is loaded
from the path given in config['demographic_model'] (relative to the config
file) using the demes library.

Usage:
    python sim/simulate_db.py [--config sim/config.yml] [--out-dir raw/]
"""
import argparse
import io
import subprocess
from pathlib import Path

import demes
import msprime
import numpy as np
import polars as pl
import yaml


def load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def simulate(config: dict, config_dir: Path, out_dir: Path, prefix: str = "") -> None:
    seq = config["sequence"]
    seed = seq.get("random_seed", 42)
    rng = np.random.default_rng(seed)

    # ── Demography ────────────────────────────────────────────────────────────
    model_path = config_dir / config["demographic_model"]
    graph = demes.load(model_path)
    demography = msprime.Demography.from_demes(graph)

    pop_cfg = config["populations"]
    samples = {pop: cfg["simulate"] for pop, cfg in pop_cfg.items()}
    query_pops = {pop for pop, cfg in pop_cfg.items() if cfg.get("query_cases", False)}
    db_targets = {pop: cfg["db_target"] for pop, cfg in pop_cfg.items()}

    # ── Ancestry + mutations ──────────────────────────────────────────────────
    print("Running ancestry simulation...")
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=seq["length_bp"],
        recombination_rate=seq["recombination_rate"],
        random_seed=seed,
    )
    print(f"  {ts.num_individuals} individuals, {ts.num_trees} trees")

    print("Adding mutations...")
    ts = msprime.sim_mutations(ts, rate=seq["mutation_rate"], random_seed=seed)
    print(f"  {ts.num_sites} sites before filtering")

    # ── MAF + biallelic filter ────────────────────────────────────────────────
    maf = seq["maf_threshold"]
    n_haps = ts.num_samples
    remove = []
    for var in ts.variants():
        if len(var.alleles) != 2:
            remove.append(var.site.id)
            continue
        af = var.genotypes.sum() / n_haps
        if min(af, 1 - af) < maf:
            remove.append(var.site.id)
    ts = ts.delete_sites(remove)
    n_sites = ts.num_sites
    print(f"  {n_sites} sites after biallelic + MAF>={maf} filter")

    # ── Demographics ──────────────────────────────────────────────────────────
    pop_id_to_name = {i: pop.name for i, pop in enumerate(demography.populations)}
    n_ind = ts.num_individuals
    pop_labels = np.array([
        pop_id_to_name[ts.node(ind.nodes[0]).population]
        for ind in ts.individuals()
    ])

    age_cfg = config["age"]
    age = np.clip(
        rng.lognormal(np.log(age_cfg["median"]), age_cfg["sigma"], n_ind),
        age_cfg["min"], age_cfg["max"],
    )
    sex = rng.integers(0, 2, n_ind)  # 0=female, 1=male

    # ── Phenotypes ────────────────────────────────────────────────────────────
    pheno_results: list[tuple[str, np.ndarray, np.ndarray]] = []
    for pheno in config["phenotypes"]:
        n_causal = min(pheno["n_causal_snps"], n_sites)
        causal_ids = set(rng.choice(n_sites, size=n_causal, replace=False).tolist())

        G = np.zeros(n_ind)
        for var in ts.variants():
            if var.site.id in causal_ids:
                effect = rng.standard_normal()
                G += effect * var.genotypes.reshape(n_ind, 2).sum(axis=1)

        h2 = pheno["heritability"]
        g_std = G.std()
        if g_std > 0:
            G = G / g_std * np.sqrt(h2)

        D = pheno["beta_sex"] * sex + pheno["beta_age"] * (age - age.mean())
        var_res = max(1.0 - h2 - D.var(), 0.05)
        liability = G + D + rng.standard_normal(n_ind) * np.sqrt(var_res)

        threshold = np.quantile(liability, 1.0 - pheno["prevalence"])
        case_mask = liability > threshold
        print(f"  {pheno['name']}: {case_mask.sum()} cases ({case_mask.mean():.1%} prevalence)")
        pheno_results.append((pheno["name"], case_mask, liability))

    # ── Query / DB split (primary phenotype drives case routing) ──────────────
    query_fraction = config.get("query_fraction", 0.97)
    primary_case_mask = pheno_results[0][1]

    latam_case_idx = np.where(
        np.isin(pop_labels, list(query_pops)) & primary_case_mask
    )[0]
    n_query = int(len(latam_case_idx) * query_fraction)
    query_idx = rng.choice(latam_case_idx, size=n_query, replace=False)

    db_pool_idx = np.where(~np.isin(np.arange(n_ind), query_idx))[0]
    db_pool_pops = pop_labels[db_pool_idx]

    db_idx_parts = []
    for pop, count in db_targets.items():
        available = db_pool_idx[db_pool_pops == pop]
        n_draw = min(count, len(available))
        if n_draw < count:
            print(f"  Warning: only {n_draw} available for {pop} (requested {count})")
        db_idx_parts.append(rng.choice(available, size=n_draw, replace=False))
    db_idx = np.concatenate(db_idx_parts)

    n_hidden = primary_case_mask[db_idx].sum()
    print(f"  Query: {len(query_idx)} cases")
    print(f"  DB: {len(db_idx)} individuals ({n_hidden} hidden cases, {n_hidden/len(db_idx):.2%})")

    # ── VCF output ────────────────────────────────────────────────────────────
    out_dir.mkdir(parents=True, exist_ok=True)

    def write_vcf_gz(individuals, path):
        sorted_idx = sorted(individuals)
        names = [f"tsk_{i}" for i in sorted_idx]
        with open(path, "wb") as f_out:
            with subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=f_out) as proc:
                wrapper = io.TextIOWrapper(proc.stdin)
                ts.write_vcf(wrapper, individuals=sorted_idx, individual_names=names)
                wrapper.flush()
                proc.stdin.close()

    print("Writing VCFs...")
    write_vcf_gz(query_idx, out_dir / f"{prefix}query.vcf.gz")
    write_vcf_gz(db_idx,    out_dir / f"{prefix}db.vcf.gz")

    # ── Metadata ──────────────────────────────────────────────────────────────
    def build_meta(indices) -> pl.DataFrame:
        cols: dict = {
            "sample_id":  [f"tsk_{i}" for i in indices],
            "population": pop_labels[indices].tolist(),
            "sex":        sex[indices].tolist(),
            "age":        age[indices].tolist(),
        }
        for name, case_mask, liability in pheno_results:
            cols[f"case_{name}"] = case_mask[indices].tolist()
            cols[f"liability_{name}"] = liability[indices].tolist()
        return pl.DataFrame(cols)

    query_meta = build_meta(query_idx)
    db_meta    = build_meta(db_idx)

    query_meta.write_parquet(out_dir / f"{prefix}query_meta.parquet")
    db_meta.write_parquet(out_dir / f"{prefix}db_meta.parquet")

    # TSV of sample_id/sex/age for cloistr-encode --sample-meta
    (
        query_meta
        .select(["sample_id", "sex", "age"])
        .write_csv(out_dir / f"{prefix}query_sample_meta.tsv", separator="\t")
    )

    print("\nQuery population counts:")
    print(query_meta.group_by("population").len().sort("population"))
    print("\nDB population counts:")
    print(db_meta.group_by("population").len().sort("population"))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).parent / "config.yml",
        help="Simulation config YAML (default: sim/config.yml)",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("raw"),
        help="Output directory for VCFs and metadata (default: raw/)",
    )
    p.add_argument(
        "--prefix",
        default="",
        help="Filename prefix for all outputs, e.g. 'toy_' (default: none)",
    )
    args = p.parse_args()
    config = load_config(args.config)
    simulate(config, args.config.parent, args.out_dir, args.prefix)


if __name__ == "__main__":
    main()
