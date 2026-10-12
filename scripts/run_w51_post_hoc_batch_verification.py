#!/usr/bin/env python3
"""W51 post hoc check: does the Pineda sample-ID date prefix act as a technical batch?

Pineda S01 (donor_information, sample_mapping) and the Cell 2024 STAR Methods
give no batch field, so the YYMMDD prefix of Sample_ID is tested empirically.
Within a fixed diagnosis group (and within fixed ALS genotype), the variance of
the top five astrocyte pseudobulk PCs explained by prefix is compared with a
label-permutation null. Prefix structure inside one diagnosis group cannot be
disease biology, so a strong effect indicates technical batch.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from p01_paths import project_root

OUT = project_root() / "results/W51_Pineda_BA4_discovery_matched_limma"
N_PERM = 5000
N_PC = 5
N_TOP_GENES = 2000
SEED = 20261010


def prefix_r2(counts: pd.DataFrame, sub: pd.DataFrame, rng: np.random.Generator) -> tuple[float, float]:
    x = counts[sub["sample_folder"]]
    x = x[(x > 0).mean(axis=1) >= 0.5]
    logcpm = np.log2(x / x.sum(axis=0) * 1e6 + 1)
    top = logcpm.loc[logcpm.var(axis=1).sort_values().index[-N_TOP_GENES:]]
    z = (top.T - top.T.mean()).to_numpy()
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    pcs, weights = u[:, :N_PC] * s[:N_PC], s[:N_PC] ** 2

    def r2(labels: np.ndarray) -> float:
        total = 0.0
        for k in range(N_PC):
            y = pcs[:, k]
            within = sum(((y[labels == g] - y[labels == g].mean()) ** 2).sum() for g in np.unique(labels))
            total += weights[k] * (1 - within / ((y - y.mean()) ** 2).sum())
        return total / weights.sum()

    labels = sub["prefix"].to_numpy()
    observed = r2(labels)
    null = np.array([r2(rng.permutation(labels)) for _ in range(N_PERM)])
    return observed, float((np.sum(null >= observed) + 1) / (N_PERM + 1))


def main() -> None:
    rng = np.random.default_rng(SEED)
    counts = pd.read_csv(OUT / "w51_pineda_astrocyte_pseudobulk_counts.tsv.gz", sep="\t", index_col=0)
    samples = pd.read_csv(OUT / "w51_pineda_astrocyte_pseudobulk_samples.tsv", sep="\t")
    samples["prefix"] = samples["sample_folder"].str.split("_").str[0]
    samples = samples[samples["lib_size"] > 0]
    region, cond, group = samples["Region"], samples["Condition"], samples["Group"]
    tests = [
        ("BA4 ALS, all prefixes", samples[(region == "BA4") & (cond == "ALS")]),
        ("BA4 C9ALS, 191112 vs 210910", samples[(region == "BA4") & (group == "C9ALS") & samples["prefix"].isin(["191112", "210910"])]),
        ("BA4 SALS, 191112 vs 201019", samples[(region == "BA4") & (group == "SALS") & samples["prefix"].isin(["191112", "201019"])]),
        ("BA4 FTLD, all prefixes", samples[(region == "BA4") & (cond == "FTLD")]),
        ("BA9 ALS, 210526 vs 210910", samples[(region == "BA9") & (cond == "ALS")]),
    ]
    rows = []
    for name, sub in tests:
        observed, p_value = prefix_r2(counts, sub, rng)
        counts_by_prefix = ";".join(f"{k}:{v}" for k, v in sub["prefix"].value_counts().sort_index().items())
        rows.append({"test": name, "n_samples": len(sub), "prefix_counts": counts_by_prefix,
                     "prefix_R2_top5PC": round(observed, 3), "perm_p": round(p_value, 4)})
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "w51_post_hoc_batch_verification.tsv", sep="\t", index=False)

    layout = samples[samples["Region"].isin(["BA4", "BA9"])]
    pd.crosstab([layout["Region"], layout["prefix"]], layout["Condition"]).to_csv(
        OUT / "w51_pineda_prefix_by_condition.tsv", sep="\t"
    )
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
