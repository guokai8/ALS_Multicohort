#!/usr/bin/env python3
"""W53 summary-statistic lookup of the final W34 astrocyte modules in Masrori 2025.

Source: Masrori et al. 2025 Nat Neurosci (PMC12586167), Supplementary Tables
24, 25, 27 and 28: astrocyte DE healthy vs C9-ALS or sALS in motor cortex (MC)
and spinal cord (SC); 5 donors per group. Raw data are restricted (EGA
EGAS00001006711). The tables are Seurat-style cell-level results, so p values
treat nuclei as replicates and are not donor-level evidence; only effect
directions are used.

Sign convention: positive avg_log2FC is higher in healthy. Homeostatic genes
(SLC1A2, NTRK2, GJA1) are positive and spinal reactive markers (SERPINA3,
CHI3L1, VIM) are negative, so the ALS effect is -avg_log2FC.
"""

from __future__ import annotations

import pandas as pd
from scipy.stats import binomtest, spearmanr

from p01_paths import project_root, source_data_root

WORKBOOK = source_data_root() / "Other_Datasets/Masrori2025_NatNeuro/files/41593_2025_2075_MOESM2_ESM.xlsx"
SEEDS = project_root() / "results/W34_synapse_resolved_regulatory_reanalysis/inputs/rna_seed_lists"
OUT = project_root() / "results/W53_Masrori2025_astrocyte_DE_lookup"
TABLES = {
    "MC_sALS": "STable25",
    "SC_sALS": "STable28",
    "MC_C9ALS": "STable24",
    "SC_C9ALS": "STable27",
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    modules = pd.concat(
        pd.read_csv(SEEDS / f"GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_{d}.tsv", sep="\t")[["common_name", "logFC"]]
        .rename(columns={"common_name": "gene", "logFC": "discovery_logFC"}).assign(module=d)
        for d in ["up", "down"]
    )
    workbook = pd.ExcelFile(WORKBOOK)
    long_rows, summary = [], []
    for comparison, sheet in TABLES.items():
        de = pd.read_excel(workbook, sheet).rename(columns={"Unnamed: 0": "gene"})
        de["als_logFC"] = -de["avg_log2FC"]
        merged = modules.merge(de[["gene", "als_logFC", "p_val_adj"]], on="gene", how="left").assign(comparison=comparison)
        long_rows.append(merged)
        for module, sub in merged.groupby("module"):
            tested = sub.dropna(subset=["als_logFC"])
            concordant = int((tested["als_logFC"].gt(0) if module == "up" else tested["als_logFC"].lt(0)).sum())
            summary.append({
                "comparison": comparison, "module": module, "n_reported": len(tested),
                "n_concordant": concordant,
                "sign_test_p_one_sided": binomtest(concordant, len(tested), 0.5, alternative="greater").pvalue if len(tested) else None,
                "spearman_vs_discovery": spearmanr(tested["als_logFC"], tested["discovery_logFC"]).statistic if len(tested) > 2 else None,
                "median_als_logFC": tested["als_logFC"].median(),
            })
    pd.concat(long_rows).to_csv(OUT / "w53_module_gene_lookup.tsv", sep="\t", index=False)
    result = pd.DataFrame(summary)
    result.to_csv(OUT / "w53_module_direction_summary.tsv", sep="\t", index=False)
    print(result.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
