#!/usr/bin/env python3
"""W56 E1 (Amendment 2): cell-type specificity of the frozen module genes.

Per lineage, the mean CPM across sample (Pineda) or donor x chip (Ulm) strata
with lib_size > 0. Astrocyte share = astrocyte mean CPM / sum of lineage mean
CPMs. Up-module genes with Pineda astrocyte share >= 0.5 form the
astrocyte-enriched subset used in E3b.
"""

from __future__ import annotations

import pandas as pd

from p01_paths import project_root

OUT = project_root() / "results/W56_cell_source_check"
SEEDS = project_root() / "results/W34_synapse_resolved_regulatory_reanalysis/inputs/rna_seed_lists"
ASTRO_SHARE_MIN = 0.5


def lineage_mean_cpm(counts: pd.DataFrame, strata: pd.DataFrame, lineage_col: str) -> pd.DataFrame:
    strata = strata[strata["lib_size"] > 0]
    cpm = counts[strata["stratum"]] / strata.set_index("stratum")["lib_size"] * 1e6
    return pd.DataFrame({lin: cpm[g["stratum"]].mean(axis=1) for lin, g in strata.groupby(lineage_col)})


def main() -> None:
    modules = pd.concat(
        pd.read_csv(SEEDS / f"GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_{d}.tsv", sep="\t")[["gene_id", "common_name"]].assign(module=d)
        for d in ["up", "down"]
    )
    p_counts = pd.read_csv(OUT / "w56_pineda_lineage_counts.tsv.gz", sep="\t", index_col=0)
    p_strata = pd.read_csv(OUT / "w56_pineda_lineage_strata.tsv", sep="\t")
    u_counts = pd.read_csv(OUT / "w56_ulm_lineage_counts.tsv.gz", sep="\t", index_col=0)
    u_strata = pd.read_csv(OUT / "w56_ulm_lineage_strata.tsv", sep="\t")

    p_mean = lineage_mean_cpm(p_counts, p_strata, "lineage")
    u_mean = lineage_mean_cpm(u_counts, u_strata, "lineage")
    p_share = p_mean.div(p_mean.sum(axis=1), axis=0)
    u_share = u_mean.div(u_mean.sum(axis=1), axis=0)

    table = modules.copy()
    for lin in p_share.columns:
        table[f"pineda_share_{lin}"] = table["gene_id"].map(p_share[lin])
    for lin in u_share.columns:
        table[f"ulm_share_{lin}"] = table["common_name"].map(u_share[lin])
    table["pineda_top_lineage"] = table["gene_id"].map(p_share.dropna(how="all").idxmax(axis=1))
    table["ulm_top_lineage"] = table["common_name"].map(u_share.dropna(how="all").idxmax(axis=1))
    table["astrocyte_enriched"] = table["pineda_share_Astrocyte"].ge(ASTRO_SHARE_MIN)
    table.to_csv(OUT / "w56_e1_module_gene_specificity.tsv", sep="\t", index=False)

    subset = table[(table["module"] == "up") & table["astrocyte_enriched"]]
    subset[["gene_id", "common_name"]].to_csv(OUT / "w56_e1_up_astrocyte_enriched_subset.tsv", sep="\t", index=False)
    cols = ["common_name", "module", "pineda_share_Astrocyte", "pineda_share_Microglia", "pineda_share_Endothelial",
            "pineda_share_Vascular/Meningeal", "pineda_top_lineage", "ulm_share_Astrocytes", "ulm_share_Microglia",
            "astrocyte_enriched"]
    print(table[cols].round(3).sort_values(["module", "pineda_share_Astrocyte"], ascending=[False, False]).to_string(index=False))
    print(f"\nAstrocyte-enriched up subset ({len(subset)}): {', '.join(subset['common_name'])}")


if __name__ == "__main__":
    main()
