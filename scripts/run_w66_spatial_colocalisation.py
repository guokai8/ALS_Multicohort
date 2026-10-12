#!/usr/bin/env python3
"""W66 A and D (Amendment 9; exploratory): where the frozen ALS-up genes sit in GSE330224 Visium
tissue, and whether candidate ligand-receptor genes from W65 are detected in the same spots.

A1. Region share: for each module gene, mean CPM per anatomical region group (donor-level
    pseudobulk), divided by the sum across region groups; median over independent donors.
A2. Spot co-localisation: within each section, Spearman correlation across spots (total counts
    >= 500) between the ALS-up module score (and each ALS-up gene, log1p CPM x 1e4) and each
    lineage marker score; summarised over sections (median rho, fraction of sections with rho > 0).
D.  Ligand-receptor co-detection: within each section, observed / expected fraction of spots in
    which both ligand and receptor are detected (expected = product of the marginal detection
    fractions); donor mean of log2((observed + 1e-3) / (expected + 1e-3)); ALS vs control
    Wilcoxon on independent donors (exploratory, nominal).
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_w66_visium_pseudobulk import MARKERS, MIN_COUNTS, REGION_GROUPS  # noqa: E402
from p01_paths import project_root  # noqa: E402

logger = logging.getLogger(__name__)
ROOT = project_root()
W66 = ROOT / "results" / "W66_spatial_context"
INP = W66 / "inputs"
LINEAGES = [*MARKERS, "vascular"]
LR_PAIRS = [("C3", "C3AR1"), ("C3", "CR1"), ("C3", "LRP1"), ("ANGPT1", "TEK"), ("KITLG", "KIT"),
            ("EFNA5", "EPHA4"), ("EGF", "ERBB4"), ("ANXA1", "EGFR")]


def region_share(tissue: str, units: pd.DataFrame, counts: pd.DataFrame, up_ids: list[str]) -> pd.DataFrame:
    u = units[(units["unit_type"] == "region") & ~units["gse330130_donor"] & (units["lib_size"] > 0)]
    cpm = counts[u["unit"]].div(u["lib_size"].to_numpy(), axis=1) * 1e6
    rows = []
    for donor, d in u.groupby("donor"):
        sub = cpm.loc[up_ids, d["unit"]]
        sub.columns = d["stratum"].to_numpy()
        share = sub.div(sub.sum(axis=1).replace(0, np.nan), axis=0)
        for gid, r in share.iterrows():
            for region, v in r.items():
                rows.append({"tissue": tissue, "donor": donor, "gene_id": gid, "region": region, "share": v})
    out = pd.DataFrame(rows)
    return out.groupby(["tissue", "gene_id", "region"], as_index=False)["share"].median()


def colocalisation(tissue: str, spots: pd.DataFrame, up_symbols: list[str]) -> pd.DataFrame:
    s = spots[spots["total_counts"] >= MIN_COUNTS].copy()
    for g in up_symbols:
        s[f"gene_{g}"] = np.log1p(1e4 * s[f"count_{g}"] / s["total_counts"])
    features = ["module_up"] + [f"gene_{g}" for g in up_symbols]
    rows = []
    for section, d in s.groupby("section"):
        for f in features:
            for lin in LINEAGES:
                rho = stats.spearmanr(d[f], d[f"score_{lin}"]).statistic if d[f].nunique() > 1 else np.nan
                rows.append({"tissue": tissue, "section": section, "donor": d["donor"].iloc[0], "feature": f,
                             "lineage": lin, "rho": rho})
    per = pd.DataFrame(rows)
    per.to_csv(W66 / f"w66_colocalisation_by_section_{tissue}.tsv.gz", sep="\t", index=False)
    return per.groupby(["tissue", "feature", "lineage"]).agg(
        n_sections=("rho", "count"), median_rho=("rho", "median"), frac_positive=("rho", lambda r: float((r > 0).mean()))
    ).reset_index()


def codetection(tissue: str, spots: pd.DataFrame, units: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    s = spots[spots["total_counts"] >= MIN_COUNTS]
    rows = []
    for section, d in s.groupby("section"):
        for lig, rec in LR_PAIRS:
            pl, pr = (d[f"count_{lig}"] > 0).mean(), (d[f"count_{rec}"] > 0).mean()
            both = ((d[f"count_{lig}"] > 0) & (d[f"count_{rec}"] > 0)).mean()
            rows.append({"tissue": tissue, "section": section, "donor": d["donor"].iloc[0], "pair": f"{lig}-{rec}",
                         "p_ligand": pl, "p_receptor": pr, "p_both": both,
                         "log2_obs_exp": np.log2((both + 1e-3) / (pl * pr + 1e-3))})
    per = pd.DataFrame(rows)
    donor_meta = units.drop_duplicates("donor")[["donor", "group", "gse330130_donor"]]
    don = per.groupby(["tissue", "donor", "pair"], as_index=False)[["p_both", "log2_obs_exp"]].mean().merge(donor_meta, on="donor")
    summ = []
    for pair, d in don[~don["gse330130_donor"]].groupby("pair"):
        a, c = d.loc[d["group"] == "ALS", "log2_obs_exp"], d.loc[d["group"] == "Control", "log2_obs_exp"]
        summ.append({"tissue": tissue, "pair": pair, "n_als": len(a), "n_control": len(c),
                     "median_log2_obs_exp_all": d["log2_obs_exp"].median(), "median_ALS": a.median(), "median_control": c.median(),
                     "wilcoxon_p": stats.mannwhitneyu(a, c).pvalue if len(a) and len(c) else np.nan})
    return don, pd.DataFrame(summ)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    seed = ROOT / "results" / "W34_synapse_resolved_regulatory_reanalysis" / "inputs" / "rna_seed_lists"
    up = pd.read_csv(seed / "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv", sep="\t")
    shares, colocs, dons, summs = [], [], [], []
    for tissue in sys.argv[1:] or ["LSC", "MTC"]:
        units = pd.read_csv(INP / f"{tissue}_units.tsv", sep="\t")
        units = units[units["group"].isin(["ALS", "Control"])]
        spots = pd.read_csv(INP / f"{tissue}_spot_table.tsv.gz", sep="\t")
        spots = spots[spots["donor"].isin(units["donor"])]
        counts = pd.read_csv(INP / f"{tissue}_pseudobulk_counts.tsv.gz", sep="\t", index_col=0)
        up_ids = [g for g in up["gene_id"] if g in counts.index]
        up_symbols = [s for s in up["common_name"] if f"count_{s}" in spots.columns]
        if REGION_GROUPS[tissue]:
            shares.append(region_share(tissue, units, counts, up_ids))
        colocs.append(colocalisation(tissue, spots[~spots["donor"].isin(units.loc[units["gse330130_donor"], "donor"])], up_symbols))
        don, summ = codetection(tissue, spots, units)
        dons.append(don); summs.append(summ)
        logger.info("%s done", tissue)
    sym = dict(zip(up["gene_id"], up["common_name"]))
    if shares:
        sh = pd.concat(shares); sh.insert(2, "symbol", sh["gene_id"].map(sym))
        sh.to_csv(W66 / "w66_region_share.tsv", sep="\t", index=False)
    pd.concat(colocs).to_csv(W66 / "w66_colocalisation_summary.tsv", sep="\t", index=False)
    pd.concat(dons).to_csv(W66 / "w66_lr_codetection_by_donor.tsv", sep="\t", index=False)
    pd.concat(summs).to_csv(W66 / "w66_lr_codetection_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
