#!/usr/bin/env python3
"""W66 post hoc sensitivity analyses (Amendment 10 A; review-prompted, already observed).

A1. Empirical null for ligand-receptor co-detection: for each W66 D pair, 200 random gene pairs
    matched to the ligand and the receptor on overall detection-fraction decile (genes detected in
    >= 1% of eligible spots); same statistic as W66 D; empirical percentile and one-sided P.
A2. Depth-partial co-localisation: W66 A2 with ranks residualised on the rank of log total counts
    within each section.
Independent donors only (not in GSE330130); spots with total counts >= 500.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_w66_visium_pseudobulk import H5AD, MIN_COUNTS  # noqa: E402
from p01_paths import project_root  # noqa: E402
from run_w66_spatial_colocalisation import LR_PAIRS  # noqa: E402

logger = logging.getLogger(__name__)
ROOT = project_root()
W66 = ROOT / "results" / "W66_spatial_context"
OUT = W66 / "post_hoc"
N_RANDOM = 200
SEED = 20261013
LINEAGES = ["astrocyte", "microglia", "endothelial", "mural_fibroblast", "oligodendrocyte", "neuron", "vascular"]


def independent_donors(tissue: str) -> set[str]:
    u = pd.read_csv(W66 / "inputs" / f"{tissue}_units.tsv", sep="\t")
    u = u[(u["unit_type"] == "whole") & u["group"].isin(["ALS", "Control"])]
    return set(u.loc[~u["gse330130_donor"].astype(bool), "donor"])


def codetection_null(tissue: str, rng: np.random.Generator) -> pd.DataFrame:
    a = ad.read_h5ad(H5AD[tissue])
    obs = a.obs.reset_index(drop=True)
    X = a.X.tocsr()
    total = np.asarray(X.sum(axis=1)).ravel()
    keep = (total >= MIN_COUNTS) & obs["Subject ID"].isin(independent_donors(tissue)).to_numpy()
    Xk = X[keep]
    section = obs.loc[keep, "array"].astype(str).to_numpy()
    donor = obs.loc[keep, "Subject ID"].astype(str).to_numpy()
    det = np.asarray((Xk > 0).mean(axis=0)).ravel()
    sym = pd.Series(np.arange(a.n_vars), index=a.var["common_name"].astype(str).to_numpy())
    sym = sym[~sym.index.duplicated()]
    pool = np.flatnonzero(det >= 0.01)
    edges = np.quantile(det[pool], np.linspace(0, 1, 11))
    decile = np.clip(np.searchsorted(edges, det, side="right") - 1, 0, 9)

    pairs_a, pairs_b, labels = [], [], []
    for lig, rec in LR_PAIRS:
        ia, ib = int(sym[lig]), int(sym[rec])
        ra = rng.choice(pool[decile[pool] == decile[ia]], N_RANDOM)
        rb = rng.choice(pool[decile[pool] == decile[ib]], N_RANDOM)
        pairs_a += [ia, *ra]; pairs_b += [ib, *rb]
        labels += [(f"{lig}-{rec}", "observed")] + [(f"{lig}-{rec}", f"random_{k}") for k in range(N_RANDOM)]
    genes = np.unique(np.r_[pairs_a, pairs_b])
    col = {g: i for i, g in enumerate(genes)}
    B = (Xk[:, genes] > 0).toarray()
    ia, ib = np.array([col[g] for g in pairs_a]), np.array([col[g] for g in pairs_b])

    rows = []
    for sec in np.unique(section):
        m = section == sec
        Bs = B[m]
        p = Bs.mean(axis=0)
        both = (Bs[:, ia] & Bs[:, ib]).mean(axis=0)
        rows.append(pd.DataFrame({"section": sec, "donor": donor[m][0], "k": np.arange(len(ia)),
                                  "value": np.log2((both + 1e-3) / (p[ia] * p[ib] + 1e-3))}))
    per = pd.concat(rows)
    med = per.groupby(["donor", "k"])["value"].mean().groupby("k").median()
    lab = pd.DataFrame(labels, columns=["pair", "kind"])
    lab["median_log2_obs_exp"] = med.to_numpy()
    out = []
    for pair, d in lab.groupby("pair", sort=False):
        obs_v = d.loc[d["kind"] == "observed", "median_log2_obs_exp"].iloc[0]
        rnd = d.loc[d["kind"] != "observed", "median_log2_obs_exp"].to_numpy()
        out.append({"tissue": tissue, "pair": pair, "observed": obs_v, "random_median": np.median(rnd),
                    "random_q95": np.quantile(rnd, 0.95), "percentile": float((rnd < obs_v).mean()),
                    "empirical_p": (1 + int((rnd >= obs_v).sum())) / (N_RANDOM + 1), "n_random": N_RANDOM})
    return pd.DataFrame(out)


def depth_partial(tissue: str) -> pd.DataFrame:
    s = pd.read_csv(W66 / "inputs" / f"{tissue}_spot_table.tsv.gz", sep="\t")
    s = s[(s["total_counts"] >= MIN_COUNTS) & s["donor"].isin(independent_donors(tissue))]
    rows = []
    for sec, d in s.groupby("section"):
        rd = np.log(d["total_counts"]).rank().to_numpy()
        def resid(v: np.ndarray) -> np.ndarray:
            return v - np.polyval(np.polyfit(rd, v, 1), rd)
        y = d["module_up"].rank().to_numpy()
        for lin in LINEAGES:
            x = d[f"score_{lin}"].rank().to_numpy()
            rows.append({"tissue": tissue, "section": sec, "lineage": lin, "rho": np.corrcoef(y, x)[0, 1],
                         "rho_depth_partial": np.corrcoef(resid(y), resid(x))[0, 1]})
    per = pd.DataFrame(rows)
    return per.groupby(["tissue", "lineage"], as_index=False)[["rho", "rho_depth_partial"]].median()


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    nulls, partials = [], []
    for tissue in ("LSC", "MTC"):
        nulls.append(codetection_null(tissue, rng))
        partials.append(depth_partial(tissue))
        logger.info("%s done", tissue)
    pd.concat(nulls).to_csv(OUT / "w66_lr_codetection_empirical_null.tsv", sep="\t", index=False)
    pd.concat(partials).to_csv(OUT / "w66_colocalisation_depth_partial.tsv", sep="\t", index=False)
    print(pd.concat(nulls).round(3).to_string(index=False))
    print(pd.concat(partials).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
