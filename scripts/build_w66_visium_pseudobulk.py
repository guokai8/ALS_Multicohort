#!/usr/bin/env python3
"""Build W66 inputs from GSE330224 Visium (Amendment 9): spot table and donor-level pseudobulk.

For each tissue (LSC lumbar spinal cord, MTC motor cortex) this writes:
- a spot table: section, donor, anatomical annotation, total counts, lineage marker scores,
  frozen-module scores and the counts of the ligand-receptor genes used in W66 D;
- pseudobulk raw counts (genes x units) for three unit types: donor x whole section,
  donor x anatomical region group, and donor x lineage niche (defined per section from
  marker scores, Amendment 9);
- unit metadata (donor, diagnosis, sex, age, PMI, number of spots, discovery overlap).

No disease comparison is made here.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from p01_paths import project_root, source_data_root  # noqa: E402

logger = logging.getLogger(__name__)

ROOT = project_root()
OUT = ROOT / "results" / "W66_spatial_context" / "inputs"
H5AD = {t: source_data_root() / "GSE330224" / f"GSE330224_{t}_VISIUM_COUNTS.h5ad" for t in ("LSC", "MTC")}

# Lineage marker sets (W55 composition markers, W56 E3a activated-myeloid and vascular markers,
# plus fibroblast/mural markers). SPP1 is left out because it is a consensus gene. No module gene is a marker.
MARKERS = {
    "astrocyte": ["AQP4", "GJA1", "SLC1A3", "SLC1A2", "ALDH1L1", "SLC4A4"],
    "microglia": ["P2RY12", "CX3CR1", "CSF1R", "PTPRC", "CD68", "TREM2", "TYROBP", "CD163", "ITGAX"],
    "endothelial": ["CLDN5", "PECAM1", "VWF", "FLT1"],
    "mural_fibroblast": ["PDGFRB", "DCN", "COL1A2", "LUM"],
    "oligodendrocyte": ["MBP", "PLP1", "MOBP", "MOG", "ST18"],
    "neuron": ["SNAP25", "RBFOX3", "SYT1", "NRGN", "SLC17A7", "GAD1"],
}
MIN_COUNTS = 500  # spot total counts for niche definition and co-localisation (Amendment 9)
LR_GENES = ["C3", "C3AR1", "CR1", "LRP1", "ANGPT1", "TEK", "KITLG", "KIT", "EFNA5", "EPHA4", "EGF", "ERBB4", "ANXA1", "EGFR"]

# Anatomical region groups fixed in Amendment 9.
REGION_GROUPS = {
    "LSC": {"Vent_Horn": "grey", "Dors_Horn": "grey", "Med_Grey": "grey",
            "Dors_Med_White": "white", "Vent_Med_White": "white", "Med_Lat_White": "white", "Vent_Lat_White": "white",
            "Vent_Edge": "edge", "Dors_Edge": "edge", "Lat_Edge": "edge"},
    "MTC": {"Layer_1": "L1", "Layer_2": "L2_3", "Layer_3": "L2_3", "Layer_5": "L5_6", "Layer_6": "L5_6",
            "White_matter": "WM"},
}


def module_genes() -> dict[str, list[str]]:
    seed = ROOT / "results" / "W34_synapse_resolved_regulatory_reanalysis" / "inputs" / "rna_seed_lists"
    mods = {m: pd.read_csv(seed / f"GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_{m}.tsv", sep="\t")["gene_id"].tolist()
            for m in ("up", "down")}
    sig = pd.read_csv(ROOT / "results" / "W59_NYGC_consensus_signature_holdout" / "w59_frozen_consensus_signature.tsv", sep="\t")
    mods["up_consensus"] = sig.loc[sig["direction"] == "up", "gene"].tolist()
    mods["down_consensus"] = sig.loc[sig["direction"] == "down", "gene"].tolist()
    return mods


def discovery_donors() -> set[str]:
    """All GSE330130 motor-cortex donors (used in discovery or not); excluded from the primary W66 tests."""
    audit = pd.read_csv(ROOT / "results" / "W33_synapse_resolved_diagnosis_sensitivity" / "MTC_Astrocyte_RNA.library_label_audit.tsv", sep="\t")
    return set(audit["donor"])


def log_score(counts: np.ndarray, total: np.ndarray) -> np.ndarray:
    return np.log1p(1e4 * counts / np.maximum(total, 1))


def indicator(labels: pd.Series) -> tuple[sp.csr_matrix, list[str]]:
    """Units x spots 0/1 matrix for non-missing labels."""
    keep = labels.notna().to_numpy()
    codes, uniques = pd.factorize(labels[keep])
    rows = codes
    cols = np.flatnonzero(keep)
    m = sp.csr_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)), shape=(len(uniques), len(labels)))
    return m, list(uniques)


def build(tissue: str) -> None:
    a = ad.read_h5ad(H5AD[tissue])
    obs = a.obs.reset_index(drop=True)
    sym = a.var["common_name"].astype(str).to_numpy()
    X = a.X.tocsr()
    total = np.asarray(X.sum(axis=1)).ravel()
    sym_idx = pd.Series(np.arange(len(sym)), index=sym)
    sym_idx = sym_idx[~sym_idx.index.duplicated()]
    id_idx = pd.Series(np.arange(a.n_vars), index=a.var_names)

    spots = pd.DataFrame({"section": obs["array"].astype(str), "donor": obs["Subject ID"].astype(str),
                          "annotation": obs["annotation"].astype(str), "array_x": obs["x"], "array_y": obs["y"],
                          "total_counts": total})
    for lin, genes in MARKERS.items():
        cols = sym_idx.reindex(genes).dropna().astype(int).to_numpy()
        spots[f"score_{lin}"] = log_score(np.asarray(X[:, cols].sum(axis=1)).ravel(), total)
    spots["score_vascular"] = log_score(
        np.asarray(X[:, sym_idx.reindex(MARKERS["endothelial"] + MARKERS["mural_fibroblast"]).dropna().astype(int).to_numpy()].sum(axis=1)).ravel(), total)
    mods = module_genes()
    for m, genes in mods.items():
        cols = id_idx.reindex(genes).dropna().astype(int).to_numpy()
        spots[f"module_{m}"] = log_score(np.asarray(X[:, cols].sum(axis=1)).ravel(), total)
        spots[f"module_{m}_n_genes"] = len(cols)
    up_symbols = [sym[id_idx[g]] for g in mods["up"] if g in id_idx.index]
    for g in dict.fromkeys(LR_GENES + up_symbols):
        spots[f"count_{g}"] = np.asarray(X[:, int(sym_idx[g])].todense()).ravel()

    # Lineage niches per section (Amendment 9): top quartile of the lineage score within the section.
    # Deviation (Supplementary Table S3; fixed before any disease test): markers are zero in most low-count
    # spots, so niches use spots with total counts >= MIN_COUNTS, "rich" also requires a score > 0,
    # and "below the median" is "at or below the median".
    lin_cols = ["score_astrocyte", "score_vascular", "score_microglia"]
    ok = spots["total_counts"] >= MIN_COUNTS
    q = spots[ok].groupby("section")[lin_cols]
    q75 = q.quantile(0.75).reindex(spots["section"]).set_axis(spots.index)
    q50 = q.quantile(0.50).reindex(spots["section"]).set_axis(spots.index)
    rich = {c: ok & (spots[c] >= q75[c]) & (spots[c] > 0) for c in lin_cols}
    niche = {
        "astrocyte_rich": rich["score_astrocyte"],
        "vascular_rich": rich["score_vascular"],
        "microglia_rich": rich["score_microglia"],
        "astrocyte_only": rich["score_astrocyte"] & (spots["score_vascular"] <= q50["score_vascular"])
                          & (spots["score_microglia"] <= q50["score_microglia"]),
    }
    for k, v in niche.items():
        spots[f"niche_{k}"] = v.to_numpy()
    spots.to_csv(OUT / f"{tissue}_spot_table.tsv.gz", sep="\t", index=False)

    region = spots["annotation"].map(REGION_GROUPS[tissue])
    unit_labels = {"whole": spots["donor"] + "|whole", "region": spots["donor"] + "|" + region}
    for k in niche:
        unit_labels[f"niche_{k}"] = (spots["donor"] + "|" + k).where(spots[f"niche_{k}"])
    meta_cols = {"Diagnosis": "diagnosis", "Level 1": "group", "Sex": "sex", "Age": "age", "PMI": "pmi"}
    donor_meta = obs[["Subject ID", *meta_cols]].drop_duplicates("Subject ID").rename(columns={"Subject ID": "donor", **meta_cols})
    disc = discovery_donors()
    blocks, metas = [], []
    for unit_type, labels in unit_labels.items():
        ind, names = indicator(labels.astype(object).where(labels.notna(), None))
        blocks.append(ind @ X)
        m = pd.DataFrame({"unit": [f"{unit_type}:{n}" for n in names], "unit_type": unit_type,
                          "donor": [n.split("|")[0] for n in names], "stratum": [n.split("|")[1] for n in names],
                          "n_spots": np.asarray(ind.sum(axis=1)).ravel()})
        metas.append(m)
    pb = sp.vstack(blocks).T.tocsr()
    meta = pd.concat(metas, ignore_index=True).merge(donor_meta, on="donor", how="left")
    meta["gse330130_donor"] = meta["donor"].isin(disc)
    meta["lib_size"] = np.asarray(pb.sum(axis=0)).ravel()
    counts = pd.DataFrame(pb.toarray(), index=a.var_names, columns=meta["unit"])
    counts.insert(0, "symbol", sym)
    counts.to_csv(OUT / f"{tissue}_pseudobulk_counts.tsv.gz", sep="\t")
    meta.to_csv(OUT / f"{tissue}_units.tsv", sep="\t", index=False)
    logger.info("%s: %d spots, %d sections, %d units", tissue, len(spots), spots["section"].nunique(), len(meta))


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    for tissue in sys.argv[1:] or ["LSC", "MTC"]:
        build(tissue)


if __name__ == "__main__":
    main()
