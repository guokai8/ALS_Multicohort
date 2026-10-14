#!/usr/bin/env python3
"""W72 (Amendment 14): disease-state lineage share of the frozen ALS-up genes in human spinal cord.

GSE287257 (Zelic et al.; 8 ALS, 4 control cervical cords), raw counts. Every nucleus is assigned with
the W54 marker-dominance rule generalised to all six lineages: the top-scoring lineage, if it exceeds
the second by >= 0.5 and the nucleus has >= 500 UMIs; others are unassigned. Counts are summed per
donor x lineage; the share of each gene is the mean CPM over donors of a group in a lineage divided by
the sum over lineages. The within-lineage ALS vs control model is fitted by run_w72_lineage_de.R.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from p01_paths import project_root, source_data_root
from w54_marker_annotation import MARGIN, MIN_UMI, lineage_scores

logger = logging.getLogger(__name__)

H5AD = source_data_root() / "GSE287257/files/pre/processed/GSE287257_als_spinal.h5ad"
OUT = project_root() / "results" / "W72_zelic_disease_cell_source"
SEED = project_root() / "results/W34_synapse_resolved_regulatory_reanalysis/inputs/rna_seed_lists"
FOCUS = ["GPNMB", "C3", "ANXA1", "LTF", "EMP3"]


def call_lineages(scores: np.ndarray, lineages: list[str], totals: np.ndarray) -> np.ndarray:
    ordered = np.sort(scores, axis=1)
    ok = ((ordered[:, -1] - ordered[:, -2]) >= MARGIN) & (totals >= MIN_UMI)
    return np.where(ok, np.array(lineages, dtype=object)[scores.argmax(axis=1)], "unassigned")


def share_reading(row: pd.Series) -> str:
    if row["Microglia"] >= 0.5:
        return "myeloid in ALS spinal cord"
    if row["Astrocyte"] >= 0.5 or row["Neuron"] >= 0.5:
        return "astrocytic or neuronal (consistent with Tanaka et al.)"
    return "mixed"


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    (OUT / "inputs").mkdir(parents=True, exist_ok=True)
    up = pd.read_csv(SEED / "GSE330130_MTC_RNA__Astrocyte.FDR05_logFC025_up.tsv", sep="\t")["common_name"].tolist()
    a = ad.read_h5ad(H5AD)
    x = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    if not np.allclose(x.data[:100000], np.round(x.data[:100000])):
        raise ValueError("GSE287257 X is not raw integer counts")
    genes = a.var_names.astype(str).tolist()
    totals = np.asarray(x.sum(axis=1)).ravel()
    scores, lineages = lineage_scores(x, genes, totals)
    lin = call_lineages(scores, lineages, totals)
    donors = a.obs["donor_id"].astype(str).to_numpy()
    disease = a.obs.groupby("donor_id", observed=True)["disease"].first().astype(str)

    sums, rows = {}, []
    for d in sorted(set(donors)):
        for l in lineages:
            sel = (donors == d) & (lin == l)
            sums[f"{d}__{l}"] = np.asarray(x[sel].sum(axis=0)).ravel()
            rows.append({"unit": f"{d}__{l}", "donor": d, "disease": disease[d], "lineage": l, "n_nuclei": int(sel.sum())})
    counts = pd.DataFrame(sums, index=pd.Index(genes, name="gene")).groupby(level=0).sum().astype(np.int64)
    units = pd.DataFrame(rows)
    units["lib_size"] = units["unit"].map(counts.sum(axis=0))
    counts.to_csv(OUT / "inputs" / "w72_lineage_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")
    units.to_csv(OUT / "inputs" / "w72_lineage_units.tsv", sep="\t", index=False)
    calls = pd.crosstab(pd.Series(donors, name="donor"), pd.Series(lin, name="lineage"))

    # disease-state lineage share (units with >= 20 nuclei)
    keep = units[units["n_nuclei"] >= 20]
    cpm = counts[keep["unit"]].div(counts[keep["unit"]].sum(axis=0), axis=1) * 1e6
    g = [s for s in up if s in cpm.index]
    long = cpm.loc[g].T.join(keep.set_index("unit")[["disease", "lineage"]])
    mean = long.groupby(["disease", "lineage"]).mean().T  # genes x (disease, lineage)
    share = {}
    for grp in ("ALS", "Control"):
        m = mean[grp]
        share[grp] = m.div(m.sum(axis=1), axis=0)
    als = share["ALS"].copy()
    als["top_lineage"] = als.idxmax(axis=1)
    als["reading"] = als.apply(share_reading, axis=1)
    shift = (share["ALS"] - share["Control"]).add_suffix("_shift")
    table = als.join(share["Control"].add_suffix("_control")).join(shift).rename_axis("gene").reset_index()

    table.to_csv(OUT / "w72_lineage_share.tsv", sep="\t", index=False)
    calls.to_csv(OUT / "w72_lineage_calls_by_donor.tsv", sep="\t")
    logger.info("lineage calls:\n%s", pd.Series(lin).value_counts().to_string())
    cols = ["gene"] + lineages + ["top_lineage", "reading"]
    logger.info("ALS-donor lineage share:\n%s", table[cols].round(2).to_string(index=False))
    logger.info("shift ALS - control (focus genes):\n%s",
                table.set_index("gene").loc[[f for f in FOCUS if f in table["gene"].values], [f"{l}_shift" for l in lineages]].round(2).to_string())


if __name__ == "__main__":
    main()
