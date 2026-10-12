#!/usr/bin/env python3
"""Build the W58 Zelic 2025 (GSE287257) astrocyte pseudobulk (Amendment 4 A4.1).

Applies the fixed W54 marker-dominance astrocyte rule to all nuclei of
GSE287257_als_spinal.h5ad (raw counts, gene symbols) and sums counts per donor.
Covariates come from the GEO-workbook subject bridge.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from p01_paths import project_root, source_data_root
from w54_marker_annotation import call_astrocytes, lineage_scores

H5AD = source_data_root() / "GSE287257/files/pre/processed/GSE287257_als_spinal.h5ad"
BRIDGE = project_root() / "results/Supplemental_QC_covariate_audit/gse287257_geo_sample_workbook_subject_bridge_2026-08-12.tsv"
OUT = project_root() / "results/W58_multicohort_meta"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(H5AD)
    x = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    if not np.allclose(x.data[:100000], np.round(x.data[:100000])):
        raise ValueError("GSE287257 X is not raw integer counts")
    genes = adata.var_names.astype(str).tolist()
    totals = np.asarray(x.sum(axis=1)).ravel()
    scores, lineages = lineage_scores(x, genes, totals)
    astro = call_astrocytes(scores, lineages, totals)
    donors = adata.obs["donor_id"].astype(str).to_numpy()
    sums, rows = {}, []
    for donor in sorted(set(donors)):
        sel = (donors == donor) & astro
        sums[donor] = np.asarray(x[sel].sum(axis=0)).ravel()
        rows.append({"donor_id": donor, "disease": adata.obs.loc[donors == donor, "disease"].iloc[0],
                     "n_nuclei": int((donors == donor).sum()), "n_astrocytes": int(sel.sum())})
    counts = pd.DataFrame(sums, index=pd.Index(genes, name="gene")).astype(np.int64)
    counts = counts.groupby(level=0).sum()
    counts.to_csv(OUT / "w58_zelic_astrocyte_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")

    bridge = pd.read_csv(BRIDGE, sep="\t", dtype=str)[["local_h5ad_donor_id", "Age", "RIN", "Subject Sex", "PMI (hours)"]]
    bridge.columns = ["donor_id", "age", "RIN", "sex", "PMI"]
    samples = pd.DataFrame(rows).merge(bridge, on="donor_id", how="left")
    samples["lib_size"] = samples["donor_id"].map(counts.sum(axis=0))
    samples.to_csv(OUT / "w58_zelic_astrocyte_pseudobulk_samples.tsv", sep="\t", index=False)
    print(samples.to_string(index=False))


if __name__ == "__main__":
    main()
