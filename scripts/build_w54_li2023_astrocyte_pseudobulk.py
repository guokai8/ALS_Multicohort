#!/usr/bin/env python3
"""Build W54 Li et al. 2023 motor-cortex astrocyte pseudobulk counts.

Source: GSE219280 (Li et al. 2023 Nat Commun, PMID 37714849), per-sample
CellBender-corrected filtered 10x h5 files for primary motor cortex (BA4):
6 C9-ALS, 5 C9-FTD and 6 controls. GEO provides no cell annotation, so
astrocytes are called with the fixed marker-dominance rule in
`w54_marker_annotation.py` (precision 0.94 against Ulm WNN labels; no module
gene is used as a marker). Counts are summed over called astrocytes per sample
(one sample per donor), with ENSG IDs as rows.
"""

from __future__ import annotations

import io
import zipfile

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from p01_paths import project_root, source_data_root
from w54_marker_annotation import MARKERS, call_astrocytes, lineage_scores

SOURCE = source_data_root() / "GSE219280"
H5_DIR = SOURCE / "motor_cortex_h5"
META_ZIP = SOURCE / "files/ResearchD_root/GSE219280.zip"
OUT = project_root() / "results/W54_Li2023_motor_cortex_validation"


def read_cellbender_h5(path) -> tuple[sp.csr_matrix, list[str], list[str]]:
    """Return a cells x genes CSR matrix with ENSG IDs and gene symbols."""
    with h5py.File(path, "r") as h5:
        m = h5["matrix"]
        n_genes, n_cells = (int(x) for x in m["shape"][:])
        csc = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=(n_genes, n_cells))
        ids = [x.decode() for x in m["features/id"][:]]
        names = [x.decode() for x in m["features/name"][:]]
    return csc.T.tocsr(), ids, names


def load_metadata() -> pd.DataFrame:
    with zipfile.ZipFile(META_ZIP) as z:
        sub = pd.read_csv(io.BytesIO(z.read("GSE219280.subject_metadata.csv")))
        bio = pd.read_csv(io.BytesIO(z.read("GSE219280.biospecimen_metadata.csv")))
    bio = bio[bio["CNSRegion"].eq("primary motor cortex")].drop_duplicates("GSM")[["GSM", "originalSubjectId", "PMI"]]
    meta = bio.merge(sub[["originalSubjectId", "age", "sex", "disease"]].drop_duplicates(), on="originalSubjectId")
    meta["PMI"] = pd.to_numeric(meta["PMI"], errors="coerce")
    return meta.rename(columns={"originalSubjectId": "donor"})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    meta = load_metadata().set_index("GSM")
    files = sorted(H5_DIR.glob("*_MotorCortex_cellBender_corrected_filtered.h5"))
    if len(files) != len(meta):
        raise ValueError(f"Expected {len(meta)} motor-cortex h5 files, found {len(files)}")

    sums, rows, gene_ids, symbols = {}, [], None, None
    for path in files:
        gsm = path.name.split("_")[0]
        x, ids, names = read_cellbender_h5(path)
        if gene_ids is None:
            gene_ids, symbols = ids, names
        elif ids != gene_ids:
            raise ValueError(f"Gene order differs in {path.name}")
        totals = np.asarray(x.sum(axis=1)).ravel()
        scores, lineages = lineage_scores(x, names, totals)
        astro = call_astrocytes(scores, lineages, totals)
        best = np.array(lineages)[scores.argmax(axis=1)]
        sums[gsm] = np.asarray(x[astro].sum(axis=0)).ravel()
        rows.append({
            "GSM": gsm, "n_nuclei": x.shape[0], "n_nuclei_min_umi": int((totals >= 500).sum()),
            "n_astrocytes": int(astro.sum()), "astro_median_umi": float(np.median(totals[astro])) if astro.any() else np.nan,
            **{f"frac_best_{lin}": float(np.mean(best[totals >= 500] == lin)) for lin in MARKERS},
        })
        print(f"[W54] {path.name}: nuclei={x.shape[0]} astrocytes={int(astro.sum())}", flush=True)

    counts = pd.DataFrame(sums, index=pd.Index(gene_ids, name="gene_id")).round().astype(np.int64)
    counts.to_csv(OUT / "w54_li2023_astrocyte_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")
    pd.DataFrame({"gene_id": gene_ids, "symbol": symbols}).to_csv(OUT / "w54_li2023_gene_annotation.tsv", sep="\t", index=False)
    samples = pd.DataFrame(rows).set_index("GSM").join(meta).reset_index()
    samples["lib_size"] = samples["GSM"].map(counts.sum(axis=0))
    samples.to_csv(OUT / "w54_li2023_astrocyte_pseudobulk_samples.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
