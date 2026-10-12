#!/usr/bin/env python3
"""Build W52 Ulm/DZNE motor-cortex astrocyte pseudobulk counts.

Source: Zenodo 10.5281/zenodo.18370646 (Ruf et al. 2026, Nat Commun), human
primary motor cortex 10x multiome; RNA raw counts are genes x cells with gene
symbols. Astrocytes are WNN_L2 == "Astrocytes" from the released metadata.

Raw counts are summed per donor x chip stratum (donors were run on up to two
chips), so the model can block on chip and account for repeated donor profiles.
No minimum-cell threshold is applied, matching discovery. Per-stratum astrocyte
QC proxies (mean percent_mito, median nCount_RNA) are recorded because the
clinical table (age, PMI, RIN) uses an ID scheme that does not map to the
released metadata IDs.
"""

from __future__ import annotations


import numpy as np
import pandas as pd

from p01_paths import project_root, source_data_root

SOURCE = source_data_root() / "Other_Datasets/Multiome_Dataset/files"
OUT = project_root() / "results/W52_Ulm_motor_cortex_validation"
CHUNK_ROWS = 25_000_000


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(SOURCE / "Multiome_Dataset_RNA_features.tsv", sep="\t", header=None)[0].astype(str)
    barcodes = pd.read_csv(SOURCE / "Multiome_Dataset_RNA_barcodes.tsv", sep="\t", header=None)[0]
    meta = pd.read_csv(
        SOURCE / "Multiome_Dataset_Metadata.txt", sep="\t",
        usecols=["CellId", "Well", "ID", "Sex", "Case", "C9ORF72", "Case_Type", "percent_mito", "nCount_RNA", "WNN_L2"],
    )
    if not (meta["CellId"].to_numpy() == barcodes.to_numpy()).all():
        raise ValueError("Metadata CellId order does not match RNA barcodes")

    meta["Chip"] = meta["Well"].str.extract(r"^(Chip\d+)")[0]
    astro = meta["WNN_L2"].eq("Astrocytes").to_numpy()
    meta["stratum"] = meta["ID"] + "|" + meta["Chip"]
    strata = sorted(meta.loc[astro, "stratum"].unique())
    stratum_index = {s: i for i, s in enumerate(strata)}
    # One-based column index -> stratum index (-1 for non-astrocytes).
    col_to_stratum = np.full(len(meta) + 1, -1, dtype=np.int64)
    col_to_stratum[np.flatnonzero(astro) + 1] = meta.loc[astro, "stratum"].map(stratum_index).to_numpy()

    matrix = SOURCE / "Multiome_Dataset_RNA_counts_raw.mtx"
    n_header = 0
    with open(matrix) as handle:
        for line in handle:
            n_header += 1
            if not line.startswith("%"):
                n_rows, n_cols = map(int, line.split()[:2])
                break
    if n_rows != len(genes) or n_cols != len(barcodes):
        raise ValueError("MatrixMarket dimensions do not match features/barcodes")

    n_strata = len(strata)
    sums = np.zeros(n_rows * n_strata, dtype=np.float64)
    reader = pd.read_csv(
        matrix, sep=" ", header=None, skiprows=n_header, names=["row", "col", "value"],
        dtype={"row": np.int64, "col": np.int64, "value": np.float64}, chunksize=CHUNK_ROWS, engine="c",
    )
    for i, chunk in enumerate(reader, 1):
        s = col_to_stratum[chunk["col"].to_numpy()]
        keep = s >= 0
        flat = (chunk["row"].to_numpy()[keep] - 1) * n_strata + s[keep]
        sums += np.bincount(flat, weights=chunk["value"].to_numpy()[keep], minlength=sums.size)
        print(f"[W52] chunk {i} rows={len(chunk)}", flush=True)

    counts = pd.DataFrame(sums.reshape(n_rows, n_strata), index=genes, columns=strata).round().astype(np.int64)
    if counts.index.duplicated().any():
        counts = counts.groupby(level=0).sum()
    counts.index.name = "gene"
    counts.to_csv(OUT / "w52_ulm_astrocyte_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")

    a = meta[astro]
    samples = a.groupby("stratum").agg(
        ID=("ID", "first"), Chip=("Chip", "first"), Sex=("Sex", "first"), Case=("Case", "first"),
        Case_Type=("Case_Type", "first"), C9ORF72=("C9ORF72", "first"), n_cells=("CellId", "size"),
        mean_percent_mito=("percent_mito", "mean"), median_nCount_RNA=("nCount_RNA", "median"),
    ).reset_index()
    samples["lib_size"] = samples["stratum"].map(counts.sum(axis=0))
    samples.to_csv(OUT / "w52_ulm_astrocyte_pseudobulk_samples.tsv", sep="\t", index=False)
    print(f"[W52] wrote {n_strata} donor x chip strata, {len(counts)} genes", flush=True)


if __name__ == "__main__":
    main()
