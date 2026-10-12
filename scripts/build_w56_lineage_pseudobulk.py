#!/usr/bin/env python3
"""Build W56 all-lineage pseudobulk counts for Pineda and Ulm (Amendment 2, E1/E2).

Pineda: raw counts summed per human sample x broad lineage (W08 mapping; T
cells dropped; mouse BL6 samples excluded).
Ulm: raw counts summed per donor x chip x WNN_L2 lineage.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from p01_paths import project_root, source_data_root
from run_w08_pineda_validation import broad_celltype, resolve_source_path

OUT = project_root() / "results/W56_cell_source_check"
PINEDA_MANIFEST = project_root() / "metadata/pineda_sample_manifest.tsv"
ULM = source_data_root() / "Other_Datasets/Multiome_Dataset/files"
PINEDA_LINEAGES = ["Astrocyte", "Microglia", "Endothelial", "Vascular/Meningeal", "Oligodendrocyte", "OPC",
                   "Excitatory neuron", "Inhibitory neuron"]
CHUNK_ROWS = 25_000_000


def matrix_header(path: Path) -> tuple[int, int, int]:
    n_header = 0
    with open(path) as handle:
        for line in handle:
            n_header += 1
            if not line.startswith("%"):
                n_rows, n_cols = map(int, line.split()[:2])
                return n_header, n_rows, n_cols
    raise ValueError(f"No dimension line in {path}")


def stratum_sums(matrix: Path, col_to_stratum: np.ndarray, n_strata: int) -> np.ndarray:
    """Sum a genes x cells MatrixMarket file into genes x strata (one-based columns)."""
    n_header, n_rows, _ = matrix_header(matrix)
    sums = np.zeros(n_rows * n_strata, dtype=np.float64)
    reader = pd.read_csv(matrix, sep=" ", header=None, skiprows=n_header, names=["row", "col", "value"],
                         dtype={"row": np.int64, "col": np.int64, "value": np.float64}, chunksize=CHUNK_ROWS, engine="c")
    for chunk in reader:
        s = col_to_stratum[chunk["col"].to_numpy()]
        keep = s >= 0
        flat = (chunk["row"].to_numpy()[keep] - 1) * n_strata + s[keep]
        sums += np.bincount(flat, weights=chunk["value"].to_numpy()[keep], minlength=sums.size)
    return sums.reshape(n_rows, n_strata)


def pineda_sample(sample_dir: str) -> tuple[str, list[str], np.ndarray, dict[str, int]]:
    path = Path(sample_dir)
    genes = pd.read_csv(path / "row_metadata.tsv", sep="\t")["ENSEMBL"].astype(str).tolist()
    col = pd.read_csv(path / "col_metadata.tsv", sep="\t")
    lineage = col.apply(broad_celltype, axis=1)
    idx = {lin: i for i, lin in enumerate(PINEDA_LINEAGES)}
    col_to_stratum = np.full(len(col) + 1, -1, dtype=np.int64)
    col_to_stratum[1:] = lineage.map(idx).fillna(-1).astype(np.int64).to_numpy()
    sums = stratum_sums(path / "counts_fil.mtx", col_to_stratum, len(PINEDA_LINEAGES))
    return path.name, genes, sums, lineage.value_counts().to_dict()


def build_pineda(workers: int) -> None:
    manifest = pd.read_csv(PINEDA_MANIFEST, sep="\t")
    manifest = manifest[manifest["counts_fil_mtx_exists"].astype(str).eq("True") & ~manifest["sample_folder"].str.contains("_BL6_")]
    jobs = [str(resolve_source_path(p)) for p in manifest["path"]]
    columns, rows, genes_ref = {}, [], None
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for name, genes, sums, n_cells in pool.map(pineda_sample, jobs):
            if genes_ref is None:
                genes_ref = genes
            elif genes != genes_ref:
                raise ValueError(f"Gene order differs in {name}")
            for j, lin in enumerate(PINEDA_LINEAGES):
                columns[f"{name}|{lin}"] = sums[:, j]
                rows.append({"stratum": f"{name}|{lin}", "sample_folder": name, "lineage": lin, "n_cells": int(n_cells.get(lin, 0))})
    counts = pd.DataFrame(columns, index=pd.Index(genes_ref, name="gene_id")).round().astype(np.int64)
    counts.to_csv(OUT / "w56_pineda_lineage_counts.tsv.gz", sep="\t", compression="gzip")
    strata = pd.DataFrame(rows).merge(manifest[["sample_folder", "Donor", "Region", "Sex", "Condition", "Group"]], on="sample_folder")
    strata["lib_size"] = strata["stratum"].map(counts.sum(axis=0))
    strata.to_csv(OUT / "w56_pineda_lineage_strata.tsv", sep="\t", index=False)
    print(f"[W56] Pineda: {len(strata)} sample x lineage strata", flush=True)


def build_ulm() -> None:
    genes = pd.read_csv(ULM / "Multiome_Dataset_RNA_features.tsv", sep="\t", header=None)[0].astype(str)
    meta = pd.read_csv(ULM / "Multiome_Dataset_Metadata.txt", sep="\t", usecols=["Well", "ID", "Sex", "Case", "Case_Type", "WNN_L2"])
    meta["Chip"] = meta["Well"].str.extract(r"^(Chip\d+)")[0]
    meta["stratum"] = meta["ID"] + "|" + meta["Chip"] + "|" + meta["WNN_L2"]
    strata = sorted(meta["stratum"].unique())
    index = {s: i for i, s in enumerate(strata)}
    col_to_stratum = np.full(len(meta) + 1, -1, dtype=np.int64)
    col_to_stratum[1:] = meta["stratum"].map(index).to_numpy()
    sums = stratum_sums(ULM / "Multiome_Dataset_RNA_counts_raw.mtx", col_to_stratum, len(strata))
    counts = pd.DataFrame(sums, index=genes, columns=strata).round().astype(np.int64)
    counts = counts.groupby(level=0).sum()
    counts.index.name = "gene"
    counts.to_csv(OUT / "w56_ulm_lineage_counts.tsv.gz", sep="\t", compression="gzip")
    info = meta.groupby("stratum").agg(ID=("ID", "first"), Chip=("Chip", "first"), lineage=("WNN_L2", "first"),
                                       Sex=("Sex", "first"), Case=("Case", "first"), Case_Type=("Case_Type", "first"),
                                       n_cells=("ID", "size")).reset_index()
    info["lib_size"] = info["stratum"].map(counts.sum(axis=0))
    info.to_csv(OUT / "w56_ulm_lineage_strata.tsv", sep="\t", index=False)
    print(f"[W56] Ulm: {len(info)} donor x chip x lineage strata", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    build_pineda(args.workers)
    build_ulm()


if __name__ == "__main__":
    main()
