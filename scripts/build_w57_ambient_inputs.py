#!/usr/bin/env python3
"""Build W57 inputs for the Amendment 3 contamination test.

For each library/channel, writes astrocyte pseudobulk counts and the pooled
all-nuclei counts (ambient proxy):
- GSE330130 MTC (discovery): library = library_name; astrocytes = Cell_Type "Astrocyte".
- Ulm: channel = Well; astrocyte counts per donor x well; ambient per well.
- Pineda BA9 batch 210526: library = sample; ambient = all nuclei in the sample.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from build_w56_lineage_pseudobulk import stratum_sums
from p01_paths import project_root, source_data_root

ROOT = project_root()
OUT = ROOT / "results/W57_ambient_contamination_test"
DISCOVERY_H5AD = source_data_root().parent / "GSE330130_MTC_SNRNASEQ_COUNTS.h5ad"
ULM = source_data_root() / "Other_Datasets/Multiome_Dataset/files"
PINEDA_DIR = source_data_root() / "Other_Datasets/Pineda_snRNAseq/synapse_download/Processed Data"
W56 = ROOT / "results/W56_cell_source_check"


def build_discovery() -> None:
    adata = ad.read_h5ad(DISCOVERY_H5AD)
    x = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    obs = adata.obs
    astro, ambient = {}, {}
    for lib, idx in obs.groupby("library_name", observed=True).indices.items():
        idx = np.asarray(idx)
        is_astro = obs["Cell_Type"].to_numpy()[idx] == "Astrocyte"
        ambient[lib] = np.asarray(x[idx].sum(axis=0)).ravel()
        astro[lib] = np.asarray(x[idx[is_astro]].sum(axis=0)).ravel()
    genes = pd.Index(adata.var_names, name="gene_id")
    pd.DataFrame(astro, index=genes).astype(np.int64).to_csv(OUT / "w57_discovery_astro_counts.tsv.gz", sep="\t", compression="gzip")
    pd.DataFrame(ambient, index=genes).astype(np.int64).to_csv(OUT / "w57_discovery_ambient_counts.tsv.gz", sep="\t", compression="gzip")
    adata.var[["common_name"]].to_csv(OUT / "w57_discovery_gene_annotation.tsv", sep="\t")
    print(f"[W57] discovery: {len(astro)} libraries", flush=True)


def build_ulm() -> None:
    genes = pd.read_csv(ULM / "Multiome_Dataset_RNA_features.tsv", sep="\t", header=None)[0].astype(str)
    meta = pd.read_csv(ULM / "Multiome_Dataset_Metadata.txt", sep="\t", usecols=["Well", "ID", "WNN_L2"])
    astro = meta["WNN_L2"].eq("Astrocytes")
    meta["astro_stratum"] = np.where(astro, meta["ID"] + "|" + meta["Well"], "")
    astro_strata = sorted(meta.loc[astro, "astro_stratum"].unique())
    wells = sorted(meta["Well"].unique())
    # One pass: astrocyte donor x well strata first, then one ambient stratum per well.
    strata = astro_strata + [f"ambient|{w}" for w in wells]
    index = {s: i for i, s in enumerate(strata)}
    matrix = ULM / "Multiome_Dataset_RNA_counts_raw.mtx"
    astro_map = np.full(len(meta) + 1, -1, dtype=np.int64)
    astro_map[1:] = meta["astro_stratum"].map(index).fillna(-1).astype(np.int64).to_numpy()
    well_map = np.full(len(meta) + 1, -1, dtype=np.int64)
    well_map[1:] = ("ambient|" + meta["Well"]).map(index).to_numpy()
    astro_sums = stratum_sums(matrix, astro_map, len(strata))
    well_sums = stratum_sums(matrix, well_map, len(strata))
    counts = pd.DataFrame(astro_sums + well_sums, index=genes, columns=strata).round().astype(np.int64)
    counts = counts.groupby(level=0).sum()
    counts.index.name = "gene"
    counts.to_csv(OUT / "w57_ulm_astro_and_ambient_counts.tsv.gz", sep="\t", compression="gzip")
    print(f"[W57] Ulm: {len(astro_strata)} astro donor x well strata, {len(wells)} wells", flush=True)


def build_pineda() -> None:
    strata = pd.read_csv(W56 / "w56_pineda_lineage_strata.tsv", sep="\t")
    samples = sorted(strata.loc[strata["Region"].eq("BA9") & strata["sample_folder"].str.startswith("210526"), "sample_folder"].unique())
    ambient, genes_ref = {}, None
    for name in samples:
        path = Path(PINEDA_DIR) / name
        genes = pd.read_csv(path / "row_metadata.tsv", sep="\t")["ENSEMBL"].astype(str).tolist()
        genes_ref = genes_ref or genes
        n_cells = len(pd.read_csv(path / "col_metadata.tsv", sep="\t", usecols=[0]))
        col_map = np.zeros(n_cells + 1, dtype=np.int64)
        col_map[0] = -1
        ambient[name] = stratum_sums(path / "counts_fil.mtx", col_map, 1)[:, 0]
    pd.DataFrame(ambient, index=pd.Index(genes_ref, name="gene_id")).round().astype(np.int64).to_csv(
        OUT / "w57_pineda_ba9_210526_ambient_counts.tsv.gz", sep="\t", compression="gzip")
    print(f"[W57] Pineda BA9 210526: {len(samples)} samples", flush=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    build_discovery()
    build_pineda()
    build_ulm()


if __name__ == "__main__":
    main()
