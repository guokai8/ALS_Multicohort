#!/usr/bin/env python3
"""Build W51 Pineda astrocyte pseudobulk counts for discovery-matched limma-voom.

Sums raw counts over all genes for Pineda astrocyte nuclei (CellType == "Astro")
per sample, matching the discovery donor/library/region/cell-type pseudobulk
unit. No minimum-cell threshold is applied here, matching discovery; n_cells is
recorded so the model script can run a >=20-cell sensitivity. Mouse M1 samples
(220616_BL6_MCX_*) are excluded because they do not share the human gene space.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from p01_paths import project_root
from run_w08_pineda_validation import broad_celltype, resolve_source_path

ROOT = project_root()
MANIFEST = ROOT / "metadata/pineda_sample_manifest.tsv"
CLINICAL = ROOT / "results/W22_Pineda_S01_clinical_pathology_context/w22_pineda_lead_module_donor_context.tsv"
OUT = ROOT / "results/W51_Pineda_BA4_discovery_matched_limma"
CHUNK_ROWS = 20_000_000


def astrocyte_gene_sums(sample_dir: str) -> dict[str, object]:
    sample_path = Path(sample_dir)
    row_meta = pd.read_csv(sample_path / "row_metadata.tsv", sep="\t")
    col_meta = pd.read_csv(sample_path / "col_metadata.tsv", sep="\t")
    col_meta["broad_cell_type"] = col_meta.apply(broad_celltype, axis=1)
    meta_cols = ["Sample_ID", "Donor", "Region", "Sex", "Condition", "Group", "C9_pos"]
    meta = {c: col_meta[c].iloc[0] if c in col_meta.columns else np.nan for c in meta_cols}

    is_astro = np.zeros(col_meta.shape[0] + 1, dtype=bool)
    is_astro[np.flatnonzero(col_meta["broad_cell_type"].eq("Astrocyte").to_numpy()) + 1] = True

    n_header = 0
    with open(sample_path / "counts_fil.mtx") as handle:
        for line in handle:
            n_header += 1
            if not line.startswith("%"):
                n_rows, n_cols = map(int, line.split()[:2])
                break
    if n_rows != row_meta.shape[0] or n_cols != col_meta.shape[0]:
        raise ValueError(f"Dimension mismatch in {sample_path}")

    gene_sums = np.zeros(n_rows, dtype=np.float64)
    cell_totals = np.zeros(n_cols + 1, dtype=np.float64)
    reader = pd.read_csv(
        sample_path / "counts_fil.mtx", sep=" ", header=None, skiprows=n_header,
        names=["row", "col", "value"], dtype={"row": np.int64, "col": np.int64, "value": np.float64},
        chunksize=CHUNK_ROWS, engine="c",
    )
    for chunk in reader:
        cols = chunk["col"].to_numpy()
        keep = is_astro[cols]
        vals = chunk["value"].to_numpy()[keep]
        gene_sums += np.bincount(chunk["row"].to_numpy()[keep] - 1, weights=vals, minlength=n_rows)
        cell_totals += np.bincount(cols[keep], weights=vals, minlength=n_cols + 1)
    return {
        "sample_folder": sample_path.name, **meta,
        "n_cells": int((cell_totals[is_astro] > 0).sum()),
        "lib_size": float(gene_sums.sum()),
        "genes": row_meta["ENSEMBL"].astype(str).tolist(),
        "gene_sums": gene_sums,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(MANIFEST, sep="\t")
    manifest = manifest[manifest["counts_fil_mtx_exists"].astype(str).eq("True")]
    manifest = manifest[~manifest["sample_folder"].str.contains("_BL6_")]
    jobs = [str(resolve_source_path(p)) for p in manifest["path"]]

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        records = []
        for i, rec in enumerate(pool.map(astrocyte_gene_sums, jobs), 1):
            records.append(rec)
            print(f"[W51] {i}/{len(jobs)} {rec['sample_folder']} astro_cells={rec['n_cells']}", flush=True)

    genes = records[0]["genes"]
    if any(rec["genes"] != genes for rec in records):
        raise ValueError("Pineda samples do not share one gene order")
    counts = pd.DataFrame({rec["sample_folder"]: rec["gene_sums"] for rec in records}, index=genes)
    counts.index.name = "gene_id"
    counts.round().astype(np.int64).to_csv(OUT / "w51_pineda_astrocyte_pseudobulk_counts.tsv.gz", sep="\t", compression="gzip")

    samples = pd.DataFrame([{k: v for k, v in rec.items() if k not in {"genes", "gene_sums"}} for rec in records])
    clinical = pd.read_csv(CLINICAL, sep="\t").drop_duplicates("Donor")[["Donor", "PMI_h", "Age of death (Y)"]]
    clinical = clinical.rename(columns={"Age of death (Y)": "age", "PMI_h": "PMI"})
    samples["Donor"] = samples["Donor"].astype(str)
    clinical["Donor"] = clinical["Donor"].astype(str)
    samples = samples.merge(clinical, on="Donor", how="left")
    samples.to_csv(OUT / "w51_pineda_astrocyte_pseudobulk_samples.tsv", sep="\t", index=False)
    gene_names = pd.read_csv(Path(jobs[0]) / "row_metadata.tsv", sep="\t")[["ENSEMBL", "Gene"]]
    gene_names.to_csv(OUT / "w51_pineda_gene_annotation.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
