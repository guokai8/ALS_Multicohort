#!/usr/bin/env python3
"""Export W03 ATAC pseudobulk inputs for donor-level DA modeling."""

from __future__ import annotations

import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from p01_paths import project_root


P01 = project_root()
PSEUDOBULK = P01 / "pseudobulk"
OUT = P01 / "intermediate/W03_ATAC_DA"


INPUTS = {
    "GSE330130_MTC_ATAC": PSEUDOBULK / "GSE330130_MTC_ATAC.main_candidates.pseudobulk.h5ad",
    "GSE330130_LSC_ATAC": PSEUDOBULK / "GSE330130_LSC_ATAC.main_candidates.pseudobulk.h5ad",
}


def project_rel(path: Path) -> str:
    try:
        return path.relative_to(P01).as_posix()
    except ValueError:
        return path.name


def safe_name(value: str) -> str:
    value = str(value).strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "NA"


def dense_counts_for_rows(adata: ad.AnnData, row_mask: np.ndarray) -> np.ndarray:
    x = adata.X[row_mask, :]
    if sp.issparse(x):
        x = x.toarray()
    x = np.asarray(x)
    if not np.issubdtype(x.dtype, np.integer):
        x = np.rint(x).astype(np.int64)
    return x


def feature_table(adata: ad.AnnData) -> pd.DataFrame:
    features = pd.DataFrame(index=adata.var_names.copy())
    features.index.name = "peak_id"
    features = features.reset_index()
    for col in ["id", "genome", "interval", "feature_type", "chrom", "start", "end", "common_name", "n_counts"]:
        if col in adata.var.columns:
            features[col] = adata.var[col].astype("string").fillna("").astype(str).to_numpy()
    if "interval" not in features.columns:
        features["interval"] = features["peak_id"]
    return features


def export_dataset(dataset_id: str, h5ad_path: Path, manifest_rows: list[dict[str, str]]) -> None:
    adata = ad.read_h5ad(h5ad_path)
    obs = adata.obs.copy()
    obs["pseudobulk_id"] = obs.index.astype(str)

    if "primary_model_als_control" in obs.columns:
        primary = obs["primary_model_als_control"].astype(str).str.lower().isin(["true", "1", "yes"])
    else:
        primary = pd.Series(True, index=obs.index)
    disease_ok = obs["disease_binary"].astype(str).isin(["ALS", "Control"])
    export_obs = obs[primary & disease_ok].copy()

    excluded = obs[~disease_ok | ~primary].copy()
    if not excluded.empty:
        excluded_out = OUT / dataset_id / "excluded_non_primary_disease.tsv"
        excluded_out.parent.mkdir(parents=True, exist_ok=True)
        excluded_cols = [
            "dataset_id",
            "donor",
            "library",
            "analysis_region",
            "tissue_subregion",
            "disease_raw",
            "disease_binary",
            "cell_type",
            "pseudobulk_id",
            "n_cells",
        ]
        excluded[[c for c in excluded_cols if c in excluded.columns]].to_csv(
            excluded_out, sep="\t", index=False
        )

    features = feature_table(adata)
    dataset_dir = OUT / dataset_id
    dataset_dir.mkdir(parents=True, exist_ok=True)
    features.to_csv(dataset_dir / "peaks.tsv", sep="\t", index=False)

    for cell_type, cell_obs in export_obs.groupby("cell_type", sort=True, observed=True):
        cell_obs = cell_obs.copy()
        cell_obs = cell_obs[cell_obs["disease_binary"].astype(str).isin(["ALS", "Control"])]
        if cell_obs.empty:
            continue

        row_mask = np.asarray(obs.index.isin(cell_obs.index))
        counts = dense_counts_for_rows(adata, row_mask)

        sample_cols = [f"S{i + 1:03d}" for i in range(counts.shape[0])]
        cell_obs["sample_col"] = sample_cols
        cell_obs["sample_id"] = cell_obs["pseudobulk_id"].astype(str)

        sample_cols_keep = [
            "sample_col",
            "sample_id",
            "dataset_id",
            "study",
            "modality",
            "donor",
            "sample",
            "library",
            "analysis_region",
            "tissue_subregion",
            "disease_raw",
            "disease_binary",
            "clinical_subgroup",
            "mutation",
            "sex",
            "age",
            "PMI",
            "cell_type",
            "n_cells",
        ]
        for col in sample_cols_keep:
            if col not in cell_obs.columns:
                cell_obs[col] = pd.NA
        samples = cell_obs[sample_cols_keep].copy()

        counts_df = pd.DataFrame(counts.T, columns=sample_cols)
        counts_df.insert(0, "peak_id", features["peak_id"].to_numpy())

        cell_dir = dataset_dir / safe_name(cell_type)
        cell_dir.mkdir(parents=True, exist_ok=True)
        counts_path = cell_dir / "counts.tsv.gz"
        samples_path = cell_dir / "samples.tsv"

        counts_df.to_csv(counts_path, sep="\t", index=False, compression="gzip")
        samples.to_csv(samples_path, sep="\t", index=False)

        n_by_group = samples.groupby("disease_binary", observed=True)["donor"].nunique().to_dict()
        manifest_rows.append(
            {
                "dataset_id": dataset_id,
                "cell_type": str(cell_type),
                "cell_type_safe": safe_name(cell_type),
                "counts_path": project_rel(counts_path),
                "samples_path": project_rel(samples_path),
                "peaks_path": project_rel(dataset_dir / "peaks.tsv"),
                "n_samples": str(samples.shape[0]),
                "n_peaks": str(counts.shape[1]),
                "n_als_donors": str(n_by_group.get("ALS", 0)),
                "n_control_donors": str(n_by_group.get("Control", 0)),
            }
        )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, str]] = []
    for dataset_id, h5ad_path in INPUTS.items():
        if not h5ad_path.exists():
            raise FileNotFoundError(h5ad_path)
        export_dataset(dataset_id, h5ad_path, manifest_rows)

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(OUT / "w03_export_manifest.tsv", sep="\t", index=False)
    manifest[
        ["dataset_id", "cell_type", "n_samples", "n_peaks", "n_als_donors", "n_control_donors"]
    ].sort_values(["dataset_id", "cell_type"]).to_csv(OUT / "w03_export_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
