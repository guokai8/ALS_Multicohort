#!/usr/bin/env python3
"""Construct donor/library-level pseudobulk count objects for P01.

This script only creates pseudobulk-ready count-sum AnnData objects for frozen
main-candidate strata. It does not run DE, DA, motif, normalization, or result
interpretation.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from p01_paths import project_root, source_data_root


PROJECT_ROOT = project_root()
SOURCE_DATA_ROOT = source_data_root()
META = PROJECT_ROOT / "metadata"
OUT = PROJECT_ROOT / "pseudobulk"


OBJECTS = [
    (
        "GSE330130_MTC_RNA",
        "GSE330130",
        "RNA",
        SOURCE_DATA_ROOT / "GSE330130/GSE330130_MTC_SNRNASEQ_COUNTS.h5ad",
    ),
    (
        "GSE330130_MTC_ATAC",
        "GSE330130",
        "ATAC",
        SOURCE_DATA_ROOT / "GSE330130/GSE330130_MTC_SNATACSEQ_COUNTS.h5ad",
    ),
    (
        "GSE330130_LSC_RNA",
        "GSE330130",
        "RNA",
        SOURCE_DATA_ROOT / "GSE330130/GSE330130_LSC_SNRNASEQ_COUNTS.h5ad",
    ),
    (
        "GSE330130_LSC_ATAC",
        "GSE330130",
        "ATAC",
        SOURCE_DATA_ROOT / "GSE330130/GSE330130_LSC_SNATACSEQ_COUNTS.h5ad",
    ),
    (
        "GSE290359_RNA",
        "GSE290359",
        "RNA",
        SOURCE_DATA_ROOT / "GSE290359/GSE290359_BA44_BA46_COUNTS.h5ad",
    ),
    (
        "GSE290359_ATAC",
        "GSE290359",
        "ATAC",
        SOURCE_DATA_ROOT / "GSE290359/GSE290359_BA44_BA46_ATAC_COUNTS.h5ad",
    ),
]


def source_data_id(path: Path) -> str:
    try:
        return path.relative_to(SOURCE_DATA_ROOT).as_posix()
    except ValueError:
        return path.name


def project_rel(path: Path) -> str:
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.name


def first_existing(obs: pd.DataFrame, names: list[str]) -> pd.Series:
    for name in names:
        if name in obs.columns:
            return obs[name]
    return pd.Series([pd.NA] * len(obs), index=obs.index)


def disease_series(obs: pd.DataFrame) -> pd.Series:
    if "Diagnosis" in obs.columns:
        return obs["Diagnosis"].astype("string")
    if "Level 1" in obs.columns:
        return obs["Level 1"].astype("string")
    if "Disease?" in obs.columns and "Classifcation" in obs.columns:
        return obs["Classifcation"].astype("string")
    if "disease" in obs.columns:
        return obs["disease"].astype("string")
    return pd.Series([pd.NA] * len(obs), index=obs.index, dtype="string")


def canonical_obs(dataset_id: str, study: str, modality: str, path: Path) -> pd.DataFrame:
    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs.copy()

    can = pd.DataFrame(index=obs.index)
    can["cell_index"] = np.arange(adata.n_obs)
    can["dataset_id"] = dataset_id
    can["study"] = study
    can["modality"] = modality
    can["source_path"] = source_data_id(path)
    can["library"] = first_existing(obs, ["library_name", "library_id", "library"]).astype("string")
    can["sample"] = first_existing(obs, ["Sample ID", "sample"]).astype("string")
    if study == "GSE290359":
        parsed_sample = can["library"].astype("string").str.replace(r"_(BA44|BA46)$", "", regex=True)
        can["sample"] = can["sample"].fillna(parsed_sample)
    can["donor"] = first_existing(obs, ["Subject ID", "Edinburgh ID", "donor_id", "Donor"]).astype("string")

    region = first_existing(obs, ["Tissue Type", "Tissue", "tissue", "Region"]).astype("string")
    can["region_raw"] = region
    can["analysis_region"] = region.replace(
        {
            "Cortex_Motor_BA4": "MTC_all",
            "Cortex_Motor_Lateral": "MTC_all",
            "Cortex_Motor_Medial": "MTC_all",
            "Lumbar Spinal Cord": "LSC",
        }
    )
    can["tissue_subregion"] = can["region_raw"]

    disease = disease_series(obs)
    can["disease_raw"] = disease
    can["disease_binary"] = disease.replace(
        {
            "Non-neurological control": "Control",
            "Control": "Control",
            "control": "Control",
            "ALS": "ALS",
            "ALSci": "ALS",
            "ALS - borderline": "ALS_borderline",
            "ALS-FTD": "ALS_FTD",
            "ALS Spectrum MND, Other Neurological Disorders": "ALS_spectrum_other_neuro",
        }
    )
    can["clinical_subgroup"] = first_existing(obs, ["Level 2", "Classifcation"]).astype("string")
    can["mutation"] = first_existing(obs, ["Mutation", "C9_pos"]).astype("string")
    can["sex"] = first_existing(obs, ["Sex", "sex"]).astype("string")
    can["age"] = pd.to_numeric(first_existing(obs, ["Age at Death", "Age", "age"]), errors="coerce")
    can["PMI"] = pd.to_numeric(first_existing(obs, ["PMI"]), errors="coerce")

    cell_type = first_existing(obs, ["Cell_Type", "cell_type", "PrimaryAnnotation", "CellClass"]).astype("string")
    can["cell_type_raw"] = cell_type
    can["cell_type"] = cell_type.replace(
        {
            "Oligodendrocytes": "Oligodendrocyte",
            "OPCs": "OPC",
            "Oligodendrocyte Progenitor Cell": "OPC",
            "Oligodendrocyte Producing Cell": "OPC",
            "Astrocytes": "Astrocyte",
            "Inhibitory Neuron": "Inhibitory neuron",
            "Inhibitory neurons": "Inhibitory neuron",
            "Excitatory  Neuron": "Excitatory neuron",
            "Motor Neurons": "Motor neuron",
            "Motor neurons": "Motor neuron",
        }
    )
    can["cell_subtype"] = first_existing(obs, ["Cell_Subtype", "ThirdAnnotation", "SubType"]).astype("string")

    for col in [
        "Language impairment",
        "Verbal fluency impairment",
        "Executive impairment",
        "ECAS_Total Score (CUTOFF=105)",
        "Language (CUTOFF=26)",
        "Verbal_Fluency (CUTOFF=14)",
        "Executive (CUTOFF=33)",
    ]:
        can[col] = obs[col] if col in obs.columns else pd.NA
    can["good_cell"] = obs["good_cell"] if "good_cell" in obs.columns else pd.NA

    adata.file.close()
    return can.reset_index(drop=True)


def add_gse290359_atac_clinical_mapping(all_obs: pd.DataFrame) -> pd.DataFrame:
    rna = all_obs[all_obs["dataset_id"].eq("GSE290359_RNA")].drop_duplicates("library")
    clinical_cols = [
        "library",
        "donor",
        "sample",
        "disease_raw",
        "disease_binary",
        "clinical_subgroup",
        "mutation",
        "sex",
        "age",
        "PMI",
        "Language impairment",
        "Verbal fluency impairment",
        "Executive impairment",
        "ECAS_Total Score (CUTOFF=105)",
        "Language (CUTOFF=26)",
        "Verbal_Fluency (CUTOFF=14)",
        "Executive (CUTOFF=33)",
    ]
    libmap = rna[clinical_cols]
    mask = all_obs["dataset_id"].eq("GSE290359_ATAC")
    if not mask.any():
        return all_obs
    atac = all_obs.loc[mask].drop(columns=[c for c in clinical_cols if c != "library"])
    atac = atac.merge(libmap, on="library", how="left")
    return pd.concat([all_obs.loc[~mask], atac[all_obs.columns]], ignore_index=True)


def eligible_cells(obs: pd.DataFrame, policy: pd.DataFrame, dataset_id: str) -> pd.DataFrame:
    main = policy[
        policy["dataset_id"].eq(dataset_id)
        & policy["final_use_tier"].eq("main_model_candidate")
    ][["analysis_region", "cell_type"]].drop_duplicates()
    filtered = obs.merge(main, on=["analysis_region", "cell_type"], how="inner")

    if dataset_id == "GSE290359_RNA":
        filtered = filtered[~filtered["cell_type_raw"].isin(["Low-Quality", "Doublet"])]
        filtered = filtered[filtered["good_cell"].astype("string").fillna("True").eq("True")]
    if dataset_id == "GSE290359_ATAC":
        filtered = filtered[filtered["good_cell"].astype("string").fillna("True").eq("True")]
        filtered = filtered[filtered["donor"].notna() & filtered["disease_binary"].notna()]

    return filtered.reset_index(drop=True)


def sanitize_dataframe_for_h5ad(df: pd.DataFrame) -> pd.DataFrame:
    """Convert mixed object/string columns into h5ad-safe string columns."""
    safe = df.copy()
    safe.index = safe.index.astype(str)
    safe.columns = safe.columns.astype(str)
    for col in safe.columns:
        if pd.api.types.is_numeric_dtype(safe[col]):
            continue
        if pd.api.types.is_bool_dtype(safe[col]):
            safe[col] = safe[col].astype(str)
            continue
        safe[col] = safe[col].astype("string").fillna("").astype(str)
    return safe


def aggregate_dataset(
    dataset_id: str,
    study: str,
    modality: str,
    path: Path,
    obs_all: pd.DataFrame,
    policy: pd.DataFrame,
) -> tuple[ad.AnnData, pd.DataFrame]:
    adata = ad.read_h5ad(path, backed="r")
    obs = eligible_cells(obs_all[obs_all["dataset_id"].eq(dataset_id)].copy(), policy, dataset_id)

    group_cols = [
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
    ]
    clinical_cols = [
        "Language impairment",
        "Verbal fluency impairment",
        "Executive impairment",
        "ECAS_Total Score (CUTOFF=105)",
        "Language (CUTOFF=26)",
        "Verbal_Fluency (CUTOFF=14)",
        "Executive (CUTOFF=33)",
    ]
    group_cols = [c for c in group_cols if c in obs.columns]

    matrices: list[sp.csr_matrix] = []
    group_obs_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []

    for group_id, (keys, sub) in enumerate(obs.groupby(group_cols, dropna=False), start=1):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys, strict=False))
        indices = np.sort(sub["cell_index"].to_numpy(dtype=np.int64))
        x_sum = adata.X[indices, :].sum(axis=0)
        x_sum = sp.csr_matrix(np.asarray(x_sum))
        matrices.append(x_sum)

        row["pseudobulk_id"] = (
            f"{dataset_id}|{row.get('analysis_region')}|{row.get('disease_binary')}|"
            f"{row.get('donor')}|{row.get('library')}|{row.get('cell_type')}"
        )
        row["n_cells"] = int(len(indices))
        row["cell_index_min"] = int(indices.min()) if len(indices) else -1
        row["cell_index_max"] = int(indices.max()) if len(indices) else -1
        row["formal_analysis_started"] = "no"
        row["count_transform"] = "raw_count_sum"
        row["biological_replicate_unit"] = "donor"
        row["aggregation_unit"] = "donor/library/analysis_region/cell_type"
        row["primary_model_als_control"] = row.get("disease_binary") in {"ALS", "Control"}
        for c in clinical_cols:
            if c in sub.columns:
                values = sorted(sub[c].dropna().astype(str).unique())
                row[c] = ";".join(values[:10])
        group_obs_rows.append(row)
        manifest_rows.append({**row, "source_h5ad": source_data_id(path), "n_vars": int(adata.n_vars)})

    if not matrices:
        raise RuntimeError(f"No eligible cells for {dataset_id}")

    x = sp.vstack(matrices, format="csr")
    pb_obs = sanitize_dataframe_for_h5ad(pd.DataFrame(group_obs_rows))
    pb_obs.index = pb_obs["pseudobulk_id"].astype(str)
    pb_var = sanitize_dataframe_for_h5ad(adata.var.copy())
    pb = ad.AnnData(X=x, obs=pb_obs, var=pb_var)
    pb.uns["p01_note"] = (
        "Pseudobulk-ready raw count sums for frozen main-candidate strata. "
        "No DE/DA/motif analysis has been run."
    )
    pb.uns["source_h5ad"] = source_data_id(path)
    pb.uns["dataset_id"] = dataset_id
    pb.uns["modality"] = modality
    adata.file.close()
    return pb, pd.DataFrame(manifest_rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    policy = pd.read_csv(META / "celltype_analysis_policy.tsv", sep="\t")

    all_obs = pd.concat([canonical_obs(*obj) for obj in OBJECTS], ignore_index=True)
    all_obs = add_gse290359_atac_clinical_mapping(all_obs)

    all_manifests = []
    summary_rows = []
    for dataset_id, study, modality, path in OBJECTS:
        pb, manifest = aggregate_dataset(dataset_id, study, modality, path, all_obs, policy)
        output_path = OUT / f"{dataset_id}.main_candidates.pseudobulk.h5ad"
        manifest_path = OUT / f"{dataset_id}.main_candidates.pseudobulk_manifest.tsv"
        pb.write_h5ad(output_path, compression="gzip")
        manifest.to_csv(manifest_path, sep="\t", index=False)
        all_manifests.append(manifest)
        summary_rows.append(
            {
                "dataset_id": dataset_id,
                "modality": modality,
                "output_h5ad": project_rel(output_path),
                "manifest": project_rel(manifest_path),
                "n_pseudobulk_groups": int(pb.n_obs),
                "n_features": int(pb.n_vars),
                "n_cells_aggregated": int(manifest["n_cells"].sum()),
                "formal_analysis_started": "no",
            }
        )
        print(
            f"{dataset_id}: groups={pb.n_obs} features={pb.n_vars} "
            f"cells={int(manifest['n_cells'].sum())}"
        )

    all_manifest = pd.concat(all_manifests, ignore_index=True)
    all_manifest.to_csv(OUT / "all_main_candidates.pseudobulk_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(OUT / "pseudobulk_build_summary.tsv", sep="\t", index=False)

    print(f"wrote summary: {OUT / 'pseudobulk_build_summary.tsv'}")
    print(f"wrote combined manifest: {OUT / 'all_main_candidates.pseudobulk_manifest.tsv'}")


if __name__ == "__main__":
    main()
