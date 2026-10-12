#!/usr/bin/env python3
"""Inventory public GSE330130 libraries and flag the Takeuchi subset gap."""

import anndata as ad
import pandas as pd

from p01_paths import project_root, source_data_root

GSE = source_data_root() / "GSE330130"
OUT = project_root() / "results" / "W32_takeuchi_gse330130_library_reconciliation"

FILES = [
    ("MTC_RNA", GSE / "GSE330130_MTC_SNRNASEQ_COUNTS.h5ad", "library_id", "Diagnosis", "MTC"),
    ("MTC_ATAC", GSE / "GSE330130_MTC_SNATACSEQ_COUNTS.h5ad", "library", "Diagnosis", "MTC"),
    ("LSC_RNA", GSE / "GSE330130_LSC_SNRNASEQ_COUNTS.h5ad", "library_id", "Level 1", "LSC"),
    ("LSC_ATAC", GSE / "GSE330130_LSC_SNATACSEQ_COUNTS.h5ad", "library", "Level 1", "LSC"),
]


def inventory() -> pd.DataFrame:
    rows = []
    for assay, path, library_column, diagnosis_column, region in FILES:
        obj = ad.read_h5ad(path, backed="r")
        obs = obj.obs
        age_column = "Age at Death" if "Age at Death" in obs.columns else "Age"
        frame = pd.DataFrame(
            {
                "assay": assay,
                "region": region,
                "library_name": obs["library_name"].astype(str).to_numpy(),
                "subject_id": obs["Subject ID"].astype(str).to_numpy(),
                "tissue_type": obs["Tissue Type"].astype(str).to_numpy(),
                "diagnosis_public": obs[diagnosis_column].astype(str).to_numpy(),
                "pmi": obs["PMI"].astype(str).to_numpy(),
                "sex": obs["Sex"].astype(str).to_numpy(),
                "age": obs[age_column].astype(str).to_numpy(),
            }
        ).drop_duplicates()
        rows.append(frame)
        obj.file.close()

    merged = pd.concat(rows, ignore_index=True)
    merged = (
        merged.pivot_table(
            index=[
                "region",
                "library_name",
                "subject_id",
                "tissue_type",
                "diagnosis_public",
                "pmi",
                "sex",
                "age",
            ],
            columns="assay",
            values="region",
            aggfunc="size",
            fill_value=0,
        )
        .reset_index()
    )
    for column in ("MTC_RNA", "MTC_ATAC", "LSC_RNA", "LSC_ATAC"):
        if column not in merged.columns:
            merged[column] = 0
    merged["rna_available"] = (
        (merged["MTC_RNA"] > 0) | (merged["LSC_RNA"] > 0)
    ).astype(int)
    merged["atac_available"] = (
        (merged["MTC_ATAC"] > 0) | (merged["LSC_ATAC"] > 0)
    ).astype(int)
    merged["paper_inclusion_status"] = "unresolved_public_library_vs_paper_subset"
    merged["paper_sample_count_note"] = (
        "Takeuchi main text/Table 1 reports 6 ALS + 6 control per tissue; "
        "the exact public-library subset is not recovered from the current full text."
    )
    return merged


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    table = inventory()
    table.to_csv(OUT / "w32_gse330130_unique_public_library_inventory.tsv", sep="\t", index=False)
    summary = table.assign(
        disease_binary=table["diagnosis_public"].str.contains("ALS", case=False, na=False).map(
            {True: "ALS", False: "Control"}
        )
    ).groupby(["region", "disease_binary"]).agg(
        libraries=("library_name", "nunique"), subjects=("subject_id", "nunique")
    ).reset_index()
    summary.to_csv(OUT / "w32_public_library_summary.tsv", sep="\t", index=False)
    (OUT / "README.md").write_text(
        "# W32 Takeuchi/GSE330130 library reconciliation\n\n"
        "Takeuchi main text/Table 1 reports 6 ALS and 6 controls per tissue, "
        "whereas the public GSE330130 h5ad objects contain additional libraries. "
        "The exact paper-included subset remains unresolved until the supplementary "
        "sample table or original inclusion manifest is recovered. Do not describe "
        "all public libraries as independent samples from the published analysis "
        "without resolving this gate.\n"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
