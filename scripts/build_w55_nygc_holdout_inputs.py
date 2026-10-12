#!/usr/bin/env python3
"""Build W55 NYGC bulk hold-out inputs (plan Section 4, Amendment 1 A1-A3).

Writes the discovery-donor exclusion table, the hold-out sample sheet and the
count matrix for retained motor-cortex and spinal-cord samples. No module or
clinical outcome is computed here.
"""

from __future__ import annotations

import pandas as pd

from p01_paths import project_root, source_data_root

NYGC = source_data_root() / "Other_Datasets/NYGC_ALS_Consortium/files"
W32 = project_root() / "results/W32_takeuchi_gse330130_library_reconciliation/w32_gse330130_unique_public_library_inventory.tsv"
OUT = project_root() / "results/W55_NYGC_bulk_holdout"

MOTOR = {"Unspecified_Motor_Cortex", "Lateral_Motor_Cortex", "Medial_Motor_Cortex"}
SPINAL = {"Cervical_Spinal_Cord", "Lumbar_Spinal_Cord", "Thoracic_Spinal_Cord", "Unspecified_Spinal_Cord"}
GROUPS = {"ALS", "Control", "ALS-FTD"}


def numeric_age(value: object) -> float:
    text = str(value).strip().lower()
    if text.startswith("90 or"):
        return 90.0
    return pd.to_numeric(text, errors="coerce")


def diagnosis_class(value: object) -> str:
    text = str(value)
    if "control" in text.lower():
        return "Control"
    if "ALS" in text:
        return "ALS"
    return "Other"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    workbook = NYGC / "NYGC_Flagship_ST_1_3_Metadata.xlsx"
    st1 = pd.read_excel(workbook, "ST1", header=1, dtype=str)
    st3 = pd.read_excel(workbook, "ST3", header=1, dtype=str)
    disc = pd.read_csv(W32, sep="\t", dtype=str).drop_duplicates("subject_id")

    # A1: exact donor-ID exclusion, then a conservative fingerprint for unmatched subjects.
    nygc_ids = set(st3["donor_id"]) | set(st1["donor_id"])
    rows = [{"nygc_donor_id": s, "gse330130_subject": s, "rule": "exact_donor_id"}
            for s in disc["subject_id"] if s in nygc_ids]
    st1["age_num"] = st1["age_at_death"].map(numeric_age)
    st1["dx_class"] = st1["disease_group"].map(diagnosis_class)
    for _, d in disc[~disc["subject_id"].isin(nygc_ids)].iterrows():
        age = float(d["age"])
        dx = diagnosis_class(d["diagnosis_public"])
        age_ok = (st1["age_num"] - age).abs().le(1) | (st1["age_num"].eq(90) & (age >= 89))
        hit = st1[age_ok & st1["dx_class"].eq(dx) & st1["sex_call"].eq(d["sex"])]
        rows += [{"nygc_donor_id": h, "gse330130_subject": d["subject_id"], "rule": "fingerprint_dx_sex_age"}
                 for h in hit["donor_id"]]
    exclusions = pd.DataFrame(rows).drop_duplicates()
    exclusions.to_csv(OUT / "w55_discovery_overlap_exclusions.tsv", sep="\t", index=False)
    excluded = set(exclusions["nygc_donor_id"])

    # A2/A3: samples, groups, regions and covariates.
    s = st3[st3["used_in_analysis"].eq("TRUE") & st3["disease_group"].isin(GROUPS)].copy()
    s["region"] = pd.Series(pd.NA, index=s.index, dtype="object")
    s.loc[s["tissue_clean"].isin(SPINAL), "region"] = "spinal_cord"
    s.loc[s["tissue_clean"].isin(MOTOR) | s["tissue_full"].eq("Cortex_Motor_BA4"), "region"] = "motor_cortex"
    s = s[s["region"].notna()]
    n_before = s.groupby("region")["donor_id"].nunique()
    s = s[~s["donor_id"].isin(excluded)]
    donor_cols = ["donor_id", "sex_call", "age_at_death", "age_at_symptom_onset", "site_of_motor_onset",
                  "disease_duration_in_months", "c9orf72_rp_pcr", "c9orf72_max_call"]
    s = s.merge(st1[donor_cols].drop_duplicates("donor_id"), on="donor_id", how="left")
    s["age"] = s["age_at_death"].map(numeric_age)
    s["pmi"] = pd.to_numeric(s["pmi"], errors="coerce")
    s["rin"] = pd.to_numeric(s["rin"], errors="coerce")
    keep_cols = ["rna_id", "donor_id", "region", "tissue_clean", "disease_group", "sex_call", "age", "pmi", "rin",
                 "rna_platform", "age_at_symptom_onset", "site_of_motor_onset", "disease_duration_in_months",
                 "c9orf72_rp_pcr", "c9orf72_max_call"]
    s[keep_cols].to_csv(OUT / "w55_holdout_samples.tsv", sep="\t", index=False)

    counts = pd.read_csv(NYGC / "count_data_NYGC_ALS_Consortium.csv.gz", index_col=0, usecols=lambda c: c == "Unnamed: 0" or c in set(s["rna_id"]))
    missing = sorted(set(s["rna_id"]) - set(counts.columns))
    counts.index = counts.index.str.replace(r"\.\d+$", "", regex=True)
    counts = counts.groupby(level=0).sum()
    counts.index.name = "gene_id"
    counts.to_csv(OUT / "w55_holdout_counts.tsv.gz", sep="\t", compression="gzip")

    summary = s.groupby(["region", "disease_group"]).agg(samples=("rna_id", "size"), donors=("donor_id", "nunique"))
    summary.to_csv(OUT / "w55_holdout_sample_summary.tsv", sep="\t")
    print("excluded NYGC donors:", len(excluded), exclusions["rule"].value_counts().to_dict())
    print("donors per region before exclusion:", n_before.to_dict())
    print(summary.to_string())
    print("samples without counts:", len(missing))


if __name__ == "__main__":
    main()
