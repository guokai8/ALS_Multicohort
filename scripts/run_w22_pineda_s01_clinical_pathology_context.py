#!/usr/bin/env python3
"""Run W22: bounded Pineda S01 donor-level clinical/pathology context audit.

W22 tests whether the already selected MTC astrocyte ALS-up module has a
relationship with clinical or neuropathology-context variables *within the
independent Pineda cohort*.  It does not add metadata to GSE330130 and does
not establish disease progression, survival, quantitative pathology burden,
or a causal mechanism.

The workflow uses one BA4/BA9 aggregate score per donor.  S01 is the Synapse
metadata source of record.  The tests and their correction family are fixed
in ANALYSIS_REGISTRY to prevent post-hoc endpoint selection.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from p01_paths import project_root, source_data_root


ROOT = project_root()
SOURCE_DATA_ROOT = source_data_root()
W08_SCORES = ROOT / "results/W08_Pineda_validation/w08_donor_mean_module_scores.tsv"
S01 = SOURCE_DATA_ROOT / "Other_Datasets/Pineda_snRNAseq/files/pineda/S01_sample_information.xlsx"
OUT = ROOT / "results/W22_Pineda_S01_clinical_pathology_context"

LEAD_MODULE = "GSE330130_MTC_RNA__Astrocyte__up"

# These seven analyses are one bounded, exploratory correction family.  They
# address Pineda-only clinical/pathology context and are not P01 endpoints.
ANALYSIS_REGISTRY = [
    {
        "analysis_id": "duration_pooled_als_ftld",
        "endpoint": "Disease duration (Y)",
        "population": "ALS_or_FTLD",
        "covariates": "disease_family;C9_pos;Sex",
        "interpretation_boundary": "external disease-stage context only; ALS and FTLD duration are not interchangeable severity measures",
    },
    {
        "analysis_id": "duration_als",
        "endpoint": "Disease duration (Y)",
        "population": "ALS",
        "covariates": "C9_pos;Sex",
        "interpretation_boundary": "ALS cortical context only; not GSE330130 progression",
    },
    {
        "analysis_id": "duration_ftld",
        "endpoint": "Disease duration (Y)",
        "population": "FTLD",
        "covariates": "C9_pos;Sex",
        "interpretation_boundary": "FTLD cortical context only; not an ALS endpoint",
    },
    {
        "analysis_id": "motor_onset_age_als",
        "endpoint": "Age of onset motor (Y)",
        "population": "ALS",
        "covariates": "C9_pos;Sex",
        "interpretation_boundary": "age-at-motor-onset context only; not disease progression",
    },
    {
        "analysis_id": "cognitive_onset_age_ftld",
        "endpoint": "Age of onset cognitive (Y)",
        "population": "FTLD",
        "covariates": "C9_pos;Sex",
        "interpretation_boundary": "FTLD cognitive-onset context only; not ALS cognition",
    },
    {
        "analysis_id": "braak_context_all_groups",
        "endpoint": "Braak_stage",
        "population": "ALL",
        "covariates": "disease_family;C9_pos;Sex",
        "interpretation_boundary": "co-pathology context only; Braak is not quantitative ALS pTDP-43 pathology",
    },
    {
        "analysis_id": "thal_context_all_groups",
        "endpoint": "Thal_phase",
        "population": "ALL",
        "covariates": "disease_family;C9_pos;Sex",
        "interpretation_boundary": "co-pathology context only; Thal is not quantitative ALS pTDP-43 pathology",
    },
]


def clean_donor(value: object) -> str | None:
    """Normalize integer-like donor identifiers without numeric-only cross-cohort matching."""
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in {"nan", "na", "nr", "none"}:
        return None
    try:
        numeric = float(text)
    except ValueError:
        return text
    if math.isfinite(numeric) and numeric.is_integer():
        return str(int(numeric))
    return text


def numeric_column(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str).str.strip().replace({"NR": np.nan, "NA": np.nan, "": np.nan}), errors="coerce")


def prepare_s01() -> pd.DataFrame:
    if not S01.exists():
        raise FileNotFoundError(f"Missing Pineda Synapse S01 workbook: {S01}")
    raw = pd.read_excel(S01, sheet_name="donor_information")
    raw["Donor"] = raw["Donor"].map(clean_donor)
    raw = raw[raw["Donor"].notna()].copy()
    valid_groups = {"ALS", "C9ALS", "FTLD", "C9FTLD", "PN"}
    raw = raw[raw["Group"].astype(str).isin(valid_groups)].copy()
    if raw["Donor"].duplicated().any():
        raise RuntimeError("S01 has duplicate valid donor IDs after cleaning")
    raw["disease_family"] = np.select(
        [raw["Group"].isin(["ALS", "C9ALS"]), raw["Group"].isin(["FTLD", "C9FTLD"]), raw["Group"].eq("PN")],
        ["ALS", "FTLD", "PN"],
        default="unknown",
    )
    raw["C9_pos"] = raw["Group"].isin(["C9ALS", "C9FTLD"]).astype(int)
    for col in [
        "PMI_h",
        "Braak_stage",
        "Thal_phase",
        "Age of onset cognitive (Y)",
        "Age of onset motor (Y)",
        "Age of death (Y)",
        "Disease duration (Y)",
    ]:
        raw[col] = numeric_column(raw[col])
    return raw.sort_values("Donor").reset_index(drop=True)


def prepare_scores() -> pd.DataFrame:
    if not W08_SCORES.exists():
        raise FileNotFoundError(f"Missing W08 donor scores: {W08_SCORES}")
    scores = pd.read_csv(W08_SCORES, sep="\t")
    scores["Donor"] = scores["Donor"].map(clean_donor)
    scores = scores[
        scores["module_id"].eq(LEAD_MODULE) & scores["target_cell_type"].eq("Astrocyte") & scores["region_scope"].eq("BA4_BA9_mean")
    ].copy()
    if scores.empty:
        raise RuntimeError(f"No BA4/BA9 astrocyte scores found for {LEAD_MODULE}")
    if scores["Donor"].duplicated().any():
        raise RuntimeError("W08 lead-module input has more than one score per donor")
    return scores.sort_values("Donor").reset_index(drop=True)


def choose_population(df: pd.DataFrame, population: str) -> pd.DataFrame:
    if population == "ALL":
        return df.copy()
    if population == "ALS_or_FTLD":
        return df[df["disease_family"].isin(["ALS", "FTLD"])].copy()
    return df[df["disease_family"].eq(population)].copy()


def categorical_design(df: pd.DataFrame, column: str) -> tuple[pd.DataFrame, list[str]]:
    if column not in df.columns or df[column].isna().any() or df[column].astype(str).nunique() < 2:
        return pd.DataFrame(index=df.index), []
    dummies = pd.get_dummies(df[column].astype(str), prefix=column, drop_first=True, dtype=float)
    return dummies, dummies.columns.tolist()


def fit_registered_model(merged: pd.DataFrame, spec: dict[str, str]) -> dict[str, object]:
    endpoint = spec["endpoint"]
    model_df = choose_population(merged, spec["population"])
    result: dict[str, object] = {
        **spec,
        "module_id": LEAD_MODULE,
        "target_cell_type": "Astrocyte",
        "status": "not_run",
        "n_donors_candidate": int(model_df["Donor"].nunique()),
        "n_donors_model": 0,
        "n_als": 0,
        "n_ftld": 0,
        "n_pn": 0,
        "model_formula": "",
        "covariates_used": "",
        "covariates_dropped": "",
        "beta_per_year_or_stage": np.nan,
        "se": np.nan,
        "ci95_low": np.nan,
        "ci95_high": np.nan,
        "p_value": np.nan,
        "adj_p_global_BH": np.nan,
        "message": "",
    }
    if endpoint not in model_df.columns:
        result.update(status="skipped", message="endpoint_column_missing")
        return result
    model_df = model_df.dropna(subset=["module_score", endpoint]).copy()
    result["n_donors_model"] = int(model_df["Donor"].nunique())
    result["n_als"] = int(model_df.loc[model_df["disease_family"].eq("ALS"), "Donor"].nunique())
    result["n_ftld"] = int(model_df.loc[model_df["disease_family"].eq("FTLD"), "Donor"].nunique())
    result["n_pn"] = int(model_df.loc[model_df["disease_family"].eq("PN"), "Donor"].nunique())
    if result["n_donors_model"] < 12:
        result.update(status="skipped", message="fewer_than_12_donors_with_score_and_endpoint")
        return result

    x_parts = [pd.DataFrame({"endpoint_value": model_df[endpoint].astype(float)}, index=model_df.index)]
    used = []
    dropped = []
    requested = spec["covariates"].split(";") if spec["covariates"] else []
    for covariate in requested:
        if covariate == "C9_pos":
            if model_df[covariate].nunique(dropna=True) >= 2:
                x_parts.append(pd.DataFrame({"C9_pos": model_df[covariate].astype(float)}, index=model_df.index))
                used.append(covariate)
            else:
                dropped.append(f"{covariate}:not_estimable")
        elif covariate == "Sex":
            design, columns = categorical_design(model_df, covariate)
            if columns:
                x_parts.append(design)
                used.append(covariate)
            else:
                dropped.append(f"{covariate}:not_estimable")
        elif covariate == "disease_family":
            design, columns = categorical_design(model_df, covariate)
            if columns:
                x_parts.append(design)
                used.append(covariate)
            else:
                dropped.append(f"{covariate}:not_estimable")
        else:
            dropped.append(f"{covariate}:unsupported")

    x = pd.concat(x_parts, axis=1).astype(float)
    x = sm.add_constant(x, has_constant="add")
    # Retain the endpoint while dropping only non-essential rank-deficient covariates.
    while np.linalg.matrix_rank(x.to_numpy()) < x.shape[1]:
        removable = [col for col in x.columns[::-1] if col not in {"const", "endpoint_value"}]
        if not removable:
            result.update(status="skipped", message="endpoint_rank_deficient")
            return result
        removed = removable[0]
        x = x.drop(columns=removed)
        dropped.append(f"{removed}:rank_deficient")

    residual_df = x.shape[0] - x.shape[1]
    if residual_df < 8:
        result.update(status="skipped", message="fewer_than_8_residual_degrees_of_freedom")
        return result
    try:
        fit = sm.OLS(model_df["module_score"].astype(float), x).fit()
    except Exception as exc:  # noqa: BLE001
        result.update(status="error", message=str(exc))
        return result
    ci = fit.conf_int().loc["endpoint_value"]
    result.update(
        status="completed",
        model_formula="module_score ~ endpoint_value + " + " + ".join(used or ["1"]),
        covariates_used=";".join(used),
        covariates_dropped="|".join(dropped),
        beta_per_year_or_stage=float(fit.params["endpoint_value"]),
        se=float(fit.bse["endpoint_value"]),
        ci95_low=float(ci.iloc[0]),
        ci95_high=float(ci.iloc[1]),
        p_value=float(fit.pvalues["endpoint_value"]),
        message="completed_OLS_one_BA4_BA9_score_per_donor",
    )
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    s01 = prepare_s01()
    scores = prepare_scores()
    merged = scores.merge(s01, on="Donor", how="left", suffixes=("_w08", "_synapse"), validate="one_to_one", indicator=True)

    linkage = merged[["Donor", "module_score", "Group_w08", "Sex_w08", "Group_synapse", "Sex_synapse", "_merge"]].copy()
    linkage = linkage.rename(columns={"Group_synapse": "Group", "Sex_synapse": "Sex"})
    linkage["group_conflict_note"] = np.where(
        linkage["_merge"].eq("both") & linkage["Group_w08"].astype(str).ne(linkage["Group"].astype(str)),
        "Synapse_S01_retained_as_source_of_record",
        "none",
    )
    linkage["sex_conflict_note"] = np.where(
        linkage["_merge"].eq("both") & linkage["Sex_w08"].astype(str).ne(linkage["Sex"].astype(str)),
        "Synapse_S01_retained_as_source_of_record",
        "none",
    )
    linkage.to_csv(OUT / "w22_synapse_w08_linkage_audit.tsv", sep="\t", index=False)
    if not merged["_merge"].eq("both").all():
        raise RuntimeError("Not all W08 lead-module donors map to a valid Synapse S01 donor record")

    # W08 fields are retained only for the source-record conflict audit above.
    # Downstream clinical models use the Synapse S01 fields by their canonical
    # names, consistent with the project metadata conflict policy.
    merged = merged.rename(
        columns={
            "Group_synapse": "Group",
            "Sex_synapse": "Sex",
            "C9_pos_synapse": "C9_pos",
        }
    ).drop(columns=["_merge"])
    endpoint_cols = [item["endpoint"] for item in ANALYSIS_REGISTRY]
    coverage = []
    for endpoint in sorted(set(endpoint_cols)):
        coverage.append(
            {
                "endpoint": endpoint,
                "n_donors_with_score": int(merged["Donor"].nunique()),
                "n_non_missing": int(merged[endpoint].notna().sum()),
                "n_non_missing_als": int(merged.loc[merged["disease_family"].eq("ALS"), endpoint].notna().sum()),
                "n_non_missing_ftld": int(merged.loc[merged["disease_family"].eq("FTLD"), endpoint].notna().sum()),
                "n_non_missing_pn": int(merged.loc[merged["disease_family"].eq("PN"), endpoint].notna().sum()),
                "metadata_provenance": "Pineda_S01_Synapse",
                "p01_claim_boundary": "Pineda external context only; not transferred to GSE330130",
            }
        )
    pd.DataFrame(coverage).to_csv(OUT / "w22_endpoint_coverage.tsv", sep="\t", index=False)
    merged.to_csv(OUT / "w22_pineda_lead_module_donor_context.tsv", sep="\t", index=False)

    registry = pd.DataFrame(ANALYSIS_REGISTRY)
    registry["module_id"] = LEAD_MODULE
    registry["correction_family"] = "all_completed_W22_tests"
    registry["analysis_unit"] = "one_BA4_BA9_donor_score"
    registry.to_csv(OUT / "w22_analysis_registry.tsv", sep="\t", index=False)

    results = pd.DataFrame([fit_registered_model(merged, spec) for spec in ANALYSIS_REGISTRY])
    completed = results["status"].eq("completed")
    if completed.any():
        results.loc[completed, "adj_p_global_BH"] = multipletests(
            results.loc[completed, "p_value"].astype(float), method="fdr_bh"
        )[1]
    results["result_tier"] = np.select(
        [
            results["status"].eq("completed") & results["adj_p_global_BH"].lt(0.05),
            results["status"].eq("completed") & results["p_value"].lt(0.05),
            results["status"].eq("completed"),
        ],
        ["exploratory_BH_supported_pineda_context", "nominal_only_not_claim_eligible", "no_association_supported"],
        default="not_tested",
    )
    results["manuscript_eligibility"] = np.where(
        results["result_tier"].eq("exploratory_BH_supported_pineda_context"),
        "supplementary_external_context_only",
        "not_a_main_claim",
    )
    results.to_csv(OUT / "w22_model_results.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "metric": "lead_module",
                "value": LEAD_MODULE,
            },
            {
                "metric": "linked_donors",
                "value": int(merged["Donor"].nunique()),
            },
            {
                "metric": "synapse_linkage_conflicts",
                "value": int((linkage["group_conflict_note"].ne("none") | linkage["sex_conflict_note"].ne("none")).sum()),
            },
            {
                "metric": "completed_models",
                "value": int(completed.sum()),
            },
            {
                "metric": "BH05_supported_models",
                "value": int((results["result_tier"] == "exploratory_BH_supported_pineda_context").sum()),
            },
            {
                "metric": "main_claim_eligible_models",
                "value": 0,
            },
            {
                "metric": "interpretation_boundary",
                "value": "No W22 result supplies GSE330130 clinical/pathology endpoint, survival, progression rate, quantitative pTDP-43 burden, or causal mechanism.",
            },
        ]
    )
    summary.to_csv(OUT / "w22_run_summary.tsv", sep="\t", index=False)
    print(f"[W22] wrote {OUT}")
    print(results[["analysis_id", "n_donors_model", "beta_per_year_or_stage", "p_value", "adj_p_global_BH", "result_tier"]].to_string(index=False))


if __name__ == "__main__":
    main()
