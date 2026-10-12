#!/usr/bin/env python3
"""Build a Synapse/NYGC-prioritized GSE330130 clinical bridge.

This script uses the latest NYGC ALS Consortium Zenodo workbook, where ST3 is
populated with RNA sample metadata. The exact `CGND-HRA-*` RNA IDs are used as
the bridge to GSE330130 GEO sample IDs. For exact matches, resolved clinical
and diagnosis fields prioritize Synapse/NYGC curated metadata; GEO fields are
retained as raw sequencing provenance and conflict notes.

Inputs can be overridden by environment variables:
- GSE330130_SERIES_MATRIX
- NYGC_ZENODO_ST1_ST3_XLSX
"""

from __future__ import annotations

import hashlib
import os
import re
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from build_gse330130_nygc_clinical_linkage_reconciliation import (
    SERIES_DURATION_ALIAS,
    broad_region,
    parse_series_libraries,
)
from p01_paths import plan_root, project_root, source_data_root


DATE_TAG = "2026-08-12"
ROOT = project_root()
PLAN_ROOT = plan_root()
OUT = ROOT / "results" / "Supplemental_QC_covariate_audit"
RESCUE_ROOT = ROOT / "metadata_endpoint_rescue_2026-08-12"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_input(env_name: str, candidates: Iterable[Path]) -> Path:
    value = os.environ.get(env_name)
    if value:
        return Path(value).expanduser().resolve()
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError(
        f"Set {env_name}; checked " + "; ".join(str(x) for x in candidates)
    )


def clean_text(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    return "" if text.lower() in {"", "nan", "unknown", "not applicable"} else text


def first_clean(*values: object) -> str:
    for value in values:
        text = clean_text(value)
        if text:
            return text
    return ""


def numeric_value(value: object) -> float | None:
    text = clean_text(value)
    if not text:
        return None
    match = re.search(r"[-+]?\d+(?:\.\d+)?", text)
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def numeric_close(a: object, b: object, tol: float = 0.51) -> bool:
    av = numeric_value(a)
    bv = numeric_value(b)
    if av is None or bv is None:
        return True
    return abs(av - bv) <= tol


def norm_disease(value: object) -> str:
    text = clean_text(value).lower().replace(" ", "_")
    if text in {"non-neurological_control", "control"}:
        return "Control"
    if text == "als-ftd":
        return "ALS-FTD"
    if text == "als":
        return "ALS"
    return clean_text(value)


def norm_onset(value: object) -> str:
    text = clean_text(value).lower()
    if "bulbar" in text and "limb" in text:
        return "Bulbar and Limb"
    if "bulbar" in text:
        return "Bulbar"
    if "limb" in text or "spinal" in text:
        return "Limb"
    if "axial" in text:
        return "Axial"
    return clean_text(value)


def add_resolved_fields(
    row: dict[str, object],
    hard_conflicts: list[str],
    diagnosis_flags: list[str],
    warnings: list[str],
) -> None:
    """Add Synapse/NYGC-prioritized fields without overwriting raw inputs."""
    direct_match = row["nygc_latest_direct_rna_id_match"] == "yes"
    if direct_match:
        row["resolved_source_policy"] = "synapse_nygc_curated_priority"
        row["resolved_endpoint_source"] = (
            "NYGC_latest_Zenodo_ST3_exact_rna_id_match_plus_ST1_donor_metadata"
        )
        row["resolved_donor_id"] = first_clean(row.get("nygc_st1_donor_id"), row.get("nygc_st3_donor_id"))
        row["resolved_rna_id"] = first_clean(row.get("nygc_st3_rna_id"), row.get("series_cgnd_id"))
        row["resolved_tissue_clean"] = first_clean(row.get("nygc_st3_tissue_clean"), row.get("series_tissue"))
        row["resolved_disease_group"] = first_clean(
            row.get("nygc_st1_disease_group"),
            row.get("nygc_st3_disease_group"),
            norm_disease(row.get("series_genotype")),
        )
        row["resolved_als_ftd_detail"] = first_clean(
            row.get("nygc_st1_subject_group_subcategory"),
            row.get("nygc_st1_subject_group"),
            row.get("resolved_disease_group"),
        )
        row["resolved_site_of_motor_onset"] = first_clean(
            row.get("nygc_st1_site_of_motor_onset"),
            row.get("series_site_of_motor_onset"),
        )
        row["resolved_site_of_motor_onset_detail"] = first_clean(
            row.get("nygc_st1_site_of_motor_onset_detail"),
            row.get("resolved_site_of_motor_onset"),
        )
        row["resolved_age_at_symptom_onset"] = numeric_value(
            first_clean(row.get("nygc_st1_age_at_symptom_onset"), row.get("series_age_at_symptom_onset"))
        )
        row["resolved_age_at_death"] = numeric_value(
            first_clean(row.get("nygc_st1_age_at_death"), row.get("series_age_at_death"))
        )
        row["resolved_disease_duration_months"] = numeric_value(
            first_clean(row.get("nygc_st1_disease_duration_in_months"), row.get(SERIES_DURATION_ALIAS))
        )
        row["resolved_cause_of_death_group"] = first_clean(row.get("nygc_st1_cause_of_death_group"))
        row["resolved_cause_of_death_related_to_als_or_mnd"] = first_clean(
            row.get("nygc_st1_cause_of_death_related_to_als_or_mnd")
        )
        row["resolved_current_alsfrs_r"] = numeric_value(row.get("nygc_st1_current_alsfrs_r"))
        row["resolved_fvc_percent_predicted"] = numeric_value(row.get("nygc_st1_fvc_percent_predicted"))
    else:
        row["resolved_source_policy"] = "geo_series_matrix_fallback_no_synapse_match"
        row["resolved_endpoint_source"] = "GSE330130_GEO_series_matrix"
        row["resolved_donor_id"] = first_clean(row.get("mo_id"))
        row["resolved_rna_id"] = first_clean(row.get("series_cgnd_id"), row.get("series_u54_id"))
        row["resolved_tissue_clean"] = first_clean(row.get("series_tissue"))
        row["resolved_disease_group"] = norm_disease(row.get("series_genotype"))
        row["resolved_als_ftd_detail"] = first_clean(row.get("resolved_disease_group"))
        row["resolved_site_of_motor_onset"] = first_clean(row.get("series_site_of_motor_onset"))
        row["resolved_site_of_motor_onset_detail"] = first_clean(row.get("series_site_of_motor_onset"))
        row["resolved_age_at_symptom_onset"] = numeric_value(row.get("series_age_at_symptom_onset"))
        row["resolved_age_at_death"] = numeric_value(row.get("series_age_at_death"))
        row["resolved_disease_duration_months"] = numeric_value(row.get(SERIES_DURATION_ALIAS))
        row["resolved_cause_of_death_group"] = ""
        row["resolved_cause_of_death_related_to_als_or_mnd"] = ""
        row["resolved_current_alsfrs_r"] = None
        row["resolved_fvc_percent_predicted"] = None

    row["resolved_broad_region"] = broad_region(row.get("resolved_tissue_clean"))
    note_parts = [
        "Synapse/NYGC curated metadata is authoritative for resolved clinical/diagnosis fields."
        if direct_match
        else "No exact Synapse/NYGC ST3 RNA match; resolved fields fall back to GEO series metadata."
    ]
    if hard_conflicts:
        note_parts.append("GEO-vs-Synapse field differences noted: " + ";".join(hard_conflicts) + ".")
    if diagnosis_flags:
        note_parts.append("GEO-vs-Synapse diagnosis label differences noted: " + ";".join(diagnosis_flags) + ".")
    if warnings:
        note_parts.append("Additional label warnings: " + ";".join(warnings) + ".")
    row["resolved_conflict_note"] = " ".join(note_parts)


def read_nygc_latest(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    st1 = pd.read_excel(path, sheet_name="ST1", header=1)
    st3 = pd.read_excel(path, sheet_name="ST3", header=1)
    st1 = st1[st1["donor_id"].notna()].copy()
    st3 = st3[st3["rna_id"].notna()].copy()
    st1["_st1_source_md5"] = md5(path)
    st1["_st1_source_sha256"] = sha256(path)
    st3["_st3_source_md5"] = md5(path)
    st3["_st3_source_sha256"] = sha256(path)
    return st1, st3


def build_bridge(series: pd.DataFrame, st1: pd.DataFrame, st3: pd.DataFrame) -> pd.DataFrame:
    st3_keep = [
        "rna_id",
        "donor_id",
        "tissue_clean",
        "region",
        "disease_group",
        "rin",
        "pmi",
        "used_in_analysis",
        "sra_run",
        "_st3_source_md5",
        "_st3_source_sha256",
    ]
    st1_keep = [
        "donor_id",
        "disease_group",
        "subject_group",
        "subject_group_subcategory",
        "sex_call",
        "site_of_motor_onset",
        "site_of_motor_onset_detail",
        "age_at_symptom_onset",
        "age_at_death",
        "disease_duration_in_months",
        "age_at_diagnosis",
        "cause_of_death",
        "cause_of_death_group",
        "cause_of_death_related_to_als_or_mnd",
        "current_alsfrs_r",
        "fvc_percent_predicted",
        "family_history_of_als_ftd",
        "family_members_with_ftd",
        "c9orf72_genotype",
        "c9orf72_max_call",
        "atxn2_genotype",
        "_st1_source_md5",
        "_st1_source_sha256",
    ]
    st3_small = st3[st3_keep].copy().add_prefix("nygc_st3_")
    st1_small = st1[st1_keep].drop_duplicates("donor_id").copy().add_prefix("nygc_st1_")

    merged = series.merge(
        st3_small,
        left_on="series_cgnd_id",
        right_on="nygc_st3_rna_id",
        how="left",
    )
    merged = merged.merge(
        st1_small,
        left_on="nygc_st3_donor_id",
        right_on="nygc_st1_donor_id",
        how="left",
    )
    merged["nygc_latest_direct_rna_id_match"] = merged["nygc_st3_rna_id"].notna().map(
        {True: "yes", False: "no"}
    )

    records: list[dict[str, object]] = []
    for _, raw in merged.iterrows():
        row = raw.to_dict()
        hard_conflicts: list[str] = []
        diagnosis_flags: list[str] = []
        warnings: list[str] = []

        if row["nygc_latest_direct_rna_id_match"] == "yes":
            series_disease = norm_disease(row.get("series_genotype"))
            st3_disease = norm_disease(row.get("nygc_st3_disease_group"))
            st1_disease = norm_disease(row.get("nygc_st1_disease_group"))
            if series_disease and st3_disease and series_disease != st3_disease:
                diagnosis_flags.append("series_vs_st3_disease_group")
            if series_disease and st1_disease and series_disease != st1_disease:
                diagnosis_flags.append("series_vs_st1_disease_group")

            series_region = broad_region(row.get("series_tissue"))
            nygc_region = broad_region(row.get("nygc_st3_tissue_clean"))
            if series_region != "unknown" and nygc_region != "unknown" and series_region != nygc_region:
                hard_conflicts.append("broad_region")
            elif clean_text(row.get("series_tissue")) != clean_text(row.get("nygc_st3_tissue_clean")):
                warnings.append("fine_tissue_label")

            if clean_text(row.get("series_sex")) and clean_text(row.get("nygc_st1_sex_call")):
                if clean_text(row.get("series_sex")) != clean_text(row.get("nygc_st1_sex_call")):
                    hard_conflicts.append("sex")

            for series_col, nygc_col, label in [
                ("series_pmi", "nygc_st3_pmi", "pmi"),
                ("series_age_at_death", "nygc_st1_age_at_death", "age_at_death"),
                ("series_age_at_symptom_onset", "nygc_st1_age_at_symptom_onset", "age_at_symptom_onset"),
                (
                    SERIES_DURATION_ALIAS,
                    "nygc_st1_disease_duration_in_months",
                    "disease_duration_months",
                ),
            ]:
                if not numeric_close(row.get(series_col), row.get(nygc_col)):
                    hard_conflicts.append(label)

            series_onset = norm_onset(row.get("series_site_of_motor_onset"))
            st1_onset = norm_onset(row.get("nygc_st1_site_of_motor_onset"))
            if series_onset and st1_onset and series_onset != st1_onset:
                hard_conflicts.append("site_of_motor_onset")

        row["nygc_latest_hard_conflict_fields"] = ";".join(hard_conflicts) if hard_conflicts else "none"
        row["nygc_latest_diagnosis_detail_flags"] = ";".join(diagnosis_flags) if diagnosis_flags else "none"
        row["nygc_latest_warning_fields"] = ";".join(warnings) if warnings else "none"
        row["nygc_st1_current_alsfrs_r_numeric"] = numeric_value(row.get("nygc_st1_current_alsfrs_r"))
        row["nygc_st1_disease_duration_months_numeric"] = numeric_value(
            row.get("nygc_st1_disease_duration_in_months")
        )
        row["nygc_st1_age_at_symptom_onset_numeric"] = numeric_value(
            row.get("nygc_st1_age_at_symptom_onset")
        )
        row["nygc_st1_age_at_death_numeric"] = numeric_value(row.get("nygc_st1_age_at_death"))
        row["nygc_st1_fvc_percent_predicted_numeric"] = numeric_value(
            row.get("nygc_st1_fvc_percent_predicted")
        )
        add_resolved_fields(row, hard_conflicts, diagnosis_flags, warnings)

        if row["nygc_latest_direct_rna_id_match"] == "no":
            decision = "no_direct_nygc_st3_match"
        elif diagnosis_flags:
            decision = "direct_match_synapse_resolved_with_diagnosis_note"
        elif hard_conflicts:
            decision = "direct_match_synapse_resolved_with_geo_conflict_note"
        else:
            decision = "direct_match_synapse_resolved_candidate_for_endpoint_sensitivity"
        row["nygc_latest_endpoint_use_decision"] = decision
        records.append(row)

    return pd.DataFrame(records)


def endpoint_coverage(bridge: pd.DataFrame) -> pd.DataFrame:
    groups = {
        "all_gse330130_mo_libraries": bridge,
        "direct_nygc_st3_matches": bridge[bridge["nygc_latest_direct_rna_id_match"].eq("yes")],
        "direct_matches_without_geo_conflict_note": bridge[
            bridge["nygc_latest_direct_rna_id_match"].eq("yes")
            & bridge["nygc_latest_hard_conflict_fields"].eq("none")
            & bridge["nygc_latest_diagnosis_detail_flags"].eq("none")
        ],
        "candidate_endpoint_sensitivity": bridge[
            bridge["nygc_latest_endpoint_use_decision"].isin(
                [
                    "direct_match_synapse_resolved_candidate_for_endpoint_sensitivity",
                    "direct_match_synapse_resolved_with_diagnosis_note",
                    "direct_match_synapse_resolved_with_geo_conflict_note",
                ]
            )
        ],
    }
    fields = [
        ("disease_group", "resolved_disease_group"),
        ("als_ftd_detail", "resolved_als_ftd_detail"),
        ("disease_duration_months", "resolved_disease_duration_months"),
        ("alsfrs_r_current", "resolved_current_alsfrs_r"),
        ("age_at_symptom_onset", "resolved_age_at_symptom_onset"),
        ("site_of_motor_onset", "resolved_site_of_motor_onset"),
        ("age_at_death", "resolved_age_at_death"),
        ("cause_of_death_group", "resolved_cause_of_death_group"),
        ("als_or_mnd_related_death", "resolved_cause_of_death_related_to_als_or_mnd"),
        ("fvc_percent_predicted", "resolved_fvc_percent_predicted"),
    ]
    rows: list[dict[str, object]] = []
    for group_name, df in groups.items():
        for endpoint, col in fields:
            if col not in df.columns:
                usable = 0
            elif pd.api.types.is_numeric_dtype(df[col]):
                usable = int(df[col].notna().sum())
            else:
                usable = int(df[col].map(clean_text).astype(bool).sum())
            rows.append(
                {
                    "group": group_name,
                    "endpoint": endpoint,
                    "usable_libraries": usable,
                    "total_libraries": len(df),
                    "coverage_fraction": round(usable / len(df), 4) if len(df) else 0,
                }
            )
    return pd.DataFrame(rows)


def conflict_rows(bridge: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    field_pairs = {
        "broad_region": ("series_tissue", "nygc_st3_tissue_clean"),
        "sex": ("series_sex", "nygc_st1_sex_call"),
        "pmi": ("series_pmi", "nygc_st3_pmi"),
        "age_at_death": ("series_age_at_death", "nygc_st1_age_at_death"),
        "age_at_symptom_onset": (
            "series_age_at_symptom_onset",
            "nygc_st1_age_at_symptom_onset",
        ),
        "site_of_motor_onset": ("series_site_of_motor_onset", "nygc_st1_site_of_motor_onset"),
        "disease_duration_months": (
            SERIES_DURATION_ALIAS,
            "nygc_st1_disease_duration_in_months",
        ),
        "series_vs_st3_disease_group": ("series_genotype", "nygc_st3_disease_group"),
        "series_vs_st1_disease_group": ("series_genotype", "nygc_st1_disease_group"),
    }
    for _, row in bridge.iterrows():
        for field_type, fields in [
            ("hard_conflict", str(row["nygc_latest_hard_conflict_fields"])),
            ("diagnosis_detail_flag", str(row["nygc_latest_diagnosis_detail_flags"])),
        ]:
            if fields == "none":
                continue
            for field in fields.split(";"):
                series_col, nygc_col = field_pairs.get(field, ("", ""))
                records.append(
                    {
                        "mo_id": row.get("mo_id"),
                        "series_cgnd_id": row.get("series_cgnd_id"),
                        "field_type": field_type,
                        "field": field,
                        "series_value": row.get(series_col, ""),
                        "nygc_value": row.get(nygc_col, ""),
                        "decision": "use_synapse_resolved_value_with_note",
                        "resolved_conflict_note": row.get("resolved_conflict_note", ""),
                    }
                )
    return pd.DataFrame(records)


def summary_rows(bridge: pd.DataFrame, st1: pd.DataFrame, st3: pd.DataFrame, workbook: Path) -> pd.DataFrame:
    direct = bridge[bridge["nygc_latest_direct_rna_id_match"].eq("yes")]
    no_note = direct[
        direct["nygc_latest_hard_conflict_fields"].eq("none")
        & direct["nygc_latest_diagnosis_detail_flags"].eq("none")
    ]
    candidates = bridge[
        bridge["nygc_latest_endpoint_use_decision"].isin(
            [
                "direct_match_synapse_resolved_candidate_for_endpoint_sensitivity",
                "direct_match_synapse_resolved_with_diagnosis_note",
                "direct_match_synapse_resolved_with_geo_conflict_note",
            ]
        )
    ]
    rows = [
        (
            "nygc_latest_zenodo_record",
            "10.5281/zenodo.20090702",
            "Latest version found via Zenodo versions API; ST3 is populated.",
        ),
        ("nygc_latest_workbook_md5", md5(workbook), "Checksum for downloaded ST1/ST3 workbook."),
        ("nygc_st1_donor_rows", len(st1), "Donor-level clinical rows in latest ST1."),
        ("nygc_st3_rna_rows", len(st3), "RNA sample metadata rows in latest ST3."),
        ("gse330130_mo_libraries", len(bridge), "Collapsed MO libraries from GSE330130 GEO matrix."),
        (
            "gse330130_mo_libraries_with_exact_st3_rna_id_match",
            len(direct),
            "Exact match between series CGND-HRA ID and NYGC ST3 rna_id.",
        ),
        (
            "exact_st3_matches_without_geo_conflict_or_diagnosis_note",
            len(no_note),
            "Direct matches where GEO and Synapse/NYGC labels are concordant under the current audit checks.",
        ),
        (
            "endpoint_sensitivity_candidate_libraries",
            len(candidates),
            "Direct ST3 matches with resolved clinical fields from Synapse/NYGC; conflict/diagnosis notes stay visible.",
        ),
        (
            "direct_matches_with_diagnosis_detail_flag",
            int(direct["nygc_latest_diagnosis_detail_flags"].ne("none").sum()),
            "Mostly GEO ALS versus Synapse ALS-FTD/Other label differences; use Synapse labels and keep notes.",
        ),
        (
            "direct_matches_with_hard_conflict",
            int(direct["nygc_latest_hard_conflict_fields"].ne("none").sum()),
            "GEO-vs-Synapse clinical differences; Synapse resolved values are used with notes.",
        ),
        (
            "candidate_libraries_with_duration_months",
            int(candidates["resolved_disease_duration_months"].notna().sum()),
            "Disease duration can support descriptive/sensitivity analysis, not progression-rate claim.",
        ),
        (
            "candidate_libraries_with_current_alsfrs_r",
            int(candidates["resolved_current_alsfrs_r"].notna().sum()),
            "Sparse current ALSFRS-R; not enough for Brain main endpoint unless effect is strong.",
        ),
    ]
    return pd.DataFrame(rows, columns=["metric", "value", "interpretation"])


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.copy().where(pd.notna(df), "NA").to_csv(path, sep="\t", index=False)


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    if df.empty:
        return "No rows."
    safe = df.copy().where(pd.notna(df), "NA").astype(str)
    lines = [
        "| " + " | ".join(safe.columns) + " |",
        "| " + " | ".join(["---"] * len(safe.columns)) + " |",
    ]
    for _, row in safe.iterrows():
        lines.append("| " + " | ".join(row[col].replace("|", "/") for col in safe.columns) + " |")
    return "\n".join(lines)


def write_report(summary: pd.DataFrame, coverage: pd.DataFrame, conflicts: pd.DataFrame) -> Path:
    report = PLAN_ROOT / f"gse330130_nygc_latest_zenodo_bridge_audit_{DATE_TAG}.md"
    text = f"""# GSE330130 / NYGC latest Zenodo bridge audit ({DATE_TAG})

## Bottom line

The latest NYGC ALS Consortium Zenodo version (`10.5281/zenodo.20090702`) changes the endpoint-rescue picture. Unlike the older v1 workbook, latest ST3 is populated with 2574 RNA sample rows. This allows exact `CGND-HRA-*` RNA-ID matching for a larger subset of GSE330130/Takeuchi MO libraries.

This improves P01 clinical endpoint context, but it still does not turn P01 into a clinical progression, survival, or pathology-burden paper. For conflicting IDs, Synapse/NYGC curated metadata is treated as the authoritative source for resolved clinical and diagnosis fields. GEO labels remain raw sequencing provenance and are preserved as conflict notes.

## Summary

{markdown_table(summary)}

## Endpoint coverage

{markdown_table(coverage)}

## Conflict preview

{markdown_table(conflicts, max_rows=30)}

## P01 use boundary

- Use `resolved_*` fields for clinical/diagnosis interpretation. For exact ST3 RNA-ID matches, these fields prioritize Synapse/NYGC ST3/ST1 curated metadata.
- Keep raw `series_*` GEO fields as provenance. If GEO and Synapse differ, use the Synapse resolved value and report the conflict in `resolved_conflict_note`.
- Use latest NYGC ST1 duration, ALSFRS-R, FVC, cause-of-death and ALS/FTD detail initially only as descriptive or sensitivity context.
- Treat ALS vs ALS-FTD/Other/Other differences as biological diagnosis-detail notes; do not silently collapse them in phenotype figures.
- Do not claim progression rate, survival/time-to-event, pTDP-43 burden, motor neuron loss score, CD68/gliosis score, spinal atrophy score, STMN2 cryptic exon, UNC13A cryptic exon, ECAS/language/executive, or pathology endpoint rescue from this bridge alone.

## Machine-readable outputs

- `results/Supplemental_QC_covariate_audit/gse330130_nygc_latest_zenodo_bridge_{DATE_TAG}.tsv`
- `results/Supplemental_QC_covariate_audit/gse330130_nygc_latest_zenodo_bridge_conflicts_{DATE_TAG}.tsv`
- `results/Supplemental_QC_covariate_audit/gse330130_nygc_latest_zenodo_endpoint_coverage_{DATE_TAG}.tsv`
- `results/Supplemental_QC_covariate_audit/gse330130_nygc_latest_zenodo_bridge_summary_{DATE_TAG}.tsv`
"""
    report.write_text(text, encoding="utf-8")
    return report


def main() -> None:
    series_path = resolve_input(
        "GSE330130_SERIES_MATRIX",
        [
            source_data_root() / "GSE330130" / "GSE330130_series_matrix",
            source_data_root() / "GSE330130" / "GSE330130_series_matrix.txt",
        ],
    )
    workbook = resolve_input(
        "NYGC_ZENODO_ST1_ST3_XLSX",
        [
            RESCUE_ROOT / "zenodo_nygc_20090702" / "NYGC_Flagship_ST_1_3_Metadata.xlsx",
            # Same Zenodo 10.5281/zenodo.20090702 file at the standard source-data location.
            source_data_root() / "Other_Datasets" / "NYGC_ALS_Consortium" / "files" / "NYGC_Flagship_ST_1_3_Metadata.xlsx",
        ],
    )

    series = parse_series_libraries(series_path, "{GSE330130_SERIES_MATRIX}")
    st1, st3 = read_nygc_latest(workbook)
    bridge = build_bridge(series, st1, st3)
    coverage = endpoint_coverage(bridge)
    conflicts = conflict_rows(bridge)
    summary = summary_rows(bridge, st1, st3, workbook)

    OUT.mkdir(parents=True, exist_ok=True)
    write_tsv(bridge, OUT / f"gse330130_nygc_latest_zenodo_bridge_{DATE_TAG}.tsv")
    write_tsv(conflicts, OUT / f"gse330130_nygc_latest_zenodo_bridge_conflicts_{DATE_TAG}.tsv")
    write_tsv(coverage, OUT / f"gse330130_nygc_latest_zenodo_endpoint_coverage_{DATE_TAG}.tsv")
    write_tsv(summary, OUT / f"gse330130_nygc_latest_zenodo_bridge_summary_{DATE_TAG}.tsv")
    report = write_report(summary, coverage, conflicts)

    print(f"bridge rows={len(bridge)}")
    print(f"direct matches={bridge['nygc_latest_direct_rna_id_match'].eq('yes').sum()}")
    print(f"hard conflicts={bridge['nygc_latest_hard_conflict_fields'].ne('none').sum()}")
    print(f"conflict rows={len(conflicts)}")
    print(f"report={report}")


if __name__ == "__main__":
    main()
